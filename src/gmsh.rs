//! Import deliberately limited, validated subsets of Gmsh MSH 4.1 and 2.2 ASCII.
//!
//! Supported: planar nodes with z=0, non-parametric node blocks, point, line
//! and linear triangle elements, geometric entities, physical tags and names.
//! Boundary line elements must correspond to external mesh edges. Partitions,
//! periodic constraints, fields, binary files, 3D and other element types are
//! rejected rather than silently ignored. Triangles must be counterclockwise.
//! [`parse`] and [`read`] support 4.1; [`parse_msh22`] and [`read_msh22`]
//! support a separate strict 2.2 subset and preserve at most one physical tag
//! per element. See their documentation for the differences.
use std::collections::{BTreeMap, BTreeSet};
use std::io::Read;

use crate::base::Point2;
use crate::grid::{BoundaryEdgeId, CellId, Mesh2, VertexId};
use crate::{error, Error, Result};

#[path = "gmsh22.rs"]
mod legacy;
pub use legacy::{parse_msh22, read_msh22};

/// A mesh and physical-group metadata aligned with its vertex, cell and
/// boundary-edge IDs. IDs must be used with the mesh in this object only.
#[derive(Debug, Clone)]
pub struct ImportedMesh {
    mesh: Mesh2,
    node_tags: Vec<u64>,
    vertex_groups: Vec<Vec<u64>>,
    cell_groups: Vec<Vec<u64>>,
    cell_element_tags: Vec<u64>,
    boundary_groups: Vec<Vec<u64>>,
    boundary_element_tags: Vec<Option<u64>>,
    physical_names: BTreeMap<(u8, u64), String>,
}

impl ImportedMesh {
    pub fn mesh(&self) -> &Mesh2 {
        &self.mesh
    }
    /// Original Gmsh node tag for this vertex; tags need not be contiguous.
    pub fn node_tag(&self, id: VertexId) -> Result<u64> {
        self.node_tags
            .get(id.0)
            .copied()
            .ok_or_else(|| error("vertex index out of bounds"))
    }
    /// Physical tags from a point entity (4.1) or point element at this vertex.
    pub fn vertex_physical_tags(&self, id: VertexId) -> Result<&[u64]> {
        self.vertex_groups
            .get(id.0)
            .map(Vec::as_slice)
            .ok_or_else(|| error("vertex index out of bounds"))
    }
    /// Physical tags of the triangle: its entity groups in 4.1, or its
    /// element physical tag in the supported 2.2 subset.
    pub fn cell_physical_tags(&self, id: CellId) -> Result<&[u64]> {
        self.cell_groups
            .get(id.0)
            .map(Vec::as_slice)
            .ok_or_else(|| error("cell index out of bounds"))
    }
    /// Original Gmsh element tag for a triangle.
    pub fn cell_element_tag(&self, id: CellId) -> Result<u64> {
        self.cell_element_tags
            .get(id.0)
            .copied()
            .ok_or_else(|| error("cell index out of bounds"))
    }
    /// Original Gmsh line-element tag, or `None` if the boundary edge was not listed.
    pub fn boundary_element_tag(&self, id: BoundaryEdgeId) -> Result<Option<u64>> {
        self.boundary_element_tags
            .get(id.0)
            .copied()
            .ok_or_else(|| error("boundary edge index out of bounds"))
    }
    /// Physical tags of the corresponding boundary line element. An empty
    /// list can mean either no line element or an untagged line element; use
    /// `boundary_element_tag` to distinguish those cases.
    pub fn boundary_physical_tags(&self, id: BoundaryEdgeId) -> Result<&[u64]> {
        self.boundary_groups
            .get(id.0)
            .map(Vec::as_slice)
            .ok_or_else(|| error("boundary edge index out of bounds"))
    }
    /// Find boundary edges belonging to a physical curve group.
    pub fn boundary_edges_with_tag(&self, tag: u64) -> Vec<BoundaryEdgeId> {
        self.boundary_groups
            .iter()
            .enumerate()
            .filter_map(|(i, tags)| {
                if tags.contains(&tag) {
                    Some(BoundaryEdgeId(i))
                } else {
                    None
                }
            })
            .collect()
    }
    /// Find triangles belonging to a physical surface group.
    pub fn cells_with_tag(&self, tag: u64) -> Vec<CellId> {
        self.cell_groups
            .iter()
            .enumerate()
            .filter_map(|(i, groups)| {
                if groups.contains(&tag) {
                    Some(CellId(i))
                } else {
                    None
                }
            })
            .collect()
    }
    /// Unique vertices of boundary line elements belonging to a physical group.
    pub fn boundary_vertices_with_tag(&self, tag: u64) -> Vec<VertexId> {
        self.boundary_edges_with_tag(tag)
            .into_iter()
            .flat_map(|id| self.mesh.boundary_edges()[id.0].iter().copied())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect()
    }
    /// Name from `$PhysicalNames`, if one was declared for dimension and tag.
    pub fn physical_name(&self, dimension: u8, tag: u64) -> Option<&str> {
        self.physical_names
            .get(&(dimension, tag))
            .map(String::as_str)
    }
}

/// Read MSH 4.1 ASCII text. Files with unsupported sections, cell types or
/// geometry are rejected; see the module documentation for the subset.
pub fn read<R: Read>(mut input: R) -> Result<ImportedMesh> {
    let mut content = String::new();
    input
        .read_to_string(&mut content)
        .map_err(|e| Error(format!("Gmsh read error: {}", e)))?;
    parse(&content)
}

/// Parse MSH 4.1 ASCII text.
pub fn parse(input: &str) -> Result<ImportedMesh> {
    let sections = split_sections(input)?;
    let format = required(&sections, "MeshFormat")?;
    let mut format_tokens = format.split_whitespace();
    if format_tokens.next() != Some("4.1")
        || format_tokens.next() != Some("0")
        || format_tokens.next() != Some("8")
        || format_tokens.next().is_some()
    {
        return Err(error(
            "only Gmsh MSH 4.1 ASCII with data size 8 is supported",
        ));
    }
    let names = match sections.get("PhysicalNames") {
        Some(text) => physical_names(text)?,
        None => BTreeMap::new(),
    };
    let entities = parse_entities(required(&sections, "Entities")?)?;
    let nodes = parse_nodes(required(&sections, "Nodes")?, &entities)?;
    let elements = parse_elements(required(&sections, "Elements")?, &entities, &nodes)?;
    let mesh = Mesh2::new(nodes.points, elements.cells)
        .map_err(|e| Error(format!("invalid Gmsh triangle mesh: {}", e)))?;
    let mut boundary_groups = vec![Vec::new(); mesh.boundary_edges().len()];
    let mut boundary_element_tags = vec![None; mesh.boundary_edges().len()];
    let boundary_lookup: BTreeMap<_, _> = mesh
        .boundary_edges()
        .iter()
        .enumerate()
        .map(|(i, edge)| ((edge[0], edge[1]), i))
        .collect();
    for (edge, (tag, groups)) in elements.lines {
        let index = boundary_lookup
            .get(&edge)
            .ok_or_else(|| error("Gmsh line element is not an exterior mesh edge"))?;
        boundary_groups[*index] = groups;
        boundary_element_tags[*index] = Some(tag);
    }
    Ok(ImportedMesh {
        mesh,
        node_tags: nodes.tags,
        vertex_groups: elements.vertex_groups,
        cell_groups: elements.cell_groups,
        cell_element_tags: elements.cell_element_tags,
        boundary_groups,
        boundary_element_tags,
        physical_names: names,
    })
}

fn split_sections(input: &str) -> Result<BTreeMap<String, String>> {
    let mut lines = input.lines();
    let mut result = BTreeMap::new();
    while let Some(line) = lines.next() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let name = line
            .strip_prefix('$')
            .ok_or_else(|| error("expected Gmsh section header"))?;
        if ![
            "MeshFormat",
            "PhysicalNames",
            "Entities",
            "Nodes",
            "Elements",
        ]
        .contains(&name)
        {
            return Err(Error(format!("unsupported Gmsh section ${}", name)));
        }
        if result.contains_key(name) {
            return Err(Error(format!("duplicate Gmsh section ${}", name)));
        }
        let end = format!("$End{}", name);
        let mut body = String::new();
        let mut found_end = false;
        for item in lines.by_ref() {
            let text = item.trim();
            if text == end {
                found_end = true;
                break;
            }
            if text.starts_with('$') {
                return Err(Error(format!("missing {} before {}", end, text)));
            }
            body.push_str(item);
            body.push('\n');
        }
        if !found_end {
            return Err(Error(format!("missing {}", end)));
        }
        result.insert(name.to_string(), body);
    }
    Ok(result)
}
fn required<'a>(sections: &'a BTreeMap<String, String>, name: &str) -> Result<&'a str> {
    sections
        .get(name)
        .map(String::as_str)
        .ok_or_else(|| Error(format!("missing Gmsh section ${}", name)))
}

struct Tokens<'a> {
    values: std::str::SplitWhitespace<'a>,
    section: &'static str,
}
impl<'a> Tokens<'a> {
    fn new(input: &'a str, section: &'static str) -> Self {
        Self {
            values: input.split_whitespace(),
            section,
        }
    }
    fn take<T: std::str::FromStr>(&mut self, what: &str) -> Result<T> {
        let text = self
            .values
            .next()
            .ok_or_else(|| Error(format!("{}: missing {}", self.section, what)))?;
        text.parse()
            .map_err(|_| Error(format!("{}: invalid {}: {}", self.section, what, text)))
    }
    fn finite(&mut self, what: &str) -> Result<f64> {
        let value: f64 = self.take(what)?;
        if !value.is_finite() {
            return Err(Error(format!("{}: nonfinite {}", self.section, what)));
        }
        Ok(value)
    }
    fn end(&mut self) -> Result<()> {
        if self.values.next().is_some() {
            return Err(Error(format!("{}: unexpected extra data", self.section)));
        }
        Ok(())
    }
}

fn positive_tag(tokens: &mut Tokens<'_>, label: &str) -> Result<u64> {
    let tag: u64 = tokens.take(label)?;
    if tag == 0 {
        return Err(Error(format!(
            "{}: {} must be positive",
            tokens.section, label
        )));
    }
    Ok(tag)
}

fn physical_names(text: &str) -> Result<BTreeMap<(u8, u64), String>> {
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let count: usize = lines
        .next()
        .ok_or_else(|| error("PhysicalNames: missing count"))?
        .trim()
        .parse()
        .map_err(|_| error("PhysicalNames: invalid count"))?;
    let mut result = BTreeMap::new();
    for _ in 0..count {
        let line = lines
            .next()
            .ok_or_else(|| error("PhysicalNames: truncated entries"))?
            .trim();
        let quote = line
            .find('"')
            .ok_or_else(|| error("PhysicalNames: expected quoted name"))?;
        let mut prefix = Tokens::new(&line[..quote], "PhysicalNames");
        let dim: u8 = prefix.take("dimension")?;
        let tag: u64 = prefix.take("tag")?;
        prefix.end()?;
        if dim > 2 {
            return Err(error(
                "PhysicalNames: only dimensions 0, 1 and 2 are supported",
            ));
        }
        let quoted = &line[quote..];
        if quoted.len() < 2 || !quoted.ends_with('"') {
            return Err(error("PhysicalNames: unterminated name"));
        }
        let mut decoded = String::new();
        let mut chars = quoted[1..quoted.len() - 1].chars();
        while let Some(ch) = chars.next() {
            match ch {
                '"' => return Err(error("PhysicalNames: unescaped quote")),
                '\\' => match chars.next() {
                    Some('"') => decoded.push('"'),
                    Some('\\') => decoded.push('\\'),
                    _ => return Err(error("PhysicalNames: unsupported escape")),
                },
                _ => decoded.push(ch),
            }
        }
        if result.insert((dim, tag), decoded).is_some() {
            return Err(error("PhysicalNames: duplicate dimension/tag"));
        }
    }
    if lines.next().is_some() {
        return Err(error("PhysicalNames: unexpected extra entries"));
    }
    Ok(result)
}

type Entities = BTreeMap<(u8, u64), Vec<u64>>;
fn parse_entities(text: &str) -> Result<Entities> {
    let mut t = Tokens::new(text, "Entities");
    let counts = [
        t.take::<usize>("point count")?,
        t.take::<usize>("curve count")?,
        t.take::<usize>("surface count")?,
        t.take::<usize>("volume count")?,
    ];
    if counts[3] != 0 {
        return Err(error("Entities: 3D volumes are unsupported"));
    }
    let mut entities = Entities::new();
    let mut boundaries = Vec::new();
    for dim in 0..=2u8 {
        for _ in 0..counts[dim as usize] {
            let tag = positive_tag(&mut t, "entity tag")?;
            let ncoords = if dim == 0 { 3 } else { 6 };
            for _ in 0..ncoords {
                t.finite("entity coordinate")?;
            }
            let n: usize = t.take("physical tag count")?;
            let mut groups = Vec::new();
            for _ in 0..n {
                let group: u64 = t.take("physical tag")?;
                if groups.contains(&group) {
                    return Err(error("Entities: duplicate physical tag"));
                }
                groups.push(group);
            }
            if entities.insert((dim, tag), groups).is_some() {
                return Err(error("Entities: duplicate dimension/entity tag"));
            }
            if dim > 0 {
                let n: usize = t.take("bounding entity count")?;
                for _ in 0..n {
                    let reference: i64 = t.take("bounding entity tag")?;
                    if reference == 0 {
                        return Err(error("Entities: zero bounding entity tag"));
                    }
                    boundaries.push((dim - 1, reference.unsigned_abs()));
                }
            }
        }
    }
    t.end()?;
    for (dim, tag) in boundaries {
        if !entities.contains_key(&(dim, tag)) {
            return Err(error("Entities: unknown bounding entity"));
        }
    }
    Ok(entities)
}

struct Nodes {
    points: Vec<Point2>,
    tags: Vec<u64>,
    lookup: BTreeMap<u64, VertexId>,
    vertex_groups: Vec<Vec<u64>>,
    node_entities: Vec<(u8, u64)>,
}
fn check_range(
    label: &str,
    min: u64,
    max: u64,
    tags: &BTreeSet<u64>,
    expected: usize,
) -> Result<()> {
    if tags.len() != expected {
        return Err(Error(format!(
            "{}: declared count does not match entries",
            label
        )));
    }
    let actual_min = tags.iter().next().copied().unwrap_or(0);
    let actual_max = tags.iter().next_back().copied().unwrap_or(0);
    if min != actual_min || max != actual_max {
        return Err(Error(format!(
            "{}: declared tag range does not match entries",
            label
        )));
    }
    Ok(())
}
fn parse_nodes(text: &str, entities: &Entities) -> Result<Nodes> {
    let mut t = Tokens::new(text, "Nodes");
    let blocks: usize = t.take("block count")?;
    let expected: usize = t.take("node count")?;
    let min: u64 = t.take("minimum node tag")?;
    let max: u64 = t.take("maximum node tag")?;
    let mut points = Vec::new();
    let mut tags = Vec::new();
    let mut lookup = BTreeMap::new();
    let mut all = BTreeSet::new();
    let mut vertex_groups = Vec::new();
    let mut node_entities = Vec::new();
    for _ in 0..blocks {
        let dim: u8 = t.take("entity dimension")?;
        let entity = positive_tag(&mut t, "entity tag")?;
        let parametric: u8 = t.take("parametric flag")?;
        let count: usize = t.take("block node count")?;
        if parametric != 0 {
            return Err(error("Nodes: parametric nodes are unsupported"));
        }
        let physical = entities
            .get(&(dim, entity))
            .ok_or_else(|| error("Nodes: unknown entity"))?;
        if count > expected.saturating_sub(tags.len()) {
            return Err(error("Nodes: block exceeds declared node count"));
        }
        let mut block_tags = Vec::new();
        for _ in 0..count {
            let tag = positive_tag(&mut t, "node tag")?;
            if !all.insert(tag) {
                return Err(error("Nodes: duplicate node tag"));
            }
            block_tags.push(tag);
        }
        for tag in block_tags {
            let x = t.finite("x coordinate")?;
            let y = t.finite("y coordinate")?;
            let z = t.finite("z coordinate")?;
            if z != 0.0 {
                return Err(error("Nodes: only planar z=0 geometry is supported"));
            }
            lookup.insert(tag, VertexId(points.len()));
            points.push(Point2::new(x, y)?);
            tags.push(tag);
            vertex_groups.push(if dim == 0 {
                physical.clone()
            } else {
                Vec::new()
            });
            node_entities.push((dim, entity));
        }
    }
    t.end()?;
    check_range("Nodes", min, max, &all, expected)?;
    Ok(Nodes {
        points,
        tags,
        lookup,
        vertex_groups,
        node_entities,
    })
}

struct Elements {
    cells: Vec<[VertexId; 3]>,
    cell_groups: Vec<Vec<u64>>,
    cell_element_tags: Vec<u64>,
    lines: BTreeMap<(VertexId, VertexId), (u64, Vec<u64>)>,
    vertex_groups: Vec<Vec<u64>>,
}
fn parse_elements(text: &str, entities: &Entities, nodes: &Nodes) -> Result<Elements> {
    let mut t = Tokens::new(text, "Elements");
    let blocks: usize = t.take("block count")?;
    let expected: usize = t.take("element count")?;
    let min: u64 = t.take("minimum element tag")?;
    let max: u64 = t.take("maximum element tag")?;
    let mut all = BTreeSet::new();
    let mut cells = Vec::new();
    let mut cell_groups = Vec::new();
    let mut cell_element_tags = Vec::new();
    let mut lines = BTreeMap::new();
    let mut vertex_groups = nodes.vertex_groups.clone();
    let mut point_elements = BTreeSet::new();
    for _ in 0..blocks {
        let dim: u8 = t.take("entity dimension")?;
        let entity = positive_tag(&mut t, "entity tag")?;
        let kind: u32 = t.take("element type")?;
        let count: usize = t.take("block element count")?;
        let groups = entities
            .get(&(dim, entity))
            .ok_or_else(|| error("Elements: unknown entity"))?;
        let nnodes = match (dim, kind) {
            (0, 15) => 1,
            (1, 1) => 2,
            (2, 2) => 3,
            _ => {
                return Err(Error(format!(
                    "Elements: unsupported dimension {} / element type {}",
                    dim, kind
                )))
            }
        };
        if count > expected.saturating_sub(all.len()) {
            return Err(error("Elements: block exceeds declared element count"));
        }
        for _ in 0..count {
            let tag = positive_tag(&mut t, "element tag")?;
            if !all.insert(tag) {
                return Err(error("Elements: duplicate element tag"));
            }
            let mut vertices = Vec::new();
            for _ in 0..nnodes {
                let node = positive_tag(&mut t, "element node tag")?;
                vertices.push(
                    *nodes
                        .lookup
                        .get(&node)
                        .ok_or_else(|| Error(format!("Elements: unknown node {}", node)))?,
                );
            }
            match (dim, kind) {
                (0, 15) => {
                    let vertex = vertices[0];
                    if nodes.node_entities[vertex.0] != (dim, entity) {
                        return Err(error("Elements: point node belongs to a different entity"));
                    }
                    if !point_elements.insert(vertex) {
                        return Err(error("Elements: duplicate point element"));
                    }
                    for &group in groups {
                        if !vertex_groups[vertex.0].contains(&group) {
                            vertex_groups[vertex.0].push(group);
                        }
                    }
                }
                (1, 1) => {
                    let a = vertices[0];
                    let b = vertices[1];
                    if a == b {
                        return Err(error("Elements: degenerate line element"));
                    }
                    let key = if a < b { (a, b) } else { (b, a) };
                    if lines.insert(key, (tag, groups.clone())).is_some() {
                        return Err(error("Elements: duplicate line element"));
                    }
                }
                (2, 2) => {
                    cells.push([vertices[0], vertices[1], vertices[2]]);
                    cell_groups.push(groups.clone());
                    cell_element_tags.push(tag);
                }
                _ => unreachable!(),
            }
        }
    }
    t.end()?;
    check_range("Elements", min, max, &all, expected)?;
    Ok(Elements {
        cells,
        cell_groups,
        cell_element_tags,
        lines,
        vertex_groups,
    })
}
