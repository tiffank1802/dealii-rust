//! Strict Gmsh 2.2 ASCII reader for planar first-order triangle meshes.
use std::collections::{BTreeMap, BTreeSet};
use std::io::Read;

use super::{physical_names, positive_tag, required, split_sections, ImportedMesh, Tokens};
use crate::base::Point2;
use crate::grid::{BoundaryEdgeId, Mesh2, VertexId};
use crate::{error, Error, Result};

/// Read the supported subset of Gmsh MSH 2.2 ASCII text from a stream.
/// See [`parse_msh22`] for supported elements, tags and limitations.
pub fn read_msh22<R: Read>(mut input: R) -> Result<ImportedMesh> {
    let mut content = String::new();
    input
        .read_to_string(&mut content)
        .map_err(|e| Error(format!("Gmsh read error: {}", e)))?;
    parse_msh22(&content)
}

/// Parse a strict subset of Gmsh MSH 2.2 ASCII (`2.2 0 8`). Supports planar
/// non-parametric nodes (`z=0`) and first-order point (15), line (1) and
/// triangle (2) elements. Triangles must be counterclockwise, and listed
/// lines must be exterior mesh edges. `PhysicalNames` is optional.
///
/// Element records may contain 0, 1 or 2 tags: the first is the single
/// physical tag (0 means none), the second is the elementary entity tag
/// (parsed but not retained). More tags, including partitions, are rejected.
/// Physical tags of point elements are associated to their node. Point and
/// line elements may be omitted; repeated elements, unknown nodes, 3D,
/// binary files, other elements and extra sections are rejected. The format
/// does not preserve multiple physical memberships per element. `node_tag`,
/// cell and boundary element tags retain their original values; geometric
/// entity tags and point-element IDs are not exposed. All input is read into
/// memory; global mesh conformity and absence of intersections are assumed.
pub fn parse_msh22(input: &str) -> Result<ImportedMesh> {
    let sections = split_sections(input)?;
    if sections.contains_key("Entities") {
        return Err(error("Entities: unsupported in MSH 2.2"));
    }
    let mut format = required(&sections, "MeshFormat")?.split_whitespace();
    if format.next() != Some("2.2")
        || format.next() != Some("0")
        || format.next() != Some("8")
        || format.next().is_some()
    {
        return Err(error(
            "only Gmsh MSH 2.2 ASCII with data size 8 is supported",
        ));
    }
    let physical_names = match sections.get("PhysicalNames") {
        Some(text) => physical_names(text)?,
        None => BTreeMap::new(),
    };
    let mut nodes = required(&sections, "Nodes")?
        .lines()
        .filter(|line| !line.trim().is_empty());
    let node_count: usize = nodes
        .next()
        .ok_or_else(|| error("Nodes: missing count"))?
        .trim()
        .parse()
        .map_err(|_| error("Nodes: invalid count"))?;
    let mut points = Vec::new();
    let mut node_tags = Vec::new();
    let mut node_index = BTreeMap::new();
    for _ in 0..node_count {
        let line = nodes
            .next()
            .ok_or_else(|| error("Nodes: truncated entries"))?;
        let mut t = Tokens::new(line, "Nodes");
        let tag = positive_tag(&mut t, "node tag")?;
        let x = t.finite("x coordinate")?;
        let y = t.finite("y coordinate")?;
        let z = t.finite("z coordinate")?;
        t.end()?;
        if z != 0.0 {
            return Err(error("Nodes: only planar z=0 geometry is supported"));
        }
        if node_index.insert(tag, VertexId(points.len())).is_some() {
            return Err(error("Nodes: duplicate node tag"));
        }
        node_tags.push(tag);
        points.push(Point2::new(x, y)?);
    }
    if nodes.next().is_some() {
        return Err(error("Nodes: unexpected extra entries"));
    }
    let mut entries = required(&sections, "Elements")?
        .lines()
        .filter(|line| !line.trim().is_empty());
    let count: usize = entries
        .next()
        .ok_or_else(|| error("Elements: missing count"))?
        .trim()
        .parse()
        .map_err(|_| error("Elements: invalid count"))?;
    let mut element_ids = BTreeSet::new();
    let mut point_elements = BTreeSet::new();
    let mut vertex_groups = vec![Vec::new(); points.len()];
    let mut cells = Vec::new();
    let mut cell_groups = Vec::new();
    let mut cell_element_tags = Vec::new();
    let mut lines = BTreeMap::new();
    for _ in 0..count {
        let line = entries
            .next()
            .ok_or_else(|| error("Elements: truncated entries"))?;
        let mut t = Tokens::new(line, "Elements");
        let id = positive_tag(&mut t, "element tag")?;
        if !element_ids.insert(id) {
            return Err(error("Elements: duplicate element tag"));
        }
        let kind: u32 = t.take("element type")?;
        let ntags: usize = t.take("number of tags")?;
        if ntags > 2 {
            return Err(error(
                "Elements: partitions and more than two tags are unsupported",
            ));
        }
        let physical = if ntags > 0 {
            t.take::<u64>("physical tag")?
        } else {
            0
        };
        if ntags == 2 {
            let _: u64 = t.take("elementary entity tag")?;
        }
        let groups = if physical == 0 {
            Vec::new()
        } else {
            vec![physical]
        };
        let n = match kind {
            15 => 1,
            1 => 2,
            2 => 3,
            _ => {
                return Err(Error(format!(
                    "Elements: unsupported element type {}",
                    kind
                )))
            }
        };
        let mut vertex_ids = Vec::new();
        for _ in 0..n {
            let tag = positive_tag(&mut t, "element node tag")?;
            vertex_ids.push(
                *node_index
                    .get(&tag)
                    .ok_or_else(|| Error(format!("Elements: unknown node {}", tag)))?,
            );
        }
        t.end()?;
        match kind {
            15 => {
                let v = vertex_ids[0];
                if !point_elements.insert(v) {
                    return Err(error("Elements: duplicate point element"));
                }
                vertex_groups[v.0] = groups;
            }
            1 => {
                let a = vertex_ids[0];
                let b = vertex_ids[1];
                if a == b {
                    return Err(error("Elements: degenerate line element"));
                }
                let edge = if a < b { (a, b) } else { (b, a) };
                if lines.insert(edge, (id, groups)).is_some() {
                    return Err(error("Elements: duplicate line element"));
                }
            }
            2 => {
                cells.push([vertex_ids[0], vertex_ids[1], vertex_ids[2]]);
                cell_groups.push(groups);
                cell_element_tags.push(id);
            }
            _ => unreachable!(),
        }
    }
    if entries.next().is_some() {
        return Err(error("Elements: unexpected extra entries"));
    }
    let mesh = Mesh2::new(points, cells)
        .map_err(|e| Error(format!("invalid Gmsh triangle mesh: {}", e)))?;
    let mut boundary_groups = vec![Vec::new(); mesh.boundary_edges().len()];
    let mut boundary_element_tags = vec![None; mesh.boundary_edges().len()];
    let lookup: BTreeMap<_, _> = mesh
        .boundary_edges()
        .iter()
        .enumerate()
        .map(|(i, e)| ((e[0], e[1]), BoundaryEdgeId(i)))
        .collect();
    for (edge, (id, groups)) in lines {
        let index = lookup
            .get(&edge)
            .ok_or_else(|| error("Gmsh line element is not an exterior mesh edge"))?
            .0;
        boundary_element_tags[index] = Some(id);
        boundary_groups[index] = groups;
    }
    Ok(ImportedMesh {
        mesh,
        node_tags,
        vertex_groups,
        cell_groups,
        cell_element_tags,
        boundary_groups,
        boundary_element_tags,
        physical_names,
    })
}
