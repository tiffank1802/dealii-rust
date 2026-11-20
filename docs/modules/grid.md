# grid — Topologie, géométrie et cycle de vie des maillages

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/grid)

**Statut :** maillage triangulaire expérimental avec normales et raffinement uniforme, et maillages distincts de quadrilatères Q1 et tétraèdres P1 ; voir le [README](../../README.md). Les autres extensions et noms ci-dessous restent des propositions.

## Relations et contrats

- **Prérequis :** base : points et erreurs.
- **Consommateurs :** fe, dofs, numerics, hp, distributed, particles.
- **Interface et invariants :** Mesh2 avec VertexId/CellId, connectivité validée, marqueurs et emprunts en lecture ; invalider ou revalider les IDs après modifications.
- **Dépendances externes :** P0 : standard. `mshio` candidat pour importer, via `gmsh` ; aucun moteur géométrique obligatoire. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P0 | `tria.h, reference_cell.h, grid_generator.h` | Triangle orienté, cellules/sommets, arêtes et génération structurée du carré ; annoncer rejet des orientations inversées et dégénérescences. | Effectifs, aire totale 1, frontières, indices hors bornes et triangles dégénérés/inversés. |
| P1 | `tria_accessor.h, tria_orientation.h, grid_tools_topology.h` | Voisinage, faces et orientations ; étendre aux quadrangles et tétraèdres avec conventions locales explicites. | Deux cellules voisines, permutations et normales opposées ; aire/volume analytiques. |
| P1 | `grid_in.h, grid_out.h` | Marqueurs de frontière/matériau et interface vers I/O sans faire dépendre le noyau du parseur. | Préservation des marqueurs et rejet des cellules ou dimensions non supportées. |
| P2 | `grid_refinement.h, cell_id.h, persistent_tria.h` | Hiérarchie et raffinement/déraffinement ; versionner la topologie, invalider caches et IDs transitoires. | Relations parent/enfant, aire conservée, voisins non conformes et cycle raffinement/déraffinement. |
| P4 | `manifold*.h, grid_tools_cache.h` | Projection courbe et recherche spatiale, caches avec règles d’invalidation après mutation. | Projection sur cercle/sphère, invalidation après raffinement, point hors domaine. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.

## Détail d’implémentation P0

- `mesh` stocke les points et les triplets d’indices, sans références internes autoréférentielles. `CellId`/`VertexId` ne sont pas interchangeables. Rejeter sommets manquants, répétitions, triangles d’aire nulle et orientation négative ; préciser la tolérance à l’échelle du domaine.
- `generator` : pour `n ≥ 1` divisions par côté, générer `(n+1)²` sommets et `2n²` triangles en orientation positive ; `n=0` est une erreur. Marquer les arêtes externes et en déduire les sommets de bord, sans classer les diagonales comme frontières.
- `topology` : compter les occurrences d’arêtes non orientées ; une arête vue une fois est au bord, deux fois est intérieure. Rejeter une arête vue plus de deux fois pour le maillage conforme P0.
- **Vérifications :** `n=1`, `n=2`, aire totale 1, connectivité, numéros hors plage, triangle inversé et maillage non manifold. Reporter les tests d’adaptation à P2 ; aucune stabilité des IDs à travers mutation n’est promise en P0.

**État :** les groupes physiques du sous-ensemble Gmsh importé sont conservés dans `gmsh::ImportedMesh`, pas dans `Mesh2` ; voir le [README](../../README.md#import-gmsh-et-résolution-sur-maillage-étiqueté). Les marqueurs généraux modifiables restent à concevoir.

Le raffinement uniforme `UniformRefinement` conserve les indices des anciens sommets, expose les parents des cellules et arêtes de frontière et permet le prolongement P1. Les marqueurs ne sont pas propagés automatiquement. Pas de raffinement local, déraffinement ou hiérarchie persistante.

**Quadrilatères :** `quad::QuadMesh` valide cellules convexes orientées et topologie d’arêtes, mais ne détecte pas les intersections et nœuds pendants. Les API Gmsh et de raffinement triangulaires ne s’appliquent pas.

**Tétraèdres :** `tetra::TetMesh` valide l’orientation, la non-dégénérescence et les faces partagées ; il fournit faces de bord orientées et normales sortantes. Les intersections, le raffinement et l’import Gmsh 3D restent hors périmètre.
