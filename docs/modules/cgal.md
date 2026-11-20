# cgal — Algorithmes de géométrie et maillages de surface

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/cgal)

**Statut :** conception uniquement, aucune API ni fonctionnalité implémentée. Les noms Rust proposés sont indicatifs, pas des signatures stabilisées.

## Relations et contrats

- **Prérequis :** base et grid ; non_matching selon usage.
- **Consommateurs :** Géométries complexes et adaptateurs CAD.
- **Interface et invariants :** Frontière d’intégration isolée : conversions validées, erreurs natives et mémoire gérée par RAII.
- **Dépendances externes :** CGAL backend natif candidat ; étudier alternatives Rust, licences des composants et coût de FFI avant décision. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P4 | `point_conversion.h, polygon.h` | Convertir points/polygones et normaliser les orientations. | Aller-retour, polygone dégénéré et frontières. |
| P4 | `intersections.h, triangulation.h, surface_mesh.h` | Prototype intersections et triangulation de surface. | Cas analytiques, intersections tangentes et validité du maillage ; erreurs du backend. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.
