# non_matching — Géométrie immergée et couplage de maillages

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/non_matching)

**Statut :** conception uniquement, aucune API ni fonctionnalité implémentée. Les noms Rust proposés sont indicatifs, pas des signatures stabilisées.

## Relations et contrats

- **Prérequis :** grid, fe, dofs, lac et recherche spatiale.
- **Consommateurs :** Applications multiphysiques, numerics et distributed.
- **Interface et invariants :** Classifier chaque cellule, produire quadratures immergées et couples de DoFs ; documenter dimensions et comportement hors domaine.
- **Dépendances externes :** `rstar` candidat pour recherche ; géométrie CGAL optionnelle, non choisie pour le cœur. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P4 | `mesh_classifier.h, quadrature_generator.h` | Classer intérieur/extérieur/coupé et calculer quadrature sur cellules coupées. | Aires/volumes analytiques, coupe tangentielle et cas quasi dégénérés. |
| P4 | `fe_immersed_values.h, mapping_info.h` | Évaluer fonctions et gradients sur quadratures immergées. | Intégrales de constantes/linéaires et normes de surface connues. |
| P4 | `coupling.h, closest_surface_point.h` | Détecter recouvrement et assembler couplage entre domaines et composantes. | Conservation, point hors domaine et comparaison de deux maillages analytiques. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.
