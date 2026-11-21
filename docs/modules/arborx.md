# arborx — Recherche géométrique spatiale

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/arborx)

**Statut :** conception uniquement, aucune API ni fonctionnalité implémentée. Les noms Rust proposés sont indicatifs, pas des signatures stabilisées.

## Relations et contrats

- **Prérequis :** base et grid : points/boîtes ; MPI pour recherche distribuée.
- **Consommateurs :** grid, non_matching et particles.
- **Interface et invariants :** Trait de requêtes spatiales et résultats stables, avec opérations exactes de validation après filtrage ; ne pas exposer un backend particulier.
- **Dépendances externes :** `rstar` candidat Rust pour R-tree ; ArborX/Kokkos en FFI seulement après étude, non engagé. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P4 | `bvh.h, access_traits.h` | Indexation, intersections et plus proche voisin sur un backend abstrait. | Comparer à recherche exhaustive, objets vides, doublons et frontières. |
| P4 | `distributed_tree.h` | Prototype de recherche par rang avec résultats identifiés globalement. | Comparaison à recherche exhaustive globale, erreurs MPI et partitions vides. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.
