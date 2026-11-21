# multigrid — Transferts, lisseurs et cycles multiniveaux

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/multigrid)

**Statut :** conception uniquement, aucune API ni fonctionnalité implémentée. Les noms Rust proposés sont indicatifs, pas des signatures stabilisées.

## Relations et contrats

- **Prérequis :** grid hiérarchique, dofs de niveau, lac et contraintes.
- **Consommateurs :** Préconditionnement des solveurs, matrix_free.
- **Interface et invariants :** Level, Transfer, Smoother et CoarseSolver ; cycle V explicite, état de chaque niveau et contraintes.
- **Dépendances externes :** Aucune nouvelle crate obligatoire ; backend creux de lac et MPI optionnel. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P2 | `mg_transfer.h, mg_constrained_dofs.h` | Construire prolongation/restriction et contraintes par niveau. | Reproduction de constantes/affines et propriétés de restriction annoncées. |
| P2 | `mg_smoother.h, mg_coarse.h, multigrid.h` | Lisseur Jacobi, solveur grossier et cycle V comme préconditionneur. | Résidu en baisse ; solution comparée à un solveur direct et itérations sous raffinement. |
| P3 | `mg_transfer_global_coarsening.h, mg_transfer_matrix_free.h` | Transferts distribués et opérateurs sans matrice après validation locale. | Résultats par rang comparés au séquentiel et transfert avec cellules fantômes. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.
