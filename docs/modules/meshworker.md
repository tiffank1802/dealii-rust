# meshworker — Parcours de maillage et assemblage

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/meshworker)

**Statut :** conception uniquement, aucune API ni fonctionnalité implémentée. Les noms Rust proposés sont indicatifs, pas des signatures stabilisées.

## Relations et contrats

- **Prérequis :** grid, fe, dofs, lac ; intégrateurs indépendants.
- **Consommateurs :** numerics, méthodes DG et assembleurs utilisateurs.
- **Interface et invariants :** Boucle cell/face avec ScratchData par tâche, CopyData local et insertion contrôlée ; erreurs propagées.
- **Dépendances externes :** P1 : standard ; `rayon` candidat pour P3 après mesure. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P1 | `mesh_loop.h, scratch_data.h, copy_data.h` | Exécuter contributions locales et insertion séquentielle avec ordre défini. | Accord avec boucle manuelle et propagation d’erreurs à la cellule concernée. |
| P2 | `assemble_flags.h, integration_info.h` | Traiter faces intérieures/frontières et sous-faces une seule fois selon convention. | Deux cellules DG, normales opposées, contributions et absence de double comptage. |
| P3 | `assembler.h, loop.h` | Construire contributions en parallèle puis réduire explicitement ; politique mémoire et ordre. | Comparaison séquentiel/parallèle avec tolérance, tests de courses et réduction stable. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.
