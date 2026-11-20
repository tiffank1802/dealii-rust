# distributed — Maillages, transferts et données distribués

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/distributed)

**Statut :** conception uniquement, aucune API ni fonctionnalité implémentée. Les noms Rust proposés sont indicatifs, pas des signatures stabilisées.

## Relations et contrats

- **Prérequis :** grid, dofs et lac distribuables ; MPI sélectionné séparément.
- **Consommateurs :** matrix_free, multigrid, numerics et particles.
- **Interface et invariants :** Communicateur explicitement fourni, cellules owned/ghost/artificial, synchronisation et époques de modification documentées.
- **Dépendances externes :** `mpi` candidat pour MPI ; p4est backend natif optionnel à évaluer, pas de crate de binding choisie. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P3 | `tria_base.h, shared_tria.h, tria.h` | Propriété, fantômes et échanges pour maillage partagé puis partitionné. | Propriétaire unique, voisins entre rangs, accord avec séquentiel sur 1/2/N rangs. |
| P3 | `grid_refinement.h, repartitioning_policy_tools.h` | Marquage global, raffinement et repartitionnement ; invalidation des caches. | Conservation aire/volume, voisins et indices après répartition déséquilibrée. |
| P3 | `cell_data_transfer.h, field_transfer.h` | Transférer données de cellule et champs entre partitions/maillages. | Conservation des constantes et contrôle d’identités sur migration. |
| P4 | `fully_distributed_tria.h, p4est_wrappers.h` | Prototype d’un backend de forêt et construction totalement distribuée. | Scaling mémoire/temps et échecs d’initialisation ; tests multi-rangs reproductibles. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.

## Protocole de données distribuées

1. Définir au niveau d’un communicateur la répartition des cellules et DoFs : possédés localement, fantômes reçus et données non locales non disponibles. Identifier les opérations collectives et leurs préconditions.
2. Planifier échanges de valeurs et d’additions séparément : les contributions hors propriétaire sont réduites chez le propriétaire, puis les fantômes sont mis à jour avant lecture. Empêcher la lecture de fantômes périmés dans l’API ou la détecter.
3. Après raffinement ou repartitionnement, transférer champs et propriétés, reconstruire les dépendances et invalider les anciens caches/indices. Documenter les erreurs collectives pour éviter les blocages silencieux.
4. **Matrice de tests :** 1, 2 et au moins 4 rangs disponibles, partition déséquilibrée et rang vide, conservation d’une constante, solution Poisson comparée au séquentiel, transferts aller-retour. Activer ce chantier seulement dans une CI MPI dédiée.
