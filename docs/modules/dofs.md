# dofs — Numérotation, accès et construction des contraintes FEM

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/dofs)

**Statut :** première implémentation P0 partielle et expérimentale ; voir le [README](../../README.md) pour le périmètre et la validation. Les extensions et noms ci-dessous restent des propositions.

## Relations et contrats

- **Prérequis :** grid, fe, base ; lac uniquement pour la représentation algébrique des contraintes.
- **Consommateurs :** numerics, meshworker, hp, distributed, matrix_free, multigrid.
- **Interface et invariants :** DoFHandler lié à une génération du maillage, DofId typé et table cellule→DoFs ; états périmés refusés.
- **Dépendances externes :** P0 : standard ; pas de crate externe propre au module. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P0 | `dof_handler.h, dof_accessor.h` | Associer un DoF par sommet P1, calculer la table locale→globale et l’ensemble de bord. | Deux triangles partagent deux DoFs ; cardinalité, indices locaux et frontières corrects. |
| P1 | `dof_faces.h, block_info.h, dof_renumbering.h` | Numéroter entités, composantes/blocs et appliquer permutations d’orientation ; renumérotation explicite. | Unicité et continuité des composantes ; solution invariante après renumérotation. |
| P2 | `dof_tools.h, dof_levels.h` | Construire contraintes h/hp depuis le maillage et les éléments ; produire des lignes pour lac. | Hanging node milieu = moyenne pour P1 ; conformité et rejet des cycles par lac. |
| P3 | `number_cache.h, dof_handler_policy.h` | Définir DoFs possédés/partagés/fantômes et échanges après repartitionnement. | Un seul propriétaire, indices cohérents et solution égale sur 1/2/N rangs. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.

## Détail d’implémentation P0

- Construire `DofHandler` à partir d’un maillage valide P1 : numéro global correspondant initialement à chaque sommet, et trois numéros par cellule dans l’ordre FE. Le nombre de DoFs est le nombre de sommets dans ce périmètre précis.
- Extraire les DoFs de bord à partir des arêtes marquées par `grid` ; ne pas les déduire d’un simple test numérique `x==0` afin de préparer les maillages importés.
- État : conserver une référence empruntée ou une génération de maillage ; les tables ne doivent pas rester valides silencieusement après raffinement. Le choix exact dépendra de l’API de mutation.
- **Tests ciblés :** carré à deux triangles : quatre DoFs globaux, trois par cellule et deux en commun ; bord complet, ordre local, indices hors plage et changement de maillage détecté dès l’introduction de la mutation.
