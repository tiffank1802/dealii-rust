# particles — Particules, propriétés et déplacements

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/particles)

**Statut :** conception uniquement, aucune API ni fonctionnalité implémentée. Les noms Rust proposés sont indicatifs, pas des signatures stabilisées.

## Relations et contrats

- **Prérequis :** grid, recherche spatiale, fe/dofs pour interpolation ; distributed pour MPI.
- **Consommateurs :** Applications de transport et post-traitement.
- **Interface et invariants :** ParticleId stable, position et propriétés ; ParticleHandler propriétaire des particules et de leur association aux cellules.
- **Dépendances externes :** `rstar` candidat pour recherche ; `rayon`/`mpi` optionnels selon profil. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P4 | `particle.h, property_pool.h, generators.h` | Génération, stockage typé et association cellule ; contrôler dimensions et identifiants. | Cardinalité, propriétés et position sur frontière/hors domaine. |
| P4 | `particle_handler.h, utilities.h` | Déplacement, localisation et évaluation des champs ; migration après changement de cellule. | Trajectoire connue et interpolation d’un champ affine. |
| P4 | `partitioner.h, data_out.h` | Migration entre rangs, transfert après raffinement et export. | Conserver IDs/nombre/propriétés après migration et sortie lisible. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.
