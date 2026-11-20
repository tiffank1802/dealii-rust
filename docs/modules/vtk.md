# vtk — Échange et visualisation VTK

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/vtk)

**Statut :** première implémentation P0 partielle et expérimentale ; voir le [README](../../README.md) pour le périmètre et la validation. Les extensions et noms ci-dessous restent des propositions.

## Relations et contrats

- **Prérequis :** grid et champs de numerics.
- **Consommateurs :** Exemples, post-traitement et visualisation.
- **Interface et invariants :** Writer vers flux fourni par l’utilisateur ; positions, connectivité et données point/cellule à tailles vérifiées.
- **Dépendances externes :** P0 : export VTK Legacy ASCII minimal sans crate. `vtkio` candidat pour formats étendus. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P0 | `base/data_out_base.h, numerics/data_out.h` | Écrire unstructured grid triangulaire et scalaire nodal ; préciser types et indexation. | Petit fichier de référence, nombre de sommets/cellules, valeurs et lecture indépendante. |
| P1 | `vtk/utilities.h, numerics/data_out_faces.h` | Étendre champs vectoriels, données cellule et autres cellules prises en charge. | Association, composantes et valeurs lues ; erreurs de dimensions. |
| P3 | `numerics/data_out*.h` | Sortie parallèle partitionnée et métadonnées si nécessaire. | Chaque partition et ensemble global lisibles, données fantômes non dupliquées selon convention. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.
