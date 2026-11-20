# boost_adaptors — Pont entre représentations géométriques

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/boost_adaptors)

**Statut :** conception uniquement, aucune API ni fonctionnalité implémentée. Les noms Rust proposés sont indicatifs, pas des signatures stabilisées.

## Relations et contrats

- **Prérequis :** base et grid.
- **Consommateurs :** Adaptateurs géométriques externes ; pas de dépendance inverse du cœur.
- **Interface et invariants :** Conversions explicites point/segment/boîte ; définir coordonnées, précision et propriété des données.
- **Dépendances externes :** Privilégier un backend Rust si pertinent ; pas de dépendance Boost choisie. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P4 | `point.h, segment.h, bounding_box.h` | Définir conversions vers le backend retenu et vérification de dimensions. | Aller-retour, points limites, boîte dégénérée et coordonnées non finies. |
| P4 | `bounding_box.h` | Tester filtres de recherche et ordre des bornes indépendamment du backend. | Intersections et séparations comparées à un calcul direct. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.
