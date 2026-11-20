# physics — Outils tensoriels et modèles physiques

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/physics)

**Statut :** conception uniquement, aucune API ni fonctionnalité implémentée. Les noms Rust proposés sont indicatifs, pas des signatures stabilisées.

## Relations et contrats

- **Prérequis :** base tenseurs, fe et integrators ; differentiation facultatif.
- **Consommateurs :** numerics et applications mécaniques.
- **Interface et invariants :** Tenseurs et conventions de déformation explicites ; séparer cinématique, loi constitutive et assemblage.
- **Dépendances externes :** `nalgebra` candidat pour calcul tensoriel ; aucune dépendance propre requise pour le MVP. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P1 | `notation.h, transformations.h` | Conversions vectorielles/tensorielles, conventions de composantes et invariants. | Identités connues et transformations de coordonnées en 2D/3D. |
| P2 | `elasticity/kinematics.h, elasticity/standard_tensors.h` | Cinématique et lois d’élasticité linéaire avec hypothèses plane stress/strain explicites. | Tenseur analytique, symétrie et patch test en mécanique. |
| P4 | `elasticity/kinematics.h` | Déformations finies et dérivées constitutives pour matériaux non linéaires. | Objectivité sous rotation, déformation identité, Jacobien tangent comparé aux différences finies. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.
