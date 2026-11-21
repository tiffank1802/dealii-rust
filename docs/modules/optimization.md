# optimization — Minimisation et recherche de pas

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/optimization)

**Statut :** conception uniquement, aucune API ni fonctionnalité implémentée. Les noms Rust proposés sont indicatifs, pas des signatures stabilisées.

## Relations et contrats

- **Prérequis :** lac vecteurs et diagnostics ; differentiation optionnelle.
- **Consommateurs :** Applications de calibration ; trilinos pour ROL.
- **Interface et invariants :** Objectif et gradient, recherche linéaire et rapport d’arrêt, politique explicite de gestion d’échec.
- **Dépendances externes :** P2 : standard. ROL via adaptateur natif éventuel, non choisi. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P2 | `line_minimization.h` | Recherche de pas avec budget, domaines invalides et critère documenté. | Quadratique analytique, rejet NaN et absence de pas acceptable. |
| P2 | `solver_bfgs.h` | BFGS ou L-BFGS avec gestion de courbure et tolérance. | Minimum de fonctions tests, comparaison du gradient, non-convergence. |
| P4 | `rol/vector_adaptor.h` | Adaptateur de vecteur pour ROL si besoin avéré. | Conversions et objectif/gradient identiques au backend local. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.
