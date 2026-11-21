# trilinos — Adaptateurs NOX et ROL de Trilinos

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/trilinos)

**Statut :** conception uniquement, aucune API ni fonctionnalité implémentée. Les noms Rust proposés sont indicatifs, pas des signatures stabilisées.

## Relations et contrats

- **Prérequis :** lac backends Trilinos, algorithms et optimization ; MPI éventuel.
- **Consommateurs :** Applications non linéaires et optimisation.
- **Interface et invariants :** Wrappers spécifiques versionnés, API Rust sûre et erreurs natives ; garder les solveurs linéaires Trilinos dans lac.
- **Dépendances externes :** Trilinos natif candidat, aucun binding sélectionné ; installation et licences par composant à vérifier. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P4 | `nox.h` | Adapter résidu/Jacobien et diagnostics à NOX si besoin identifié. | Solution non linéaire analytique, échec solveur et cleanup. |
| P4 | `rol_adaptor.h` | Adapter vecteurs et objectif ROL, sans imposer ROL au cœur. | Minimum quadratique, cohérence des vecteurs et gestion d’erreurs. |
| P4 | `lac/trilinos_*.h` | Documenter compatibilité avec backend linéaire et vecteurs distribués. | Combinaisons de versions annoncées, 1/2 rangs et dimensions incompatibles. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.
