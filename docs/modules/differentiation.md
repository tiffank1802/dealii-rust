# differentiation — Différentiation automatique et symbolique

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/differentiation)

**Statut :** conception uniquement, aucune API ni fonctionnalité implémentée. Les noms Rust proposés sont indicatifs, pas des signatures stabilisées.

## Relations et contrats

- **Prérequis :** base : scalaires/tenseurs ; algorithms pour validation des résidus.
- **Consommateurs :** physics, optimization et applications non linéaires.
- **Interface et invariants :** API de Jacobien/Hessien locale : différencier le résidu et renvoyer dimensions et erreurs ; backends séparés.
- **Dépendances externes :** Aucune crate AD choisie ; comparer candidates Rust à ADOL-C/Sacado/SymEngine via FFI seulement après prototype et contrôle des licences. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P2 | `ad/ad_helpers.h, ad/ad_drivers.h` | Prototype AD mode avant sur fonctions locales, sélectionner backend après mesures. | Gradients/Jacobiens analytiques, comparaison différences finies et cas multi-variables. |
| P4 | `ad/ad_number_types.h` | Évaluer mode inverse et Hessiennes si volume mémoire justifie. | Hessien symétrique d’un potentiel, précision et ressources mesurées. |
| P4 | `sd/symengine_optimizer.h, sd/symengine_tensor_operations.h` | Expressions symboliques et optimisation si cas utilisateur identifié. | Expressions analytiques, erreurs de domaine et équivalence avec évaluation numérique. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.
