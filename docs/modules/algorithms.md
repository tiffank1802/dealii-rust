# algorithms — Orchestration non linéaire et temporelle

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/algorithms)

**Statut :** conception uniquement, aucune API ni fonctionnalité implémentée. Les noms Rust proposés sont indicatifs, pas des signatures stabilisées.

## Relations et contrats

- **Prérequis :** lac et contrats de résidus/états ; numerics pour applications.
- **Consommateurs :** Applications, optimisation et adaptateurs sundials/trilinos.
- **Interface et invariants :** Résidu/Jacobien et TimeStepper ; rapport d’état, budgets et erreurs explicites ; pas d’abstraction AnyData imposée.
- **Dépendances externes :** P2 : standard. `sundials` externe uniquement via module adaptateur séparé. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P2 | `newton.h, operator.h` | Newton : résidu, Jacobien, solveur linéaire, tolérance et nombre maximal. | Racine analytique, Jacobien incorrect, non-convergence et singularité. |
| P2 | `theta_timestepping.h, timestep_control.h` | Euler explicite/implicite et Crank–Nicolson ; gestion du temps et pas. | ODE analytique : ordre 1/1/2 ; état, pas invalide et conditions initiales. |
| P2 | `any_data.h, general_data_storage.h` | Données typées et callbacks ; définir la propriété et durée de vie des états. | Transitions, reprise après échec et comportement des callbacks. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.
