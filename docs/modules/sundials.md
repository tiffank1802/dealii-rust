# sundials — Intégration des solveurs SUNDIALS

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/sundials)

**Statut :** conception uniquement, aucune API ni fonctionnalité implémentée. Les noms Rust proposés sont indicatifs, pas des signatures stabilisées.

## Relations et contrats

- **Prérequis :** algorithms, lac et conventions FFI.
- **Consommateurs :** Applications ODE, DAE et non linéaires ; backend optionnel.
- **Interface et invariants :** Wrapper RAII isolant pointeurs, callbacks, durées de vie et erreurs ; configuration facultative par backend et version.
- **Dépendances externes :** SUNDIALS natif candidat ; aucun binding Rust choisi avant audit des versions, licences et maintenance. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P4 | `arkode.h, arkode_stepper.h` | Prototype ARKODE avec ODE, callbacks et critères d’arrêt. | Solution analytique et ordre, callback en erreur et libération des ressources. |
| P4 | `ida.h` | Prototype IDA sur DAE avec conditions initiales cohérentes. | DAE analytique, initialisation impossible et rapport d’échec. |
| P4 | `kinsol.h, n_vector.h, sunlinsol_wrapper.h` | KINSOL et adaptation vecteurs/solveurs avec ownership défini. | Racine connue, dimensions incompatibles et destruction après erreur. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.
