# integrators — Formes bilinéaires et linéaires locales

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/integrators)

**Statut :** intégrateur local P1 expérimental de diffusion symétrique, réaction et source dans `src/integrators.rs`, assemblé par `elliptic` ; voir le [README](../../README.md#diffusion-anisotrope-et-réaction). Les autres noms et extensions restent des propositions.

## Relations et contrats

- **Prérequis :** fe pour valeurs/gradients, base pour quadrature.
- **Consommateurs :** numerics et meshworker.
- **Interface et invariants :** Fonctions pures recevant géométrie/évaluations et coefficients, retournant matrices/vecteurs locaux dimensionnés.
- **Dépendances externes :** P0/P1 : standard ; backend dense éventuel issu de lac, pas de dépendance directe. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P1 | `l2.h, laplace.h` | Calculer masse, rigidité et source indépendamment du stockage global. | Matrice P1 analytique ; symétrie, somme des lignes et intégrale constante. |
| P2 | `advection.h, divergence.h, grad_div.h` | Ajouter advection et couplages mixtes ainsi que contributions de faces. | Signes de flux, dimensions de blocs et conservation sur deux cellules. |
| P2 | `elasticity.h` | Former élasticité linéaire avec tenseur constitutif et composantes définies. | Mouvements rigides, symétrie et patch test sur maillage conforme. |
| P4 | `maxwell.h, patches.h` | Noyaux Maxwell et patches après prise en charge H(curl). | Convergence et contrôle de l’orientation tangentielle par cellule. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.

**Validation partielle P1 :** matrices analytiques à coefficients constants et variables, permutation locale, coefficients invalides et comparaison avec l’assemblage Poisson. Le contrôle de positivité porte sur les points de quadrature, pas sur tout le domaine.
