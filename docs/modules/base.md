# base — Fondations géométriques, numériques et techniques

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/base)

**Statut :** première implémentation P0 partielle et expérimentale ; voir le [README](../../README.md) pour le périmètre et la validation. Les extensions et noms ci-dessous restent des propositions.

## Relations et contrats

- **Prérequis :** Aucun module FEM.
- **Consommateurs :** Tous les modules ; `grid`, `fe` et `lac` en premier.
- **Interface et invariants :** Points et vecteurs finis, erreurs structurées, règles de quadrature immuables ; préciser repères, unités, dimensions et ordre des composantes.
- **Dépendances externes :** P0 : bibliothèque standard. Candidats ultérieurs : `nalgebra` pour matrices/tenseurs, `serde` pour paramètres, `tracing` pour diagnostics ; comparer leur coût et leur MSRV avant adoption. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P0 | `point.h, tensor.h, exceptions.h` | Point2, opérations et déterminant 2×2 ; rejeter données non finies en entrée des constructeurs validés. | Normes et produits connus, cas limites, NaN/infini et erreurs documentées. |
| P0 | `quadrature.h, quadrature_lib.h` | Quadrature du triangle de référence : couples point/poids, somme des poids égale à l’aire ; règles suffisantes pour assemblage et mesure d’erreur. | Monômes jusqu’au degré garanti et vérification de l’aire ; règles plus précises pour la solution sinusoïdale. |
| P1 | `polynomial*.h, tensor*.h` | Bases polynomiales, tenseurs et gradients dimensionnés ; éviter de figer un trait numérique universel. | Identités de dérivation, symétries et résultats analytiques en 2D/3D. |
| P1 | `function*.h, parameter_handler.h` | Fonctions scalaires/vectorielles et paramètres typés ; validation explicite sans parser obligatoire. | Fonctions constantes, temporelles et vectorielles ; paramètres manquants ou invalides. |
| P3 | `index_set.h, mpi*.h, work_stream.h` | Ensembles d’indices, répartition et orchestration : états, frontières et communications explicites. | Ensembles vides/disjoints ; accord multi-rangs et tests de concurrence. |
| P4 | `hdf5.h, symbolic_function.h, std_cxx*` | Isoler services optionnels ; ne pas reproduire les adaptations C++ de la bibliothèque standard. | Tests activés seulement pour les services implémentés, erreurs d’environnement explicites. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.

## Détail d’implémentation P0

- `geometry` : `Point2([f64; 2])`, soustraction, produit scalaire et déterminant. Vérifier la finitude à la frontière d’entrée ; les opérations internes signalent un résultat non fini là où il est utilisé. Les règles d’échelle/tolérance de dégénérescence doivent être explicites et testées.
- `quadrature` : triangle de référence de sommets `(0,0)`, `(1,0)`, `(0,1)` ; points et poids sur ce triangle, pas normalisés à 1. Exposer une règle exacte pour les polynômes d’ordre annoncé et une règle suffisamment riche pour les fonctions non polynomiales du test de convergence.
- `error` : variantes pour dimension, indice, donnée non finie et géométrie invalide ; inclure contexte sans `panic` pour les entrées utilisateur. Un test par variante et un test de propagation au niveau assemblage.
- **Dépendances de livraison :** geometry puis quadrature ; le module ne doit pas connaître `grid`, `fe` ou `lac`. Commencer avec `std`, réévaluer `nalgebra` uniquement si le prototype démontre un bénéfice.
