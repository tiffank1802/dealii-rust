# matrix_free — Application d’opérateurs sans matrice assemblée

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/matrix_free)

**Statut :** conception uniquement, aucune API ni fonctionnalité implémentée. Les noms Rust proposés sont indicatifs, pas des signatures stabilisées.

## Relations et contrats

- **Prérequis :** grid, fe, dofs, lac et contraintes ; P3 MPI selon backend.
- **Consommateurs :** multigrid et solveurs itératifs.
- **Interface et invariants :** LinearOperator appliqué cellule par cellule ; caches de géométrie liés à une version de maillage, vecteurs locaux et fantômes explicites.
- **Dépendances externes :** `rayon` candidat CPU ; backend accélérateur non choisi, aucune API GPU promise. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P3 | `matrix_free.h, fe_evaluation.h, shape_info.h` | Opérateur Laplace tensorisé, lecture/écriture et contraintes ; définir le sous-ensemble FE. | Comparer action à CSR sur vecteurs aléatoires, bords et contraintes inclus. |
| P3 | `tensor_product_kernels.h, mapping_info.h` | Somme factorisée et cache géométrique ; invalider lors de modification. | Concordance aux noyaux directs, invalidation après raffinement et benchmark documenté. |
| P4 | `evaluation_kernels_face.h, hanging_nodes_internal.h` | Faces DG, nœuds contraints et maillages adaptés. | Opérateur DG assemblé vs sans matrice et orientation des interfaces. |
| P4 | `portable_matrix_free.h, vector_data_exchange.h` | Prototype backend portable et échanges distribués derrière interface. | Comparaison CPU/backend, 1/2/N rangs, mémoire et coût de transfert. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.

## Protocole opérateur et validation

1. Définir `apply(src, dst)` avec dimensions, propriété des vecteurs et convention claire pour `dst` ; traiter séparément l’échange des fantômes et les contraintes.
2. D’abord implémenter un chemin de référence naïf sans matrice, puis la somme factorisée pour les cellules tensorisées ; éviter de l’annoncer pour les triangles tant qu’elle n’est pas prise en charge.
3. Mettre en cache géométrie et fonctions de forme en les associant à une génération de maillage/DoFs. Refuser ou reconstruire après modification ; mesurer le compromis mémoire/temps.
4. Sur plusieurs vecteurs déterministes, comparer `A_assemblée x` et `A_sans_matrice x`, avant puis après contraintes ; comparer résidu et solution du problème complet.
5. Mesurer les performances CPU avant d’introduire parallélisme, SIMD ou GPU ; si MPI est activé, tester partitions vides et états fantômes périmés.
