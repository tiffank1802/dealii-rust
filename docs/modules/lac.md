# lac — Algèbre linéaire, contraintes et solveurs

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/lac)

**Statut :** première implémentation P0 partielle et expérimentale ; voir le [README](../../README.md) pour le périmètre et la validation. Les extensions et noms ci-dessous restent des propositions.

## Relations et contrats

- **Prérequis :** base : erreurs et indices ; aucune dépendance vers dofs.
- **Consommateurs :** numerics, algorithms, multigrid, matrix_free et adaptateurs.
- **Interface et invariants :** CooBuilder→CsrMatrix, Vector, LinearOperator et SolverReport ; dimensions, tri, doublons, tolérances et erreurs explicites.
- **Dépendances externes :** P0 : standard, COO/CSR et CG minimaux. `faer` ou `sprs` candidats P1 ; `mpi` et backends natifs seulement en option. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P0 | `sparse_matrix.h, sparsity_pattern.h` | Assembler triplets, trier et sommer doublons, produire CSR avec dimensions validées. | Doublons, lignes vides, indices invalides, produit matrice-vecteur connu. |
| P0 | `solver_cg.h, solver_control.h` | CG sur SPD, critères absolu/relatif et budget d’itérations ; signaler rupture, NaN et non-convergence. | Système SPD connu, résidu explicite, dimensions incompatibles et échec volontaire. |
| P0 | `affine_constraints.h` | Élimination Dirichlet symétrique : colonnes, lignes et correction du second membre ; homogène dans exemple MVP. | Valeur non nulle, symétrie et solution réduite comparée à un calcul direct. |
| P1 | `precondition.h, solver_gmres.h, full_matrix.h` | Jacobi, GMRES et matrices locales ; choisir backend creux après prototype comparatif. | Solutions connues, préconditionneur singulier, résidus et non-convergence. |
| P2 | `affine_constraints.h, block_linear_operator.h` | Contraintes x=Cy+g : fermeture, cycles, condensation CᵀAC et reconstruction ; blocs/Schur. | Chaînes, cycles, inhomogénéités, contraintes dépendantes et comparaison au système réduit. |
| P3 | `la_parallel_vector.h, petsc_*.h, trilinos_*.h` | Protocole possédé/fantôme et adaptateurs optionnels, sans exposer FFI dans le cœur. | Synchronisation multi-rangs, tests de dimension et comparaison avec séquentiel. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.

## Détail d’implémentation P0 et contraintes

1. `CooBuilder` vérifie `row < nrows`, `col < ncols` et finitude. Conversion : trier `(row,col)`, sommer les doublons dans un ordre défini, construire `row_ptr`, `col_idx`, `values` ; lignes vides admises.
2. `CsrMatrix::mul` vérifie la dimension du vecteur et renvoie une erreur en cas de données invalides. Ne pas imposer qu’une matrice générale soit symétrique.
3. Pour une valeur Dirichlet `gᵢ`, soustraire `A[j,i]gᵢ` à chaque second membre libre avant de mettre à zéro ligne et colonne, puis poser `A[i,i]=1`, `b[i]=gᵢ`. Traiter simultanément plusieurs indices contraints et tester les répétitions.
4. CG : initialiser résidu et direction ; contrôler `max_iter`, tolérances, produit de courbure et valeurs finies. Retourner itérations, résidu final et statut ; l’atteinte du budget et une rupture ne sont pas des succès.
5. En P2, fermer les contraintes affines, détecter les cycles et calculer la forme réduite ; conserver cette mécanique algébrique indépendante de `dofs`.

**Tests ciblés :** CSR avec doublons/lignes vides, multiplication connue, Dirichlet non homogène et contraintes multiples, petit système SPD, second membre nul, dimension erronée, budget épuisé et rupture. Ne pas présenter un test de rupture comme une détection complète de non-SPD.

**CG initialisé :** `conjugate_gradient_with_initial_guess` accepte un vecteur initial, contrôle le résidu réel et conserve le critère basé sur `||b||`. Une estimation prolongée doit être adaptée aux valeurs de Dirichlet du système éliminé ; aucun gain en itérations n’est garanti.

**Jacobi :** `JacobiPreconditioner` copie la diagonale CSR et valide sa positivité ; `conjugate_gradient_jacobi` et sa variante avec estimation initiale conservent le contrôle du résidu non préconditionné. Construire le préconditionneur après les contraintes et le reconstruire si la matrice change. La positivité de la diagonale ne prouve pas que la matrice est SPD, et le gain n’est pas garanti. GMRES et autres préconditionneurs restent à faire.
