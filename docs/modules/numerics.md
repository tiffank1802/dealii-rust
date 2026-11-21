# numerics — Assemblage, outils de solution et post-traitement

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/numerics)

**Statut :** assemblage Poisson et diffusion-réaction anisotrope P1 et export expérimental, avec conditions de Dirichlet, Neumann mixte et Neumann pur normalisé, interpolation et projection P1, évaluation ponctuelle et indicateur résiduel Poisson limité ; voir le [README](../../README.md). Les autres extensions et noms ci-dessous restent des propositions.

## Relations et contrats

- **Prérequis :** base, grid, fe, dofs et lac.
- **Consommateurs :** Exemples utilisateur ; hp, particles et visualisation.
- **Interface et invariants :** Assembler Poisson à partir de données locales ; distinguer problème, opérateur, résolution et diagnostic ; retourner erreurs mesurées.
- **Dépendances externes :** P0 : export VTK minimal maison. `vtkio` candidat pour I/O étendu ; `rayon` candidat P3. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P0 | `matrix_creator.h, vector_tools_rhs.h` | Assembler Kᵢⱼ=∫∇φᵢ·∇φⱼ et bᵢ=∫fφᵢ, puis imposer Dirichlet et résoudre. | Matrice locale analytique, somme globale sur deux triangles et solution Poisson. |
| P0 | `vector_tools_integrate_difference.h, data_out.h` | Calculer erreurs L2/H1 avec quadrature adaptée ; exporter valeurs nodales et cellules. | u=sin(πx)sin(πy) : pentes proches de 2/1 ; VTK lu par outil indépendant. |
| P1 | `vector_tools_interpolate.h, vector_tools_project.h` | Interpolation, projection, champs de bord, Neumann et évaluation ponctuelle. | Reproduction polynomiale, charge de frontière, point hors domaine et Neumann pur singulier. |
| P2 | `error_estimator.h, solution_transfer.h` | Estimation a posteriori, marquage et transfert après adaptation. | Indicateurs et erreur sur cas de référence ; conservation du champ constant après raffinement. |
| P3 | `data_out*.h` | Sorties distribuées et post-traitement de champs vectoriels. | Associations points/cellules, partitions et lecture des fichiers produits. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.

## Détail d’implémentation du parcours Poisson

1. Construire le maillage, la règle, l’élément et les DoFs ; pour chaque cellule, évaluer `∇φᵢ·∇φⱼ JxW` et `f(xq)φᵢ JxW` et sommer localement.
2. Insérer chaque contribution dans COO et le vecteur global, consolider CSR, puis imposer Dirichlet avec le contrat de `lac`. Vérifier la symétrie sur les petits exemples avant résolution.
3. Calculer les erreurs par intégration de `(u−u_h)²` et `|∇u−∇u_h|²` à l’aide d’une quadrature de précision suffisante ; pas seulement aux sommets. Afficher taille, itérations, résidu et normes.
4. Exporter en VTK Legacy ASCII un champ nodal. L’exemple doit créer le répertoire ou signaler une erreur I/O, et documenter la commande et le chemin de sortie.
5. Raffiner uniformément une série prédéfinie de maillages. Tester la décroissance et les pentes attendues, avec seuils et quadrature calibrés par essais réels ; ne pas inventer des tolérances non vérifiées.

**État P2 partiel :** `estimator::poisson_dirichlet` calcule les termes de volume et de saut intérieur pour la diffusion unitaire et Dirichlet homogène. `UniformRefinement` permet le transfert P1 sur un niveau uniforme. Le marquage, le raffinement local et les conditions Neumann restent hors périmètre de cet indicateur.

**État P1 partiel :** `field` propose l’interpolation nodale et de bord, la masse consistante, la projection L2 via CG et la recherche/évaluation ponctuelle par balayage. Les tolérances, la quadrature et les limites de classification sont détaillées dans le README ; champs vectoriels, spatialisation et formulation générale restent à faire.

**Diffusion-réaction :** `elliptic::assemble` sur P1 triangulaire, `elliptic_quad::assemble` sur Q1 et `elliptic_tetra::assemble` sur P1 tétraédrique acceptent des coefficients dépendant de la cellule. Utiliser respectivement `numerics::assemble_neumann`, `quad::assemble_neumann` ou `tetra::assemble_neumann` pour le flux conormal. Voir le README pour les limites de quadrature, de coercivité et du solveur Neumann pur.

**Conditions mixtes 3D :** `tetra::assemble_neumann` intègre le flux sur les faces triangulaires orientées et `tetra::interpolate_boundary` construit les valeurs de Dirichlet. `neumann3d::solve` traite le Neumann pur tétraédrique avec contrôle par composante et normalisation intégrale.

**Conditions mixtes Q1 :** `quad::assemble_neumann` intègre le flux et `quad::interpolate_boundary` produit les valeurs de Dirichlet. `neumann_quad::solve` traite le Neumann pur Q1, avec contrôle de compatibilité par composante et moyenne intégrale calculée à partir des poids Q1.

**Neumann pur anisotrope :** `neumann::solve_diffusion`, `neumann_quad::solve_diffusion` et `neumann3d::solve_diffusion` traitent la diffusion scalaire sans réaction ; le flux conormal est fourni par l’appelant et la compatibilité assemblée est vérifiée par composante.

**Champs Q1 :** `quad_field` fournit interpolation nodale, matrice de masse, projection L² et évaluation en coordonnées physiques sur quadrilatères ; la recherche est linéaire et les correspondances par tolérance sont projetées sur le bord. `quad::interpolate_boundary` traite les valeurs de bord. Voir le README pour les limites de quadrature et de conditionnement.

**Champs P1 tétraédriques :** `tetra_field` fournit interpolation nodale, matrice de masse consistante et projection L² ; `tetra::interpolate_boundary` traite les faces de bord. La recherche et l’évaluation ponctuelle 3D sont disponibles par balayage : la tolérance barycentrique peut conduire à une extrapolation. Voir le README pour la quadrature et le conditionnement.

**Multi-level non-conforming Poisson:** `grid::NonconformingRefinement::assemble_poisson`
assembles the P1 volume terms of `-Δu=f` on any active level, including meshes
built by successive `refine` calls or stored in `NonconformingRefinementHistory`.
`hanging_constraints` emits midpoint rules for every hanging node, including
those at different levels; `AffineConstraints::close` composes the chains before
`condense` and `reconstruct`. Dirichlet values come from `interpolate_boundary`.
`examples/poisson_nonconforming_adaptive.rs` demonstrates a two-level
mark/refine/solve/coarsen cycle. The usual `Mesh2` boundary routines, residual
estimators and Neumann integration still do not support this seam.
