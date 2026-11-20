# dealii-rust

Bibliothèque d’éléments finis en Rust inspirée de deal.II. **Le projet est expérimental** : il propose un premier parcours Poisson 2D, mais ne vise pas encore la parité avec deal.II.

## Premier exemple

Le programme résout `−Δu = 2π² sin(πx) sin(πy)` sur le carré unité, avec `u = 0` au bord et `u = sin(πx) sin(πy)` comme solution de référence. Il utilise des triangles conformes P1, une quadrature de degré 5, une matrice creuse CSR et le gradient conjugué.

```sh
cargo test
cargo run --example poisson_2d -- 16
```

L’argument facultatif est le nombre de subdivisions par côté (16 par défaut). L’exemple écrit `poisson_2d.vtk` dans le répertoire courant, au format VTK Legacy ASCII, visualisable par exemple avec ParaView. Il affiche le nombre d’itérations, le résidu et les erreurs L2 et en semi-norme H1. Les fichiers générés ne sont pas versionnés.

Le gradient conjugué démarre à zéro et utilise par défaut `max(1e-12, 1e-10 × ||b||)` comme seuil de norme du résidu, avec 10 000 itérations au maximum. La matrice doit être symétrique définie positive ; une rupture ou une non-convergence est signalée par une erreur. Le résidu annoncé est recalculé à partir de `b − Ax`.

## Poisson sur quadrilatères Q1

```sh
cargo run --example poisson_q1 -- 16
```

Ce second parcours résout le même problème manufacturé sur un maillage de quadrilatères Q1, et écrit `poisson_q1.vtk` (VTK Legacy ASCII, cellules de type quadrilatère). `quad::QuadMesh` et `quad::BilinearMap` sont séparés de l’API triangulaire ; les IDs de sommets et cellules ne doivent pas être échangés entre maillages. Chaque cellule est ordonnée dans le sens antihoraire selon les sommets de référence `(0,0), (1,0), (1,1), (0,1)`. Le mapping isoparamétrique est bilinéaire et le Jacobien est évalué à chaque point. Les quadrilatères non convexes, inversés ou dégénérés sont rejetés ; le déterminant aux quatre coins doit dépasser `64 ε` fois le carré de la plus grande arête. Les intersections, recouvrements et nœuds pendants entre cellules ne sont pas détectés : le maillage fourni doit être conforme.

L’assemblage et les erreurs utilisent la règle de Gauss tensorisée 3×3, exacte sur le carré de référence jusqu’au degré 5 **dans chaque coordonnée**. Pour un mapping non affine, l’intégrande de rigidité est généralement rationnel et l’intégration est approchée. Les tests comprennent un patch affine sur maillage déformé, des matrices analytiques et la convergence sur grilles régulières et déformées. Le parcours Q1 couvre Poisson avec Dirichlet, conditions mixtes ou Neumann pur ; les API triangulaires de diffusion-réaction, import Gmsh, raffinements et indicateurs ne s’appliquent pas aux quadrilatères. Les tétraèdres 3D disposent d’un parcours séparé, décrit ci-dessous.

### Conditions mixtes sur quadrilatères Q1

```sh
cargo run --example poisson_q1_mixed -- 16
```

L’exemple résout `−Δu = −4` avec `u = 1 + x + x² + 2y + y²`, Dirichlet sur `x = 0` et `g = ∂u/∂n` ailleurs. Il écrit `poisson_q1_mixed.vtk` ; la sélection par coordonnées est propre au carré structuré. `QuadMesh::boundary_edge_geometry` fournit les arêtes dirigées selon la cellule et leur normale sortante. `quad::interpolate_boundary` construit les contraintes nodales, et `quad::assemble_neumann` calcule `∫g φᵢ ds` à ajouter au second membre **avant** l’élimination de Dirichlet. Les arêtes Q1 bilinéaires sont droites et les traces des fonctions de forme sont linéaires ; la quadrature Gauss de degré 5 est exacte pour les flux polynomiaux de degré ≤ 4, à l’arrondi près. Les IDs invalides ou répétés et les valeurs non finies sont rejetés.

L’assemblage du flux ne résout pas à lui seul Neumann pur : la matrice a un noyau constant et requiert compatibilité et normalisation. Utiliser `neumann_quad::solve` pour ce cas ; ne pas appliquer directement CG au système singulier.

### Neumann pur sur quadrilatères Q1

```sh
cargo run --example poisson_q1_neumann -- 16
```

L’exemple utilise `−Δu = −4`, `g = 2(x nₓ + y nᵧ)` et la solution de moyenne nulle `u = x² + y² − 2/3` sur le carré unité. Il écrit `poisson_q1_neumann.vtk`. `neumann_quad::solve` vérifie, sur chaque composante du graphe des degrés de liberté, la compatibilité **assemblée** `|∑bᵢ| ≤ tol_abs + tol_rel ∑|bᵢ|` (défauts `1e-12` et `1e-10`). Il rejette les charges incompatibles sans projection ; la quadrature peut influencer ce contrôle.

Le solveur ancre temporairement un sommet par composante, résout par CG puis retranche la moyenne intégrale Q1. Les poids `∫φᵢ` sont calculés avec la quadrature 3×3 : sur un quadrilatère bilinéaire non affine, ils ne sont en général pas égaux, contrairement à une moyenne nodale. La règle intègre ces poids exactement à l’arrondi près. Le rapport indique l’aire, l’ancre, le déséquilibre et les résidus du système ancré et du système original ; ils peuvent différer si un déséquilibre est accepté dans la tolérance. Ce solveur traite uniquement le Laplacien scalaire unitaire avec Neumann sur tout le bord d’un maillage conforme, sans coefficients variables ni réaction.

## Poisson P1 sur tétraèdres 3D

```sh
cargo run --example poisson_3d -- 8
```

Cet exemple résout `−Δu = 3π²u` sur le cube unité avec Dirichlet homogène et `u = sin(πx)sin(πy)sin(πz)`. `tetra::TetMesh::unit_cube` découpe chaque cube en six tétraèdres conformes ; `tetra::assemble_poisson` assemble le système P1 et `tetra::write_vtk` écrit `poisson_3d.vtk` (cellules VTK de type tétraèdre). L’API est distincte des maillages 2D ; les IDs sont propres à chaque maillage.

Les tétraèdres utilisent les sommets de référence `(0,0,0), (1,0,0), (0,1,0), (0,0,1)` et un mapping affine d’orientation positive. Un déterminant inférieur ou égal à `64 ε` fois le cube de la longueur de la plus grande arête est rejeté. Les faces partagées doivent avoir des orientations opposées ; les faces de frontière sont orientées vers l’extérieur et leurs normales et aires sont accessibles. Comme pour les maillages 2D, les intersections, recouvrements et nœuds pendants ne sont pas détectés : la conformité géométrique reste une précondition.

La source et les erreurs utilisent une quadrature Gauss-Duffy à 64 points, exacte pour les polynômes de degré total ≤ 5 sur le tétraèdre de référence, à l’arrondi près. Les tests vérifient les monômes, la matrice locale analytique, les normales, un patch affine et la convergence L2/H1. Ce parcours 3D est limité à Poisson scalaire avec Dirichlet, conditions mixtes ou Neumann pur ; la diffusion-réaction, l’import Gmsh 3D et le raffinement des tétraèdres ne sont pas implémentés.

### Conditions mixtes sur tétraèdres

```sh
cargo run --example poisson_3d_mixed -- 8
```

L’exemple utilise `u = 1 + x + x² + 2y + y² + 3z + z²`, avec Dirichlet sur `x = 0` et le flux sortant `g = ∂u/∂n` sur les autres faces. La sélection par coordonnées est propre au cube structuré ; il n’y a pas encore de marqueurs Gmsh 3D. `tetra::interpolate_boundary` déduplique les sommets des faces de Dirichlet ; `tetra::assemble_neumann` intègre `∫g φᵢ dS` sur des `BoundaryFaceId` uniques, avec une quadrature triangulaire de degré 5. Le vecteur obtenu doit être ajouté à la charge volumique **avant** l’élimination de Dirichlet. Les IDs invalides ou répétés et les données non finies sont rejetés.

L’assemblage de charge seul ne résout pas le cas Neumann pur : la matrice de diffusion possède un noyau constant par composante et exige compatibilité et normalisation. Utiliser `neumann3d::solve` pour ce cas ; ne pas appliquer directement CG sur la matrice singulière.

### Neumann pur 3D

```sh
cargo run --example poisson_3d_neumann -- 8
```

Cet exemple résout `−Δu = −6`, avec `g = 2(x nₓ + y nᵧ + z n_z)` et la solution de moyenne nulle `u = x² + y² + z² − 1`. `neumann3d::solve` vérifie la compatibilité **discrète** séparément pour chaque composante du graphe P1 : `|∑bᵢ| ≤ tol_abs + tol_rel ∑|bᵢ|`, avec les mêmes options et défauts (`1e-12`, `1e-10`) que `neumann::solve` en 2D. Il rejette une charge incompatible sans la projeter, fixe temporairement un sommet par composante, puis retire sa moyenne intégrale FEM, calculée à partir des volumes des tétraèdres. Le rapport fournit le volume, l’ancre, le déséquilibre et les résidus du système fixé et du système original. Si une incompatibilité est acceptée dans la tolérance, ces deux résidus peuvent différer. La compatibilité est évaluée après quadrature et non sur les données continues exactes.

Le solveur couvre uniquement le Laplacien scalaire avec Neumann sur tout le bord d’un maillage conforme. Les coefficients variables et les termes de réaction ne sont pas pris en charge.

## Conditions aux limites mixtes

```sh
cargo run --example poisson_mixed -- 16
```

Cet exemple résout `−Δu = −4` avec `u = 1 + x + x² + 2y + y²` sur le carré unité : Dirichlet non homogène à gauche et flux de Neumann `g = ∂u/∂n` sur les trois autres côtés. Il écrit `poisson_mixed.vtk` et affiche les erreurs intégrées. Le choix des côtés par coordonnées est propre au maillage structuré de l’exemple ; l’exemple Gmsh ci-dessous utilise plutôt des groupes physiques importés.

L’API `numerics::assemble_neumann` prend des `BoundaryEdgeId` de ce maillage et une fonction de flux recevant la position et la normale **sortante**. Elle retourne le vecteur `∫g v ds` à ajouter au second membre volumique **avant** l’élimination de Dirichlet. Des arêtes répétées ou invalides et un flux non fini sont rejetés.

### Neumann pur

```sh
cargo run --example poisson_neumann -- 16
```

Cet exemple résout `−Δu = −4` avec `∂u/∂n = 2x nₓ + 2y nᵧ`. La solution de référence, de moyenne nulle, est `u = x² + y² − 2/3`. Il écrit `poisson_neumann.vtk` et affiche les erreurs et les diagnostics du solveur.

Un problème de Neumann pur a un noyau constant et exige la compatibilité `∫Ω f + ∫∂Ω g = 0` sur chaque composante connexe. La fonction `neumann::solve` vérifie séparément chaque composante du graphe des DoFs : `|∑bᵢ| ≤ tol_abs + tol_rel ∑|bᵢ|`. Les tolérances par défaut sont `1e-12` et `1e-10`. Le test porte sur le **second membre assemblé** et inclut donc l’erreur de quadrature ; des données non compatibles sont rejetées.

Pour résoudre, on fixe temporairement un DoF par composante, puis on retire la **moyenne intégrale FEM** de chaque composante. Aucune projection du second membre n’est effectuée. Si une petite incompatibilité est admise par la tolérance, le résidu du système original peut dépasser celui du système fixé : les deux sont rapportés. Le CG ordinaire ne doit pas être appliqué directement à la matrice singulière.

## Diffusion anisotrope et réaction

```sh
cargo run --example elliptic_mixed -- 16
```

Cet exemple résout `−div(A∇u) + cu = f` sur le carré unité avec `A = [[2+x, 0,5], [0,5, 3+y]]`, `c = 1+x+y` et la solution manufacturée `u = 1+x+x²+2y+y²`. Il impose Dirichlet à gauche et un flux conormal sur les trois autres côtés, et écrit `elliptic_mixed.vtk`.

`integrators::diffusion_reaction_p1` calcule une matrice et un vecteur locaux ; `elliptic::assemble` les assemble en CSR sans imposer de conditions aux limites. Ses fonctions de coefficients reçoivent l’ID de cellule et le point physique : il est donc possible d’utiliser des matériaux constants par cellule à partir des groupes physiques Gmsh. Le tenseur `DiffusionTensor` est symétrique par construction. Il doit être défini positif, la réaction doit être non négative et la source finie ; ces propriétés sont contrôlées **uniquement aux points de quadrature**, pas entre eux. L’utilisateur doit garantir l’ellipticité uniforme et les hypothèses du problème continu. Des coefficients très mal conditionnés ou des plages numériques extrêmes peuvent rendre le système peu fiable même après ces contrôles.

La quadrature triangulaire de degré 5 intègre exactement, à l’arrondi près, un tenseur polynomial de degré ≤ 5, une réaction de degré ≤ 3 et une source de degré ≤ 4 dans ces formes P1 ; pour les autres données, elle introduit une erreur de quadrature. Ajouter `∫g v` au second membre avec `numerics::assemble_neumann`, où **`g = n · A∇u`** est le flux conormal sortant, puis imposer Dirichlet si nécessaire. L’assemblage ne garantit pas que le système soit défini positif : avec `c = 0` et Neumann partout, il existe un noyau constant par composante. `neumann::solve` reste réservé au Laplacien unitaire ; il ne résout pas ce problème elliptique généralisé. Une réaction positive peut supprimer le noyau. L’indicateur `estimator::poisson_dirichlet` reste lui aussi réservé au Laplacien unitaire sans réaction.

## Import Gmsh et résolution sur maillage étiqueté

```sh
cargo run --example poisson_gmsh
# Ou : cargo run --example poisson_gmsh -- chemin/maillage.msh
```

L’exemple lit `tests/fixtures/square41.msh` par défaut, résout un problème mixte de solution `1 + x + 2y`, et écrit `poisson_gmsh.vtk`. Il exige le groupe physique de courbe 40 pour Dirichlet et les groupes 10, 20 ou 30 pour Neumann ; toutes les arêtes extérieures doivent appartenir à exactement un type de condition.

L’API `gmsh::read`/`gmsh::parse` retourne un `ImportedMesh` contenant un `Mesh2` et des métadonnées séparées, consultables par `VertexId`, `CellId` et `BoundaryEdgeId` du maillage importé. Les tags de nœuds et d’éléments triangle/ligne d’origine sont conservés ; les groupes physiques peuvent être multiples par entité. `cells_with_tag`, `boundary_edges_with_tag` et `boundary_vertices_with_tag` permettent de sélectionner les régions. Les noms optionnels sont indexés par dimension et tag. Une arête sans groupe peut être absente du fichier ou présente sans groupe : `boundary_element_tag` distingue les deux.

**Sous-ensemble strict :** MSH 4.1 ASCII, taille de données 8, sections `MeshFormat`, `Entities`, `Nodes`, `Elements` obligatoires, `PhysicalNames` optionnelle ; coordonnées planes `z = 0`, nœuds non paramétriques, éléments point (15), segment linéaire (1) et triangle linéaire (2) uniquement. Les triangles doivent être orientés positivement. Chaque segment doit être une arête extérieure ; les interfaces internes ne sont pas prises en charge. Les points physiques sont associés aux nœuds classés sur l’entité ponctuelle. Un maillage peut omettre des segments de bord, mais ces arêtes seront alors sans marqueurs. Les sections supplémentaires (notamment partitions, périodicité et données), le binaire, les quadrangles, les éléments d’ordre supérieur et la 3D sont rejetés. Les tags de point-élément et les identifiants géométriques ne sont pas exposés. Les contrôles géométriques restent ceux de `Mesh2` : le lecteur ne prouve pas la conformité globale ni l’absence d’intersections. Les fichiers sont lus intégralement en mémoire ; ne pas les traiter comme une interface durcie pour des entrées non fiables.

## Interpolation, projection et évaluation P1

```sh
cargo run --example field_tools
```

Le module `field` travaille sur des champs scalaires P1 continus : `interpolate` échantillonne une fonction aux sommets ; `interpolate_boundary` produit les contraintes de Dirichlet sur des arêtes sélectionnées en dédupliquant les sommets communs. `mass_matrix` assemble la matrice de masse consistante et `projection_system` le système `M u = b`, avec `Mᵢⱼ = ∫φᵢφⱼ` et `bᵢ = ∫fφᵢ`. `project_l2` résout ce système par CG **sans condition de bord** ; pour une projection contrainte, appliquer `with_dirichlet` au système avant résolution. La quadrature du second membre est de degré 5 : elle intègre exactement `fφᵢ` lorsque `f` est polynomiale de degré ≤ 4, à l’arrondi près. Pour d’autres fonctions, la projection est approximée par cette quadrature. Les tolérances CG, en particulier la tolérance absolue, doivent être adaptées à l’échelle des données ; la non-convergence est signalée.

`locate` et `evaluate` recherchent un point par parcours linéaire des cellules. Une tolérance barycentrique explicite, sans dimension, autorise de petits dépassements ; un point réellement contenu numériquement est préféré à un point accepté par tolérance. Sinon la première cellule compatible est choisie. Les poids ne sont pas bornés : une correspondance par tolérance **extrapole** légèrement, elle ne projette pas sur le bord. Un point hors domaine retourne `None` ; les données non finies et débordements détectés retournent une erreur. La classification des points de bord dépend des arrondis et de la qualité géométrique ; l’API suppose un maillage conforme sans recouvrement, comme le reste de `Mesh2`. Elle ne fournit ni index spatial ni évaluation de gradient sur les interfaces.

## Raffinement uniforme et indicateur résiduel

```sh
cargo run --example poisson_refinement
```

L’exemple résout successivement le problème Poisson homogène sur quatre niveaux uniformes, puis affiche l’erreur H1 et l’indicateur. `grid::UniformRefinement::new` partage les milieux des arêtes et remplace chaque triangle par quatre triangles conformes. Le nouvel objet expose les relations parent/enfants pour les cellules, ainsi que le parent de chaque arête de frontière. Les anciens sommets conservent leurs indices, mais les IDs de cellules et d’arêtes appartiennent à leur maillage respectif. `prolongate_p1` interpole le champ nodal aux milieux et préserve la fonction P1 par morceaux. Les métadonnées Gmsh ne sont **pas** transférées automatiquement ; les relations parentales permettent de les recopier explicitement.

`estimator::poisson_dirichlet` fournit des indicateurs par cellule et leur norme globale pour `−Δu = f`, diffusion constante unitaire, P1 triangulaire et Dirichlet homogène sur tout le bord, supposé imposé par l’appelant. Il calcule

```text
η_K² = h_K² ∫_K f² + ½ Σ_(e intérieur de K) h_e ∫_e [∂u_h/∂n]².
```

`h_K` est la longueur de la plus grande arête, `h_e` celle de l’arête ; la quadrature triangulaire de degré 5 approxime l’intégrale de `f²`. La contribution de saut de gradient est constante sur chaque arête P1 et répartie à parts égales entre les deux cellules. C’est un **indicateur**, non une borne garantie ni une reproduction de l’estimateur deal.II. Il ne traite pas les résidus de bord Neumann, les coefficients variables ou les maillages non conformes. La décroissance observée sur l’exemple manufacturé ne prouve pas sa fiabilité générale. Le raffinement est global uniquement : aucun marquage, déraffinement, hanging node ou adaptation locale n’est implémenté.

## Périmètre et validation

- Implémentés : points 2D/3D, triangles et tétraèdres P1, quadrilatères Q1, maillages structurés, mappings affine et bilinéaire, quadratures, assemblage Poisson et diffusion-réaction anisotrope P1, Dirichlet, Neumann mixte ou pur normalisé, COO/CSR, CG, export VTK, import d’un sous-ensemble Gmsh 4.1 ASCII, raffinement uniforme, indicateur résiduel Poisson limité et opérations sur champs P1.
- Limités au calcul séquentiel scalaire ; pas encore de raffinement local adaptatif, MPI, Trilinos, autres familles d’éléments que P1/Q1 ou prise en charge complète de Gmsh. Le parcours 3D est restreint à Poisson sur tétraèdres P1 avec Dirichlet, conditions mixtes ou Neumann pur. La validation de `Mesh2::new` ne détecte pas les intersections géométriques ni tous les cas non conformes : le maillage fourni doit être conforme.
- Aucune dépendance Cargo externe. Rust 1.70 est la version minimale **cible**, configurée dans la CI avec la version stable ; ces deux configurations restent à confirmer par leur exécution en CI.
- Les tests couvrent les matrices locales et globales, les coefficients elliptiques et données invalides, le solveur, les normales, les conditions mixtes et Neumann pur 2D/3D, l’import et les marqueurs Gmsh, le raffinement et le transfert P1, les indicateurs analytiques, l’interpolation, la projection et l’évaluation ponctuelle, les mappings et patch tests Q1 et tétraédriques, ainsi que la convergence des solutions manufacturées. Le fichier VTK a aussi été lu avec `meshio` lors d’une vérification locale.

## Documentation et plans

- [Feuille de route initiale](PLAN_IMPLEMENTATION.md) : objectif Poisson 2D et critères de validation.
- [Plan des modules deal.II](docs/PLAN_MODULES_DEALII.md) et [fiches détaillées](docs/modules/README.md).
- [Dépendances de deal.II et stratégie Rust](docs/DEPENDANCES_DEALII.md) : Trilinos, PETSc, MPI et autres intégrations.
- [Crates candidates et stratégie de test](docs/modules/DEPENDANCES_ET_TESTS.md) et [inventaire des en-têtes amont](docs/INVENTAIRE_ENTETES_DEALII.md).

Le clone deal.II de référence, non versionné et potentiellement temporaire, se trouve dans `.cache/dealii`. La procédure pour le recréer est indiquée dans le plan détaillé.

## Licence

Le code original du projet est proposé au choix sous [MIT](LICENSE-MIT) ou [Apache-2.0](LICENSE-APACHE). Les dépendances et les éléments éventuellement repris de tiers conservent leurs propres licences.
