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

- Implémentés : points et triangles 2D, maillage structuré, DoFs P1 par sommet, mapping affine, quadratures, assemblage Poisson et diffusion-réaction anisotrope P1, Dirichlet, Neumann mixte ou pur normalisé, COO/CSR, CG, export VTK, import d’un sous-ensemble Gmsh 4.1 ASCII, raffinement uniforme, indicateur résiduel Poisson limité et opérations sur champs P1.
- Limités au calcul séquentiel scalaire ; pas encore de 3D, raffinement local adaptatif, MPI, Trilinos, autres familles d’éléments ou prise en charge complète de Gmsh. La validation de `Mesh2::new` ne détecte pas les intersections géométriques ni tous les cas non conformes : le maillage fourni doit être conforme.
- Aucune dépendance Cargo externe. Rust 1.70 est la version minimale **cible**, configurée dans la CI avec la version stable ; ces deux configurations restent à confirmer par leur exécution en CI.
- Les tests couvrent les matrices locales et globales, les coefficients elliptiques et données invalides, le solveur, les normales, les conditions mixtes et de Neumann pur, l’import et les marqueurs Gmsh, le raffinement et le transfert P1, les indicateurs analytiques, l’interpolation, la projection et l’évaluation ponctuelle, ainsi que la convergence des solutions manufacturées. Le fichier VTK a aussi été lu avec `meshio` lors d’une vérification locale.

## Documentation et plans

- [Feuille de route initiale](PLAN_IMPLEMENTATION.md) : objectif Poisson 2D et critères de validation.
- [Plan des modules deal.II](docs/PLAN_MODULES_DEALII.md) et [fiches détaillées](docs/modules/README.md).
- [Dépendances de deal.II et stratégie Rust](docs/DEPENDANCES_DEALII.md) : Trilinos, PETSc, MPI et autres intégrations.
- [Crates candidates et stratégie de test](docs/modules/DEPENDANCES_ET_TESTS.md) et [inventaire des en-têtes amont](docs/INVENTAIRE_ENTETES_DEALII.md).

Le clone deal.II de référence, non versionné et potentiellement temporaire, se trouve dans `.cache/dealii`. La procédure pour le recréer est indiquée dans le plan détaillé.

## Licence

Le code original du projet est proposé au choix sous [MIT](LICENSE-MIT) ou [Apache-2.0](LICENSE-APACHE). Les dépendances et les éléments éventuellement repris de tiers conservent leurs propres licences.
