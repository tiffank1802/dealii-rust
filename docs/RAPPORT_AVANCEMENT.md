# Rapport d’avancement — deal.II en Rust

État de la branche `arena/1c071f35-dealii-rust`, incluant les modifications locales non commitées. La bibliothèque reste **expérimentale et séquentielle** : elle ne vise pas encore la parité avec deal.II. La [feuille de route](../PLAN_IMPLEMENTATION.md) décrit les critères complets ; le [plan des modules](PLAN_MODULES_DEALII.md) répertorie 26 fiches et 91 sous-chantiers. Aucun pourcentage global d’avancement n’est déduit de ces chiffres : les sous-chantiers n’ont pas tous la même taille.

## Bilan par jalon

| Jalon | État et éléments concrets | Reste à faire |
| --- | --- | --- |
| J0 — Projet | Crate Cargo sans dépendance externe, licence `MIT OR Apache-2.0`, CI stable et MSRV 1.70. | Compléter la gouvernance et la politique de contribution ; vérifier la CI après intégration des changements locaux. |
| J1–J4 — Parcours Poisson 2D | Maillage triangulaire, P1 et quadrature, assemblage COO/CSR, Dirichlet, CG, exemple, erreurs et export VTK ; tests analytiques et de convergence. | Comparaison indépendante avec deal.II et validation élargie aux plateformes et cas réels. |
| J5 — Robustesse et interopérabilité | Dirichlet non homogène, conditions mixtes, Neumann pur et import Gmsh 4.1/2.2 limité ; tests de données invalides et de convergence. | Comparer avec un solveur de référence ; étendre formats, variantes et maillages pris en charge. |
| J6 — Autres cellules | Q1 quadrilatère 2D et P1 tétraédrique 3D : Poisson, conditions mixtes, Neumann pur, diffusion-réaction anisotrope et exemples. Champs : interpolation, projection L² et évaluation P1/Q1. | Import 3D/quadrilatères, autres formulations et champs vectoriels. |
| J7 — Adaptativité | Raffinement global pour les trois cellules ; marquage de Dörfler, red-green P1 triangulaire et pile de niveaux avec retour complet et injection nodale ; localement, contraintes affines et un niveau non conforme triangulaire avec détection des nœuds pendants. | Déraffinement sélectif, adaptation locale Q1/3D, topologie non conforme et contraintes multiniveau ; éléments d’ordre supérieur et `hp`. |
| J8 — Performance/distribution | CG avec estimation initiale, préconditionnement Jacobi, exemple comparatif et premières mesures séquentielles. | Benchmarks plus approfondis et mesure mémoire réelle/maximale, parallélisme, MPI et étude de solveurs/backends optionnels ; aucun backend PETSc ou Trilinos intégré. |

## Champs tétraédriques

- `tetra_field::{interpolate,mass_matrix,projection_system,project_l2}` : interpolation nodale, masse consistante et projection L² par CG sur éléments P1 3D. La quadrature de degré 5 intègre la charge d’un polynôme physique de degré ≤ 4 à l’arrondi près. Une projection contrainte peut être construite avec `tetra::interpolate_boundary` et `with_dirichlet`.
- `tetra_field::{locate,evaluate}` : recherche linéaire des cellules et évaluation barycentrique. Une cellule contenant numériquement le point est préférée à une cellule acceptée seulement par tolérance ; sinon la première admissible est choisie. La tolérance barycentrique peut entraîner une **extrapolation**, sans projection. Pour Q1 quadrilatère, la tolérance conduit au contraire à une projection sur la frontière.
- `examples/tetra_field_tools.rs` illustre projection, évaluation et prolongation ; `tests/tetra_field.rs` et `tests/tetra_locate.rs` couvrent notamment résultats analytiques, convergence, faces partagées, maillage disjoint, tolérance, géométrie déformée, transfert et entrées invalides.

## Avancées récentes sur les priorités

1. **Référence deal.II :** [protocole Q1](../comparison/dealii_q1/README.md) préparé pour comparer rigidité, charge, projection L² et erreur ; le comparateur Rust passe un jeu analytique. deal.II et CMake ne sont pas installés ici : le code C++ n’a **pas été compilé ni exécuté**, donc aucune comparaison indépendante avec deal.II n’est encore revendiquée.
2. **Interopérabilité :** import strict MSH 2.2 ASCII triangulaire, testé notamment avec un fichier généré par `meshio` ; export VTK Legacy d’un scalaire nodal et d’un scalaire par cellule. Le résultat adaptatif a été relu par `meshio` 5.3.5 : 287 points, 535 triangles, associations et valeurs finies vérifiées. Les groupes Gmsh et entiers arbitraires ne sont pas exportés par cette API.
3. **Raffinement local :** `LocalRefinement` produit des triangles conformes red-green et transfère les valeurs P1 ; `mark_bulk` fournit le marquage de Dörfler. Un exemple Poisson relie résolution, estimation, marquage et raffinement. Tests des 256 ensembles de cellules marquées sur un petit maillage, de la généalogie, de la continuité et de la résolution. `RefinementHistory` conserve les niveaux et retire le dernier niveau complet ; la restriction nodale par injection ne préserve pas les détails introduits aux nouveaux sommets. Ce raffinement red-green conforme n’a pas de hanging nodes, de déraffinement sélectif, ni de garantie de qualité des éléments.
4. **Performance :** un banc séquentiel mesure maillage, assemblage, contraintes et CG ; protocole et résultats dans [BENCHMARK_POISSON](BENCHMARK_POISSON.md). Il expose les octets logiques des tableaux CSR après Dirichlet, mais ne mesure pas la mémoire réelle ni le pic et ne démontre pas de gain. Pas de nouveaux backends externes.

## Dernier travail local : contraintes affines

`constraints::AffineConstraints` décrit des relations entre degrés de liberté et leurs éventuels décalages. La fermeture développe les dépendances sans récursion et détecte les cycles ; `condense` construit `Cᵀ A C` et `Cᵀ(b − A g)`, puis `reconstruct` restitue `x = C y + g`. Le cas entièrement contraint produit un système 0×0. Des tests couvrent des systèmes analytiques, non symétriques et SPD, les contraintes de Dirichlet, les cycles et erreurs numériques, ainsi qu’une chaîne de 20 000 dépendances.

L’exemple `affine_constraints` assemble explicitement un petit patch non conforme avec un nœud pendant et une solution affine. Il ne constitue pas une prise en charge générale des maillages non conformes : la topologie, le bord et les contraintes sont fournis manuellement, sans découverte automatique. L’expansion et le produit creux peuvent augmenter fortement le remplissage mémoire ; SPD et conditionnement ne sont pas vérifiés.

## Dernier travail local : nœuds pendants sur un niveau

`grid::NonconformingRefinement` raffine uniquement les cellules marquées d’un parent triangulaire conforme. Les interfaces grossière/fine sont identifiées par connectivité parentale et `hanging_constraints` génère les relations de milieu P1. La frontière est dérivée du parent ; les interfaces internes ne sont pas prises pour des arêtes extérieures. Ce type spécialisé ne se convertit pas en `Mesh2` et propose l’assemblage volumique Poisson, la prolongation et les relations parentales. L’exemple `poisson_hanging` applique les contraintes et Dirichlet pour une solution quadratique.

Les tests couvrent les 256 sélections d’un petit maillage, les frontières et normales, la continuité, des composantes disjointes, des données invalides, une solution affine et la convergence d’une solution quadratique. La prise en charge reste limitée à un niveau depuis un parent conforme ; pas de raffinement multiniveau non conforme, d’estimateur associé, de Neumann, de déraffinement sélectif ou de détection géométrique générale.

## Validation observée

- Le dernier commit poussé est `e2e8734` ; [sa CI GitHub](https://github.com/tiffank1802/dealii-rust/actions/runs/44892253304) a réussi pour Rust stable et Rust 1.70, ainsi que l’étape de style. **Cette CI ne couvre pas les changements locaux décrits ci-dessus.**
- Sur l’état local courant, `cargo test --offline --all-targets` a réussi : **273 tests**, dont 11 au raffinement non conforme à un niveau, 11 aux contraintes affines, 19 aux champs tétraédriques, 6 au raffinement local, 6 à la pile de niveaux, 8 à l’import Gmsh 2.2, 2 au marquage et 2 à l’export VTK. `cargo fmt --check`, Clippy strict, doctests et `cargo doc --no-deps` ont également réussi avec la toolchain AOSP Rust `1.83.0-dev`. Cette validation locale ne remplace pas une exécution stable ou 1.70 du nouvel état.
- Exemple tétraédrique sur le champ `x²+y²+z²` : erreurs L² de projection `3.120036e-2`, `8.008928e-3`, `2.013494e-3` pour `n=2,4,8` respectivement. Ces résultats illustrent ce cas précis, pas une garantie générale de performance.

## Limites et risques importants

- Les maillages fournis à `Mesh2` et comme parents de raffinements doivent être conformes et non recouvrants : les intersections, recouvrements et nœuds pendants ne sont pas tous détectés. Le type non conforme spécialisé n’accepte qu’un niveau issu d’un tel parent. La qualité du raffinement tétraédrique n’est pas garantie.
- Les classifications de points sur les frontières et les calculs sur cellules très minces ou fortement translatées sont sensibles aux arrondis. Les recherches de points sont linéaires, sans index spatial ; aucun gradient d’interface n’est fourni.
- La positivité des coefficients de diffusion-réaction n’est vérifiée qu’aux points de quadrature ; CG suppose un système SPD et la positivité diagonale seule ne le prouve pas. La convergence et le conditionnement ne sont pas garantis.
- L’import Gmsh est limité à des sous-ensembles MSH 4.1 et 2.2 ASCII, triangles/segments/points plans ; en 2.2, un seul groupe physique par élément, au plus deux tags et pas de partitions ; pas de quadrilatère ou 3D. L’export est VTK Legacy ASCII.
- Pas encore de comparaison indépendante exécutée avec deal.II sur les mêmes données, de mesure mémoire réelle ou maximale, de raffinement local Q1/3D, de champs vectoriels généraux ou de calcul distribué.

## Prochaines priorités proposées

1. Exécuter et valider le protocole de comparaison Q1 avec deal.II, puis l’étendre à des problèmes et maillages plus riches.
2. Étendre et tester l’interopérabilité au-delà du sous-ensemble Gmsh actuel.
3. Étendre la pile de niveaux au déraffinement sélectif et généraliser les contraintes et la topologie non conforme au multiniveau, avant l’adaptation `hp`.
4. Mesurer les performances et évaluer les backends et préconditionneurs complémentaires en fonction de besoins concrets.
