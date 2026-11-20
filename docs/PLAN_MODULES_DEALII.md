# Cartographie détaillée et plan des modules deal.II → Rust

Voir aussi les [26 fiches détaillées par module et sous-chantier](modules/README.md) et leur [catalogue de dépendances, relations et tests](modules/DEPENDANCES_ET_TESTS.md).

## 1. Référence et méthode

Analyse de [deal.II](https://github.com/dealii/dealii), branche amont `master`, au commit **`1d230ad0be5ab2a3913cf4de8187301430721ed5`** (2026-11-20, version `9.9.0-pre`). C’est un instantané de développement, pas une version stable. Le clone superficiel est dans `.cache/dealii` ; ce cache est ignoré par Git et peut ne pas être conservé entre sessions. Pour le reconstruire depuis la racine :

```sh
git clone --depth 1 --single-branch https://github.com/dealii/dealii.git .cache/dealii
git -C .cache/dealii fetch --depth 1 origin 1d230ad0be5ab2a3913cf4de8187301430721ed5
git -C .cache/dealii checkout --detach 1d230ad0be5ab2a3913cf4de8187301430721ed5
```

L’inventaire porte sur les **26 répertoires de premier niveau** de [`include/deal.II`](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II) et sur leurs en-têtes `*.h` récursifs. Les regroupements et noms ci-dessous résultent de l’examen des en-têtes, de `source/`, des catégories de tests et de la documentation des exemples ; ils ne prétendent pas reproduire exactement les groupes Doxygen. La [liste exhaustive des 689 en-têtes](INVENTAIRE_ENTETES_DEALII.md) inclut internes et fichiers templates et permet de retrouver chaque fichier. Les fichiers générés par configuration et les extensions autres que `.h` ne sont pas comptés.

À ne pas confondre avec ces 26 modules : `module/` contient l’infrastructure des modules C++ et leurs partitions d’intégration ; `cmake/` et `bundled/` concernent le build et les dépendances, `doc/`, `examples/` et `tests/` la documentation et la validation. On prévoit leur **fonction**, pas une traduction de CMake, des mécanismes de templates C++ ou des copies de dépendances tierces.

La bibliothèque amont est sous double licence **Apache-2.0 avec LLVM-exception / LGPL-2.1-or-later** ; certains fichiers, notamment dans `bundled/`, `examples/` et `tests/`, ont des licences distinctes. Vérifier leurs en-têtes SPDX et les conditions applicables avant toute réutilisation. Voir [LICENSE.md amont](https://github.com/dealii/dealii/blob/1d230ad0be5ab2a3913cf4de8187301430721ed5/LICENSE.md). Une inspiration d’architecture n’implique pas de copier du code.

## 2. Niveaux, dépendances et principes

**Niveaux de livraison** — ils indiquent le premier objectif, pas une promesse de parité :

- **P0** : J0–J4 du [plan initial](../PLAN_IMPLEMENTATION.md) : Poisson scalaire 2D P1 sur triangles, séquentiel ; une crate Cargo avec modules internes.
- **P1** : généralisation séquentielle : 2D/3D, Q1, géométrie et frontières, éléments vectoriels, I/O et solveurs usuels. Chaque famille est livrée avec son propre exemple et ses tests.
- **P2** : adaptation locale et `hp`, DG, interfaces, contraintes, multigrille, non-linéaire et temps ; ces chantiers sont indépendants quand leurs prérequis sont satisfaits.
- **P3** : parallélisme partagé et distribué, matrix-free, accélérateurs et passage à l’échelle ; une partie peut être prototypée avant P2 mais ne doit pas bloquer P0.
- **P4** : domaines spécialisés et parité fonctionnelle élargie (géométries non conformes, particules, familles avancées).
- **X** : intégration optionnelle d’un outil externe, à décider selon besoin, licence et plateformes. Un X dépend du niveau fonctionnel concerné et ne bloque pas le cœur.

Pour chaque section : les chemins indiqués sont **des exemples représentatifs**, non une correspondance un-pour-un. La liste complète des fichiers est dans l’inventaire. Les priorités sont au niveau des sous-chantiers ; « critères » fixe la condition minimale pour annoncer le sous-champ concerné comme pris en charge. Les détails comme les `*.templates.h`, caches, sélecteurs et helpers internes sont couverts par les sous-chantiers correspondants ou remplacés par des mécanismes Rust, pas nécessairement portés comme API.

### Dépendances fonctionnelles principales

```text
base (points/tensors/quadrature/fonctions)
  ├─ grid (cellules/géométrie/topologie) ── distributed, hp, particles
  ├─ fe (éléments/mappings/FEValues) ───────┬─ hp, non_matching
  └─ lac (vecteurs/matrices/solveurs)       │
grid + fe → dofs → contraintes(lac) ────────┤
grid + fe + dofs + lac → assemblage/numerics, integrators, meshworker
ces couches → algorithms, physics, optimization, differentiation
assemblage + lac + niveaux de grille → multigrid
fe + dofs + lac + contraintes → matrix_free → multigrid matrix-free
MPI + grid + dofs + lac → distributed → échanges, transferts et solveurs distribués
adaptateurs externes → chaque sous-système concerné ; jamais prérequis de P0
```

`AffineConstraints` est dans `lac` en amont, mais les contraintes d’interpolation dépendent des DoFs et de la géométrie ; éviter une dépendance circulaire en séparant représentation algébrique et construction FEM. De même `hp`, `distributed` et `matrix_free` traversent plusieurs couches : ne pas les construire comme une simple chaîne de modules.

**Stratégie Rust :** privilégier les types concrets et identifiants indexés pour commencer ; extraire des traits lorsque plusieurs implémentations les motivent. Séparer concepts, orchestration et backends ; limiter `unsafe`, encapsuler les FFI, documenter les invariants, les états fantômes et la validité des caches/identifiants après raffinement. Les modules ci-dessous sont une cible logique : scinder en crates seulement après stabilisation des frontières. La prise en charge de chaque dimension, type cellulaire, famille FE, nombre scalaire et backend doit être explicite.

## 3. Plan par module amont

### 3.1 `base` — fondamentaux

Références : `point.h`, `tensor.h`, `quadrature*.h`, `polynomial*.h`, `function*.h`, `index_set.h`, `mpi*.h`, `work_stream.h`, `parameter_handler.h`.

- **P0 — géométrie et erreurs :** points/vecteurs 2D, opérations géométriques, types d’indices, erreurs typées, contrôles de dimensions et de finitude.
- **P0 → P1 — quadrature et polynômes :** points/poids sur triangle puis règles de Gauss, produits tensoriels et autres cellules ; bases polynomiales et dérivées, avec conventions et degré d’exactitude.
- **P1 — fonctions et tenseurs :** fonctions scalaires/vectorielles, temps, tenseurs et tenseurs symétriques ; primitives de tableaux et vues si utiles. Les fonctions spécialisées sont ajoutées à la demande.
- **P1 → P3 — infrastructure :** ensembles d’indices, temporisation, paramètres validés, partitionnement et vues locales/fantômes ; orchestration de tâches, SIMD et MPI seulement avec un besoin mesuré.
- **P4/X — services spécialisés :** parsers, HDF5, dérivation symbolique et utilitaires de communication ; utiliser des crates/adaptateurs plutôt que traduire les couches de compatibilité C++ (`std_cxx*`, pointeurs observateurs, macros).
- **Critères :** quadrature exacte sur monômes annoncés, calculs tensoriels de référence, validation des paramètres ; tests d’indices et de communications séparés pour les fonctionnalités parallèles.

### 3.2 `grid` — maillages et géométrie

Références : `tria.h`, `tria_accessor.h`, `cell_id.h`, `reference_cell.h`, `grid_generator.h`, `grid_in.h`, `grid_refinement.h`, `manifold*.h`.

- **P0 :** sommets, triangles conformes, générateur du carré unité, arêtes et marqueurs de frontière ; accès par IDs et itérateurs, validation de connectivité et orientation.
- **P1 :** cellules de référence triangle/quadrilatère/tétraèdre puis hexaèdre, voisinage et orientations de faces, matériaux et marqueurs ; primitives d’import/export et transformations géométriques.
- **P2 :** hiérarchie parent/enfant, raffinement/déraffinement, identifiants persistants et transfert ; rechercher voisins et faces non conformes, poser des règles explicites d’invalidation.
- **P2 → P4 :** géométries/manifolds, projection et interpolation sur surfaces, caches spatiaux, géométrie courbe et codimension ; pyramides/prismes seulement après spécification complète des conventions.
- **Critères :** tests de validité, aire/volume, orientations et frontières ; puis conservation topologique, raffinements aller-retour et géométrie courbe avec tolérances documentées.

### 3.3 `fe` — éléments, mappings et évaluation

Références : `fe.h`, `fe_simplex_p.h`, `fe_q.h`, `fe_system.h`, `fe_values*.h`, `fe_interface_values.h`, `mapping*.h`, familles `fe_*`.

- **P0 :** élément de référence P1 triangle, mapping affine, valeurs et gradients physiques ; petit équivalent de `FEValues` avec données réellement demandées.
- **P1 :** Q1 et P1 tétraèdre, mapping Q1 et quadrature associée ; Lagrange d’ordre supérieur, champs vectoriels/composants, évaluation de faces et normales.
- **P2 :** éléments discontinus (`fe_dg*`), évaluation d’interfaces et sauts, éléments mixtes, interpolation et transferts. Définir explicitement conventions de faces et de sous-faces.
- **P4 :** H(div) Raviart–Thomas/BDM, H(curl) Nédélec et autres familles (`abf`, Hermite, Bernstein, enrichis, trace, pyramide/prisme) par lots distincts ; transformations de Piola, orientations et degrés de liberté associés aux entités.
- **Critères :** nodalité/partition de l’unité, gradients, mappings, quadrature et convergence par cellule/famille ; tests d’orientation, continuité tangentielle/normale et diagrammes de commutation pour H(curl)/H(div). Aucune famille n’est considérée supportée sur simple existence d’un trait.

### 3.4 `dofs` — numérotation et accès

Références : `dof_handler.h`, `dof_accessor.h`, `dof_tools.h`, `dof_renumbering.h`, `block_info.h`.

- **P0 :** numérotation P1 par sommet, local → global, DoFs de frontière et cohérence entre cellules.
- **P1 :** DoFs associés aux entités, composantes/blocs, orientation de faces et renumérotation ; masque de composantes et énumération par type d’élément.
- **P2 :** numérotation par niveau, éléments hétérogènes et `hp`, construction des contraintes de continuité et interpolation vers `lac`.
- **P3 :** propriété locale, indices fantômes, numérotation distribuée et synchronisation après adaptation.
- **Critères :** DoFs partagés, permutations/orientations, unicité et invariance de la solution après renumérotation ; cohérence multi-rangs pour le distribué.

### 3.5 `lac` — algèbre linéaire et contraintes

Références : `vector.h`, `sparse_matrix.h`, `sparsity_pattern.h`, `affine_constraints.h`, `solver_*.h`, `precondition*.h`, `linear_operator.h`, adaptateurs PETSc/Trilinos.

- **P0 :** vecteurs, assemblage COO puis CSR, fusion des doublons, produit matrice-vecteur, CG et diagnostics ; élimination symétrique des Dirichlet avec correction du second membre.
- **P1 :** matrices denses/locales, motifs creux, blocs, opérateurs et préconditionneurs Jacobi ; GMRES/BiCGStab et solveurs directs au besoin. Évaluer des backends Rust avant de les choisir.
- **P2 :** contraintes affines `x = Cy + g`, fermeture et cycles, condensation et reconstruction ; opérateurs de blocs, Schur et préconditionneurs adaptés aux systèmes mixtes.
- **P3/X :** vecteurs distribués, échanges propriétaire/fantôme, solveurs et préconditionneurs via PETSc, Trilinos ou autres backends optionnels ; GPU seulement si demandé et mesuré.
- **P4/X :** problèmes aux valeurs propres, matrices complexes et backends ARPACK/SLEPc, ScaLAPACK et autres disponibles en amont, chacun à évaluer séparément.
- **Critères :** dimensions et doublons, conservation de symétrie, systèmes de référence, résidus, échec/non-convergence ; tests distribués de propriété et communication. CG exige SPD.

### 3.6 `numerics` — opérations FEM et post-traitement

Références : `vector_tools_*`, `matrix_creator.h`, `error_estimator.h`, `solution_transfer.h`, `data_out*.h`, `adaptation_strategies.h`.

- **P0 :** assemblage Poisson, second membre, Dirichlet, erreurs L2/H1 et export VTK du champ.
- **P1 :** interpolation, projection, évaluation ponctuelle, intégration, calcul de moyennes, export des données de cellules/faces ; conditions de Neumann et Dirichlet non homogènes.
- **P2 :** estimateurs a posteriori et de régularité, marquage, transfert de solution lors de l’adaptation, dérivées et postprocesseurs ; adaptation de maillage pilotée par erreurs.
- **P3/P4 :** sortie parallèle, échantillonnage, suivi de points et cas de codimension ; définir localisations et échec des requêtes hors domaine.
- **Critères :** solution manufacturée et pentes de convergence, comparaison à deal.II ; vérification de conservation/transfert et décroissance de l’erreur pour l’adaptation.

### 3.7 `hp` — collections et adaptation en ordre

Références : `fe_collection.h`, `mapping_collection.h`, `q_collection.h`, `fe_values.h`, `refinement.h`.

- **P2 :** collections d’éléments, mappings et quadratures, choix par cellule et évaluation selon indices actifs.
- **P2 → P4 :** sélection h/p, interfaces entre ordres différents, interpolation et contraintes cohérentes ; raffinements combinés, transferts et hiérarchies.
- **Prérequis :** plusieurs ordres FE, géométrie et DoFs par entité, adaptation et contraintes affines fermées.
- **Critères :** continuité sur interfaces h/p, reproduction polynomiale, transferts et convergence sur cas manufacturés ; tests d’orientation sur chaque cellule annoncée.

### 3.8 `integrators` — noyaux d’intégration locale

Références : `laplace.h`, `l2.h`, `advection.h`, `divergence.h`, `elasticity.h`, `maxwell.h`.

- **P1 :** formes locales masse et Laplace, interfaces indépendantes du stockage global ; comparer aux assemblages de référence.
- **P2 :** advection, divergence, grad-div et élasticité ; termes de frontière/interface DG et quadrature adaptée.
- **P4 :** Maxwell et noyaux spécialisés, selon familles FE prises en charge.
- **Critères :** matrices locales analytiques, symétrie/coercivité lorsque pertinentes et conservation des flux ; tests par forme et cellule.

### 3.9 `meshworker` — parcours d’assemblage

Références : `mesh_loop.h`, `scratch_data.h`, `copy_data.h`, `assembler.h`, `integration_info.h`.

- **P1 :** parcours cellulaire séparant données temporaires locales, contribution et insertion ; gestion des erreurs et politique d’ordre définie.
- **P2 :** parcours des faces intérieures/frontières et sous-faces, assemblage DG et contraintes ; éviter de compter deux fois les interfaces.
- **P3 :** tâches indépendantes et réduction des contributions, allocation locale réutilisable, compromis reproductibilité/performance.
- **Critères :** égalité avec assemblage direct sur petits maillages ; termes de faces et résultats séquentiels/parallèles comparés avec tolérances.

### 3.10 `matrix_free` — opérateurs sans matrice

Références : `matrix_free.h`, `fe_evaluation.h`, `fe_point_evaluation.h`, `operators.h`, `tensor_product_kernels.h`, `portable_*`.

- **P3 :** opérateur Laplace sans matrice sur cellules tensorisées, évaluation par lots et somme factorisée ; indices, géométrie, contraintes et échanges cohérents.
- **P3 → P4 :** intégration de faces/DG, noeuds contraints, évaluation ponctuelle distante et cache de mapping ; élargir progressivement éléments et cellules.
- **X :** chemins portables/accélérateurs par backend séparé, après validation d’une API et de benchmarks ; ne pas supposer un équivalent direct à Kokkos en Rust.
- **Critères :** action comparée à la matrice assemblée sur vecteurs aléatoires, contraintes et frontières incluses ; tests parallèles, consommation mémoire et benchmarks reproductibles.

### 3.11 `multigrid` — solveurs multiniveaux

Références : `multigrid.h`, `mg_transfer*.h`, `mg_smoother.h`, `mg_coarse.h`, `mg_constrained_dofs.h`.

- **P2 :** niveaux géométriques, prolongation/restriction, opérateurs de niveau, lisseurs, solveur grossier et cycle V ; contraintes par niveau.
- **P3 :** transferts distribués, coarsening global, blocs et opérateurs matrix-free ; prendre en compte communications et propriété des niveaux.
- **Critères :** exactitude des transferts sur polynômes, comparaison à un solveur direct sur petits systèmes, résidus et nombre d’itérations sous raffinement ; tests de cohérence distribuée.

### 3.12 `distributed` — maillage parallèle

Références : `tria.h`, `shared_tria.h`, `fully_distributed_tria.h`, `grid_refinement.h`, `cell_data_transfer.h`, `p4est_wrappers.h`.

- **P3 :** modèle de propriété des cellules, fantômes/artificielles, partitions et identifiants globaux ; maillage partagé puis véritablement distribué selon besoins.
- **P3 :** raffinement et répartition, transferts de données et de champs, synchronisation, charge et sauvegarde ; tester au moins 1, 2 et plusieurs rangs et partitions déséquilibrées.
- **X :** évaluer p4est ou autre infrastructure par adaptateur isolé ; définir format d’échange, compilation native, licences et architectures supportées.
- **Critères :** même solution et quantité conservée que le séquentiel, cohérence des fantômes et IDs, raffinement et repartitionnement ; résultats et tolérances reproductibles.

### 3.13 `algorithms` — orchestration numérique

Références : `newton.h`, `operator.h`, `theta_timestepping.h`, `timestep_control.h`, `any_data.h`.

- **P1/P2 :** interface de résidu/Jacobien, Newton et gestion d’échec ; schéma θ et boucle temporelle avec gestion explicite du pas.
- **P2 :** contrôle des tolérances et pas, callbacks/événements et collecte des diagnostics ; interfaces de données typées plutôt que reproduire `AnyData` tel quel.
- **Critères :** résolution de systèmes non linéaires analytiques, test de divergence et de non-convergence ; ordre temporel vérifié sur ODE/équation de la chaleur.

### 3.14 `differentiation` — différentiation automatique et symbolique

Références : `ad/` (drivers, helpers, ADOL-C, Sacado) et `sd/` (SymEngine).

- **P2 :** prototyper une interface dérivée/Jacobien sur un problème local ; évaluer AD avant de choisir un modèle de nombres ou un backend.
- **P4/X :** Hessiennes, dérivation symbolique, compilation/optimisation d’expressions et adaptateurs vers bibliothèques externes si nécessaires.
- **Critères :** comparaison aux dérivées analytiques et différences finies avec tolérances ; validation du Jacobien d’un problème non linéaire complet. Documenter les limites de scalaires, dimensions et threading.

### 3.15 `physics` — outils de mécanique et transformations

Références : `elasticity/kinematics.h`, `elasticity/standard_tensors.h`, `notation.h`, `transformations.h`, `vector_relations.h`.

- **P1/P2 :** transformations de tenseurs, invariants et cinématique de l’élasticité ; conventions de dimensions et mesures explicites.
- **P2/P4 :** tenseurs constitutifs, exemples d’élasticité linéaire puis non linéaire ; ne pas confondre outils de physique et solveur d’application général.
- **Critères :** identités tensorielles et tests de référence analytique ; patch test, symétries et objectivité lorsque la formulation l’exige.

### 3.16 `optimization` — minimisation

Références : `line_minimization.h`, `solver_bfgs.h`, `rol/vector_adaptor.h`.

- **P2 :** recherche linéaire et BFGS avec critères d’arrêt, rapport d’itérations et échecs explicites.
- **P4/X :** contraintes d’optimisation, adaptation de vecteurs et backend ROL si le cas d’usage le justifie ; en amont ROL est aussi présent dans `trilinos`.
- **Critères :** minima analytiques, convergence et tests d’échec, gradients comparés à AD ou différences finies ; exemples de calibration FEM séparés.

### 3.17 `non_matching` — maillages non concordants

Références : `mesh_classifier.h`, `quadrature_generator.h`, `fe_immersed_values.h`, `coupling.h`, `closest_surface_point.h`.

- **P4 :** classifier cellules coupées, rechercher intersections et générer quadratures de domaines/surfaces immergés ; gérer les configurations dégénérées.
- **P4 :** évaluer et coupler champs sur maillages différents, recherche de points, mapping et échanges distribués éventuels.
- **Critères :** aires/volumes et intégrales sur géométries analytiques, conservation du couplage et convergence sous raffinement ; tests de coupures presque dégénérées.

### 3.18 `particles` — particules et champs

Références : `particle.h`, `particle_handler.h`, `generators.h`, `property_pool.h`, `partitioner.h`, `data_out.h`.

- **P4 :** identifiants, positions, propriétés, générateurs et association particule-cellule ; interpolation du champ et sortie.
- **P4/P3 :** déplacement, réinsertion après raffinement, migration entre rangs et gestion des particules hors domaine.
- **Critères :** conservation du nombre et des propriétés, positionnement sur cellules/frontières, migration et reproductibilité sur plusieurs rangs.

### 3.19 `arborx` — recherche spatiale accélérée

Références : `bvh.h`, `distributed_tree.h`, `access_traits.h`.

- **P4/X :** interface de recherche de boîtes englobantes, intersections et plus proches voisins ; choisir entre index Rust et backend ArborX optionnel.
- **P3/X :** recherche distribuée et gestion des données distantes si le besoin existe ; coordonner avec géométrie et particules.
- **Critères :** résultats comparés à une recherche exhaustive, cas dégénérés et tests multi-rangs ; mesures de performance et prérequis du backend.

### 3.20 `boost_adaptors` — adaptation Boost.Geometry

Références : `point.h`, `segment.h`, `bounding_box.h`.

- **P4/X :** définir les conversions nécessaires pour une bibliothèque géométrique choisie ; privilégier des traits/types Rust au portage direct de Boost.
- **Critères :** conversions aller-retour, unités, dimensions et conventions de coordonnées ; aucun backend Boost requis par P0.

### 3.21 `cgal` — géométrie algorithmique

Références : `intersections.h`, `triangulation.h`, `surface_mesh.h`, `polygon.h`, `point_conversion.h`.

- **P4/X :** inventaire des opérations nécessaires : conversions, polygones, maillages de surface, intersections et triangulation ; comparer alternatives Rust et adaptateur CGAL.
- **Critères :** cas géométriques connus et dégénérés, précision et validité des maillages ; documenter compilation native et licence du backend retenu.

### 3.22 `opencascade` — CAD et manifolds

Références : `boundary_lib.h`, `manifold_lib.h`, `utilities.h`.

- **P4/X :** importer ou interroger une géométrie CAD, projection et normales ; adapter les manifolds aux courbes et surfaces via un backend optionnel.
- **Critères :** points projetés sur surfaces connues, tolérances et erreurs de projection, raffinement géométrique ; documenter versions et dépendances natives.

### 3.23 `gmsh` — intégration Gmsh

Référence : `gmsh/utilities.h` ; l’I/O générique amont se trouve aussi dans `grid/grid_in.h` et `grid/grid_out.h`.

- **P1/X :** importer un sous-ensemble explicitement versionné de `.msh`, types cellulaires et marqueurs physiques ; évaluer une crate de lecture avant le choix.
- **P4/X :** intégration de l’API native Gmsh si génération ou fonctions non couvertes par l’I/O sont nécessaires ; la distinguer du simple parsing de fichiers.
- **Critères :** fichiers de référence, frontières et marqueurs, versions/éléments non supportés retournant une erreur informative.

### 3.24 `vtk` — intégration VTK

Référence : `vtk/utilities.h` ; la sortie FEM est également dans `numerics/data_out*.h` et `base/data_out_base.h`.

- **P0 :** export limité VTK du maillage et du champ du MVP, avec association point/cellule explicitée.
- **P1/X :** évaluer `vtkio` ou backend pour lecture/écriture, formats, données vectorielles et maillages pris en charge ; sortie parallèle en P3 si nécessaire.
- **Critères :** vérifier structure, association et valeurs sur un petit exemple ; lecture avec un outil indépendant, tests d’aller-retour lorsque supportés.

### 3.25 `sundials` — intégrateurs/solveurs externes

Références : `arkode*.h`, `ida.h`, `kinsol.h`, `n_vector.h`, wrappers `sunlinsol` et `sunnonlinsol`.

- **P4/X :** adaptateur optionnel ARKODE pour ODE, IDA pour DAE et KINSOL pour non-linéaire, chacun activable/testable séparément selon besoin.
- **Prérequis :** contrats de vecteurs, résidus, Jacobien, mémoire et durée de vie ; FFI sûre, versions natives et erreurs documentées.
- **Critères :** problèmes ODE/DAE analytiques et ordre de convergence, tests d’erreurs et de destruction des ressources ; matrice CI des versions supportées.

### 3.26 `trilinos` — adaptateurs complémentaires

Références : `nox.h`, `rol_adaptor.h` ; beaucoup de wrappers linéaires Trilinos sont en réalité dans `lac/trilinos_*`.

- **P4/X :** évaluer NOX pour systèmes non linéaires et ROL pour optimisation ; coordonner avec `algorithms`, `optimization` et les backends `lac`.
- **Critères :** convergence sur problème de référence, propagation des erreurs et nettoyage FFI ; versions et combinaisons de bibliothèques testées, sans dépendance obligatoire.

## 4. Décomposition en lots et portes de validation

| Lot | Contenu concret | Prérequis / porte de sortie |
| --- | --- | --- |
| L0 | Cargo, CI, licence du projet, conventions et test de référence | MSRV et stable ; aucun backend natif requis |
| L1 | `base` géométrie/quadrature, `grid` triangles, `fe` P1/mapping, `dofs`, `lac` CSR/CG, `numerics` Poisson/VTK | Critères P0 et convergence du [MVP](../PLAN_IMPLEMENTATION.md#2-premier-objectif-livré--poisson-2d) |
| L2 | Q1/3D, champs et familles Lagrange, faces, import/export, conditions limites, noyaux locaux et solveurs usuels | Convergence par cellule/dimension, compatibilité I/O et tests unitaires de chaque interface |
| L3 | Raffinement local et contraintes, estimateurs, transferts, puis `hp` ; DG et faces comme branche parallèle après L2 | Continuité/non-conformité, reconstruction des contraintes, convergence et conservation des transferts |
| L4 | Non-linéaire, temps, élasticité, optimisation/AD au besoin ; multigrille après niveaux et transferts | Ordre temporel, patch tests, Jacobien et diagnostics solveurs ; tests de multigrille |
| L5 | Assemblage partagé, `matrix_free`, données distribuées/MPI, multigrille distribuée | Équivalence numérique, performances mesurées et tests multi-rangs ; l’ordre entre ces travaux dépend des prototypes |
| L6 | Particules, non-matching, familles FE avancées et intégrations CAD/géométriques/solveurs | Cas d’usage approuvé et critères par sous-module ; backends optionnels isolés |

Chaque lot doit être découpé en tickets portant : résultat utilisateur, type de cellule et dimensions, API concernée, prérequis, tests/benchmark, documentation, définition de terminé et dépendances externes. Il n’y a pas d’estimation calendaire fiable sans équipe, financement et profils d’application. Le parcours Poisson P0 est en cours d’implémentation et de validation ; le [README](../README.md) précise l’état actuel. Les autres capacités restent planifiées.

## 5. Validation transversale et couverture

- **Référence figée :** conserver commit, configurations, maillages, formulations, versions d’outils et tolérances des comparaisons ; ne pas rendre l’installation de deal.II obligatoire pour les tests courants.
- **Échelle des tests :** calcul local analytique → petit assemblage → solution manufacturée et convergence → exemple complet → comparaison avec deal.II. Ajouter tests des erreurs, limites numériques et conditions aux frontières à chaque niveau.
- **Exemples amont de départ :** `step-3` (Poisson), `step-6` (adaptativité), `step-12` (DG), `step-15` (non-linéaire), `step-22` (Stokes), `step-26` (chaleur), `step-37` (matrix-free), `step-40` (MPI), `step-60` (couplage immergé). Ce sont des références de scénarios à étudier, pas des cibles automatiquement identiques ni des fichiers à recopier sans contrôle de licence.
- **Distribution et parallélisme :** comparer 1/2/N rangs, partitions, contraintes et fantômes ; les sommes flottantes peuvent varier selon l’ordre, donc définir tolérances et politique de reproductibilité.
- **Architecture/qualité :** tests des invariants d’IDs et invalidation des caches, tests d’intégration avec fonctionnalités optionnelles, documentation des préconditions et matrice des versions de backends natifs.
- **Mesures :** compiler temps d’assemblage/résolution, mémoire, taille du problème, machine et configurations ; distinguer correction, robustesse et performance.

## 6. Décisions à prendre avant d’engager les lots avancés

1. Quel domaine prime après Poisson : mécanique, fluides, électromagnétisme ou autre ? Cela détermine l’ordre des familles FE, solveurs et adaptateurs.
2. Quelles cibles (OS, CPU/GPU, nombre de rangs MPI) et quel budget de maintenance pour les dépendances natives ?
3. La licence du code original est choisie : `MIT OR Apache-2.0`. Comment gérer et vérifier les licences distinctes des dépendances, fichiers amont et données de test ?
4. Quels formats et versions d’I/O sont contractuels ? Quels backends sont suffisamment maintenus et compatibles avec la MSRV choisie ?
5. À quel niveau vise-t-on la compatibilité : équivalence numérique de scénarios, concepts d’API, ou migration assistée ? La compatibilité C++ n’est pas promise.
