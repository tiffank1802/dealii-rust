# Dépendances externes, relations et validation transversale

[Retour à l’index des fiches](README.md) · [Plan général](../PLAN_MODULES_DEALII.md)

Ce document complète les 26 fiches. **Tout est au stade de proposition : aucune crate n’est ajoutée au dépôt, aucune version n’est verrouillée et aucune intégration n’a été compilée ou validée.** « Crate » désigne une dépendance Cargo ; un backend natif est une bibliothèque externe qui peut nécessiter FFI, compilation et installation système.

## 1. Règles de choix

- **P0 :** une seule crate interne, bibliothèque standard, sans bibliothèque native obligatoire. COO/CSR, CG et l’export VTK minimal sont limités aux besoins du MVP. Cette préférence est provisoire : une évaluation peut démontrer qu’une dépendance Rust est plus sûre ou économique.
- Chaque adoption exige un mini-prototype, la vérification des licences **et dépendances transitives**, de la maintenance, de la sécurité, des fonctionnalités et formats exacts, des cibles OS/CPU et de la MSRV. Consigner ces résultats avant de modifier `Cargo.toml`.
- Une version observée ci-dessous n’est **pas** une recommandation ferme. Choisir une contrainte Cargo et une MSRV après prototypage ; utiliser `Cargo.lock` pour CI/exemples reproductibles, sans présumer qu’il contraint les utilisateurs d’une bibliothèque publiée.
- Les backends natifs doivent être séparés derrière des fonctionnalités optionnelles, avec versions système, CI, comportement en cas d’absence et règles de distribution documentés. Ne jamais faire dépendre le cœur d’un outil spécialisé.
- Licence adoptée pour le code original du projet : **MIT OR Apache-2.0**, au choix de l’utilisateur. Les fichiers éventuels issus de deal.II, de ses exemples, de ses tests ou d’autres projets gardent leurs propres conditions ; vérifier les SPDX et notices avant toute copie. Les données de référence doivent avoir une provenance et des droits clairs. Voir aussi [les dépendances de deal.II](../DEPENDANCES_DEALII.md).

## 2. Catalogue des crates candidates

Les versions stables non retirées ci-dessous ont été lues dans les métadonnées de l’index public Cargo le 20 novembre 2026 ; quand elle est indiquée, la colonne Rust est la valeur `rust_version` **déclarée par la version**, pas une preuve de compilation. Le snapshot de l’index est [efbb236b4d0ea65cbb4ac5805cce0aaee42b81b9](https://github.com/rust-lang/crates.io-index/tree/efbb236b4d0ea65cbb4ac5805cce0aaee42b81b9). Les fonctionnalités, formats, licences et compatibilité **n’ont pas encore été audités**. Les liens pointent vers l’entrée d’index, pas vers une garantie de support.

| Crate candidate | Version observée / Rust déclaré | Rôle possible et périmètre | Décision avant adoption |
| --- | --- | --- | --- |
| [`faer`](https://github.com/rust-lang/crates.io-index/blob/efbb236b4d0ea65cbb4ac5805cce0aaee42b81b9/fa/er/faer) | 0.24.4 / 1.84.0 | `lac` : algèbre et solveurs ; prototype creux | Vérifier formats, méthodes, performances, MSRV et licence ; comparer à CSR du MVP. |
| [`sprs`](https://github.com/rust-lang/crates.io-index/blob/efbb236b4d0ea65cbb4ac5805cce0aaee42b81b9/sp/rs/sprs) | 0.11.6 / 1.88.0 | `lac` : stockage et opérations creuses | Mesurer assemblage/produit et vérifier maintenance et compatibilité des solveurs. |
| [`nalgebra`](https://github.com/rust-lang/crates.io-index/blob/efbb236b4d0ea65cbb4ac5805cce0aaee42b81b9/na/lg/nalgebra) | 0.35.0 / 1.89.0 | `base`, `fe`, `physics` : petites matrices et tenseurs | Mesurer coût et ergonomie ; ne pas présumer qu’elle fournit les éléments FEM. |
| [`rayon`](https://github.com/rust-lang/crates.io-index/blob/efbb236b4d0ea65cbb4ac5805cce0aaee42b81b9/ra/yo/rayon) | 1.12.0 / 1.80 | `meshworker`, `matrix_free` : tâches CPU parallèles | Prototype de réduction, ordre flottant, surcoût et compatibilité MSRV. |
| [`mshio`](https://github.com/rust-lang/crates.io-index/blob/efbb236b4d0ea65cbb4ac5805cce0aaee42b81b9/ms/hi/mshio) | 0.4.2 / non renseigné | `gmsh` : lecture éventuelle de `.msh` | Vérifier versions, variantes ASCII/binaires, marqueurs et cellules effectivement supportés. |
| [`vtkio`](https://github.com/rust-lang/crates.io-index/blob/efbb236b4d0ea65cbb4ac5805cce0aaee42b81b9/vt/ki/vtkio) | 0.6.3 / non renseigné | `vtk` : I/O étendu | Vérifier variantes Legacy/XML, associations des données et accès depuis un lecteur indépendant. |
| [`mpi`](https://github.com/rust-lang/crates.io-index/blob/efbb236b4d0ea65cbb4ac5805cce0aaee42b81b9/3/m/mpi) | 0.8.2 / 1.78 | `distributed`, `lac` : communications | Tester implémentations MPI, initialisation, ownership, erreurs et cibles ; dépendance système optionnelle. |
| [`rstar`](https://github.com/rust-lang/crates.io-index/blob/efbb236b4d0ea65cbb4ac5805cce0aaee42b81b9/rs/ta/rstar) | 0.14.0 / 1.90 | `arborx`, `non_matching`, `particles` : recherche spatiale | Comparer à recherche exhaustive, vérifier précision, concurrence et coût de mise à jour. |
| [`thiserror`](https://github.com/rust-lang/crates.io-index/blob/efbb236b4d0ea65cbb4ac5805cce0aaee42b81b9/th/is/thiserror) | 2.0.21 / 1.77 | Erreurs typées | Optionnel ; évaluer si le bénéfice justifie une dépendance alors que `std::error::Error` suffit. |
| [`serde`](https://github.com/rust-lang/crates.io-index/blob/efbb236b4d0ea65cbb4ac5805cce0aaee42b81b9/se/rd/serde) | 1.0.229 / 1.56 | Paramètres et formats de configuration | Ne pas confondre sérialisation de paramètres et format de sauvegarde stable de maillage. |
| [`tracing`](https://github.com/rust-lang/crates.io-index/blob/efbb236b4d0ea65cbb4ac5805cce0aaee42b81b9/tr/ac/tracing) | 0.1.44 / 1.65 | Diagnostics structurés | Évaluer le bruit, les coûts et l’activation facultative. |
| [`approx`](https://github.com/rust-lang/crates.io-index/blob/efbb236b4d0ea65cbb4ac5805cce0aaee42b81b9/ap/pr/approx) | 0.5.1 / non renseigné | Assertions numériques en tests | Choisir des tolérances explicites, éventuellement avec des helpers locaux pour P0. |
| [`proptest`](https://github.com/rust-lang/crates.io-index/blob/efbb236b4d0ea65cbb4ac5805cce0aaee42b81b9/pr/op/proptest) | 1.12.0 / 1.88 | Génération de cas géométriques et d’invariants | Contrôler graines, taille des cas et MSRV ; compléter, pas remplacer, les exemples analytiques. |
| [`criterion`](https://github.com/rust-lang/crates.io-index/blob/efbb236b4d0ea65cbb4ac5805cce0aaee42b81b9/cr/it/criterion) | 0.8.2 / 1.86 | Benchmarks séparés de la correction | Décrire machine, build et variabilité ; ne pas faire dépendre P0 des performances. |

Les dépendances de test/benchmark peuvent elles aussi relever la MSRV et doivent être contrôlées dans la matrice CI. Les candidats ne sont pas nécessairement tous compatibles avec une MSRV unique ; aucune MSRV précise n’est décidée ici.

### Backends natifs envisagés, sans crate/binding retenu

| Backend | Fiches | Point de décision |
| --- | --- | --- |
| PETSc / Trilinos | `lac`, `trilinos` | ABI, indices, variantes scalaires, MPI, licences et bibliothèques disponibles ; wrappers isolés. |
| p4est | `distributed` | Construction/partitionnement, contrats de propriété, MPI et licences. |
| SUNDIALS | `sundials` | Version API, `N_Vector`, callbacks, durée de vie et linkage. |
| Gmsh / VTK | `gmsh`, `vtk` | L’API native est distincte de la lecture/écriture de fichiers ; elle n’est pas requise pour P0. |
| CGAL / Open CASCADE | `cgal`, `opencascade` | Fonctions géométriques nécessaires, robustesse numérique, composants et licences. |
| ArborX / Kokkos | `arborx`, `matrix_free` | Besoin réel par rapport aux alternatives CPU Rust, coût FFI et accélérateurs. |
| ADOL-C / Sacado / SymEngine | `differentiation` | Évaluer les alternatives Rust et les modes de différentiation avant tout binding. |

## 3. Relations entre modules : contrats d’intégration

Les flèches indiquent « dépend de » au niveau fonctionnel ; elles ne prescrivent pas de crates distinctes.

| Consommateur | Fournisseur | Données et invariants traversant la frontière |
| --- | --- | --- |
| `grid` → `base` | points, indices, erreurs | Coordonnées finies ; topologie indépendante des solveurs. |
| `fe` → `base`, `grid` | cellule de référence, quadrature, géométrie | Ordre local, orientation, Jacobien et repère fixés. |
| `dofs` → `grid`, `fe` | entités, types FE | Numérotation cohérente ; invalider après changement du maillage. |
| `dofs` → `lac` | construction de contraintes | `dofs` construit les relations ; `lac` les représente et les ferme, sans importer `dofs`. |
| `numerics` / `meshworker` → `fe`, `dofs`, `lac`, `grid` | contributions locales | Dimensions locales/globales et insertion cohérentes ; frontières une seule fois. |
| `integrators` → `fe`, `base` | évaluation locale | Fonctions locales sans ownership de la matrice globale. |
| `hp` → `grid`, `fe`, `dofs`, `lac` | collections, hiérarchie, contraintes | Indices actifs et contraintes valides sur interfaces de cellules. |
| `matrix_free` → `fe`, `dofs`, `lac`, `grid` | action opérateur et cache | Même résultat que l’opérateur assemblé ; invalider le cache après mutation. |
| `multigrid` → `grid`, `dofs`, `lac` | niveaux et transferts | Dimensions, contraintes, propriété de chaque niveau ; matrix-free optionnel. |
| `distributed` → `grid`, `dofs`, `lac`, MPI | propriétaires/fantômes | Une seule source de vérité par donnée ; échanges et ordre d’appel documentés. |
| `non_matching` / `particles` → `grid`, `fe`, recherche | coordonnées et localisation | Comportement défini sur frontières, hors domaine et lors de changements de maillage. |
| Adaptateurs → cœur | conversions, ressources | Le cœur ne dépend pas des intégrations ; erreurs et ownership traversent la frontière de façon sûre. |

`physics`, `algorithms`, `optimization` et `differentiation` s’appuient sur les contrats numériques nécessaires, sans imposer leurs backends aux modules fondamentaux. Les liens des 26 fiches donnent les dépendances propres et consommateurs.

## 4. Matrice de tests et critères mesurables

| Niveau | Cas minimal | Critère / exécution prévue |
| --- | --- | --- |
| Unitaire | Triangle de référence `(0,0),(1,0),(0,1)` | Aire 1/2 ; formes P1 nodales ; rigidité `[[1,-1/2,-1/2],[-1/2,1/2,0],[-1/2,0,1/2]]` avec tolérance flottante documentée. |
| Validation | Indices invalides, géométrie dégénérée, dimensions incompatibles | Erreur explicite, pas de résultat silencieusement invalide ; tests des NaN lorsque pertinent. |
| Intégration P0 | Carré unité, Poisson, maillages raffinés | `u=sin(πx)sin(πy)`, `f=2π²u` ; pente L2 proche de 2 et semi-norme H1 proche de 1 ; fixer maillages et tolérances après caractérisation. |
| Solveurs | Matrice SPD connue, non-SPD/rupture, itérations insuffisantes | Résidu recalculé indépendamment et diagnostics cohérents ; CG ne garantit une solution que sous ses préconditions SPD. |
| Contraintes | Dirichlet inhomogène et chaîne de contraintes | Comparer à `CᵀAC y=Cᵀ(b−Ag)` et `x=Cy+g` ; rejeter cycles et indices invalides. |
| Référence externe | Cas deal.II versionné | Comparer quantités invariantes, pas nécessairement la numérotation des DoFs ; documenter commande, commit, paramètres et tolérances. |
| I/O | Triangle, marqueurs et champ nodal | Import/export selon versions et types annoncés ; vérifier fichier avec outil indépendant, échecs informatifs. |
| Adaptation | Maillage à interface non conforme | Continuité, conservation et transfert de constantes ; vérifier invalidation des données antérieures. |
| Parallèle | 1, 2 et N rangs / plusieurs threads | Une propriété globale cohérente, conservation et tolérance numérique ; pas d’exigence bit à bit sauf décision spécifique. |
| Performance | Maillages et ordres croissants | Temps et mémoire avec CPU, compilateur, configuration et variance ; pas de seuil chiffré inventé avant baseline. |

**CI configurée :** formatage, lint, documentation, tests et exemple Poisson sur Linux, avec Rust stable et une MSRV cible 1.70 ; l’exécution de ces jobs reste à confirmer. Prévoir des jobs distincts et des versions épinglées pour les intégrations natives et MPI lorsqu’elles existeront. Les tests numériques ne remplacent pas les comparaisons de convergence ; les benchmarks restent distincts des tests de correction.

## 5. Ordre de réalisation et frontières de livraison

1. Initialisation, déclaration de la licence approuvée, politique MSRV et CI ; pas de dépendances externes imposées au cœur.
2. `base` → `grid` et `fe` → `dofs` + `lac` → `numerics`/`vtk` ; livrer le MVP complet avant d’élargir l’API.
3. Formats et types cellulaires supplémentaires, champs vectoriels, `integrators`/`meshworker`, élasticité linéaire si l’application prioritaire est confirmée.
4. Contraintes, raffinement et `hp` ; en parallèle après leurs prérequis, DG, temps/non-linéaire et multigrille.
5. Parallélisme, MPI, matrix-free ; puis géométrie spécialisée, particules et backends natifs en réponse à des cas utilisateurs et budgets d’entretien identifiés.

La liste précise, les jalons et les critères par sous-module se trouvent dans les fiches individuelles. Aucune fonctionnalité future n’est garantie avant son prototype et ses tests.
