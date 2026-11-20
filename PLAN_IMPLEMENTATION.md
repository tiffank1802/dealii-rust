# Plan d’implémentation de deal.II en Rust

Pour l’analyse détaillée des 26 modules amont et leurs sous-chantiers, voir [le plan par module](docs/PLAN_MODULES_DEALII.md) et [l’inventaire des en-têtes](docs/INVENTAIRE_ENTETES_DEALII.md).

## 1. Objectif et positionnement

Construire progressivement une bibliothèque Rust de simulation par éléments finis, inspirée des concepts et des tutoriels de [deal.II](https://www.dealii.org/), avec une API idiomatique, testable et performante. **La parité avec deal.II est un objectif de long terme, pas la définition de la première version.** Une réimplémentation complète implique notamment de nombreuses familles d’éléments, la 3D, le raffinement adaptatif, les problèmes multiphysiques, les solveurs et le calcul distribué.

Le dépôt contient désormais une première crate et un parcours Poisson 2D expérimental. Le [README](README.md) décrit le périmètre réellement implémenté ; les étapes ci-dessous restent la feuille de route et leurs critères complets, notamment la validation CI sur Rust stable et MSRV, ne sont pas encore tous confirmés.

### Principes

- Livrer d’abord un **parcours complet et vérifiable**, puis étendre les abstractions à partir de cas réels.
- Séparer géométrie, topologie, discrétisation et algèbre ; préserver l’accès aux données de bas niveau.
- Viser une API Rust sûre et documentée ; n’introduire `unsafe` que si nécessaire, avec invariants et tests explicites.
- Ne pas promettre une compatibilité binaire ou une traduction classe pour classe de l’API C++.
- Évaluer les dépendances sur leurs fonctionnalités réelles, licences, maintenance et contraintes natives avant de les adopter. Les noms de crates ci-dessous sont des **candidats**, pas des choix validés.

## 2. Premier objectif livré : Poisson 2D

Résoudre sur le carré unité le problème scalaire stationnaire

```text
−Δu = f dans Ω = (0, 1)²,  u = 0 sur ∂Ω.
```

Périmètre précis du premier parcours : maillage triangulaire structuré conforme, éléments de Lagrange **P1** continus, mapping affine, quadrature sur triangle, assemblage séquentiel en COO puis CSR, conditions de Dirichlet homogènes, résolution par gradient conjugué (CG), export VTK et exemple exécutable. P1 désigne ici les polynômes affines sur triangle ; pour les quadrilatères, l’élément bilinéaire correspondant sera **Q1**, traité plus tard.

### Critères d’acceptation du MVP

1. `cargo test` et `cargo run --example poisson_2d` réussissent sur un environnement Rust stable documenté, sans dépendance native obligatoire.
2. Le programme construit son maillage, assemble le système, impose les conditions aux limites, résout et écrit un résultat visualisable ; il signale les erreurs au lieu de retourner silencieusement un résultat invalide.
3. Avec la solution manufacturée `u = sin(πx) sin(πy)` et `f = 2π²u`, les erreurs L2 et semi-norme H1 diminuent sous raffinement uniforme. Un test de convergence vérifie des pentes proches de 2 en L2 et de 1 en H1 sur une série de maillages suffisamment fins, avec des tolérances documentées et une quadrature assez précise pour l’intégration de `f` et des erreurs.
4. Les tests couvrent la géométrie, les fonctions de forme, la quadrature, les DoFs partagés, l’assemblage, les conditions aux limites et le solveur, y compris des entrées invalides.
5. La formulation, les hypothèses, les commandes, les critères d’arrêt du solveur et les limites de l’exemple sont documentés.

## 3. Architecture cible, à faire évoluer

Commencer par **une seule crate bibliothèque** organisée en modules ; extraire des crates du workspace seulement lorsque leurs frontières, besoins de versionnement et dépendances sont stabilisés. La relation entre modules n’est pas une simple chaîne : par exemple l’assemblage utilise à la fois le maillage, les DoFs, les éléments, la quadrature et l’algèbre.

| Module | Responsabilités initiales | Extensions ultérieures |
| --- | --- | --- |
| `base` | Points/vecteurs géométriques, identifiants typés, erreurs et conventions numériques | Dimensions génériques et utilitaires supplémentaires |
| `grid` | Sommets, triangles, connectivité, frontières, génération du carré unité, validation | Voisins, faces, import, 3D, hiérarchie de raffinement |
| `quadrature` | Règles sur triangle, carré et tétraèdre de référence, degrés d’exactitude documentés | Autres cellules et règles |
| `fe` / `quad` / `tetra` | P1 sur triangles et tétraèdres, Q1 sur quadrilatères, valeurs, gradients et mappings | Ordre supérieur, autres cellules 3D et transformations H(div)/H(curl) |
| `dofs` | Numérotation globale P1 conforme par sommet, indices locaux → globaux | Éléments vectoriels, orientations, contraintes et renumérotation |
| `lac` | Triplets COO, compression CSR, produit matrice-vecteur, CG et diagnostic | Backends externes, préconditionneurs et autres solveurs |
| `integrators` / `elliptic` | Intégration locale et assemblage P1 de diffusion-réaction symétrique | Advection, systèmes vectoriels et autres intégrateurs |
| `numerics` / `field` | Assemblage Poisson, conditions aux limites, erreurs, interpolation et projection P1 | Formulations générales et adaptation |
| `io` / `gmsh` | Export VTK et import limité MSH 4.1 ASCII avec groupes physiques | Autres variantes Gmsh, formats et interopérabilité |

Les contrats initiaux doivent distinguer coordonnées de référence et coordonnées physiques, et expliciter l’ordre local des sommets, l’orientation et le type de cellule. Les identifiants `CellId`, `VertexId` et `DofId` encapsulent des indices ; leur validité et leur durée de vie après modification du maillage doivent être documentées. Utiliser des itérateurs empruntés sûrs plutôt que des pointeurs internes exposés.

Ne pas figer dès le départ des traits universels `Mesh`, `Mapping` et `FiniteElement` : commencer par des types concrets, puis extraire des traits réduits lorsque deux implémentations ou un besoin d’extension les justifient. De même, commencer explicitement en 2D ; étudier ensuite `const DIM: usize`, éventuellement avec une dimension d’espace distincte, sans imposer une généralisation qui compliquerait le MVP.

### Contrats numériques importants

- Triangle de référence, orientation et convention du Jacobien fixés et testés. Rejeter les cellules dégénérées et les connectivités invalides ; définir si les triangles d’orientation inversée sont réordonnés ou rejetés, et tester ce choix.
- Transformer les gradients par `J⁻ᵀ` et les mesures par `|det J|` ; contrôler les matrices singulières ou mal conditionnées selon une politique documentée.
- Additionner les doublons COO de manière définie lors de la compression CSR ; vérifier bornes, dimensions et cohérence de la matrice.
- Pour préserver la symétrie nécessaire à CG, éliminer les Dirichlet en modifiant **à la fois lignes, colonnes et second membre** ; tester aussi des valeurs non homogènes avant de les annoncer comme prises en charge.
- CG exige une matrice symétrique définie positive : documenter cette précondition, fixer tolérance relative/absolue et nombre maximal d’itérations, et renvoyer convergence, résidu, nombre d’itérations ou une erreur explicite en cas de rupture/non-convergence.
- Pour les contraintes affines futures, exprimer `x = C y + g` et assembler ou condenser le système réduit `CᵀAC y = Cᵀ(b − Ag)` ; prévoir détection des cycles et gestion des dépendances entre contraintes. Les hanging nodes et le raffinement `hp` sont hors MVP.

## 4. Jalons et ordre de réalisation

Chaque jalon doit inclure code, documentation, tests et exemple pertinent ; il n’est terminé que lorsque ses critères sont satisfaits. Les estimations calendaires dépendront de l’équipe et des résultats des prototypes.

### J0 — Initialiser le projet

- Créer la crate Cargo, organiser les modules, définir la version minimale de Rust prise en charge (MSRV) et les commandes locales.
- Configurer CI pour `cargo fmt --check`, `cargo clippy` et `cargo test` ; ajouter rustdoc et une politique de contribution.
- Déclarer `MIT OR Apache-2.0` dans le manifeste Cargo, conformément au choix du projet ; vérifier séparément les licences des sources d’inspiration et des dépendances. Ne pas copier du code ou des tables sans vérifier les conditions applicables.
- **Terminé si** une crate minimale compile sur la MSRV choisie et Rust stable, et que les contrôles CI passent.

### J1 — Géométrie et maillage triangulaire 2D

- Définir points, identifiants et cellules ; construire un maillage structuré du carré unité par subdivision en triangles.
- Identifier les arêtes de frontière et les sommets associés. Valider indices, absence de cellules dégénérées et convention d’orientation ; fournir des erreurs utiles.
- Tester le nombre de sommets/cellules, la couverture et l’aire totale, les frontières, ainsi que des maillages invalides.
- **Terminé si** la génération et les itérateurs sont documentés et vérifiés sur plusieurs tailles, y compris les cas limites.

### J2 — Quadrature, P1 et mapping

- Définir le triangle de référence, les trois fonctions de forme P1 et leurs gradients ; implémenter le mapping affine et la transformation des gradients.
- Fournir au moins une règle de quadrature triangulaire de degré connu, puis des règles suffisamment précises pour le membre droit et le calcul d’erreurs de l’exemple.
- Tester partition de l’unité, propriété nodale, somme nulle des gradients, intégration de monômes jusqu’au degré annoncé, aire et transformation des gradients.
- **Terminé si** les résultats de référence et les cas de géométrie invalide passent les tests unitaires.

### J3 — DoFs, algèbre et assemblage

- Numéroter les DoFs P1 par sommet, exposer les indices locaux et les DoFs de frontière ; vérifier le partage entre cellules voisines.
- Implémenter ou intégrer un stockage COO → CSR et le produit matrice-vecteur ; tester doublons, lignes vides, bornes et dimensions.
- Assembler localement puis globalement la rigidité et le second membre de Poisson ; imposer Dirichlet en préservant la symétrie.
- Tester une matrice locale connue, les propriétés globales attendues et l’élimination avec valeur de frontière non nulle.
- **Terminé si** les matrices et vecteurs d’un petit maillage sont comparés à des valeurs de référence et si les invariants sont vérifiés.

### J4 — Solveur et exemple complet

- Ajouter CG et des diagnostics explicites ; tester une petite matrice SPD connue, la non-convergence, les dimensions incompatibles et les ruptures numériques.
- Ajouter `examples/poisson_2d.rs`, calcul des erreurs, export VTK et documentation des commandes.
- **Terminé si** tous les critères du MVP de la section 2 sont satisfaits en CI.

### J5 — Robustesse et interopérabilité

**En cours :** les exemples `poisson_mixed` et `poisson_neumann` sont implémentés. Le second contrôle la compatibilité discrète par composante et impose une moyenne intégrale nulle ; l’interopérabilité reste à faire. Les validations CI sur Rust stable et MSRV restent également à confirmer.

- Ajouter Dirichlet non homogène à l’exemple et, séparément, des conditions de Neumann ; documenter les conditions de compatibilité et le noyau du cas Neumann pur.
- Import Gmsh 4.1 ASCII limité implémenté sans dépendance ; étudier les variantes et crates existantes (dont `mshio`) avant extension. Tester les allers-retours là où ils sont possibles.
- Comparer sur des cas identiques les matrices, solutions et erreurs avec un code de référence deal.II ; consigner versions, maillages, tolérances et commandes.
- **Terminé si** des exemples reproductibles et des tests de non-régression couvrent les cas annoncés.

### J6 — Maillages et éléments supplémentaires

**En cours :** les parcours Poisson-Dirichlet Q1 2D et P1 tétraédrique 3D expérimentaux sont implémentés avec exemples, tests géométriques et de convergence. Les conditions mixtes et Neumann pur à moyenne intégrale nulle sont traités sur quadrilatères et tétraèdres. La validation sur Rust stable et MSRV en CI et les extensions (import et autres formulations) restent à faire.

- Ajouter Q1 sur quadrilatères, puis tétraèdres P1 en 3D ; traiter explicitement l’ordre local, les orientations, les faces et les mappings non affines lorsque nécessaires.
- Généraliser les interfaces seulement après identification des points communs et des différences ; conserver les tests de Poisson 2D.
- **Terminé si** chaque type de cellule possède ses tests géométriques, de quadrature et de convergence, et un exemple de bout en bout.

### J7 — Adaptativité et éléments d’ordre supérieur

**Amorce expérimentale :** raffinement global conforme en quatre triangles, relations parent/enfants à un niveau, transfert P1 et indicateur résiduel pour Poisson homogène. Ce n’est pas encore une adaptation locale ni une hiérarchie persistante.

- Définir la hiérarchie de cellules, le marquage, le raffinement et le transfert de solution ; introduire les contraintes affines et les hanging nodes avec tests de continuité.
- Ajouter estimation d’erreur, raffinement local puis éléments d’ordre supérieur ; traiter ensuite l’adaptation `hp`, les orientations et les contraintes entre ordres différents.
- **Terminé si** des tests de convergence et de continuité sur maillages non conformes passent, et si un exemple adaptatif est reproductible. `hp` exige des critères supplémentaires propres à chaque famille d’éléments.

### J8 — Performance et distribution

- Mesurer d’abord les coûts : construction du maillage, assemblage, mémoire et résolution, sur des cas et machines décrits.
- Prototyper l’assemblage parallèle par contributions locales puis réduction maîtrisée ; tester absence de courses et accord numérique avec la version séquentielle. La réduction flottante n’est pas nécessairement bit à bit reproductible : documenter la politique choisie.
- Évaluer préconditionneurs, solveurs externes et MPI (partitionnement, données fantômes, communication et I/O distribuées), sans rendre ces dépendances obligatoires pour le parcours de base.
- **Terminé si** des benchmarks reproductibles démontrent les gains, et que les tests de correction séquentiels/parallèles et les limites des backends sont documentés.

## 5. Dépendances : décisions à prendre

| Besoin | Candidats / décision | Vérifications avant adoption |
| --- | --- | --- |
| Algèbre creuse | Crate Rust telle que `faer`, ou implémentation minimale pour le MVP | Formats, assemblage, solveurs effectivement disponibles, performance, MSRV et licence |
| Solveurs directs | Par exemple `russell_sparse` et ses backends natifs, en option | Disponibilité réelle des backends, installation, plateformes, licence et coûts de build |
| Parallélisme | `rayon`, après profilage | Stratégie de réduction, reproductibilité et surcoût sur petits maillages |
| I/O | `mshio`, `vtkio`, ou export limité maison | Versions et variantes des formats, maintenance, erreurs et tests d’interopérabilité |
| Références architecturales | `fenris`, deal.II et autres projets FEM Rust | API et état actuels, compatibilité des licences ; inspiration plutôt que dépendance automatique |

Ne pas ajouter toutes ces dépendances à J0 : réaliser un prototype ciblé, noter la décision et son coût, puis encapsuler les backends derrière une interface seulement lorsque cela simplifie réellement l’usage. Les solveurs natifs et MPI doivent être optionnels et documenter leurs prérequis système.

## 6. Stratégie de validation continue

- **Unitaires** : invariants de chaque couche, cas limites, erreurs de dimension, géométries dégénérées et matrices connues.
- **Intégration** : petits maillages et problèmes manufacturés, avec convergence en norme L2 et semi-norme H1 ; quadrature et tolérances adaptées à la quantité mesurée.
- **Référence croisée** : cas équivalents exécutés avec deal.II, sans exiger son installation pour les tests courants ; enregistrer la provenance des valeurs de référence.
- **Régression** : `fmt`, `clippy`, tests, exemples et documentation en CI ; élargir progressivement aux plateformes et à la MSRV annoncées.
- **Performance** : benchmarks séparés des tests de correction, avec tailles de problèmes, matériel, compilateur et critères de comparaison publiés.

## 7. Risques et décisions ouvertes

1. **Portée** : la parité deal.II est pluriannuelle et dépend des moyens. Prioriser les besoins applicatifs après le MVP plutôt que promettre toutes les fonctionnalités.
2. **Topologie et adaptativité** : orientations, voisinage, identifiants et contraintes deviennent difficiles en 3D et sur maillages non conformes. Prototyper sur des cas minimaux avant de figer l’API.
3. **Qualité numérique** : mauvais Jacobien, quadrature insuffisante ou conditions aux limites incorrectes peuvent produire une solution plausible mais fausse. Multiplier tests analytiques et comparaisons indépendantes.
4. **Dépendances natives** : elles compliquent installation et portabilité. Garder un chemin de base sans dépendance native obligatoire.
5. **Abstractions prématurées** : la généricité dimensionnelle et les traits peuvent contraindre les évolutions. Les introduire sur la base d’au moins deux cas concrets.
6. **Licence et gouvernance** : appliquer la licence choisie `MIT OR Apache-2.0`, définir les règles de contribution et vérifier les licences de deal.II, des exemples et de toutes les dépendances avant toute réutilisation.

## 8. Prochaines tâches concrètes

1. Valider le périmètre du MVP et les cibles de plateforme avec les contributeurs.
2. Réaliser J0 : crate, licence, MSRV, CI et documentation de développement.
3. Ouvrir des tâches pour J1 à J4 avec les critères d’acceptation ci-dessus et leurs dépendances ; commencer par J1 et J2, puis intégrer J3 et J4.
4. À la fin de J4, mesurer la qualité du parcours complet et réviser l’ordre des jalons suivants selon les besoins utilisateurs.
