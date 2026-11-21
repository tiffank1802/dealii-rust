# Dépendances de deal.II et stratégie pour dealii-rust

## Périmètre et sources

Cette analyse porte sur l’instantané deal.II **`1d230ad0be5ab2a3913cf4de8187301430721ed5`** (`9.9.0-pre`), pas sur toutes ses versions. Sources principales : [configuration CMake](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/cmake/configure), [dépendances intégrées](https://github.com/dealii/dealii/blob/1d230ad0be5ab2a3913cf4de8187301430721ed5/bundled/README.md), [configuration Trilinos](https://github.com/dealii/dealii/blob/1d230ad0be5ab2a3913cf4de8187301430721ed5/cmake/configure/configure_20_trilinos.cmake), [en-têtes `lac`](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/lac) et [licence deal.II](https://github.com/dealii/dealii/blob/1d230ad0be5ab2a3913cf4de8187301430721ed5/LICENSE.md).

On distingue **la configuration de la bibliothèque C++**, les bibliothèques que l’on peut choisir pour **l’implémentation Rust**, et celles éventuellement utilisées pour **comparer les résultats**. Une dépendance obligatoire de deal.II n’est pas nécessairement obligatoire pour sa réimplémentation Rust.

## Ce que deal.II utilise

Pour ce commit, CMake impose Boost et Kokkos et recherche un support de threads. Boost et Kokkos peuvent être fournis par le répertoire `bundled/` lorsque la configuration l’autorise. Le code C++ et son système de construction requièrent aussi leurs outils de compilation ; la présence des sources intégrées ne signifie pas qu’elles sont toutes activées dans chaque installation.

Les éléments ci-dessous sont des **intégrations configurables** : disponibilité et options dépendent de la configuration choisie. Le tableau indique leurs rôles principaux et notre traitement proposé, pas une équivalence fonctionnelle déjà démontrée.

| Dépendance amont | Rôle dans l’écosystème deal.II | Proposition pour Rust |
| --- | --- | --- |
| **Trilinos** | Suite de matrices/vecteurs, solveurs, préconditionneurs et autres composants scientifiques. | Adaptateur natif **optionnel** pour usages HPC ; prototyper Tpetra/Belos/Ifpack2/MueLu et exposer les seuls composants réellement validés. Voir ci-dessous. |
| **PETSc** | Algèbre et solveurs, notamment distribués. | Autre backend optionnel ; même contrat de haut niveau que possible, sans supposer la parité des fonctions ou coûts. |
| **MPI** | Communications entre processus. | Backends distribués en P3 ; évaluer la crate `mpi`, soumise à une installation MPI compatible. |
| **p4est** | Forêts d’arbres pour maillages distribués adaptatifs. | Prototype d’adaptateur optionnel lorsque `distributed`/AMR sont prêts ; ne pas réimplémenter tout p4est avant analyse. |
| **METIS**, **Zoltan** | Partitionnement / répartition ; Zoltan est un composant Trilinos. | Étudier lorsque la qualité de partition et la redistribution deviennent nécessaires ; backend natif facultatif. |
| **BLAS/LAPACK** | Noyaux algébriques et algèbre dense. | Pas de prérequis natif pour P0 ; évaluer un backend optimisé après mesures. |
| **UMFPACK**, **MUMPS** | Solveurs directs creux. | Prévoir des adaptateurs séparés, licences et linkage à auditer ; ne pas promettre qu’un seul backend résout tous les cas. |
| **ARPACK**, **SLEPc** | Problèmes aux valeurs propres ; SLEPc s’appuie sur PETSc. | Plus tard, selon application ; tests spectres et prérequis séparés. |
| **ScaLAPACK**, **PSBLAS**, **AMG4PSBLAS** | Algèbre distribuée et préconditionnement. | Backends optionnels avancés, à sélectionner selon besoin et maintenance. |
| **Ginkgo** | Solveurs et préconditionneurs, avec différents backends d’exécution. | Étudier pour les solveurs performants/accélérateurs ; aucun engagement GPU sans prototype. |
| **Kokkos**, **ArborX** | Portabilité d’exécution et recherche spatiale ; ArborX utilise Kokkos. | Pas obligatoires pour le cœur Rust ; étudier adaptateurs natifs et alternatives Rust lorsque recherche spatiale ou accélérateurs sont nécessaires. |
| **Taskflow**, **TBB** | Exécution de tâches et parallélisme partagé. | En Rust, `rayon` est candidat après mesures ; aucune équivalence complète n’est supposée. |
| **SUNDIALS** | Intégration temporelle ODE/DAE et systèmes non linéaires. | Adaptateur optionnel pour ARKODE/IDA/KINSOL après stabilisation des vecteurs et callbacks. |
| **ADOL-C**, **Sacado**, **SymEngine** | Différentiation automatique / symbolique ; Sacado est un composant Trilinos. | Comparer d’abord des solutions Rust aux bindings natifs, suivant le mode de différentiation requis. |
| **Gmsh**, **VTK** | Génération/interopérabilité des maillages et données. | Export VTK limité en P0 ; évaluer `mshio`/`vtkio` pour les formats, indépendamment d’éventuelles API natives. |
| **CGAL**, **Open CASCADE**, **Assimp** | Algorithmes géométriques, CAD et import de géométries/objets. | Intégrations spécialisées optionnelles uniquement après cas d’usage et audit des composants. |
| **HDF5**, **netCDF** | Stockage et échange de données scientifiques. | Pas requis pour P0 ; choisir formats, parallélisme, versions et bibliothèques avant intégration. |
| **muParser**, **GSL**, **zlib**, **magic_enum** | Expressions, calcul scientifique, compression et utilitaires, suivant les fonctions activées. | Utiliser `std` ou une crate appropriée selon besoin ; ne pas porter mécaniquement les dépendances C++. |

Dans `bundled/`, l’instantané inclut notamment Boost, Kokkos, Taskflow, magic_enum, muParser et une version d’UMFPACK avec AMD/UFCONFIG. Leurs licences et versions sont propres à ces composants ; **le clone de référence ne les ajoute pas comme dépendances du projet Rust**. Le [README bundled amont](https://github.com/dealii/dealii/blob/1d230ad0be5ab2a3913cf4de8187301430721ed5/bundled/README.md) signale notamment des licences Boost, BSD, MIT et LGPL selon les composants.

### Exemples de contraintes de configuration amont

- La configuration deal.II de Trilinos et celle de PETSc déclarent une dépendance à **MPI** ; p4est y dépend aussi de MPI et ses vérifications exigent une construction p4est avec MPI et zlib.
- La configuration **SLEPc → PETSc** vérifie la compatibilité des versions majeure et mineure ; **UMFPACK → LAPACK** et **ScaLAPACK → MPI + LAPACK**. Ces relations sont celles du commit amont, pas des obligations générales pour toutes les utilisations possibles de ces bibliothèques.
- Kokkos doit être cohérent avec ceux utilisés pour construire Trilinos ou PETSc ; plusieurs variantes incompatibles de MPI, de types scalaires, d’indices et de compilateurs peuvent empêcher la compilation ou l’exécution.

## Focus Trilinos

Trilinos est une suite de composants, **pas une crate unique à substituer directement**. Dans l’instantané analysé, deal.II vérifie notamment les ensembles `Epetra`/`Amesos`/`Ifpack`/`AztecOO`/`Teuchos`/`ML` et `Tpetra`/`Amesos2`/`Belos`/`Ifpack2` selon les interfaces activées ; d’autres composants tels que `MueLu`, `NOX`, `ROL` et `Sacado` sont détectés suivant la configuration. La disponibilité de Trilinos ne prouve pas que tous ces paquets sont installés.

**Plan proposé, non implémenté :**

1. Définir dans `lac` un contrat minimal d’opérateur, de vecteur, de solveur et de rapport d’itérations ; garder les opérations/fonctions non communes comme extensions propres au backend.
2. Après le parcours séquentiel et un premier modèle MPI, réaliser un prototype isolé d’adaptateur Tpetra + Belos, puis Ifpack2/MueLu si disponibles. Comparer au système CSR/CG sur un petit cas déterministe.
3. Choisir une stratégie de liaison après audit : binding maintenu si adapté, ou mince couche C/C++ avec façade Rust sûre. Ne pas supposer une ABI C stable pour toutes les interfaces C++ de Trilinos. Encapsuler ownership, exceptions C++, erreurs, buffers et durées de vie.
4. Documenter la combinaison testée : version Trilinos et paquets activés, compilateur/runtime C++, MPI, Kokkos, type scalaire et largeur d’indices. Vérifier l’initialisation/finalisation MPI et l’ordre de destruction sans prendre silencieusement possession du communicateur utilisateur.
5. Tester en CI dédiée : absence de Trilinos avec backend désactivé, erreurs de découverte/configuration, petit système connu, convergence/non-convergence, puis 1/2/plusieurs rangs avec comparaison au séquentiel. Tester préconditionneurs et résolution directe séparément si ces capacités sont activées.

Le nom d’une future fonctionnalité Cargo du type `trilinos` serait **indicatif** : aucun binding, feature ou backend n’existe dans le dépôt aujourd’hui. Si la compatibilité HPC avec une installation Trilinos donnée est prioritaire, il serait judicieux de prototyper l’étape 2 plus tôt, sans rendre le backend obligatoire pour le MVP.

## Politique proposée de dépendances Rust

- **MVP :** `std`, assemblage COO/CSR limité, CG et export VTK limité ; pas d’installation native obligatoire. Les opérations et formats réellement couverts restent explicitement documentés.
- **Après validation :** prototypes comparatifs de `faer`/`sprs`, `rayon`, `mshio` et `vtkio` ; voir le [catalogue détaillé](modules/DEPENDANCES_ET_TESTS.md) et vérifier licence, MSRV et fonctionnalités avant de les retenir.
- **HPC et spécialisé :** MPI, PETSc, Trilinos, p4est et autres backends optionnels avec CI propre, configuration reproductible et tests de compatibilité. Une intégration native apporte aussi ses dépendances transitives : elle n’est ni gratuite à maintenir, ni automatiquement compatible avec tous les environnements.

## Licences et redistribution

Le choix accepté pour le code original de ce projet est **`MIT OR Apache-2.0`**, au choix de l’utilisateur. Les fichiers de licence sont à la racine ; lors de la création du manifeste, renseigner `license = "MIT OR Apache-2.0"`.

Cela **ne change pas** les licences des dépendances ou de deal.II. L’amont est sous `Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later`, avec des exceptions propres à certains fichiers. Copier du code, des exemples ou des tests, distribuer une bibliothèque native liée ou redistribuer un binaire peut créer des obligations différentes selon les composants et le mode de distribution. Conserver une provenance, vérifier chaque licence et les dépendances transitives, et faire valider les configurations de distribution si nécessaire ; isoler un backend derrière une feature ne règle pas à lui seul ces questions.
