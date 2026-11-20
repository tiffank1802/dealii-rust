# Plans détaillés des modules et sous-modules

[Plan global](../PLAN_MODULES_DEALII.md) · [Inventaire amont](../INVENTAIRE_ENTETES_DEALII.md) · [Dépendances, relations et tests transversaux](DEPENDANCES_ET_TESTS.md)

Chaque fiche est indépendante et regroupe les sous-modules fonctionnels du répertoire amont, les tâches priorisées, contrats, relations, crates candidates et critères de test. Les groupes sont conçus pour Rust : ce ne sont pas nécessairement des répertoires existants en C++ ni des crates Cargo séparées. Les fichiers internes/templates ne sont pas tous des fonctionnalités à reproduire. Un premier parcours FEM P0 est implémenté ; voir le [README du projet](../../README.md) pour son périmètre.

| Module | Plan séparé | Sous-chantiers documentés | Première priorité |
| --- | --- | ---: | --- |
| `base` | [Fondations géométriques, numériques et techniques](base.md) | 6 | P0 |
| `grid` | [Topologie, géométrie et cycle de vie des maillages](grid.md) | 5 | P0 |
| `fe` | [Éléments finis, mappings et évaluations](fe.md) | 5 | P0 |
| `dofs` | [Numérotation, accès et construction des contraintes FEM](dofs.md) | 4 | P0 |
| `lac` | [Algèbre linéaire, contraintes et solveurs](lac.md) | 6 | P0 |
| `numerics` | [Assemblage, outils de solution et post-traitement](numerics.md) | 5 | P0 |
| `vtk` | [Échange et visualisation VTK](vtk.md) | 3 | P0 |
| `gmsh` | [Import et intégration Gmsh](gmsh.md) | 3 | P1 |
| `integrators` | [Formes bilinéaires et linéaires locales](integrators.md) | 4 | P1 |
| `meshworker` | [Parcours de maillage et assemblage](meshworker.md) | 3 | P1 |
| `hp` | [Collections et adaptation locale h/p](hp.md) | 3 | P2 |
| `multigrid` | [Transferts, lisseurs et cycles multiniveaux](multigrid.md) | 3 | P2 |
| `algorithms` | [Orchestration non linéaire et temporelle](algorithms.md) | 3 | P2 |
| `physics` | [Outils tensoriels et modèles physiques](physics.md) | 3 | P1 |
| `differentiation` | [Différentiation automatique et symbolique](differentiation.md) | 3 | P2 |
| `optimization` | [Minimisation et recherche de pas](optimization.md) | 3 | P2 |
| `distributed` | [Maillages, transferts et données distribués](distributed.md) | 4 | P3 |
| `matrix_free` | [Application d’opérateurs sans matrice assemblée](matrix_free.md) | 4 | P3 |
| `particles` | [Particules, propriétés et déplacements](particles.md) | 3 | P4 |
| `non_matching` | [Géométrie immergée et couplage de maillages](non_matching.md) | 3 | P4 |
| `arborx` | [Recherche géométrique spatiale](arborx.md) | 2 | P4 |
| `boost_adaptors` | [Pont entre représentations géométriques](boost_adaptors.md) | 2 | P4 |
| `cgal` | [Algorithmes de géométrie et maillages de surface](cgal.md) | 2 | P4 |
| `opencascade` | [Géométrie CAD et manifolds](opencascade.md) | 3 | P4 |
| `sundials` | [Intégration des solveurs SUNDIALS](sundials.md) | 3 | P4 |
| `trilinos` | [Adaptateurs NOX et ROL de Trilinos](trilinos.md) | 3 | P4 |

## Légende et ordre

- **P0** : parcours Poisson 2D ; **P1** : extensions séquentielles ; **P2** : adaptativité, DG et algorithmes ; **P3** : parallélisme et distribution ; **P4** : fonctions spécialisées. **X** dans le plan global désigne une intégration externe optionnelle.
- L’ordre d’une fiche ne garantit pas que ses prérequis d’autres fiches soient terminés : consulter [les frontières de livraison](DEPENDANCES_ET_TESTS.md#5-ordre-de-réalisation-et-frontières-de-livraison).
- Une ligne de sous-chantier doit devenir un ou plusieurs tickets en précisant le périmètre exact, la configuration, les cas de tests, les erreurs et la documentation. La parité avec deal.II n’est pas promise à la livraison du MVP.
- La double licence du code original **MIT OR Apache-2.0** est approuvée. Les versions et backends indiqués restent des propositions ; aucune dépendance n’est installée.
