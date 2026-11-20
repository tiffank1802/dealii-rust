# gmsh — Import et intégration Gmsh

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/gmsh)

**Statut :** import expérimental d’un sous-ensemble MSH 4.1 ASCII dans `src/gmsh.rs`, sans dépendance externe ; voir le [README](../../README.md#import-gmsh-et-résolution-sur-maillage-étiqueté) pour le contrat et les limites. L’API native Gmsh n’est pas intégrée.

## Relations et contrats

- **Prérequis :** grid : connectivité, types cellulaires et marqueurs.
- **Consommateurs :** Exemples et applications important un maillage.
- **Interface et invariants :** Importer `.msh` en Mesh validé ; erreurs avec version, type cellulaire et marqueur concernés ; distinguer fichier et API native.
- **Dépendances externes :** `mshio` candidat ; vérifier formats réellement pris en charge. Gmsh natif facultatif et non retenu par défaut. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P1 | `grid/grid_in.h` | Définir et implémenter un sous-ensemble Gmsh 4.1 ASCII : triangles linéaires, sommets, groupes physiques et frontières ; vérifier sa faisabilité avec mshio. | Fixture minimaliste indépendante, tags, indices non contigus, version et élément non supportés. |
| P1 | `grid/grid_in.h` | Convertir connectivité et orientation selon invariants grid ; conserver identifiants de matériaux/bord. | Aire, orientation, voisinage et marqueurs après import ; données invalides rejetées. |
| P4 | `gmsh/utilities.h` | Évaluer l’intégration native Gmsh pour génération CAD, séparée du parseur. | Maillage créé puis importé, erreur de génération et gestion des ressources. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.

### Validation actuelle

La fixture originale `tests/fixtures/square41.msh` contient des tags non contigus, plusieurs groupes par surface et des marqueurs de bord. Elle a été relue indépendamment avec meshio 5.3.5 ; `tri.msh` du clone local deal.II a également été importé, sans le redistribuer. Les tests couvrent l’import, les erreurs de format et de connectivité, et un problème mixte de solution affine. Les limitations de conformité géométrique de `Mesh2` s’appliquent aussi ici. `mshio` reste candidat non évalué ni adopté.
