# opencascade — Géométrie CAD et manifolds

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/opencascade)

**Statut :** conception uniquement, aucune API ni fonctionnalité implémentée. Les noms Rust proposés sont indicatifs, pas des signatures stabilisées.

## Relations et contrats

- **Prérequis :** grid/manifolds ; cgal et gmsh facultatifs.
- **Consommateurs :** Applications avec frontières géométriques courbes.
- **Interface et invariants :** Objet CAD opaque propriétaire, projection/normale avec tolérance et résultat explicite ; isoler les durées de vie FFI.
- **Dépendances externes :** Open CASCADE candidat natif facultatif ; aucun binding Rust sélectionné. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P4 | `utilities.h, boundary_lib.h` | Ouvrir géométrie CAD et projeter un point sur courbe/surface. | Plan et cylindre analytiques, projection hors tolérance et fichier invalide. |
| P4 | `manifold_lib.h` | Utiliser la surface CAD pour placement des sommets raffinés. | Sommet projeté sur géométrie et erreur géométrique sous raffinement. |
| P4 | `manifold_lib.h` | Contrôler ownership, fermeture et propagation d’erreurs à travers FFI. | Libération après erreur, utilisation après fermeture impossible via API sûre. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.
