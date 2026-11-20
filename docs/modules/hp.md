# hp — Collections et adaptation locale h/p

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/hp)

**Statut :** conception uniquement, aucune API ni fonctionnalité implémentée. Les noms Rust proposés sont indicatifs, pas des signatures stabilisées.

## Relations et contrats

- **Prérequis :** grid raffiné, fe multi-ordre, dofs, contraintes lac.
- **Consommateurs :** numerics et meshworker ; distributed en étape parallèle.
- **Interface et invariants :** Collections d’éléments/quadratures/mappings ; indices actifs et futurs attachés aux cellules ; vérifications d’interfaces et de version du maillage.
- **Dépendances externes :** Aucune crate propre ; réutiliser les backends des modules prérequis. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P2 | `fe_collection.h, q_collection.h, mapping_collection.h` | Collections typées et choix par cellule ; dimensions et indices vérifiés. | Indices invalides, évaluation selon ordre local et quadrature compatible. |
| P2 | `fe_values.h, refinement.h` | Décision h/p, indices futurs puis application de l’adaptation et transfert. | Sur deux cellules d’ordres différents, continuité et reproduction d’un polynôme admissible. |
| P4 | `refinement.h` | Politique combinée h/p, conflits de marquage et adaptation distribuée. | Convergence en ordre, transfert et contraintes sur interfaces h/p, 1/2 rangs. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.

## Décomposition des interfaces hp

1. Stocker les collections sans confondre l’indice d’élément actif et l’indice sélectionné pour la prochaine adaptation. Vérifier que quadrature/mapping choisis sont compatibles avec la cellule.
2. Marquer les cellules via un indicateur externe, décider h ou p puis résoudre explicitement les conflits ; transférer solution et indices avant de reconstruire les DoFs.
3. Construire les contraintes de continuité entre cellules d’ordres différents ou de tailles différentes ; fermer les contraintes via `lac` et contrôler leur sens/orientation.
4. **Validation :** cas p uniforme, saut d’ordre sur deux cellules, hanging node, combinaison h/p et conservation d’un polynôme représentable. Chaque famille de FE et cellule demande son propre jeu de tests.
