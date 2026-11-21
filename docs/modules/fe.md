# fe — Éléments finis, mappings et évaluations

[Index des fiches](README.md) · [Plan global](../PLAN_MODULES_DEALII.md) · [En-têtes amont](https://github.com/dealii/dealii/tree/1d230ad0be5ab2a3913cf4de8187301430721ed5/include/deal.II/fe)

**Statut :** P1 triangulaire et tétraédrique et premier parcours Q1 quadrilatère expérimentaux ; voir le [README](../../README.md) pour le périmètre et la validation. Les extensions et noms ci-dessous restent des propositions.

## Relations et contrats

- **Prérequis :** base : quadrature ; grid : cellules et sommets.
- **Consommateurs :** dofs, integrators, numerics, hp, matrix_free.
- **Interface et invariants :** Types concrets P1 et mapping affine ; `shape(i, ξ)`, `grad_ref(i, ξ)` puis évaluation physique ; définir l’ordre local et les erreurs.
- **Dépendances externes :** P0 : standard. `nalgebra` éventuellement pour petits systèmes ; aucune hypothèse de support FE fournie par la crate. Voir la [politique des crates](DEPENDANCES_ET_TESTS.md).

## Sous-modules et tâches

Les fichiers amont cités sont des repères fonctionnels représentatifs, pas un engagement de traduction fichier par fichier. Les priorités P0–P4 et X suivent le plan global.

| Priorité | Référence amont | Travail et contrat Rust | Validation / définition de terminé |
| --- | --- | --- | --- |
| P0 | `fe_simplex_p.h, mapping_p1.h` | Fonctions nodales P1 et mapping triangle, avec gradient physique J⁻ᵀ et mesure `abs(det J)` ; rejeter Jacobien invalide. | Valeurs nodales, somme des formes=1 et gradients=0 ; triangle transformé et Jacobien singulier. |
| P0 | `fe_values.h, fe_update_flags.h` | Évaluation aux points de quadrature : valeurs, gradients et poids physiques ; données immuables par cellule. | Intégrales et matrice locale analytique ; nombre et ordre de points cohérents. |
| P1 | `fe_q.h, fe_system.h, mapping_q1.h` | Q1 quadrilatère, P1 tétraèdre et champs vectoriels ; mapping Q1 calculé au point de quadrature. | Nodalité, patch test et pentes de convergence indépendantes pour chaque cellule et dimension. |
| P2 | `fe_dg*.h, fe_interface_values.h` | Éléments DG, évaluations de faces, normales, sauts et sous-faces ; fixer l’orientation. | Sauts de fonctions constantes, conservation des flux et permutations des faces. |
| P4 | `fe_raviart_thomas.h, fe_nedelec.h, fe_bdm.h` | Familles H(div)/H(curl) et autres familles, une à la fois ; transformations de Piola et orientations. | Continuité normale/tangentielle, moments nodaux et convergence par famille. |

## Exécution et validation

- Réaliser les tâches dans l’ordre des priorités et des prérequis indiqués ; chaque ligne constitue au minimum un ticket indépendant avec code, documentation et tests.
- Pour chaque capacité livrée, préciser dimensions, cellules, scalaires, erreurs et backends réellement pris en charge ; ajouter des tests négatifs et conserver le cas de non-régression.
- Appliquer les [contrats transverses, matrices de test et politique de dépendances](DEPENDANCES_ET_TESTS.md). Les intégrations natives et capacités non encore validées ne doivent pas devenir obligatoires pour P0.

## Détail d’implémentation P0

- `reference::TriangleP1` : `φ₀=1−ξ−η`, `φ₁=ξ`, `φ₂=η`. Gradients de référence constants `(-1,-1)`, `(1,0)`, `(0,1)` et ordre aligné sur `grid`.
- `mapping::AffineTriangle` : `x=x₀+Jξ`, colonnes de `J` égales à `x₁−x₀` et `x₂−x₀`. Pour chaque point, calculer `J⁻ᵀ∇φ` et `JxW=|det J|w`; refuser Jacobien non valide selon la politique de `grid`.
- `values` : pré-calculer points physiques, valeurs, gradients et `JxW` sans supposer que les gradients seront constants pour tous les futurs mappings. Vérifier tailles entre règle et résultats ; borner les indices.
- **Tests ciblés :** identité et triangle transformé, partition de l’unité, gradients physiques exacts et intégrale de 1 égale à l’aire ; quadrature plus précise pour le membre droit non polynomial et les erreurs.

**État Q1 partiel :** `quad` fournit formes, gradients, mapping bilinéaire au point et assemblage Poisson avec quadrature 3×3. Patch affine et convergence sur maillages déformés sont testés. Les conditions mixtes, Neumann pur, l’assemblage diffusion-réaction, la prolongation sur raffinement uniforme, l’interpolation, la projection L² et l’évaluation ponctuelle Q1 sont disponibles ; les champs vectoriels restent à faire.

**État tétraèdre P1 :** `tetra` propose un mapping affine 3D avec orientation positive, des gradients constants et Poisson-Dirichlet avec quadrature de degré 5. Les conditions mixtes, Neumann pur, l’assemblage diffusion-réaction à coefficients variables et le transfert P1 sur raffinement uniforme sont disponibles.

**Champs P1 tétraédriques :** interpolation nodale, projection L² et recherche/évaluation en un point physique sont disponibles dans `tetra_field`.
