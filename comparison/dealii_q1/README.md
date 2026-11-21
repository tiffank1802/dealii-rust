# Comparaison ciblée avec deal.II : Q1 sur un carré

**État : protocole préparé, pas une comparaison deal.II exécutée dans cet environnement.** deal.II et CMake ne sont pas installés ici ; `reference.cc` et sa compilation restent à valider avec une installation réelle. `analytic_smoke.txt` contient des valeurs analytiques, **pas une sortie de deal.II**.

Le cas compare sur `[0,1]²` avec un seul quadrilatère affine Q1 :

- la matrice de rigidité locale pour le Laplacien ;
- le second membre de projection pour `f=x²+y²` ;
- la solution de la projection L² sans contrainte ;
- l’erreur L² de cette projection.

Les deux programmes utilisent une quadrature Gauss 3×3. La sortie est ordonnée selon les sommets `(0,0), (1,0), (1,1), (0,1)`, indépendamment du numérotage local deal.II. Le fichier attendu comporte quatre lignes de cinq nombres (matrice et charge), puis quatre coefficients projetés et l’erreur L², soit 25 nombres finis. Le comparateur vérifie un écart absolu maximal de `1e-10` pour chaque grandeur. Ce cas ne valide ni les maillages déformés, ni les conditions aux limites, ni la convergence globale.

## Exécution sur une machine disposant de deal.II ≥ 9.4 et CMake

Depuis la racine du dépôt :

```sh
cmake -S comparison/dealii_q1 -B comparison/dealii_q1/build -DDEAL_II_DIR=/chemin/vers/deal.II
cmake --build comparison/dealii_q1/build --config Release
mkdir -p .cache
comparison/dealii_q1/build/q1_reference > .cache/dealii_q1_reference.txt
cargo run --release --example compare_q1_reference -- .cache/dealii_q1_reference.txt
```

Enregistrer les versions exactes de deal.II, CMake, compilateur C++, Rust et système, la commande et la sortie pour rendre la comparaison traçable. En cas de différences, vérifier d’abord l’ordre des sommets, la quadrature et les options de compilation ; ne pas modifier les résultats pour masquer un écart.

Pour tester seulement le comparateur sans deal.II :

```sh
cargo run --example compare_q1_reference -- comparison/dealii_q1/analytic_smoke.txt
```
