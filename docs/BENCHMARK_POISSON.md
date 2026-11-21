# Mesure exploratoire : Poisson triangulaire séquentiel

`examples/benchmark_poisson.rs` mesure séparément la construction du maillage, l’assemblage Poisson COO/CSR, l’application de Dirichlet et CG, avec l’horloge murale monotone. Chaque répétition reconstruit le problème ; il n’y a pas de phase de chauffe séparée. Les étapes sont chronométrées indépendamment et la médiane de chaque série est affichée (médiane supérieure pour un nombre pair de répétitions). Compilation, démarrage du processus, calcul du résidu réel après résolution et I/O sont exclus des temps affichés. La sortie inclut également le nombre d’entrées CSR stockées et la taille logique des tableaux CSR après Dirichlet ; ce n’est ni un profil de mémoire réelle ou maximale, ni un benchmark statistique robuste ou une comparaison avec un autre solveur.

Cas : carré unité subdivisé en `n×n` carrés et deux triangles P1 par carré ; source `2π² sin(πx) sin(πy)`, Dirichlet homogène, CG non préconditionné à estimation initiale nulle, tolérances par défaut `rel=1e-10`, `abs=1e-12`, maximum 10 000 itérations. L’exemple vérifie la finitude du résidu réel et affiche itérations et norme `||Ax-b||₂`.

```sh
cargo run --release --example benchmark_poisson -- 16 5
cargo run --release --example benchmark_poisson -- 32 5
cargo run --release --example benchmark_poisson -- 64 5
```

Exécution exploratoire locale : Linux x86_64, Intel Xeon @ 2,60 GHz, environnement virtualisé partagé, deux processeurs visibles, code séquentiel ; Rust AOSP `1.83.0-dev (90b35a623 2024-11-26)`, compilation `--release`, cinq répétitions indépendantes par commande.

| n | Sommets | Triangles | Itérations | Maillage ms | Assemblage ms | Dirichlet ms | CG ms | Résidu réel | Entrées CSR | Octets logiques CSR |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 16 | 289 | 512 | 24 | 0,170 | 0,367 | 0,022 | 0,104 | 3,200e-11 | 1521 | 26656 |
| 32 | 1089 | 2048 | 50 | 0,869 | 1,497 | 0,087 | 0,997 | 2,858e-11 | 6609 | 114464 |
| 64 | 4225 | 8192 | 100 | 3,101 | 7,598 | 0,332 | 5,958 | 1,320e-11 | 27537 | 474400 |

Le stockage logique comptabilise uniquement les longueurs des tableaux de pointeurs de lignes, indices de colonnes et valeurs ; les zéros stockés comptent. Il exclut capacités inutilisées, en-têtes, allocateur, maillage, vecteurs, assemblage COO et pic mémoire. Les temps varient avec le matériel, la charge et le compilateur ; ils ne démontrent aucun gain de performance. Pour une étude plus solide : répéter sur une machine dédiée, avec échauffement, intervalles de dispersion et consommation mémoire ; comparer CG/Jacobi et des backends à tolérance et problème identiques, puis profiler avant optimisation.
