# Bayesian Networks — K2 Algorithm with **Mutual-Information Ordering**

Learning the topology of a Bayesian Network from data is **NP-hard** when the parent–child ordering of the nodes is unknown.
This project shows how a simple *information-theoretic* trick—ordering the variables by pair-wise **Mutual Information (MI)** before running the greedy **K2** search—shrinks the worst-case search space from `O(n! · uⁿ)` to `O(n · uⁿ)` and delivers an **≈ 1000× speed-up** compared with an exhaustive permutation search on identical hardware.

---

## Contents

| Folder / file                             | Purpose                                                |
| ----------------------------------------- | ------------------------------------------------------ |
| `part1_structure_learning_k2_algo.ipynb`  | Pure K2 implementation (exhaustive order search)       |
| `part2_structure_learning_bnstruct.ipynb` | K2 re-implemented with MI ordering (uses **bnstruct**) |
| `part3_comparison.ipynb`                  | Timing and accuracy benchmarks                         |
| `bnstruct_objects.R`                      | Re-usable R helpers for MI ordering                    |
| `man_functions.R`, `man_objects.R`        | Utility helpers (graph checks, scoring, plotting)      |
| `asia/`, `child/`, `ruiz/`, `sachs/`      | Canonical discrete data sets                           |
| `Bayesian_Networks_and_K2_algorithm.pdf`  | Short project report                                   |

---

## Quick start

```bash
git clone https://github.com/sleepy-byte/bayesian-networks-K2-algorithm
cd bayesian-networks-K2-algorithm

# Install R dependencies
R -e "install.packages(c('bnstruct','bnlearn','entropy','igraph'))"

# Run benchmarks
jupyter notebook part3_comparison.ipynb
```

> **Note** All notebooks run with the vanilla `IRkernel`; no separate Conda environment is required.

---

## Methodology

1. **Compute pair-wise MI**
   For every pair of variables $X_i, X_j$ estimate

   $$
     \mathrm{MI}(X_i,X_j)=\sum_{x_i,x_j}p(x_i,x_j)\log\frac{p(x_i,x_j)}{p(x_i)p(x_j)}.
   $$

2. **Derive a global order**

   * Start with the variable having the **highest total MI** to all others.
   * Iteratively append the variable that maximises the *average* MI to the already ordered set (greedy maximum-spanning-tree heuristic).
     The resulting permutation roughly aligns parents before children.

3. **Run K2 once** with this order instead of enumerating all `n!` permutations.

Because K2’s inner loop is unchanged, **model quality (F₁ score, structural Hamming distance) stays comparable or improves slightly**, while the combinatorial explosion is eliminated.

---

## Performance

| Data set | Variables | Exhaustive search (wall-time) | MI-ordered K2 (wall-time) | Speed-up |
| -------- | --------- | ----------------------------- | ------------------------- | -------- |
| Ruiz     | 3         | 0.24 s                        | 0.02 s                    | 12 ×     |
| Asia     | 8         | 28.3 s                        | 0.19 s                    | 149 ×    |
| Child    | 20        | 4 h 17 m                      | 15 s                      | 1030 ×   |
| Sachs    | 11        | 126 s                         | 0.12 s                    | 1050 ×   |

Benchmarks were run on an M1 MacBook Air (8 GB RAM) with R 4.3.2. Raw outputs and code are in `part3_comparison.ipynb`.

---

## Re-using the MI ordering in other projects

```r
source("bnstruct_objects.R")      # provides mi_node_order()

order <- mi_node_order(my_dataframe)        # character vector
net   <- learn.network(
           my_dataset,
           algo         = "mmhc",           # or "k2"
           layering     = order,            # <- plug it in
           scoring.func = "BDeu",
           max.parents  = ncol(my_dataframe) - 1
         )
```

---

## Dependencies

* **R ≥ 4.2** — `bnstruct`, `bnlearn`, `entropy`, `igraph`
* (Optional) **Jupyter Notebook** with `IRkernel`
* Tested on macOS 14 and Ubuntu 24.04

---

## Background reading

* G. F. Cooper & E. Herskovits — *A Bayesian Method for the Induction of Probabilistic Networks from Data* (1992)
* X.-W. Chen, G. Anantha & X. Lin — *Improving Bayesian Network Structure Learning with Mutual-Information-Based Node Ordering in the K2 Algorithm* (IEEE TKDE 2008)
* M. Scutari & J. Denis — *Bayesian Networks* (CRC Press, 2022)

---

## License

The upstream repository currently has no explicit license. Until one is added, all code and notebooks are provided **“All rights reserved”** for academic and non-commercial use. Open an issue if you need a different license.

---

## Contributing

Bug reports, pull requests, and benchmarking results on larger data sets are welcome. Please open an issue before making substantial changes.

---

© 2025 sleepy-byte — Advanced Statistics for Physics Analysis project
