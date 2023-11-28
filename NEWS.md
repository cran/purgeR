### purgeR 1.8.2

- Small fixes in documentation format

### purgeR 1.8.1

- Fix format string print warning functions.
- Update citation format.

### purgeR 1.8

- Added `ped_graph()` to convert pedigrees to igraph objects.
- Removed the restriction to use limited number of cores.

### purgeR 1.7

- Fix use of `ped_clean()` on tibbles.

### purgeR 1.6

- Enhanced portability with tibbles.
- Minor bug fixes and improvements.

### purgeR 1.5

- Fix use of `ped_clean()`.
- Performance improvements for preprocessing functions.

### purgeR 1.4

- Fix use of `exp_g()` for large numbers of generations.
- Use improvements for `ped_sort()` and internal check functions.

### purgeR 1.3

- Include genedropping options for `ip_Fij()` and `ip_op()`.
- Raw (uncorrected) estimates are always returned with `ip_op()`.
- Fix use of `ped_sort()` with tibbles as input.
- Documentation improvements.

### purgeR 1.2

- Fix computation of total opportunity of purging (*O*).
- Function `ip_op()` only computes the expressed opportunity of purging (*Oe*) by default.
- Add new function `ped_sort()` to sort pedigrees.

### purgeR 1.1

- Fixed issue in `pop_Nancestors` due to bad access to index.
- Documentation improvements.
