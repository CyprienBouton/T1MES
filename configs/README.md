# Experiment-specific parameter files

Variants of [`t1mes/T1MES_values.json`](../t1mes/T1MES_values.json) used for specific experiments.
They have the same format (see the main README, "JSON configuration") and differ only in some `T1`/`T2`
values [ms]. Pass them with `--param_json`. The `.json` files are git-ignored (`.gitignore`) until it is
decided whether they can be shared.

| File | Use with | Content (differences from the reference file) |
|------|----------|-----------------------------------------------|
| `T1MES_values_real_acquisition.json` | `t1mes_build_static` | Values measured on an acquisition of the physical phantom: T1/T2 of the post-GBCA tubes G, H and T2 of A, I changed. |
| `T1MES_values_real_acquisition_0303.json` | `t1mes_build_static` | Same, for a later acquisition (03/03): also changes T1 of A, C, I and of the 3 T native tubes J–M. |
| `T1MES_values_p18_16min.json` | `t1mes_build_static` | "p18" values 16 min after injection: T1 of the post-GBCA blood (A), short (G) and long (H) myocardium tubes set to the values at that time; all other T1 values are the earlier MOLLI-based values of the reference file (git history, commit 820046c); T2 of A and C changed. |
| `T1MES_values_p18_27min.json` | `t1mes_build_static` | Same, 27 min after injection. |
| `T1MES_values_p18_36min.json` | `t1mes_build_static` | Same, 36 min after injection. |
| `T1MES_values_p18_dynamic.json` | `t1mes_build_dynamic --times_post` (7 values) | T1 of A, G, H as 7-point time series: a pre-contrast value, then the 16, 27 and 36 min values of the files above, then three later time points. |
| `T1MES_values_simple_dynamic<N>.json` (N = 27, 36, 45, 54, 63, 72) | `t1mes_build_dynamic --times_post` (10 values) | Synthetic test: T1 of A, G and H increases linearly from its reference value by N/9 ms per frame, i.e. N ms over 10 frames. |
| `T1MES_values_simple_dynamic.json` | `t1mes_build_dynamic --times_post` (10 values) | Identical to `T1MES_values_simple_dynamic72.json`. |

To be completed by the author: what "p18" refers to, the time points of the `p18_dynamic` frames
(4th to 7th) and the acquisition protocol of the `real_acquisition` files.
