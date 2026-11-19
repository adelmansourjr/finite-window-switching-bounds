# A finite-window forecast for a nonsmooth switching system

Using only mode fields, initial conditions, step size and a known switching-width envelope, the interval for simultaneous mode occupancy covered **all 1,032 trajectories**. In the prospective 576 coupled runs its median width was **0.396**, a tightening of **0.604** against the unconditional interval. The largest matched difference between regularisations in coupled settings was **0.366**. Both marginal occupancy forecasts also covered all trajectories.

| Cohort | Covered by finite joint interval | Median width | Median absolute error, product | Median absolute error, max entropy |
|---|---:|---:|---:|---:|
| Original vector fields | 72/72 | 0.203 | 0.0426 | 0.0426 |
| New coupled fields, primary | 576/576 | 0.396 | 0.0460 | 0.0503 |
| New uncoupled fields | 384/384 | 0.336 | 0.0015 | 0.0015 |

The uncoupled models have very accurate median point forecasts. Their largest regularisation difference was still 0.334. The interval is a deterministic bound under known mode fields and a switching-width envelope. [Frozen protocol](SWITCHING_PROTOCOL.md) · [Results](switching_results.json) · [Catsis, Hall and Jeffrey author manuscript](https://nonsmooth.land/wp-content/uploads/2025/05/2024_paradigm.pdf).

## Finite-window bound

Take the strictly positive mode speeds in Fig. 2. At step t the update is x[t+1] = x[t] + h V p[t], where p[t] selects a mode or interpolates the four modes in the smooth variant. Summing over a window of length T gives V μ = (x[end] − x[start])/T, where μ is the mean mode-weight vector.

For a coordinate with minimum speed m and maximum speed M, the direct sign update pushes inward outside zero and the hysteretic rule pushes inward outside its width w. For the tanh rule, a positive coordinate beyond w log(M/m)/2 gives positive-side switch probability at least M/(M+m). Its outward drift is bounded above by (1−s)M, and its opposing drift is at least s m. So it pushes inward regardless of the other coordinate. By reflection and a maximum step hM, induction bounds the absolute coordinate by the larger of its initial magnitude and max(w, w log(M/m)/2) + hM.

Each terminal coordinate bound and the displacement identity give a linear slab for μ. Constraining occupancies to be nonnegative and sum to one yields a polytope. Exhaustive active-constraint enumeration finds the exact extrema for the *relaxed information set*. Some vectors allowed by those constraints may remain dynamically unreachable. All trials check the integrated balance and the terminal bounds, and the affine controls use analytically derived marginals.

This finite discrete-time extension is motivated by Sections 3–4 and Appendix A. It compares three Euler variants using the same mode fields. A focused follow-up for Jeffrey is to test whether observed switch timing gives tighter predictive intervals on an independently nominated kernel. [Original method and coefficients](https://nonsmooth.land/wp-content/uploads/2025/05/2024_paradigm.pdf).

## Replay

```sh
python3 switching.py --verify switching_results.json
```

Python 3.10+ standard library. The replay includes each field family, the frozen horizon and all step widths and the numerical balance checks.
