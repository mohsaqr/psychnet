# A transparent, self-certifying base-R implementation of psychometric network estimation

## Abstract

Psychometric network models are estimated almost exclusively through a small set
of R packages that wrap compiled numerical kernels. Those kernels are written in
Fortran and C, so they are not legible in the language of the analysis, they sit
on large dependency trees, their output depends on a convergence tolerance the
user does not set, and a result cannot be checked except by running a second
package. `psychnet` is a clean-room implementation of the cross-sectional
psychometric network models written entirely in base R, with `stats` as its only
import. It calls none of the reference packages and none of the compiled solvers
they depend on. Every regularized fit returns a correctness certificate: the
stationarity residual of the convex objective it solves, a number that is near
zero at the optimum and is computed without reference to any external solver.
This design delivers four properties the wrapped implementations cannot offer
together, namely transparency, reproducibility, freedom from a compiled
dependency tree, and auditability. It also makes the estimator more accurate in a
specific and reproducible sense: because it solves to the optimum and uses the
textbook information criterion, it does not reproduce the spurious edges a
loosely converged reference solver can select. We develop the graphical lasso as
a worked case, show that `psychnet` returns the textbook quantity to the digit
where the reference does not, and report agreement with the reference packages on
real questionnaire data and recovery against a known graph on synthetic data. The
cost is run time: the pure-R graphical lasso is about two orders of magnitude
slower than the Fortran kernel `qgraph` calls. The implementation is a readable
reference whose hot loops can be optimized later without changing the certified
result, and Section 5 reports the current timings.

## 1. Introduction

A psychometric network analysis is usually carried out with `qgraph`,
`IsingFit`, `mgm`, or `bootnet`, the standard tools in the field. They share a
structure that the rest of this paper builds on. Each is a wrapper that prepares
a correlation matrix or a design matrix, sets a penalty path and an information
criterion, and hands the numerical work to a compiled kernel: `qgraph`'s EBIC
graphical lasso calls the `glasso` Fortran package, and `IsingFit` and `mgm` call
`glmnet`. The wrapper is short and is described in the methods paper for each
estimator. The kernel is where the computation happens, and it is not written in
R.

This arrangement has four consequences for the analyst.

First, the estimator is **not transparent to a reader working in R**. The penalty
path and the criterion are written in R and can be read, but the solve itself is a
call into a Fortran or C kernel. That kernel is legible to a reader fluent in
Fortran or C; it is not legible in the language of the analysis, so an analyst who
wants to know what the estimator did at a given penalty, or why one edge survived
and another did not, cannot follow it in the package source without leaving R.

Second, the result is only **as reproducible as the compiled binary**. The
Fortran and C kernels are deterministic given their inputs, but they are linked
at build time, they carry their own numerical settings, and a reinstall on a new
platform can change a result at the level of the convergence tolerance. The
tolerance is a property of the kernel, not a choice the analyst makes.

Third, the estimator carries a **large dependency tree**. Installing `bootnet`
pulls more than a hundred recursive hard dependencies, most of them compiled
through `Rcpp` and `RcppArmadillo`. The tree is a cost in installation time, in
build-failure probability, and in the surface that has to be trusted and
maintained.

Fourth, a result is **not auditable from itself**. The only way to check a
network is to estimate it again with a second package and compare. No quantity
attached to the fit reports how far it is from the answer it is supposed to
compute.

`psychnet` answers these four points. It implements the same cross-sectional
models in base R, with `stats` as its only import, and it reaches the published
estimators or states precisely where it cannot. It adds one thing the wrapped
implementations do not have: every regularized fit returns the stationarity
residual of its own objective, so the distance from the optimum is printed and
can be rechecked from the fitted object. Section 3 develops this.

Two further properties follow from the design.

The implementation is **more precise than the default reference, and provably
so**. `qgraph` runs `glasso` at a default convergence threshold of `1e-4`, so its
returned precision sits about that far from the optimum of the convex objective.
`psychnet` refits the selected penalty until the stationarity residual is near
`1e-11`. Because the objective is strictly convex, the fit with the smaller
residual is closer to the unique minimizer, and the residual is printed, so the
claim is checkable rather than asserted. Section 4 shows that this difference is
not cosmetic: on a near-empty graph it is the difference between declining a noise
edge and keeping it.

The implementation is **a readable reference that can be optimized**. Pure R is
slower than Fortran, and Section 5 reports the gap. The point of a base-R
reference is that the algorithm is legible and that its output is certified.
Should a hot loop need to be moved to compiled code, the certificate is the
invariant the optimization must preserve, so speed can be recovered later without
returning the computation to a kernel the R reader cannot follow.

The paper is organized as follows. Section 2 describes the clean-room
construction and the reproducible provenance check showing that no reference
package or compiled kernel is called. Section 3 describes the certificates.
Section 4 develops the graphical lasso as a worked case, including the one
situation where the textbook implementation and the reference disagree, and
explains why the textbook implementation gives the better answer there. Section 5
reports the run-time cost. Section 6 reports the validation. Section 7 states the
scope and limitations.

## 2. Clean-room construction

A clean-room implementation here means that the estimators are written from their
published specifications, in base R, and that the package calls neither the
reference packages nor the compiled kernels those packages use. This provenance
claim was checked by a reproducible source scan. The scan covered every
executable line of the package R source and the package `DESCRIPTION`. A positive
hit was defined as a namespace-qualified call, `pkg::fn`, to a reference package
or compiled solver used in the reference stack, including `qgraph`, `IsingFit`,
`mgm`, `huge`, `bootnet`, `glasso`, `glmnet`, and `Matrix`. An occurrence inside
a roxygen or comment line was classified as non-executable documentation rather
than as a call site. The result was that no executable line invokes a reference
package or compiled kernel. The `DESCRIPTION` declares a single import, `stats`,
and one suggested package, `testthat`. The recursive hard-dependency count is
therefore the count for base R, against 71 for `qgraph`, 77 for `IsingFit` and
`mgm`, and more than 100 for `bootnet`. The core estimators install and run with
no C or Fortran toolchain.

The numerical work is carried by two kernels written in the package. The first is
the covariance block-coordinate-descent graphical lasso of Friedman, Hastie and
Tibshirani (2008), used by the Gaussian graphical model and its nonparanormal and
graph-restricted variants. It updates one covariance block at a time, reconstructs
the precision matrix from the fitted covariance, and evaluates the Gaussian
penalized likelihood on the same scale used for model selection. The second is
the penalized iteratively reweighted least squares nodewise generalized linear
model of Friedman, Hastie and Tibshirani (2010), used by the Ising and mixed
models with a logistic or Gaussian link. It solves each nodewise regression over a
penalty path, applies the estimator's EBIC or pruning rule, and symmetrizes
nodewise coefficients into an undirected edge set where the published estimator
does so. Two estimators use neither kernel: the Triangulated Maximally Filtered
Graph is a greedy planar construction, and the relative-importance network is a
Shapley decomposition of nodewise R-squared.

The eleven estimator verbs, with their objective, selection rule, and
certificate, are given in Table 1. The table is also the audit map for the
package: each row states the fitted object, the rule that selects its sparsity
pattern, and the certificate that can be recomputed from the returned object.
They are separate from the framework verbs (`centrality`, `predictability`,
`bootstrap_network`, `centrality_stability`, `nct`), which operate on an
already-fitted network.

**Table 1. Estimators.**

| verb | model | selection | certificate |
|---|---|---|---|
| `cor_network` | marginal correlation | optional significance threshold | closed form |
| `pcor_network` | partial correlation | optional significance threshold | closed form |
| `ebic_glasso` | graphical lasso | EBIC over a penalty path | `glasso_kkt` |
| `huge_network` | nonparanormal graphical lasso | EBIC over a penalty path | `glasso_kkt` |
| `ggm_modselect` | unregularized graph-restricted MLE | EBIC over candidate graphs | `ggm_support_kkt` |
| `tmfg_network` | maximally filtered graph | greedy planar construction | `tmfg_certificate` |
| `logo_network` | chordal Markov random field | filtered-graph support | `ggm_support_kkt` |
| `relimp_network` | relative importance, directed | subset enumeration | `lmg_certificate` |
| `ising_fit` | Ising, L1-penalized | per-node EBIC | `glm_lasso_kkt` |
| `ising_sampler` | Ising, unregularized | Wald pruning, optional | `glm_lasso_kkt` |
| `mgm_fit` | mixed Gaussian and binary | per-node EBIC | `glm_lasso_kkt` |

## 3. Self-certification

A from-scratch implementation cannot be trusted on the grounds that it byte-matches
the reference. An independent solver of a strictly convex objective need not
reproduce the reference Fortran path bit for bit, and doing so would make the
compiled tolerance part of the target. It has to be checked against the
mathematics it claims to compute. The device for this is a stationarity residual,
computed from the fitted object against the objective, with no reference solver in
the computation.

The graphical lasso minimizes

    f(Theta) = -log det Theta + tr(S Theta) + rho * sum_{i != j} |Theta_ij|

over positive-definite Theta, where S is the correlation matrix and rho the
penalty. The function is strictly convex, so the minimizer is unique and is
characterized exactly by its subgradient conditions. Writing W = Theta^{-1}, the
conditions are: W_ii = S_ii on the diagonal; W_ij - S_ij = rho * sign(Theta_ij)
on an edge that is present; and |W_ij - S_ij| <= rho on an edge that is absent.
The exported function `glasso_kkt(theta, S, rho)` returns the largest violation
of these conditions. A return near zero certifies that the supplied precision
matrix is the global optimum, and it certifies this independently of any solver.
Every `ebic_glasso` result stores this number; on the fits in this paper it is
near zero, below `1e-9`, and is exactly zero when the selected graph is empty.

The nodewise generalized linear model has the analogous certificate built from
the penalized-likelihood score conditions, `glm_lasso_kkt`, and the Ising and
mixed models store the worst nodewise value. The unregularized estimators that
have a closed form carry a structural certificate instead. The constrained
Gaussian Markov random field behind `ggm_modselect` and `logo_network` is checked
by `ggm_support_kkt`, which verifies that W equals S on the retained edges and
that the precision is exactly zero off them. The filtered graph is checked for
the planar edge count and chordality, and the relative-importance shares are
checked against the identity that they sum to each node's full-model R-squared.

The certificate is what makes the base-R implementation auditable. It replaces
"trust this because it matches the reference" with "trust this because its
distance from the unique optimum is printed and is near zero." For strictly
convex objectives, a smaller stationarity residual is a mathematical witness of a
fit closer to the unique minimizer. It is also the invariant for any future
optimization: a faster kernel is acceptable exactly when it returns the same
certified residual.

## 4. The textbook graphical lasso: a worked case

The graphical lasso is the most used psychometric network estimator and the one
where the difference between a textbook implementation and a wrapped solver is
sharpest. We develop it in full.

### 4.1 The textbook quantity

The estimator has two specified parts. The penalty path is a log-spaced grid of
rho from the largest off-diagonal correlation down to a small fraction of it. The
selector is the extended Bayesian information criterion of Foygel and Drton
(2010),

    EBIC = -2 L + E log n + 4 E gamma log p,

where L is the Gaussian log-likelihood of the fitted precision, E is the number
of edges, and gamma is the tuning parameter. The criterion has a closed form at
the boundary of the path. When the graph is empty and S is a correlation matrix,
the precision is the identity, the log-likelihood is -n p / 2, no parameters are
penalized, and the criterion equals n p exactly. For p = 6 and n = 120 that is
720.

`psychnet` implements the path and the criterion directly and refits the selected
penalty to machine tolerance. On the empty-graph boundary it returns 720.00, the
textbook value to the digit, and its stationarity residual at the selected
penalty is near `1e-11`.

### 4.2 What the reference does

`qgraph` calls the `glasso` Fortran kernel at its default convergence threshold
of `1e-4`, so the returned precision sits about that far from the optimum. Its
default information criterion is not the closed-form expression above; it is read
from a model-fit object, and on the empty-graph boundary it returns 728.37 rather
than 720. The constant offset between the two criteria does not change which
penalty is selected, because it is the same for every model on the same data.
What changes the selection is the solver residual at the boundary, where the
penalty path empties and the criterion curve is nearly flat.

### 4.3 The consequence: a spurious edge

The two implementations agree on the great majority of datasets. They disagree
when the criterion curve is flat, which happens when the true graph is nearly
empty and the sample is small. We examined the single dataset in a
hundred-dataset comparison where the disagreement reached a visible magnitude. It
is a six-node graph at n = 120 generated from a sparse precision whose six true
partial correlations are all below 0.18 in absolute value, below the detection
threshold at that sample size. Neither implementation recovers the true edges.
They differ on one edge, between variables 2 and 5. The true partial correlation
of that edge is exactly zero. Its sample correlation is -0.24, which is sampling
noise. `psychnet` returns the empty graph. `qgraph` returns the edge at -0.041.

On this dataset the textbook implementation declines the edge and the reference
reports a false positive. The reason is the combination developed above: the
textbook criterion and the optimum-converged solve place the boundary one penalty
step higher than the loosely converged solve does, and at that step the noise edge
is gone. This is the practical content of the precision argument. Solving to the
optimum and scoring with the textbook criterion is cleaner in principle, and on
near-empty graphs it declines edges that a loosely converged solver keeps.

This is a specific result. It concerns one borderline edge on near-empty graphs
at small n; on data with real signal the two implementations select the identical
edge set, as Section 6 shows. Where the two differ at all, the difference favors
the implementation that computes the textbook quantity exactly.

## 5. The cost: run time

A pure-R kernel is slower than a Fortran kernel. Table 2 reports median elapsed
seconds over repeated fits on an Apple-silicon laptop, from the benchmark script
in the package's local directory. Each row holds the data-generating setup, node
count, sample size, penalty grid, and selection rule fixed across the two
implementations. The reference column is the package each estimator reproduces,
and the ratio is the median `psychnet` time divided by the median reference time.

**Table 2. Run time.**

| estimator | size (p, n) | psychnet (s) | reference (s) | ratio |
|---|---|---:|---:|---:|
| EBICglasso | 10, 250 | 0.29 | 0.002 | 146 |
| EBICglasso | 20, 500 | 1.33 | 0.007 | 190 |
| EBICglasso | 30, 500 | 3.22 | 0.017 | 190 |
| Ising | 6, 500 | 0.26 | 0.013 | 20 |
| mgm | 4, 500 | 0.14 | 0.23 | 0.6 |

The graphical lasso carries the largest penalty, about two orders of magnitude,
because the pure-R block-coordinate descent competes against Fortran. The Ising
model is about twenty times slower than `IsingFit`, which calls `glmnet`. The
mixed model is faster than `mgm`, whose wrapper overhead exceeds the cost of the
base-R nodewise loop at this size.

The absolute figures matter more than the ratios. A thirty-node graphical lasso
fits in about three seconds. Psychometric networks are usually smaller than fifty
nodes, where a single fit is a few seconds and the cost in applied work is the
bootstrap, which `psychnet` runs in the same base-R path. The run time is a
property of the current implementation, not of the design: the algorithm is a
readable reference, and the certificate of Section 3 is what any future
compiled-loop optimization would have to preserve. For problems of hundreds of
nodes the Fortran kernel is the appropriate tool.

## 6. Validation

Agreement with the reference packages was checked on real questionnaire data and
on synthetic data with a known generating graph; the synthetic part is the only
one that can report recovery, because the true network there is known. The
scripts that produce the validation tables are in the package's validation
directory and load only from installed CRAN packages. For each comparison, the
preprocessed data or sufficient statistic was held fixed and passed to both
implementations. Structure agreement means equality of the zero/nonzero
off-diagonal edge pattern. Edge-weight agreement is the maximum absolute
difference between matched edge weights after both results are put on the same
scale. Synthetic recovery is reported as the F-measure against the known edge
set. When missing data required pairwise complete correlations, the effective
sample size used in the EBIC was the pairwise complete observation count attached
to the analyzed correlation matrix.

On real questionnaire instruments, the graphical lasso was compared against
`qgraph`. Both estimators received the same correlation matrix, so the comparison
isolates the penalty path, numerical solve, and selection rule from data
preprocessing. Across nineteen instruments, including the 240-item Big Five
inventory distributed with `qgraph`, the Big Five Inventory, state and trait
anxiety scales, depression and PTSD symptom sets, NEO openness, and intelligence
batteries, the two implementations selected the identical edge set on every
instrument, and the largest single edge-weight difference across all of them was
0.008. The Ising model was compared against `IsingFit` on real binary ability
and intelligence items with the same binary matrices, nodewise EBIC settings, and
symmetrization rule; the edge sets agreed exactly or within one edge. On
synthetic data drawn from a known sparse graph, `psychnet` and `qgraph` agreed
and recovered the true edge set at the same F-measure.

The contrast with Section 4 is the substance of the validation. The two
implementations diverge only on near-empty synthetic graphs at small sample size,
the regime that does not arise in real questionnaire networks. On the data the
methods are used for, they are indistinguishable in the edge set and agree on the
weights to the reference solver's own convergence tolerance.

Two convention differences were found and resolved during development. The
Gaussian nodewise EBIC in the mixed model used a profiled-variance deviance,
n log(RSS/n), while `glmnet` and `mgm` use the residual sum of squares itself;
the two differ by a logarithm that compresses the penalty and selects denser
graphs. Aligning the deviance to the residual sum of squares brought the mixed
model into edge-for-edge agreement with `mgm` on continuous data and preserved
exact ground-truth recovery on a known chain. The mixed model also gained the
Loh-Wainwright post-selection threshold that `mgm` applies by default.

One estimator was removed rather than reconciled. A non-convex graphical lasso
with the SCAD, MCP, or atan penalty has no unique solution path, so its one-step
local linear approximation depends on the warm start, the penalty grid, and the
derivative parameterization. It differed from the `GGMncv` package by about 0.2
even on identical input, while recovering structure at least as well. Because it
cannot be made reproducible across implementations by construction, and because
the reference package is little used, it was removed. Convexity is what makes
cross-package agreement and a stationarity certificate possible, and the one
estimator that lacked it was the one that could agree with neither.

## 7. Scope and limitations

`psychnet` implements the established cross-sectional models. It reconstructs
their published specifications in base R; it does not introduce a new estimator,
and the accuracy result of Section 4 follows from computing the textbook quantity
exactly rather than from any change to the model. Its scope is bounded in four
ways. The temporal models, `graphicalVAR` and `mlVAR`, are out of scope. The
mixed model handles Gaussian and binary nodes, not categorical variables with
more than two levels or counts. The network comparison test is defined for
Gaussian graphical models. The package returns a fitted network as a tidy edge
list and provides no plotting method, leaving the figure to the analyst's
preferred tool.

On the graphical lasso, byte-identity with `qgraph` is not the target. The
residual difference between the two is `qgraph`'s Fortran convergence tolerance,
not a formula, so matching it would require either the `glasso` dependency or a
deliberate reproduction of its imprecision. The target is the certified optimum
of the strictly convex objective. Section 4 shows that this target is the more
accurate quantity where the two differ: the empty-graph EBIC has the textbook
value 720, and the optimum-converged fit declines the documented noise edge.

The principal cost is run time, quantified in Section 5: the pure-R graphical
lasso is about two orders of magnitude slower than the Fortran kernel, which suits
the package to the tens-of-nodes problems of psychometrics rather than to problems
of hundreds of nodes.

The package passes `R CMD check` with no errors, warnings, or notes, and ships
180 test expectations across 54 test cases. The certificate functions are
exported, so any fit can be rechecked from the object: `glasso_kkt`,
`glm_lasso_kkt`, `ggm_support_kkt`, `tmfg_certificate`, and `lmg_certificate`.

## References

Chen, J., and Chen, Z. (2008). Extended Bayesian information criteria for model
selection. Biometrika, 95(3), 759-771.

Foygel, R., and Drton, M. (2010). Extended Bayesian information criteria for
Gaussian graphical models. Advances in Neural Information Processing Systems, 23.

Friedman, J., Hastie, T., and Tibshirani, R. (2008). Sparse inverse covariance
estimation with the graphical lasso. Biostatistics, 9(3), 432-441.

Friedman, J., Hastie, T., and Tibshirani, R. (2010). Regularization paths for
generalized linear models via coordinate descent. Journal of Statistical
Software, 33(1), 1-22.

Haslbeck, J. M. B., and Waldorp, L. J. (2020). mgm: Estimating time-varying mixed
graphical models in high-dimensional data. Journal of Statistical Software,
93(8), 1-46.

Massara, G. P., Di Matteo, T., and Aste, T. (2016). Network filtering for big
data: Triangulated Maximally Filtered Graph. Journal of Complex Networks, 5(2),
161-178.

van Borkulo, C. D., Borsboom, D., Epskamp, S., Blanken, T. F., Boschloo, L.,
Schoevers, R. A., and Waldorp, L. J. (2014). A new method for constructing
networks from binary data. Scientific Reports, 4, 5918.

Grömping, U. (2006). Relative importance for linear regression in R: the package
relaimpo. Journal of Statistical Software, 17(1), 1-27.

Liu, H., Lafferty, J., and Wasserman, L. (2009). The nonparanormal: semiparametric
estimation of high-dimensional undirected graphs. Journal of Machine Learning
Research, 10, 2295-2328.

Epskamp, S., Borsboom, D., and Fried, E. I. (2018). Estimating psychological
networks and their accuracy: a tutorial paper. Behavior Research Methods, 50(1),
195-212.
