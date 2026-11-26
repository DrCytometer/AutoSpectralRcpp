# AutoSpectralRcpp 1.3.6 (2026-11-26)

## Bug fixes

- `writeFCS()`: the TEXT delimiter is chosen per file as a character that
appears in no keyword or value, so nothing needs escaping and
`flowCore::read.FCS()` reads the files with its default `emptyValue = TRUE`.
Escaping by doubling is kept as a fallback and now covers keyword names too.
- `writeFCS()`: writes the required `$MODE` keyword, and sets `$PnB` to 32 and
`$PnE` to `0,0` for every parameter, as FCS 3.1 requires for floating-point
data. Matches AutoSpectral 1.9.0, which no longer copies `$MODE` from the raw
file.


# AutoSpectralRcpp 1.3.5 (2026-11-23)

## Bug fixes

- `readFCS()` and `writeFCS()`: FCS files with padding after the final TEXT
delimiter, or with empty keyword values (e.g. NovoCyte files), no longer have
their keyword/value pairs shifted. A doubled delimiter is read as an escaped
literal inside a value (FCS 3.1), and on writing, keywords with a blank name or
value are dropped and delimiters inside values are escaped. Matches
AutoSpectral 1.8.9.


# AutoSpectralRcpp 1.3.4 (2026-10-30)

## New features

- `fix_anchored_batch_rcpp()`: compiled, multithreaded kernel for the
`estimator = "anchored"` pair slopes in AutoSpectral's `fix.my.unmix()` (one
source against every target, warm-started, with iterative target-negativity
masking). AutoSpectral falls back to its R implementation if this is not
available.

## Improvements

- AF candidate scoring in `unmix_autospectral_joint_cpp()` uses a weighted-L2
fluorophore-leakage term in place of weighted-L1, computed for all candidates
at once by a rank-1 update rather than a per-cell pass over every candidate.
Faster, and matches the R pipeline in `unmix.autospectral.joint()`; AF
assignment can differ slightly from 1.3.3.


# AutoSpectralRcpp 1.3.3 (2026-09-22)

## Bug fixes

- A previous attempt to clamp AF candidate spectra that would produce bad results
for residual projection wasn't quite right, and produced hypernegative results
for some fluorophores when AF spectral profiles were nearly collinear with the
fluorophore and thus almost entirely in-span (not on the residual). This is now
guarded properly and matches the R pipeline in `unmix.autospectra.joint()`.


# AutoSpectralRcpp 1.3.1 (2026-09-03)

## Improvements

- More memory management improvements to the C++ pipeline to speed up processing.
Thread-local buffers, fewer heap allocations.

## Bug fixes

- Patched a bug wherein residual scoring for AF candidate selection would inflate
massively due to improperly normalized AF spectra.
- Safe solve to prevent crashes if singular values obtained.
- Reproducible random number generation for `fix.my.unmix()`.


# AutoSpectralRcpp 1.3.0 (2026-08-26)

## New features

- Hand-rolled Huber slope robust linear modelling to support adaptive unmixing.

## Improvements

- Faster FCS reading with reduced memory usage thanks to tips from Paul Heisig.
- Faster unmixing with reduced memory usage thanks to tips from Paul Heisig.
- Set "-O3" flag for C++ compiler to speed up processing.

## Bug fixes

- Fixed issue causing divergence in some events when running `unmix.poisson.fast()`.
Briefly, the weights could collapse to infinity, and are now clamped to the
approximate noise floor. Additional issues with tracking convergence and fallback
are fixed.


# AutoSpectralRcpp 1.2.1 (2026-08-04)

## Bug fixes

- Fix issues causing C++ crash when reading or writing large FCS files


# AutoSpectralRcpp 1.0.5 (2026-02-24)

## New features

- Fast kernel density estimation for gating and plotting in AutoSpectral.


# AutoSpectralRcpp 1.0.0 (2026-02-10)

## New features

- Parallelized C++ assignment and unmixing of per-cell autofluorescence

## Improvements

- Faster (~10x) per-cell optimization using residual-alignment pre-screening.
- Import functions from AutoSpectral rather than duplicating.


# AutoSpectralRcpp 0.2.0 (2025-12-07)

## New features

- Added fast Poisson–IRLS unmixing with incremental updates.
- Added OpenMP support with optimized C++ kernels.
- Implemented new SSM calculation pipeline.

## Improvements

- Faster Poisson–IRLS unmixing with fast QR decomposition
- Better handling of convergence with step halving, deviance monitoring
- Allow early exit if convergence reached
- unmix.wls updated to match AutoSpectral, ensuring non-negative weighting and
a more numerically stable solve.
- Hopefully faster compiler flags.


## Bug fixes

- Initial estimates for IRLS are no longer clamped to non-negative values
- Indentation error in optimize_unmix_rcpp_woodbury

---
