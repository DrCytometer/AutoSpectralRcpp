# AutoSpectralRcpp 1.2.2 (2026-08-26)

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
