# SLSpectra

A tool for diagnosing spectra in regional ocean simulations from Sturm-Liouville based eigenmodes.

Diagnosing spectra can be generalized to the following process

1. Create a set of basis functions, whose eigenvalues are proportional to a measure of length scale
2. Project data onto the basis functions; the modal coefficients are representative of the amount of data that can be explained by the corresponding eigenmode.

Historically, oceanographers have turned to Discrete Fourier Transforms or other Lagrangian approaches for spectra calculations.  DFTs are challenging for ocean basins, which have irregular boundaries and don't necessarily satsify constraints of periodicity, which the Fourier basis functions demand. Since the Fourier basis functions are essentially the eigenmodes of a Laplacian with periodic boundary conditions, we can generallize the generating boundary value problem to create basis functions for domains with irregular geometry and non-periodic boundary conditions.

The tools in SLSpectra can be used to 
* Diagnose model-consistent sparse Laplacian matrices using only grid information and a subroutine that calculates the action of the Laplacian. This is done using Impulse and Impulse response functions, similar to work done in (FEOTS)[https://github.com/LANL/FEOTS]

* Generate eigenvalues and eigenvectors for a sparse matrix; the eigenvectors form a basis for spectra calculation

* Calculate modal coefficients for the generated eigenvectors, collapse the modal coefficients for degenerate eigenvectors (eigenvectors that share the same eigenvalue), and generate a spectra in irregular geometry
