# ClassicalFPUT

[ClassicalFPUT](https://github.com/MF-Richter/ClassicalFPUT) is a julia project that together with its twin project [QuantumFPUT](https://github.com/MF-Richter/QuantumFPUT) is written to investigate quantum and classical information flow through FPUT chains.

The FPUT chain was firstly simulated by Enrico Fermi, John Pasta, Stanislaw Ulam and Mary Tsingou and consists of oscillating particles coupled to their nearest neighbors by the potential

$V(\Delta q) = \frac{\kappa}{2}\Delta q^2 + \frac{\alpha}{3!}\Delta q^3 + \frac{\beta}{4!}\Delta q^4$

i.e. a harmonic chain of coupling strength $\kappa$ with some anharmonic perturbation of strength $\alpha$ and $\beta$. FPUT studied this chain originally in the context of energy sharing between its normal modes to test the ergodicity hypothesis (which in fact they found not to be fulfilled for this model[1]).

[ClassicalFPUT](https://github.com/MF-Richter/ClassicalFPUT) allows to construct such FPUT chains of arbitrary length and coupling strengths and to initialize them either locally, i.e. excite a single site, or in a normal mode, i.e. an ordered collective excitation of the whole chain. Additionally, it is possible to excite the chain not just in a single state but in an ensemble, i.e. an arbitrary number of initial states distributed in phase-space according to a Gaussian distribution defined by its displacement vector and covariance matrix. The initial states are evolved by solving the Hamiltonian equations of motion using the package [DifferentialEquations](https://github.com/SciML/DifferentialEquations.jl). [ClassicalFPUT](https://github.com/MF-Richter/ClassicalFPUT) further provides modules to compute from such states and state trajectories information flow, correlations, phase-space means, energy splitting between sites and normal modes ect.

To study the information flow through a FPUT chain and the effects of anharmonic coupling to it [QuantumFPUT](https://github.com/MF-Richter/QuantumFPUT) uses the Breuer-Laine-Piilo (BLP) measure based on the trace distance between quantum density operators [2,3,4], while [ClassicalFPUT](https://github.com/MF-Richter/ClassicalFPUT) uses the classical limit of the BLP measure based on the Kolmogorov distance between phase-space distributions [5,6]. The phase-space distributions can be approximated by initializing the chain in an ensemble as described above.



## Author
- Moritz F. Richter


## Installation & Setup
Since the plotting is done by [PyPlot](https://github.com/JuliaPy/PyPlot.jl), make sure that you have a stable version of matplotlib.pyplot for python installed.

This code base is using the [Julia Language](https://julialang.org/) and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/) to make a reproducible scientific project named
To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "ClassicalFPUT"
```
which auto-activate the project and enable local path handling from DrWatson.


## Work flow
Raw data in the form of phase-space trajectories (or Monte-Carlo ensembles of them) is computed
by scripts within the folder scripts/COMPUTE and saved in corresponding subfolders in the
data directory.

In scripts/PROCESS, one finds scripts that load some raw data and process it, i.e., computing distances through time from it, etc. Such processed data might be saved as well in the data
directory. scripts/LOAD includes scripts that simply load such processed data and plot it again.

It is best to run each of the scripts in scripts/COMPUTE in order to generate some raw test data to
be loaded by some processing scripts. If you keep any parameters and names as they are now,
every file should be named such that you can also run the processing scripts directly.


## Acknowledgments
This work has been supported by the German Research Foundation (DFG) through FOR 5099.


## Data Availability
Raw data not included. Please request from author.


## References
[1] S. Ulam (and M. Tsingou) E. fermi, J. Pasta. Studies of nonlinear problems. Los Alamos Report, (LA-1940), 1955.

[2] Heinz-Peter Breuer, Elsi-Mari Laine, and Jyrki Piilo. Measure for the degree of non-Markovian behavior of quantum processes in open systems. Phys. Rev. Lett., 103:210401,2009.

[3] Elsi-Mari Laine, Jyrki Piilo, and Heinz-Peter Breuer. Measure for the non-markovianity of quantum processes. Phys. Rev. A, 81:062115, Jun 2010.

[4] E.-M. Laine, J. Piilo, and H.-P. Breuer. Witness for initial system-environment correlations in open-system dynamics. Europhysics Letters, 92(6):60010, jan 2011.

[5] Moritz Ferdinand Richter, Raphael Wiedenmann, and Heinz-Peter Breuer. Witnessing non-markovianity by quantum quasi-probability distributions. New Journal of Physics, 2022.

[6] Moritz F. Richter and Heinz-Peter Breuer. Phase-space measures of information flow in open systems: A quantum and classical perspective of non-markovianity. Phys. Rev. A, 110:062401, Dec 2024.


## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.


## Release Note
This is still a beta version. The Documentation and single commend in the code are still work in progress