"""
With this script, one can compute two ensembles of trajectories for a Fermi-Pasta-Ulam-Tsingou (FPUT)
chain, i.e., a chain of N oscillators coupled position-wise to their next neighbor via the coupling
potential

V(x) = κx²/2! + αx³/3! + β*x⁴/4!

with a harmonic term with strength κ and unharmonic couplings of strength α and β with fixed
boundary conditions (the 1st and the N-th oscillator are coupled to an imagined 0th and (N+1)th oscillator
that is forced to stay at rest).
In our setting, we assume one oscillator to be the actual system of interest (by default the 1st) and the
remaining oscillators are interpreted as the environment. The ensemble is initiated by 'n_traj' initial
points in phase-space. The coordinates of the system oscillator are randomly drawn from a Gaussian distribution
of displacement '+/-μS' and covariance matrix 'ΣS', while the coordinates of the environmental
oscillators are drawn from a distribution with displacement and covariance matrix 'μE' and 'ΣE'. If
the chain should not be initiated by local oscillators, the module 'MakeEnsemble' offers some tools
to generate a Gaussian distribution over the whole combined phase-space.
Additionally, one can attach a thermal bath to one oscillator, resulting in a dissipative term in
the Hamiltonian equations of motion and a stochastic fluctuation term. These two terms are
connected by the Fluctuation-Dissipation Theorem. If the temperature is set to zero and thus no
fluctuation is taking place, the time evolution is solved as an ordinary differential equation by the
ODE solver of DifferentialEquations.jl; otherwise, a stochastic differential equation has to be
solved with an SDE solver.
The ensemble is saved by a name containing the parameters and an optional prefix as a .jld2 file
in the data/Ensemble folder.
"""

using DrWatson
@quickactivate "ClassicalFPUT"
include(srcdir("Sorting_Coordinates.jl"))



##########################################
##########   Enter Parameters   ##########
##########################################

## Initial state of the chain:
# μS... mean value of initial distribution within the system (first) oscillators phase-space
# ΣS... covariance matrix of initial distribution within the system (first) oscillators phase-space
r = 1.5
θ = 0.5*pi
μS = r * [cos(θ), sin(θ)]
ΣS = [1.0 0.0; 0.0 1.0]

# μE... mean value of initial distribution within the environmental (2nd to Nth) oscillators phase-space
# ΣE... covariance matrix of initial distribution within the environmental (2nd to Nth) oscillators phase-space
μE = [0.0, 0.0]
ΣE = 1.0*ΣS

# choose, if the subsystems should be interpreted as normal modes of the harmonic chain
initiate_as_normalmode = false
mode = 1

# number of initial points or trajectories in the ensemble
n_traj = 10000


## Time Parameters:
# t_0...    start time
# t_end...  end time
# n_time... number of time steps
t_0    = 0.0
t_end  = 10.0
n_time = 100


## Parameters of the FPUT Dynamic:
# N... number of oscillators in the chain
# coupling potential V(x) = κ*x²/2! + α*x³/3! + β*x⁴/4!
N = 2
κ = 1.0
α = 0.0
β = 0.0


## Bath Parameter:
# γ...         coupling strength to the bath (coefficient of dissipative term)
# kbT...       mean thermal energy of the bath
# BathSites... index of oscillators where to attache the bath
γ = 0.0
kbT = 0.0
BathSites = [3]


# optional prefix for the name of the data file
Prefix = "TEST"














############################################  ###########################################################################################
##########   Generating (S)EOMs   ##########  ###########################################################################################
############################################  ###########################################################################################

using EquationsOfMotionFPUT

function Drift!(dx,x, params, t)
    N = Int64(0.5*length(x))
    p = x[1:N]
    q = x[N+1:2*N]
    dx .= HamEOM_FPUT(p,q, κ,α,β) + KineticDissipation(p,q, γ, BathSites)
end

function Diffusion!(dx,x, params, t)
    N = Int64(0.5*length(x))
    dx .= NoiseAtSites(N, γ,kbT  * 0.5/sqrt(2), BathSites)
end



#############################################
##########   Generating Ensemble   ##########
#############################################

using MakeEnsemble
import NormalModes: convertensemble__local_to_nomo

if initiate_as_normalmode==true && N > 1

    μ = SystemEnvMeanvalue(μS, N, siteS=mode)
    Σ = SystemEnvCovariances(ΣS, N, ΣE=ΣE)

    points = convertensemble__local_to_nomo(Gaussian_ensemble(μ, Σ, n_traj))
    
elseif N == 1

    μ = μS
    Σ = ΣS

    points = Gaussian_ensemble(μ, Σ, n_traj)

else

    μ = SystemEnvMeanvalue(μS, N)
    Σ = SystemEnvCovariances(ΣS, N, ΣE=ΣE)

    points = Gaussian_ensemble(μ, Σ, n_traj)

end

println("ensemble initialized")



###########################################
##########   Evolving Ensemble   ##########
###########################################

using EvolveEnsemble, OrdinaryDiffEq
tspan = (t_0, t_end * 2*pi)
times = LinRange(t_0, t_end*2*pi, n_time)
periods = LinRange(t_0, t_end, n_time)

if kbT == 0.0 || γ==0.0
    abstol = 1e-8
    reltol = 1e-8
    alg = DP8()
    states1, sorting  = @time ensembleevolutionEOM(Drift!, times, points; alg=alg, combined_sorting=true, abstol=abstol, reltol=reltol)
else
    abstol = 1e-6
    reltol = 1e-6
    alg = EulerHeun()
    states1, sorting  = @time ensembleevolutionSEOM(Drift!, Diffusion!, times, points; alg=alg, combined_sorting=true, abstol=abstol, reltol=reltol)
end

println("ensemble evolved")



#############################################
##########   Saving Trajectories   ##########
#############################################

using JLD2

if initiate_as_normalmode==true
    parameters   = @strdict N κ α β γ kbT BathSites n_traj μ Σ mode abstol reltol alg=alg
    name = savename("NoMo_"*Prefix, parameters, "jld2", accesses=["N", "α","β", "kbT","γ", "mode", "n_traj"], sort=false)
else
    parameters   = @strdict N κ α β γ kbT BathSites n_traj μ Σ abstol reltol alg=alg
    if γ==0.0
        name = savename("LocOsc_"*Prefix, parameters, "jld2", accesses=["N", "α","β", "n_traj"], sort=false)
    else
        name = savename("LocOsc_"*Prefix, parameters, "jld2", accesses=["N", "α","β", "kbT","γ", "n_traj"], sort=false)
    end
end

data = @strdict parameters periods states1 sorting

filename = datadir("Ensembles", name)
safesave(filename, data)

println("ensemble saved")