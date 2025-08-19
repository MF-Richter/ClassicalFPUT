"""
With this script one can compute a single trajectory for a Fermi-Pasta-Ulam-Tsingou (FPUT) chain,
i.e. a chain of N oscillators coupled position wise to their next neighbor via the coupling poten-
tial

V(x) = κ*x²/2! + α*x³/3! + β*x⁴/4!

i.e. by a harmonic term with strength κ and unharmonic couplings of strength α and β with fixed
boundry contitions (the 1st and the N-th oscillator are coupled to a imagined 0th and (N+1)th osc.
with is forced to stay at rest).
In our setting we assume one osc. to be the actual system of interest (by defoult the 1st) and the
remaining osc. are interpreted as environment. The trajectory is initiated by an initial point in
phase-space. The coordinates of the system osc. are given by displacemnt 'μS' while the coordinates
of the environmental osc. are displacement 'μE'. If the chain should not be initiated by local osc.
one can simply define a arbitrary vector μ = [q1,p1, ... , qN,pN]. Furthermore, we implemented a
keyword such that the first oscilator is interpreted as first normal mode of the harmonic chain ect.
Thus, the excitement of the first oscilator becomes the excitement of the first normal mode and the
chain starts in a collective oscillation instead of a local one.
Additionally one can attache a thermal bath to one oscillator, resulting in a dissipative term in
the Hamiltonian equations of motion and a stochastic Fluctuation term. This two initial terms are
connected by the Fluctuation-Dissipation Theorem. If the temperature is set to zero and thus no
fluctuation is taking place the time evolution is solved as ordinary differential equation by the
ODE solver of DifferentialEquations.jl, otherwise a stochastic differential equation has to be
solved with a SDE solver.
The trajectory is saved by a name containing the parameters and an optional prefix as a .jld2 file
in the data/Trajectory folder.
"""

using DrWatson
@quickactivate "ClassicalFPUT"

##########################################
##########   Enter Parameters   ##########
##########################################

## Initial state of the chain:
# μS... initial phase-space displacement of the system (first oscillator)
# μS... initial phase-space displacement of the environment (oscillators 2 to N)
r = 2.0
θ = +0.5*pi

δr = 0.00
δθ = 0.00

μS = (r+δr) * [cos(θ+δθ), sin(θ+δθ)]
μE = [0.0, 0.0]

# choose, if the subsystems should be interpreted as normal modes of the harmonic chain
initiate_as_normalmode = true
mode = 1


## Time Parameters:
# t_0...    start time
# t_end...  end time
# n_time... number of time steps
t_0    = 0.0
t_end  = 1000.0
n_time = 10000


## Parameters of the FPUT Dynamic:
# N... number of oscillators in the chain
# coupling potential V(x) = κ*x²/2! + α*x³/3! + β*x⁴/4!
N = 32
κ = 1.0
α = 0.5
β = 0.0


## Bath Parameter:
# γ...         coupling strength to the bath (coefficient of dissipative term)
# kbT...       mean thermal energy of the bath
# BathSites... index of oscillators where to attache the bath
γ = 0.0
kbT = 0.0
BathSites = [2]

# optinal prefix for the name of the data file
# Prefix = "LONG_μ0"
Prefix = "FPUTeffect"














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



###########################################
##########   Evolving Ensemble   ##########
###########################################

using EvolveEnsemble, OrdinaryDiffEq
import MakeEnsemble: SystemEnvMeanvalue, directsum
import NormalModes: Matrix__local_to_nomo
include(srcdir("Sorting_Coordinates.jl"))

if initiate_as_normalmode==true
    Pstc = PermMatrix__subsys_to_combined(N)
    Pcts = PermMatrix__combined_to_subsys(N)
    M = directsum(Matrix__local_to_nomo(N), Matrix__local_to_nomo(N))

    μ = Pcts*M*Pstc * SystemEnvMeanvalue(μS, N, siteS=mode)  # initial point in phase space
else
    μ = SystemEnvMeanvalue(μS, N)  # initial point in phase space
end

tspan = (t_0, t_end * 2*pi)
times = LinRange(t_0, t_end*2*pi, n_time)
periods = LinRange(t_0, t_end, n_time)

abstol = 1e-8
reltol = 1e-8

if kbT == 0.0 || γ == 0.0
    abstol = 1e-8
    reltol = 1e-8
    alg = DP8()
    states1, sorting  = @time evolutionEOM(Drift!, times, μ; combined_sorting=true, abstol=abstol, reltol=reltol, alg=alg)
else
    abstol = 1e-4
    reltol = 1e-4
    alg = EulerHeun()
    states1, sorting  = @time evolutionSEOM(Drift!, Diffusion!, times, μ; combined_sorting=true, abstol=abstol, reltol=reltol, maxiters=1e6, alg=alg)
end

println("trajectory evolved")



#############################################
##########   Saving Trajectories   ##########
#############################################

using JLD2

if initiate_as_normalmode==true
    parameters   = @strdict N κ α β γ kbT BathSites μ mode
    name = savename(Prefix, parameters, "jld2", accesses=["N", "α","β", "kbT","γ", "mode"], sort=false)
else
    parameters   = @strdict N κ α β γ kbT BathSites μ
    name = savename(Prefix, parameters, "jld2", accesses=["N", "α","β", "kbT","γ"], sort=false)
end

data = @strdict parameters periods states1 sorting

filename = datadir("Trajectories", name)
safesave(filename, data)
println("trajectory saved")