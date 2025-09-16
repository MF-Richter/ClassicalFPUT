"""
This script loads a single trajectory and computes the energy spillting between its local oscil-
lators as well as the chains normal modes over time.
"""

using DrWatson
@quickactivate "ClassicalFPUT"


############################################################
##########  Enter File Name and Plot parameters   ##########
############################################################

filename = "FPUTeffect_N=10_α=0.5_β=0.0_kbT=0.0_γ=0.0_mode=1.jld2"
folder = "Trajectories"

colorlist = ["midnightblue", "lightskyblue", "crimson", "orange", "teal", "limegreen", "purple", "plum", "maroon", "orangered"]

fs=25
ts=20

NMmax = 4
ScalingMode = 0



#####################################
##########  Loading File   ##########
#####################################

using JLD2

parameters  = load(datadir(folder, filename), "parameters")
periods     = load(datadir(folder, filename), "periods")
trajectory  = load(datadir(folder, filename), "states1")
sorting     = load(datadir(folder, filename), "sorting")
println(sorting)

t_start, t_end = first(periods), last(periods)
n_times = length(periods)




############################################
##########  Printing Parameters   ##########
############################################

N = parameters["N"]

κ_round = round(parameters["κ"], digits=2)
α_round = round(parameters["α"], digits=2)
β_round = round(parameters["β"], digits=2)

κ = parameters["κ"]
α = parameters["α"]
β = parameters["β"]

γ   = parameters["γ"]
kbT = parameters["kbT"]

BathSites = get(parameters,"BathSites", undef)

n_traj = get(parameters, "n_traj", 1)

println("#################################################################################################")
println("##################################   Parameters    ##############################################")
println()
println("Number of oscillators: N=$N")
println("Coupling Potential: V(r)= r²*$κ_round/2! + r³*$α_round/3! + r⁴*$β_round/4!")
if γ === nothing
    println("No Markovian baths attached\n")
else
    println("Baths temperature is kbT=$kbT")
    println("Baths sites are $BathSites")
    println("Baths coupling strength is γ=$γ")
end
println("Averaged over $n_traj trajectories")
println()
println("#################################################################################################")
println("#################################################################################################")



using PyPlot, LaTeXStrings, EnergySplitting
import NormalModes: Matrix__local_to_nomo, nomoEnergy, frequencyNormalMode
pygui(true)
include(srcdir("overprint.jl"))

include(srcdir("Sorting_Coordinates.jl"))
Pstc = PermMatrix__subsys_to_combined(N)
Pcts = PermMatrix__combined_to_subsys(N)

if sorting == "combined"
    combined_sorting = true
elseif sorting == "subsys"
    combinded_sorting = false
else
    error("no valid sorting info accesible")
end

if ScalingMode != 0
    println("frequency of mode $ScalingMode is Ω=",frequencyNormalMode(N, ScalingMode))
    periods *= frequencyNormalMode(N, ScalingMode)
end

energies      = @time localEnergy(trajectory, κ,α,β; combined_sorting=combined_sorting)
totalenergies = @time totalEnergy(trajectory, κ,α,β; combined_sorting=combined_sorting)
energiesNM    = @time nomoEnergy(trajectory; κ=κ,    combined_sorting=combined_sorting)


fig, ax = subplots(nrows=1, ncols=1)
ax.plot(periods, totalenergies, color="k", lw=2, ls="dashed", label=L"total energy $⟨H⟩_\mathbb{E}$")
for k in 1:NMmax
    ax.plot(periods, energiesNM[k,:], lw=2, color=colorlist[k], label="mode $k")
end
ax.grid()
ax.legend(fontsize="xx-large", loc="upper right")
ax.set_xlim([t_start-0.01, t_end+0.01])
ax.set_xlabel(L"time $t / \tau$", fontsize=fs)
ax.set_ylim([-0.01, maximum(totalenergies)*1.1])
ax.set_ylabel(L"energy $⟨\tilde{H}_k⟩_\mathbb{E}$", fontsize=fs)
ax.tick_params(labelsize=ts)
show()