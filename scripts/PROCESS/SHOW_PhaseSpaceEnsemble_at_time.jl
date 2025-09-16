"""
This script loads some ensemble of trajectories from the data/Ensembles and plays its phase space
distribution of the oscillators chosen in 'PlotSites' at a given time or time index.
"""

using DrWatson
@quickactivate "ClassicalFPUT"
include(srcdir("Sorting_Coordinates.jl"))


############################################################
##########  Enter File Name and Plot parameters   ##########
############################################################

filename = "LocOsc_TEST_N=2_α=0.0_β=0.0_n_traj=10000.jld2"
folder   = "Ensembles"

qrange = LinRange(-5, 5, 100)
prange = LinRange(-5, 5, 100)
alpha=0.05

PlotSites = [1,2]

# Time or index of plotting
T = 8.0


#####################################
##########  Loading File   ##########
#####################################

using JLD2

parameters   = load(datadir(folder, filename), "parameters")
periods      = load(datadir(folder, filename), "periods")

t_start, t_end = first(periods), last(periods)
n_time = length(periods)




############################################
##########  Printing Parameters   ##########
############################################

N = parameters["N"]

κ = round(parameters["κ"], digits=2)
α = round(parameters["α"], digits=2)
β = round(parameters["β"], digits=2)

γ   = parameters["γ"]
kbT = parameters["kbT"]

BathSites = parameters["BathSites"]

println("#################################################################################################")
println("##################################   Parameters    ##############################################")
println()
println("Number of oscillators: N=$N")
println("Coupling Potential: V(r)= r²*$κ/2! + r³*$α/3! + r⁴*$β/4!")
println("Baths attached to oscillators $BathSites with coupling strengt γ=$γ and temperature kbT=$kbT")
println()
println("#################################################################################################")
println("#################################################################################################")





using PyPlot
import DistanceEnsembles: griddensity
pygui(true)


if typeof(T)==Int64
    τ = T
elseif typeof(T)==Float64
    τ = round(Int, T/t_end * n_time)
else
    error("No valid data type for time or time index")
end
t = round(periods[τ], digits=2)

points = @time load(datadir(folder, filename), "states1")[:,:,τ]

fig, axs = subplots(ncols=length(PlotSites))
fig.suptitle("Classical ensemble at t = $t "*L"$2\pi/\sqrt{\kappa}$")

for p in eachindex(PlotSites)

    site = PlotSites[p]
    
    if length(PlotSites)==1
        ax = axs
    else
        ax = axs[p]
    end
    
    ax.plot(points[site+N,:], points[site,:], ".", ms=5, color="crimson", alpha=alpha)
    ax.set_xlim([minimum(qrange), maximum(qrange)])
    ax.set_ylim([minimum(prange), maximum(prange)])

    ax.set_xlabel(L"$q$", fontsize=12)
    if p==1
        ax.set_ylabel(L"$p$", fontsize=12)
    end
    ax.set_title("Site $site")
    ax.set_aspect(1)
    ax.tick_params(labelsize=10)
end

show()