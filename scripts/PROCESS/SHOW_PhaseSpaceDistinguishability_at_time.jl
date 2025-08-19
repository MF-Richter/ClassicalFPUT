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

filename = "Dbl_LocOsc_MidBath_N=3_α=0.0_β=0.2_kbT=0.0_γ=0.1_n_traj=10000.jld2"
folder = "DblEnsembles"

# filename1 = "LocMo_HPC+__N=2_α=0.0_β=0.2_kbT=0.0_γ=0.0_n_traj=30000000.jld2"
# filename2 = "LocMo_HPC-__N=2_α=0.0_β=0.2_kbT=0.0_γ=0.0_n_traj=30000000.jld2"
# folder = "Ensembles"

qmin, qmax = -5.0 , 5.0
pmin, pmax = -5.0 , 5.0

PlotSites = [1,2]

# Time or index of plotting
τ = 3
cellnumbers = [100, 300, 600]


#####################################
##########  Loading File   ##########
#####################################

using JLD2

parameters = load(datadir(folder, filename1), "parameters")
periods    = load(datadir(folder, filename1), "periods")
points1    = load(datadir(folder, filename1), "states1")[:,:,τ]
points2    = load(datadir(folder, filename2), "states1")[:,:,τ]

# parameters = load(datadir(folder, filename), "parameters")
# periods    = load(datadir(folder, filename), "periods")
# points1    = load(datadir(folder, filename), "states1")[:,:,τ]
# points2    = load(datadir(folder, filename), "states2")[:,:,τ]

T = periods[τ]

t_start, t_end = first(periods), last(periods)
n_time = length(periods)
println(n_time)




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



using LinearAlgebra, Random, Distributions
using CorrelationEnsemble
import DistanceEnsembles: koldis

using PyPlot, LaTeXStrings
pygui(true)

n_points = size(points1)[2]
rnd_i = rand(1:n_points)

for n_cells in cellnumbers

    qrange = LinRange(qmin,qmax, n_cells)
    prange = LinRange(pmin,pmax, n_cells)

    dkolA = round(koldis(points1[[N+PlotSites[1], PlotSites[1]] ,:], points2[[N+PlotSites[1], PlotSites[1]] ,:], qrange, prange), digits=4)
    dkolB = round(koldis(points1[[N+PlotSites[2], PlotSites[2]] ,:], points2[[N+PlotSites[2], PlotSites[2]] ,:], qrange, prange), digits=4)
    println("The distance between ensembles is dA=$dkolA and dB=$dkolB")
    println("")

    fig, axs = subplots(ncols=2, nrows=2)
    # fig_ϵ.suptitle("Time t = $T*2πω")
    fig.suptitle(L"$d_\mathrm{Kol.}(W_1^A,W_2^A)=$"*"$dkolA    "*L"$d_\mathrm{Kol.}(W_1^B,W_2^B)=$"*"$dkolB", fontsize=15)

    W1 = griddensity(points1[[N+PlotSites[1], PlotSites[1]] ,:], qrange,prange)
    W2 = griddensity(points2[[N+PlotSites[1], PlotSites[1]] ,:], qrange,prange)


    axs[1,1].pcolormesh(qrange,prange, W1, cmap="gist_heat_r")
    axs[1,1].set_xlim([minimum(qrange), maximum(qrange)])
    axs[1,1].set_ylim([minimum(prange), maximum(prange)])
    axs[1,1].set_ylabel(L"$p_A$", fontsize = 12)
    axs[1,1].set_aspect(1)
    axs[1,1].grid()

    axs[2,1].pcolormesh(qrange,prange, W2, cmap="gist_heat_r")
    axs[2,1].set_xlim([minimum(qrange), maximum(qrange)])
    axs[2,1].set_ylim([minimum(prange), maximum(prange)])
    axs[2,1].set_xlabel(L"$q_A$", fontsize=12)
    axs[2,1].set_ylabel(L"$p_A$", fontsize=12)
    axs[2,1].set_aspect(1)
    axs[2,1].grid()


    W1 = griddensity(points1[[N+PlotSites[2], PlotSites[2]] ,:], qrange,prange)
    W2 = griddensity(points2[[N+PlotSites[2], PlotSites[2]] ,:], qrange,prange)

    axs[1,2].pcolormesh(qrange,prange, W1, cmap="gist_heat_r")
    axs[1,2].set_xlim([minimum(qrange), maximum(qrange)])
    axs[1,2].set_ylim([minimum(prange), maximum(prange)])
    axs[1,2].set_ylabel(L"$p_B$", fontsize = 12)
    axs[1,2].set_aspect(1)
    axs[1,2].grid()

    axs[2,2].pcolormesh(qrange,prange, W2, cmap="gist_heat_r")
    axs[2,2].set_xlim([minimum(qrange), maximum(qrange)])
    axs[2,2].set_ylim([minimum(prange), maximum(prange)])
    axs[2,2].set_xlabel(L"$q_B$", fontsize=12)
    axs[2,2].set_ylabel(L"$p_B$", fontsize=12)
    axs[2,2].set_aspect(1)
    axs[2,2].grid()


    show()

end