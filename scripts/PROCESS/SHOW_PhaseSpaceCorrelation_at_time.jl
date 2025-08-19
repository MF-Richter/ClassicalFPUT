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

# filename = "Dbl_LocOsc_MidBath_N=3_α=0.0_β=0.2_kbT=0.0_γ=0.1_n_traj=10000.jld2"
# folder = "DblEnsembles"

filename1 = "LocMo_HPC+__N=2_α=0.0_β=0.2_kbT=0.0_γ=0.0_n_traj=10000000.jld2"
filename2 = "LocMo_HPC+__N=2_α=0.0_β=0.2_kbT=0.0_γ=0.0_n_traj=10000000.jld2"
folder = "Ensembles"

qrange = LinRange(-5, 5, 300)
prange = LinRange(-5, 5, 300)

PlotSites = [1,2]

# Time or index of plotting
τ = 9
# epsilons = [0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7]
# epsilons = [1.0, 0.5, 0.25]
epsilons = [0.25]


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

# μ = points1[[1,N+1],rnd_i]
μ = [0.0, 0.0]

for ϵ in epsilons

    cell = (μ, ϵ)
    subpoints1  = narrowPoints(points1, 1, cell)
    subensemblesize1 = size(subpoints1)[2]

    subpoints2  = narrowPoints(points2, 1, cell)
    subensemblesize2 = size(subpoints2)[2]

    dkolA = round(koldis(subpoints1[[N+PlotSites[1], PlotSites[1]] ,:], subpoints2[[N+PlotSites[1], PlotSites[1]] ,:], qrange, prange), digits=4)
    dkolB = round(koldis(subpoints1[[N+PlotSites[2], PlotSites[2]] ,:], subpoints2[[N+PlotSites[2], PlotSites[2]] ,:], qrange, prange), digits=4)
    println("The distance between subensembles of size $ϵ arround $μ is dA=$dkolA and dB=$dkolB ")
    println("The 1st subensemble contains $subensemblesize1 and the 2nd $subensemblesize2 points")
    println("")

    fig, axs = subplots(ncols=2, nrows=2)
    # fig_ϵ.suptitle("Time t = $T*2πω")
    fig.suptitle(L"$d_\mathrm{Kol.}(W_1^A,W_2^A)=$"*"$dkolA    "*L"$d_\mathrm{Kol.}(W_1^B,W_2^B)=$"*"$dkolB", fontsize=15)

    W1 = griddensity(subpoints1[[N+PlotSites[1], PlotSites[1]] ,:], qrange,prange)
    W2 = griddensity(subpoints2[[N+PlotSites[1], PlotSites[1]] ,:], qrange,prange)


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


    W1 = griddensity(subpoints1[[N+PlotSites[2], PlotSites[2]] ,:], qrange,prange)
    W2 = griddensity(subpoints2[[N+PlotSites[2], PlotSites[2]] ,:], qrange,prange)

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


# println("The distance between the ensembles is d=", koldis(points1[[N+2, 2] ,:], points2[[N+2, 2] ,:], qrange, prange))

# fig, axs = subplots(ncols=length(PlotSites), nrows=2)
# fig.suptitle("Time t = $T*2πω")

# for p in eachindex(PlotSites)
#     site = PlotSites[p]

#     W1 = griddensity(points1[[site+N, site],:], qrange,prange)
#     W2 = griddensity(points2[[site+N, site],:], qrange,prange)
    
#     if length(PlotSites)==1
#         ax1 = axs[1]
#         ax2 = axs[2]
#     else
#         ax1 = axs[1,p]
#         ax2 = axs[2,p]
#     end

#     ax1.pcolormesh(qrange, prange, W1, cmap="inferno")
#     ax1.set_xlim([minimum(qrange), maximum(qrange)])
#     ax1.set_ylim([minimum(prange), maximum(prange)])
#     ax1.set_xlabel(L"position $q$")
#     ax1.set_ylabel(L"momentum $p$")
#     ax1.set_aspect(1)
#     ax1.grid()

#     ax2.pcolormesh(qrange, prange, W2, cmap="inferno")
#     ax2.set_xlim([minimum(qrange), maximum(qrange)])
#     ax2.set_ylim([minimum(prange), maximum(prange)])
#     ax2.set_xlabel(L"position $q$")
#     ax2.set_ylabel(L"momentum $p$")
#     ax2.set_aspect(1)
#     ax2.grid()
# end

# show()
















# for ϵ in epsilons

#     Δ = ϵ*[cos(θ), sin(θ)]
#     cell1 = (μ+Δ, ϵ)
#     cell2 = (μ-Δ, ϵ)

#     subpoints1 = narrowPoints(points1, 1, cell1)
#     subpoints2 = narrowPoints(points1, 1, cell2)

#     println("The distance between subensembles of size $ϵ is d=",koldis(subpoints1[[N+2, 2] ,:], subpoints2[[N+2, 2] ,:], qrange, prange))

#     fig, axs = subplots(ncols=length(PlotSites), nrows=2)
#     fig.suptitle("Time t = $T*2πω")

#     for p in eachindex(PlotSites)

#         site = PlotSites[p]

#         if length(PlotSites)==1
#             ax1 = axs[1]
#             ax2 = axs[2]
#         else
#             ax1 = axs[1,p]
#             ax2 = axs[2,p]
#         end


#         ax1.plot(points[site+N,1:100:end], points[site,1:100:end], ".", ms=5, color="k", alpha=alpha)
#         ax1.plot(points[site+N,rnd_i], points[site,rnd_i], "x", ms=5, color="darkorange")
#         ax1.set_xlim([minimum(qrange), maximum(qrange)])
#         ax1.set_ylim([minimum(prange), maximum(prange)])

#         ax1.set_xlabel(L"position $q$")
#         ax1.set_ylabel(L"momentum $p$")
#         ax1.set_aspect(1)
#         ax1.grid()


#         ax2.plot(subpoints1[site+N,:], subpoints1[site,:], ".", alpha=0.05, ms=1, color="crimson")
#         ax2.plot(subpoints2[site+N,:], subpoints2[site,:], ".", alpha=0.05, ms=1, color="midnightblue")
#         ax2.set_xlim([minimum(qrange), maximum(qrange)])
#         ax2.set_ylim([minimum(prange), maximum(prange)])

#         ax2.set_xlabel(L"position $q$")
#         ax2.set_ylabel(L"momentum $p$")
#         ax2.set_aspect(1)
#         ax2.grid()
#     end

#     show()

# end


# Δ = epsilon*[cos(θ), sin(θ)]
# cell1 = (μ+Δ, epsilon)
# cell2 = (μ-Δ, epsilon)

# subpoints11 = narrowPoints(points1, 1, cell1)
# subpoints12 = narrowPoints(points1, 1, cell2)
# subpoints21 = narrowPoints(points2, 1, cell1)
# subpoints22 = narrowPoints(points2, 1, cell2)

# println("The distance between subensembles of size $epsilon with initial states +μ is d=",koldis(subpoints11[[N+2, 2] ,:], subpoints12[[N+2, 2] ,:], qrange, prange))
# println("The distance between subensembles of size $epsilon with initial states -μ is d=",koldis(subpoints21[[N+2, 2] ,:], subpoints22[[N+2, 2] ,:], qrange, prange))



# gridnumbers = [30, 60, 90, 120, 150, 180, 210]
# for n in gridnumbers
#     range = LinRange(-5.0, 5.0, n)
#     Corr = CorrelationMC(points, (1,2), 10, 0.5; qrangeB=range, prangeB=range, pointsteps=1000)
#     println("The correlations with $n^2 cells is C=", Corr)
# end