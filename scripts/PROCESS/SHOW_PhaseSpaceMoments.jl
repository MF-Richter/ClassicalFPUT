using DrWatson
@quickactivate "ClassicalFPUT"


############################################################
##########  Enter File Name and Plot parameters   ##########
############################################################

filename = "LocOsc_TEST_N=2_α=0.0_β=0.0_n_traj=10000.jld2"

colorlist = ["lightcoral", "midnightblue", "mediumvioletred"]

PlotSite = 1

fs=25
ts=20



#####################################
##########  Loading File   ##########
#####################################

using JLD2

parameters   = load(datadir("DblEnsembles", filename), "parameters")
periods      = load(datadir("DblEnsembles", filename), "periods")
trajectories = load(datadir("DblEnsembles", filename), "states1")
sorting      = load(datadir("DblEnsembles", filename), "sorting")
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

γ   = get(parameters,"γ", nothing)
kbT = get(parameters,"kbT", nothing)
BathSites = get(parameters, "BathSites", nothing)

println("#################################################################################################")
println("##################################   Parameters    ##############################################")
println()
println("Chain length is $N oscillators")
println("Coupling Potential: V(r)= r²*$κ/2! + r³*$α/3! + r⁴*$β/4!")
if γ === nothing
    println("No Markovian baths attached\n")
else
    println("Baths temperature is kbT=$kbT")
    println("Baths sites are $BathSites")
    println("Baths coupling strength is γ=$γ")
end
println("The phase-space coordinates ar sorted as $sorting")
println()
println("#################################################################################################")
println("#################################################################################################")


using PyPlot, LaTeXStrings, PhaseSpaceMoments, LinearAlgebra
pygui(true)

include(srcdir("Sorting_Coordinates.jl"))

means = meandisplacement(trajectories, PlotSite)
μs = meanPhSp(trajectories, PlotSite)
Vars_q, Vars_p, Covs_qp = Covariances(trajectories, PlotSite)

fig, axs = subplots()
fig.suptitle("Classical phase-space meanvalues")
axs.plot(periods, μs[1,:], color=colorlist[1], lw=3, label=L"$\langle q \rangle_W$")
axs.plot(periods, μs[2,:], color=colorlist[2], lw=3, label=L"$\langle p \rangle_W$")
axs.set_xlim([t_start-0.01, t_end+0.01])
axs.set_xlabel(L"time $t \cdot \frac{2\pi}{\sqrt{\kappa}}$", fontsize=fs)
axs.tick_params(labelsize=ts)
axs.legend(fontsize="xx-large", loc="upper right")
axs.grid()
show()

fig, axs = subplots()
fig.suptitle("Classical phase-space (co-) variances")
axs.plot(periods,  Vars_q, color=colorlist[1], lw=3, label=L"$\sigma_{qq}^W$")
axs.plot(periods,  Vars_p, color=colorlist[2], lw=3, label=L"$\sigma_{pp}^W$")
axs.plot(periods, Covs_qp, color=colorlist[3], lw=3, label=L"$\sigma_{qp}^W$")
axs.set_xlim([t_start-0.01, t_end+0.01])
axs.set_xlabel(L"time $t \cdot \frac{2\pi}{\sqrt{\kappa}}$", fontsize=fs)
axs.tick_params(labelsize=ts)
axs.legend(fontsize="xx-large", loc="upper right")
axs.grid()
show()