"""
This script loads a file of double trajectory ensembles and computes through time the
Kolmogorov distance (based on average shited histograms of each ensemble) beween the 
distributions over the phase space of the oscilators defined in PlotSites.
"""


using DrWatson
@quickactivate "ClassicalFPUT"


############################################################
##########  Enter File Name and Plot parameters   ##########
############################################################

filename = "Dbl_LocOsc_TEST_N=2_α=0.0_β=0.0_n_traj=10000.jld2"

PlotSites = [1,2]

colorlist = ["midnightblue", "lightskyblue", "crimson", "orange", "teal", "limegreen", "purple", "plum", "maroon", "orangered"]

fs=25
ts=20

## Ramnges for the average shited histrograms ov the ensembles
rngx = LinRange(-5,5,100)
rngy = LinRange(-5,5,100)




#####################################
##########  Loading File   ##########
#####################################

using JLD2
println(keys(load(datadir("DblEnsembles", filename))))
parameters  = load(datadir("DblEnsembles", filename), "parameters")
periods     = load(datadir("DblEnsembles", filename), "periods")
states1     = load(datadir("DblEnsembles", filename), "states1")
states2     = load(datadir("DblEnsembles", filename), "states2")

n_time = length(periods)
t_start, t_end = first(periods), last(periods)




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

BathSites = parameters["BathSites"]

println()
println()
println("#################################################################################################")
println("##################################   Parameters    ##############################################")
println()
println("Number of oscillators: N=$N")
println("Coupling Potential: V(r)= r²*$κ_round/2! + r³*$α_round/3! + r⁴*$β_round/4!")
println("Baths attached to oscillators $BathSites with coupling strengt γ=$γ and temperature kbT=$kbT")
println()
println("#################################################################################################")
println("#################################################################################################")
println()




using PyPlot, LaTeXStrings, EnergySplitting
using DistanceEnsembles
pygui(true)
include(srcdir("overprint.jl"))

array_distances = Matrix{Float64}(undef, n_time,length(PlotSites))
for i in eachindex(PlotSites)
    site = PlotSites[i]
    @time dis = koldis(states1[[site, N+site] ,: ,:], states2[[site, N+site] ,: ,:], rngx,rngy)
    array_distances[:,i] .= dis
end

fig, ax = subplots(nrows=1, ncols=1)
for p in eachindex(PlotSites)
    ax.plot(periods, array_distances[:,p], color=colorlist[p], label="site $p")
end
ax.plot(periods, array_distances[:,1], c=colorlist[1])
ax.grid()
ax.legend(fontsize="xx-large", loc="upper right")
ax.set_xlim([t_start-0.01, t_end+0.01])
ax.set_xlabel(L"time $t \cdot \frac{2\pi}{\sqrt{\kappa}}$", fontsize=fs)
ax.set_ylim([-0.01, 1.01])
ax.set_ylabel(L"$d_\mathrm{Kol.}(W^\mathrm{class.}_1,W^\mathrm{class.}_2)$", fontsize=fs)
ax.tick_params(labelsize=ts)
show()