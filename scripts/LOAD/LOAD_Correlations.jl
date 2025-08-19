using DrWatson
@quickactivate "QuantumFPUT"


filename = "LONG-Corr__N=4_α=0.0_β=0.0_n_mc=10_ϵ=1.0.jld2"
CorrSites = [1,2]

colorlist = ["midnightblue","orange","crimson","lightcoral","tab:brown","tab:pink"]

fs=25
ts=20


#####################################
##########  Loading File   ##########
#####################################

using JLD2
import QuantumOptics: Ket, Operator, FockBasis

parameters      = load(datadir("Correlations", filename), "parameters")
periods         = load(datadir("Correlations", filename), "periods")
initialstate    = load(datadir("Correlations", filename), "initialstate")
CorrelationDict = load(datadir("Correlations", filename), "CorrelationDict")

correlations = CorrelationDict[CorrSites]
t_start, t_end = first(periods), 20.0 # last(periods)


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

n_mc = parameters["n_mc"]
ϵ    = parameters["ϵ"]

γ   = get(parameters, "γ", nothing)
kbT = get(parameters, "kbT", nothing)

BathSites = get(parameters, "BathSites", nothing)

println()
println()
println("#################################################################################################")
println("##################################   Parameters    ##############################################")
println()
println("Number of oscillators: N=$N")
println("Coupling Potential: V(r)= r²*$κ_round/2! + r³*$α_round/3! + r⁴*$β_round/4!")
println("Correlations witnessed by $n_mc subensembles of radius $ϵ")
if γ === nothing
    println("No bath attached")
else
    println("Baths attached to oscillators $BathSites with coupling strengt γ=$γ and temperature kbT=$kbT")
end
println()
println("#################################################################################################")
println("#################################################################################################")
println()





###################################################
##################   Plot Data   ##################
###################################################

using PyPlot, LaTeXStrings
pygui(true)

# rng = LinRange(-5.0, 5.0, 100)
# fig, axs = subplots()
# axs.pcolormesh(rng, rng, W, vmin=-Wmax, vmax=Wmax, cmap="seismic")
# axs.set_aspect(1)
# axs.grid()
# show()


# fig, axs = subplots()
# fig.suptitle("Correlations in FPUT-Chain")

# axs.plot(periods, correlations, label="full corr. ρ$(CorrSites[1]) and ρ$(CorrSites[2])")
# axs.grid()
# axs.legend()
# axs.set_ylim([-0.01, 1.01])
# axs.set_xlim([t_start-0.01, t_end+0.01])
# show()


fig, axs = subplots()
# fig.suptitle("Correlations in FPUT-Chain")

axs.plot(periods, CorrelationDict[[1,2]], color=colorlist[1], lw=2,  label=L"$\Gamma^{1 + 2}$")
axs.plot(periods, CorrelationDict[[1,3]], color=colorlist[2], lw=1,  label=L"$\Gamma^{1 + 3}$", alpha=0.7)
axs.plot(periods, CorrelationDict[[1,4]], color=colorlist[3], lw=1,  label=L"$\Gamma^{1 + 4}$", alpha=0.7)
axs.grid()
axs.legend(fontsize="xx-large", loc="upper right")
axs.set_ylim([-0.01, 0.31])
axs.set_ylabel(L"$\mathcal{\tilde{C}}_\Gamma(W)$", fontsize=fs)
axs.set_xlim([t_start-0.01, t_end+0.01])
axs.set_xlabel(L"time $t \cdot \frac{2\pi}{\sqrt{\kappa}}$", fontsize=fs)
axs.tick_params(labelsize=ts)
show()