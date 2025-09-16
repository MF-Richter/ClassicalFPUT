"""
This script loads an ensemble of trajectories from the data/'folder' computes
the correlation witness between each pair of sites and save it as a dictionary
to be loaded by 'LOAD_Correlations.jl" in the scrips/LOAD folder.
"""

using DrWatson
@quickactivate "ClassicalFPUT"


############################################################
##########  Enter File Name and Plot parameters   ##########
############################################################

filename = "LocOsc_TEST_N=2_α=0.0_β=0.0_n_traj=10000.jld2"
folder   = "Ensembles"
Prefix   = "TEST"

n_mc = 3
ϵ    = 0.3
## Ranges for the average shited histograms over the ensembles
qrangeB = LinRange(-5.0, 5.0, 100)
prangeB = LinRange(-5.0, 5.0, 100)




#####################################
##########  Loading File   ##########
#####################################

using JLD2
println(keys(load(datadir(folder, filename))))
parameters = load(datadir(folder, filename), "parameters")
periods    = load(datadir(folder, filename), "periods")
states     = load(datadir(folder, filename), "states1")

initialstate = states[:,:,1]
n_time = length(periods)



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
if γ == 0
    println("No bath attached")
else
    println("Baths attached to oscillators $BathSites with coupling strengt γ=$γ and temperature kbT=$kbT")
end
println()
println("#################################################################################################")
println("#################################################################################################")
println()




######################################
##########   Compute Data   ##########
######################################
using CorrelationEnsemble

CorrelationDict = Dict{Vector{Int64},Vector{Float64}}()

for i in 1:(N-1)
    for j in (i+1):N

        correlations = @time Corr = CorrelationMC(states, (i,j), n_mc, ϵ; CondEnsSize=1000, qrangeB=qrangeB,prangeB=prangeB)
        
        CorrelationDict[[i,j]] = correlations
        CorrelationDict[[j,i]] = correlations
    end 
end

println("\n...Computation completed\n")



#################################################
##########   Saving DataTrajectories   ########## 
################################################

using JLD2

data = @strdict initialstate parameters periods CorrelationDict

parameters["n_mc"] = n_mc
parameters["ϵ"] = ϵ
filename = datadir("Correlations", savename(Prefix*"Corr_", parameters, "jld2", accesses=["N", "α","β", "n_mc", "ϵ"], sort=false))
safesave(filename, data)
println("data saved\n")

println("##################### Computations completed #####################")
println("##################### ###################### #####################\n\n")