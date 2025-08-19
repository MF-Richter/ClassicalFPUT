using DrWatson
@quickactivate "ClassicalFPUT"


filename = "LocMo_HPC+__N=2_α=0.0_β=0.2_kbT=0.0_γ=0.0_n_traj=10000000.jld2"
folder = "Ensembles"


#####################################
##########  Loading File   ##########
#####################################

using JLD2
import QuantumOptics: Ket, Operator, FockBasis
file = load(datadir(folder, filename))

println("The file has the following keys: \n", keys(file))

println(file["parameters"])