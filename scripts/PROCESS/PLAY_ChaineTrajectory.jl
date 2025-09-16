"""
This script loads some single trajectory from the data/Trajectoriers and displays the position displacement
of each oscillator, i.e. the motion of the chain.
Note, that the coordinates need to be sorted combined, i.e. first all momentum coordinates
followed by all position coordinates.
"""


using DrWatson
@quickactivate "ClassicalFPUT"
include(srcdir("Sorting_Coordinates.jl"))


############################################################
##########  Enter File Name and Plot parameters   ##########
############################################################

folder = "Trajectories"
filename = "FPUTeffect_N=10_α=0.5_β=0.0_kbT=0.0_γ=0.0_mode=1.jld2"














#####################################  ###########################################################################################
##########  Loading File   ##########  ###########################################################################################
#####################################  ###########################################################################################

using JLD2

parameters  = load(datadir(folder, filename), "parameters")
periods     = load(datadir(folder, filename), "periods")
trajectory  = load(datadir(folder, filename), "states1")




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




using PyPlot, LaTeXStrings
import NormalModes: Matrix__local_to_nomo
pygui(true)
include(srcdir("overprint.jl"))

Pcts = PermMatrix__combined_to_subsys(N)
Pstc = PermMatrix__subsys_to_combined(N)
M    = Matrix__local_to_nomo(N)

ms = 8
color = "crimson"

fig, ax = subplots(nrows=2, ncols=1)

fig.suptitle("FPUT chain at time t = 0.0000 T", fontsize=14)

pq_init = trajectory[:,1]  #* Pstc
q_init = pq_init[N+1:2*N]
ymax = 1.1*maximum(trajectory)

ax[1].plot(1:N, q_init, ".", ms=ms, color=color)
ax[1].set_xlim([0.5, N+0.5])
ax[1].set_xlabel("sites of oscillators", fontsize=10)
ax[1].set_ylim([-ymax, ymax])
ax[1].set_ylabel(L"displacement $q$ at local sites", fontsize=10)
ax[1].set_title("Local sites of the FPUT chain", fontsize=11)
ax[1].grid()



Q_init = M*q_init
ax[2].plot(1:N, Q_init, ".", ms=ms, color=color)
ax[2].set_xlim([0.5, N+0.5])
ax[2].set_xlabel("indices of normal modes", fontsize=10)
ax[2].set_ylim([-2.1, 2.1])
ax[2].set_ylabel(L"displacement $Q$ in normal modes", fontsize=10)
ax[2].set_title("Normal modes of the FPUT chain", fontsize=11)
ax[2].grid()


sleep(3)

t_0 = periods[1]
println("time t=$t_0*T")
for τ in 1:3:length(periods)

    t = round(periods[τ], digits=4)
    fig.suptitle("FPUT chain at time t = $t T and time step $τ", fontsize=13)

    pq = trajectory[:,τ] #* Pstc
    q = pq[N+1:2*N]

    ax[1].clear()
    ax[1].plot(1:N, q, ".", ms=ms, color=color)
    ax[1].set_xlim([0.5, N+0.5])
    ax[1].set_xlabel("sites of oscillators", fontsize=10)
    ax[1].set_ylim([-ymax, ymax])
    ax[1].set_ylabel(L"displacement $q$ at local sites", fontsize=10)
    ax[1].set_title("Local sites of the FPUT chain", fontsize=11)
    ax[1].grid()
    
    
    global Q = M*q
    ax[2].clear()
    ax[2].plot(1:N, Q, ".", ms=ms, color=color)
    ax[2].set_xlim([0.5, N+0.5])
    ax[2].set_xlabel("indices of normal modes", fontsize=10)
    ax[2].set_ylim([-2.01, 2.01])
    ax[2].set_ylabel(L"displacement $Q$ in normal modes", fontsize=10)
    ax[2].set_title("Normal modes of the FPUT chain", fontsize=11)
    ax[2].grid()
    

    sleep(1/48)
end
show()