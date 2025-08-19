
"""
The module 'EnergySplitting' includes functions to quantify the distribution of the chains total
energy over the local oscilators. Note that while the kenitic energy of a single oscillator is
easily determined only by its momentum the potential energy has a certain arbitrarity to it. One
can easily define the potential between two oscilators using their relative distance but how to
share this potential between thees two oscillators is gives some freedom of choice.
In here the potential between two oscillators is split eaquely.
Yet this leaves half the remaining potential between the first oscillator and the left fixed end
and another half pot. between the last oscillator and the right fixed end. Since thees two remain-
ing half potentials are understood as a potential in which the chain as a whole is put into we
distribute them equaly over all oscillators.
"""
module EnergySplitting

    export localEnergy, totalEnergy

    include("Sorting_Coordinates.jl")
    

    """
    Returning the the FPUT potential V as function of 'r', i.e. the coupling is between particle at site i and
    its neigbor at i+1 is given by V(q_i - q_i+1).
    """
    function Potential(
        κ::Float64,
        α::Float64,
        β::Float64;
        potential_factor::Float64=1.0
        )
        V(r::Float64) = (r^2*κ/2 + r^3*α/6 + r^4*β/24) * potential_factor
        return V
    end

    """
    Returns the energie of the oscillator at site 'i' for a chain of N oscillators coupled by the
    FPUT potential
    V(r) = (r²*κ/2! + r³*α/3! + r⁴*β/4!)
    with fixed ends. The state of the chain is defined by its phase-space ccordinate vector in
    combined sorting (see EvolveEnsemble.jl), i.e. pq = [p1,...,pN, q1,...,qN]. The potential en-
    ergy between two oscillators is devided between the two, the remaining half potentials to the
    fixed ends is distributed equaly over all oscillators as some outer pot. of the whole chain.
    When no index 'i' is passed the function returns a vector containg the energies for at all
    sites.
    """
    function localEnergy(
        pq::Vector{Float64}, # Vector with the local position and momentum coordinates in combiened sorting

        i::Int64,            # index of particle

        κ::Float64,          # harmonic coupling coefficient
        α::Float64,          # cubis coupling coefficient
        β::Float64;          # tetric coupling coefficient

        potential_factor::Float64 = 1.0
        )

        N = Int64(0.5*length(pq))
        p = pq[1:N]
        q = pq[N+1:2*N]

        V = Potential(κ,α,β, potential_factor=potential_factor)
        if i==1
            if length(p)==1
                return 0.5*p[i]^2 + V(-q[i]) + V(q[i])
            else
                return 0.5*p[i]^2 + 0.5*V(-q[i]) + 0.5*V(q[i]-q[i+1])  + 0.5/N*(V(-q[1]) +V(q[N]))
            end
        elseif i==N
            return 0.5*p[i]^2 + 0.5*V(q[i-1]-q[i]) + 0.5*V(q[i])  + 0.5/N*(V(-q[1]) +V(q[N]))
        elseif 1<i<N
            return 0.5*p[i]^2 + 0.5*V(q[i-1]-q[i]) + 0.5*V(q[i]-q[i+1])  + 0.5/N*(V(-q[1]) +V(q[N]))
        else
            error("unvalid particle index! i must be 1<i<$N")
        end
    end

    function localEnergy(
        trajectory::Matrix{Float64}, # Vector with the local position and momentum coordinates in combiened sorting

        i::Int64,            # index of particle

        κ::Float64,          # harmonic coupling coefficient
        α::Float64,          # cubis coupling coefficient
        β::Float64;          # tetric coupling coefficient   

        potential_factor::Float64=1.0,
        combined_sorting=true
        )

        dim = size(trajectory)[1]
        N = Int64(0.5*dim)
        n_time = size(trajectory)[2]
        energies = Vector{Float64}(undef, n_time)

        if combined_sorting==true
            M = 1.0
        else
            M = PermMatrix__subsys_to_combined(N)
        end
        
        for τ in eachindex(energies)
            pq = M*trajectory[:,τ]
            energies[τ] = localEnergy(pq, i, κ,α,β; potential_factor=potential_factor)
        end

        return energies
    end



    function localEnergy(
        pq::Vector{Float64}, # Vector with the local position and momentum coordinates in combiened sorting

        κ::Float64,          # harmonic coupling coefficient
        α::Float64,          # cubis coupling coefficient
        β::Float64;          # tetric coupling coefficient   

        potential_factor::Float64=1.0
        )

        N = Int64(0.5*length(pq))
        energies = Vector{Float64}(undef, N)

        for i in 1:N
            push!(energies, localEnergy(pq, i, κ,α,β, potential_factor=potential_factor))
        end

        return energies
    end

    function localEnergy(
        trajectory::Matrix{Float64}, # Vector with the local position and momentum coordinates in combiened sorting
        κ::Float64,α::Float64,β::Float64;  
        potential_factor::Float64=1.0,
        combined_sorting=true
        )

        dim = size(trajectory)[1]
        N = Int64(0.5*dim)
        n_time = size(trajectory)[2]
        energies = Matrix{Float64}(undef, N,n_time)
        
        for i in 1:N
            energies[i,:] = localEnergy(trajectory, i, κ,α,β; potential_factor=potential_factor, combined_sorting=combined_sorting)
        end

        return energies
    end

    function localEnergy(
        trajectories::Array{Float64,3},
        κ::Float64,α::Float64,β::Float64;
        kwargs... 
        )

        n_traj = size(trajectories)[2]
        energies = 1/n_traj * localEnergy(trajectories[:,1,:], κ,α,β; kwargs...)

        for i in 2:n_traj
            energies += 1/n_traj * localEnergy(trajectories[:,i,:], κ,α,β; kwargs...)
        end

        return energies
    end




    function KineticEnergy(
        p::Vector{Float64},
        )
        T = 0.5*p[1]^2
        for i in 2:length(p)
            T += 0.5*p[i]^2
        end
        return T
    end

    function PotentialEnergy(
        q::Vector{Float64},
        κ::Float64,α::Float64,β::Float64;
        potential_factor::Float64=1.0
        )
        Pot = Potential(κ,α,β, potential_factor=potential_factor)
        N = length(q)
        V =  Pot(-q[1])
        for i in 2:N
            V += Pot(q[i-1]-q[i])
        end
        V += Pot(q[N])
        return V
    end

    function totalEnergy(
        pq::Vector{Float64},
        κ::Float64,α::Float64,β::Float64;
        potential_factor::Float64=1.0       
        )
        N = Int64(0.5*length(pq))
        p = pq[1:N]
        q = pq[N+1:2*N]
        return KineticEnergy(p) + PotentialEnergy(q, κ,α,β, potential_factor=potential_factor)
    end

    function totalEnergy(
        trajectory::Matrix{Float64},
        κ::Float64,α::Float64,β::Float64;
        potential_factor::Float64=1.0,
        combined_sorting::Bool = true      
        )
        dim = size(trajectory)[1]
        N = Int64(0.5*dim)
        n_time = size(trajectory)[2]
        energies = Vector{Float64}(undef, n_time)

        if combined_sorting==true
            M = 1.0
        else
            M = PermMatrix__subsys_to_combined(N)
        end
        
        for τ in eachindex(energies)
            pq = M*trajectory[:,τ]
            energies[τ] = totalEnergy(pq, κ,α,β; potential_factor=potential_factor)
        end

        return energies
    end

    function totalEnergy(
        trajectories::Array{Float64,3},
        κ::Float64,α::Float64,β::Float64;
        kwargs... 
        )

        n_traj = size(trajectories)[2]
        energies = 1/n_traj * totalEnergy(trajectories[:,1,:], κ,α,β; kwargs...)

        for i in 2:n_traj
            energies += 1/n_traj * totalEnergy(trajectories[:,i,:], κ,α,β; kwargs...)
        end

        return energies
    end

end
