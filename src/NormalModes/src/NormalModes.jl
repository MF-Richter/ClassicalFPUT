"""
The module 'NormalModes' provides functions to compute normal mode quantities, such as frequency
and energy in a single mode, as well as a transformation of a vector in local mode coordinates
(i.e., position and quadratures of the local sites) into normal mode coordinates (i.e., collective
positional and momentum displacement).
"""
module NormalModes
    include("Sorting_Coordinates.jl")

    function directsum(M1::T, M2::T) where T<:AbstractMatrix
        Tentry = typeof(M1[1,1])
        n1,m1 = size(M1)
        n2,m2 = size(M2)
        return [M1 zeros(Tentry, n1,m2); zeros(Tentry, n2,m1) M2]
    end

    
    export frequencyNormalMode, coordinateNormalMode, Matrix__local_to_nomo, convertensemble__local_to_nomo, nomoEnergy
    


    """
    A function that comnputes the frequency of the normal mode 'mode' of a linear chain of 'N'
    sites assuming a harmonic coupling potential
    V(r) = κ/2 r²
    """
    function frequencyNormalMode(
        N::Int64,       # number of chain elements (i.e. chain length)
        mode::Int64,    # index of the mode
        κ::Float64=1.0  # harmonic coupling coefficient
        )
        return sqrt(κ) * 2*sin(0.5*pi*mode/(N+1))
    end


    """
    Given the vector with local coordinates 'localCs' (i.e., positions or momenta for each site
    of a harmonic chain), this function computes the corresponding coordinate 'nomoCn' of the
    normal mode with index 'mode' of this chain. The length of the chain is determined by the
    length of the given vector.
    If no mode index 'mode' is given, it returns a vector containing the coordinate for all
    normal modes of the chain.
    """
    function coordinateNormalMode(
        localCs::Vector{Float64}, # Vector with the local coordinates
        mode::Int64   # index of the mode
        )

        N = length(localCs)  # number of chain elements (i.e. chain length)
        if mode>N
            error("The chain has only $N normal modes; n=$mode is to large")
        else
            nomoCn = sqrt(2/(N+1)) * localCs[1]*sin(pi*mode*1/(N+1))
            for i in 2:N
                nomoCn += sqrt(2/(N+1)) * localCs[i]*sin(pi*mode*i/(N+1))
            end

            return nomoCn
        end

    end

    function coordinateNormalMode(
        localCs::Vector{Float64}, # Vector with the local coordinates
        )
        N = length(localCs)  # number of chain elements (i.e. chain length)
        nomoCs = Vector{Float64}()
        for mode in 1:N
            push!(nomoCs, coordinateNormalMode(localCs,mode))
        end
        return nomoCs
    end

    """
    For a linear chain of 'N' coupled oscillators, this function returns the basis change
    matrix 'M' between local phase-space coordinates localCs = [c1,...,cN] and normal mode
    phase-space coordinates nomoCs = [C1,...,CN], i.e., nomoCs = M*localCs
    """
    function Matrix__local_to_nomo(N::Int64)
        M = Array{Float64,2}(undef, N,N)
        for i in 1:N
            e = zeros(N)
            e[i] = 1.0
            M[:,i] = coordinateNormalMode(e)
        end
        return M
    end

    """
    Converts an ensemble of points given in local coordinates to nomal mode coordinates
    """
    function convertensemble__local_to_nomo(ensemble)
        N = Int64(size(ensemble)[1]/2)
        n_points = size(ensemble)[2]
        P_stc = PermMatrix__subsys_to_combined(N)
        P_cts = PermMatrix__combined_to_subsys(N)
        M = Matrix__local_to_nomo(N)
    
        ensembleNoMo = Array{Float64, 2}(undef, 2*N, n_points)
    
        for k in 1:n_points
            ensembleNoMo[:,k] .= P_cts*directsum(M,M)*P_stc*ensemble[:,k]
        end
    
        return ensembleNoMo
    end


    """
    Given the vector with local phase-space coordinates 'localPQs' = [p1,...,pN, q1,...,qN] of
    a linear chain of harmonic oscillators coupled by the potential
    V(r) = κ/2 r²,
    this function computes the energy within the normal mode with index 'mode'.
    If no mode index 'n' is given, it returns a vector containing the energies for each normal
    mode of the chain.
    """
    function nomoEnergy(
        localPQs::Vector{Float64}, # Vector with the local position and momentum coordinates in combiened sorting
        mode::Int64;                  # index of the mode
        κ::Float64=1.0            # harmonic coupling coefficient
        )

        N = Int64(0.5*length(localPQs))
        nomoP = coordinateNormalMode(localPQs[1:N],mode)
        nomoQ = coordinateNormalMode(localPQs[N+1:2*N],mode)
        Ω = frequencyNormalMode(N, mode, κ)

        return 1/2 * (nomoP^2 + Ω^2*nomoQ^2)
    end

    function nomoEnergy(
        trajectory::Array{Float64,2}, # Vector with the local position and momentum coordinates in combiened sorting
        mode::Int64;                  # index of the mode
        κ::Float64=1.0,            # harmonic coupling coefficient
        combined_sorting::Bool = true
        )

        dim = size(trajectory)[1]
        N = Int64(0.5*dim)
        n_time = size(trajectory)[2]
        nomoEnergies = Vector{Float64}(undef, n_time)

        if combined_sorting==true
            M = 1.0
        else
            M = PermMatrix__subsys_to_combined(N)
        end

        for i in eachindex(nomoEnergies)
            localPQs = M*trajectory[:,i]
            nomoEnergies[i] = nomoEnergy(localPQs, mode; κ=κ)
        end

        return nomoEnergies
    end

    function nomoEnergy(
        localPQs::Vector{Float64}; # Vector with the local position and momentum coordinates in combiened sorting
        κ::Float64=1.0             # harmonic coupling coefficient
        )
        N = Int64(0.5*length(localPQs))
        nomoEs = Vector{Float64}()
        for mode in 1:N
            push!(nomoEs, nomoEnergy(localPQs, mode, κ=κ))
        end
        return nomoEs
    end

    function nomoEnergy(
        trajectory::Array{Float64,2}; # Vector with the local position and momentum coordinates in combiened sorting
        κ::Float64=1.0,            # harmonic coupling coefficient
        combined_sorting::Bool = true
        )

        dim = size(trajectory)[1]
        N = Int64(0.5*dim)
        n_time = size(trajectory)[2]
        nomoEnergies = Matrix{Float64}(undef, N,n_time)

        if combined_sorting==true
            M = 1.0
        else
            M = PermMatrix__subsys_to_combined(N)
        end

        for i in 1:n_time
            localPQs = M*trajectory[:,i]
            nomoEnergies[:,i] = nomoEnergy(localPQs; κ=κ)
        end

        return nomoEnergies
    end

    function nomoEnergy(
        trajectories::Array{Float64,3};
        kwargs... 
        )

        n_traj = size(trajectories)[2]
        nomoEnergies = 1/n_traj * nomoEnergy(trajectories[:,1,:]; kwargs...)

        for i in 2:n_traj
            nomoEnergies += 1/n_traj * nomoEnergy(trajectories[:,i,:]; kwargs...)
        end

        return nomoEnergies
    end
    
end
