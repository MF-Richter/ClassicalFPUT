"""
The module 'CorrelationEnsemble' provides functions to compute a witness of correlations between
two sites of a bipartite phase-space system. The state is given by an ensemble of points distributed
according to a multipartite distribution function W(xA, xB) over the mutual phase-space
Γ = ΓA × ΓB.

The correlation witness is based on conditioning the ensemble locally at site A to those points
within a certain area γA ⊂ ΓA and computing the distribution function of this conditional ensem-
ble over the phase-space Γ_B at site B, i.e., its marginal distribution. This gives a conditio-
nal probability density WB(xB| γA) (the distribution at site B if only points in γA at site A
are takeninto account).

The correlation can now be witnessed by computing the Kolmogorov distance between the conditional
distribution WB(xB| γA) and the marginal distribution of the full ensemble WB(x_B). In case of no
correlations, the bipartite distribution is
W(xA, xB) = WA(xA) * WB(xB)
and any condition at site A does not change the marginal distribution of site B, i.e., the dis-
tance will be zero. Averaging this local correlation test over a number of random areas, defined
by balls of radius ϵ around some randomly drawn points x_A,i from the ensemble, increases the
accuracy of the correlation witness and minimizes the chance of missing regions with strong
correlations. This defines a Monte-Carlo correlation witness for the whole ensemble.
"""
module CorrelationEnsemble

    using LinearAlgebra, AverageShiftedHistograms, Base.Threads
    export narrowPoints, griddensity, CorrelationMC

    include("overprint.jl")


    ### ### ###   ##########################################################################   ### ### ###
    ### ### ###   Functions to reduce ensemble to subensembles of a certain cell at one site   ### ### ###
    ### ### ###   ##########################################################################   ### ### ###

    """
    Checks, if a point is within a given cell. A cell is defined by a tuple consisting of a central
    point and a radius around it.
    """
    function is_point_in_cell(
        point::Vector{Float64},
        cell::Tuple{Vector{Float64}, Float64}
        )

        μ = cell[1]
        ϵ = cell[2]

        if length(μ)!=length(point)
            error("point and cell have different dimension")
        end

        if norm(point-μ) <= ϵ
            return true
        else
            return false
        end
    end

    """
    Takes an ensemble of points in a multi-partite phase space and computes conditiones it at one
    'site' to a given 'cell' given by a tuple consisting of a central point and a radius around it.
    """
    function narrowPoints(
        points::Matrix{Float64},
        site::Int64,
        cell::Tuple{Vector{Float64}, Float64}
        )
        dim = size(points)[1]  # dimension of phase-space
        N = Int64(dim/2)  # number of subsystems
        n_points = size(points)[2]  # points in the ensemble

        indices = []

        for i in 1:n_points
            point = points[:,i]
            if is_point_in_cell(point[[site,site+N]], cell) == true
                push!(indices, i)
            end
        end

        narrowed_points = Array{Float64, 2}(undef, dim, length(indices))

        for j in eachindex(indices)
            index = indices[j]
            narrowed_points[:,j] .= points[:,index]
        end

        return narrowed_points
    end



    ### ### ###   ####################################################################################################   ### ### ###
    ### ### ###   Functions that compute for two ensembles of points the Kolmogorov distance based on their histograms   ### ### ###
    ### ### ###   ####################################################################################################   ### ### ###

    """
    The function 'griddensity()' takes an ensemble of points in a 2D space (e.g., phase-space) as
    well as ranges 'rngx' and 'rngy' over its two axes and generates an average shifted histogram
    of this ensemble over the 2D grid defined by rngx and rngy, containing the probabilities to
    find a member of the ensemble in the specific cell of the grid. This function is essentially
    a wrapper around the central functions of AverageShiftedHistograms.jl. See its documentation

    https://joshday.github.io/AverageShiftedHistograms.jl/stable/

    for more details.
    """
    function griddensity(
        points::Matrix{Float64},
        rngx::AbstractRange,
        rngy::AbstractRange
        )
        histogram = xyz( ash(points[1,:], points[2,:]; rngx=rngx, rngy=rngy) )[3]
        return histogram / sum(histogram)
    end


    function koldis_grid(
        griddensity_1::Matrix{Float64},
        griddensity_2::Matrix{Float64}
        )
        return 0.5*sum( abs.(griddensity_1 - griddensity_2) )
    end

    """
    Given two ensembles of points in a 2D space, the function 'koldis' generates for both ensembles
    an average shifted histogram over the same grid and computes their Kolmogorov distance
    as a measure of distinguishability of the two ensembles. The grid is defined by fixed ranges
    'rngx' and 'rngy'.
    """
    function koldis(
        points1::Matrix{Float64},
        points2::Matrix{Float64},
        rngx::AbstractRange, rngy::AbstractRange
        )
        density1 = griddensity(points1, rngx, rngy)
        density2 = griddensity(points2, rngx, rngy)
        return koldis_grid(density1, density2)
    end



    ### ### ###   ###############################   ### ### ###
    ### ### ###   Monte-Carlo Correlation Witness   ### ### ###
    ### ### ###   ###############################   ### ### ###

    """
    This function provides the Monte-Carlo correlation witness as described in the module docstring
    (see above). The ensemble is given by the collection 'points' of points in the mutual phase-
    space, and 'Sites' gives the indices of the two subsystems between which the correlations are to
    be computed. The correlation is computed at 'n_MC' balls of radius 'ϵ' around points randomly
    drawn from the ensemble itself.
    If an ensemble of trajectories is given instead of an ensemble of points, the correlations are
    computed for every time step.
    """
    function CorrelationMC(
        points::Matrix{Float64},
        Sites::Tuple{Int64,Int64},
        n_MC::Int64,
        ϵ::Float64;

        CondEnsSize::Int64 = 10000,      # Minimal number of points wihtin one cond. ensemble
        qrangeB = LinRange(-5, 5, 100),  # ranges defining the grid to compute the Kol. dist.s
        prangeB = LinRange(-5, 5, 100),  # ranges defining the grid to compute the Kol. dist.s
        )

        # N... number of chain elements
        # indices... indices of points within the given ensemble
        N = Int(0.5*size(points)[1])
        indices = 1:size(points)[2]

        counter = 0
        Corr = 0.0

        # println("Correlation in $counter cells computed")
        while counter < n_MC
            # randomly picking a point within the ensemble and setting cell around it
            index = rand(indices)
            μA = points[[Sites[1],Sites[1]+N], index]
            cellA = (μA, ϵ)

            # creating subensemble of points, checking if it is large enough and computing Kol. dis. at 2nd site
            conditioned_points = narrowPoints(points, Sites[1], cellA)
            if size(conditioned_points)[2] >= CondEnsSize
                counter += 1
                pointsB1 = points[[N+Sites[2], Sites[2]] , :]
                pointsB2 = conditioned_points[[N+Sites[2], Sites[2]] , :]
                Corr += koldis(pointsB1, pointsB2, qrangeB, prangeB)/n_MC
                # overprint("Correlation in $counter cells computed")
            end
        end

        return Corr
    end

    function CorrelationMC(
        trajectories::Array{Float64,3},
        Sites::Tuple{Int64,Int64},
        n_MC::Int64,
        ϵ::Float64;
        kwargs...
        )

        n_time = size(trajectories)[3]
        correlations = Vector{Float64}(undef, n_time)

        for τ in eachindex(correlations)
            points = trajectories[:,:,τ]
            Corr = CorrelationMC(points, Sites, n_MC, ϵ; kwargs...)
            correlations[τ] = Corr
        end

        return correlations
    end

end
