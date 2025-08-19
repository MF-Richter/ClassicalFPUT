module CorrelationEnsemble

    using LinearAlgebra, AverageShiftedHistograms, Base.Threads
    export narrowPoints, griddensity, CorrelationMC

    include("overprint.jl")


    ### ### ###   ##########################################################################   ### ### ###
    ### ### ###   Functions to reduce ensemble to subensembles of a certain cell at one site   ### ### ###
    ### ### ###   ##########################################################################   ### ### ###

    """
    Checks, if a point is within a given cell. A cell is defined by a tuple consisting of a central
    point and a radius arround it.
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
    Takes an ensemble of points in a multi-partite phase space and reduces it to all its points
    with coordinates at site are within a given cell.
    """
    function narrowPoints(
        points::Matrix{Float64},
        site::Int64,
        cell::Tuple{Vector{Float64}, Float64}
        )
        dim = size(points)[1]
        N = Int64(dim/2)
        n_points = size(points)[2]

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

    function narrowPoints(
        points::Matrix{Float64},
        EnvSites::Vector{Int64},
        cell::Tuple{Vector{Float64}, Float64}
        )
        dim = size(points)[1]
        N = Int64(dim/2)
        n_points = size(points)[2]
        n_sites = length(EnvSites)

        indices = []

        for i in 1:n_points
            point = points[:,i]
            if is_point_in_cell(point[[EnvSites; EnvSites+N*ones(Int, n_sites)]], cell) == true
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
    The function'griddensity()' takes an ensemble of points in a 2D space (e.g. phase-space) as
    well as ranges 'rngx' and 'rngy' over its two axis and generates an average shifted histogram
    of this ensemble over the 2D grid defined by rngx and rngy containing the probabilities to
    find a member of the ensemble in the specific cell of the grid. This function is essentially
    a wrapper arround the central functions of AverageShiftedHistograms.jl, see its documentation

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

    """
    For tow 2D histograms, i.e. prob. distributions,  over the same grid 'koldis_grid()'
    computes the Kolmogorov distance between them.
    """
    function koldis_grid(
        griddensity_1::Matrix{Float64},
        griddensity_2::Matrix{Float64}
        )
        return 0.5*sum( abs.(griddensity_1 - griddensity_2) )
    end

    """
    Given two ensembles of points in a 2D space the function 'koldis' generates for both ensem-
    bles an average shifted histogram over the same grid and computes their Kolmogorov distance
    as a measure of distinguishability of the two ensembles. The grid is defined by fixed ranges
    'rngx' and 'rngy'. This reflects an experimental approach where one would measure an actual
    ensemble of points by some aparatus which defins its fixed grid of bins such that each
    ensemble appears as a distribution over this fixed grid. Thus, the apartus can distinguish
    two points only up to the size of its grid cells.
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
    Approximates the correlations measure between the two 'Sites' based on the Kolmogorov distance between
    an ensemble of points and its decorrelation via a Monte-Carlo correlation witness. This means that the
    integral over the phase space coordinates at the first sites are approximated by picking 'n_MC' random
    points from the ensemble, creat the subensembles of all points with coordinates at the first Site with-
    in a cell of radius 'ϵ' arround the randomly picked points and summing up the Kolmogorov distances at
    the second Site between the full ensemble and the subensemble via average shifted histogramms between
    them both.
    """
    function CorrelationMC(
        points::Matrix{Float64},
        Sites::Tuple{Int64,Int64},
        n_MC::Int64,
        ϵ::Float64;

        CondEnsSize::Int64 = 10000,
        qrangeB = LinRange(-5, 5, 100),
        prangeB = LinRange(-5, 5, 100),
        pointsteps::Int64 = 1
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
                pointsB1 = points[[N+Sites[2], Sites[2]] , 1:pointsteps:end]
                pointsB2 = conditioned_points[[N+Sites[2], Sites[2]] , :]
                Corr += koldis(pointsB1, pointsB2, qrangeB, prangeB)/n_MC
                # overprint("Correlation in $counter cells computed")
            end
        end

        return Corr
    end

    function CorrelationMC(
        points::Matrix{Float64},
        SysSite::Int64,  # Site defining the system part
        n_MC::Int64,
        ϵ::Float64;

        SubEnsSize::Int64 = 10000,
        qrangeSys = LinRange(-5, 5, 100),
        prangeSys = LinRange(-5, 5, 100),
        )

        # N... number of chain elements
        # indices... indices of points within the given ensemble
        dim = size(points)[1]
        N = Int(0.5*dim)

        # Sites belonging to the environment and their coordinates in collective phase space
        EnvSites = deleteat!(collect(1:N), SysSite)
        env_coordinates = deleteat!(collect(1:dim), [SysSite, SysSite+N])

        # indices for all the points in the ensemble
        point_indices = 1:size(points)[2]

        counter = 0
        Corr = 0.0

        while counter < n_MC
            # randomly picking a point within the ensemble and setting cell around it
            point_index = rand(point_indices)
            μ_env = points[env_coordinates, point_index]
            cell_env = (μ_env, ϵ)

            # creating subensemble of points, checking if it is large enough and computing Kol. dis. at 2nd site
            subpoints = narrowPoints(points, EnvSites, cell_env)
            if size(subpoints)[2] >= SubEnsSize
                counter += 1
                pointsSys1 = points[[N+SysSite, SysSite] , :]
                pointsSys2 = subpoints[[N+SysSite, SysSite] , :]
                Corr += koldis(pointsSys1, pointsSys2, qrangeSys, prangeSys)/n_MC
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
