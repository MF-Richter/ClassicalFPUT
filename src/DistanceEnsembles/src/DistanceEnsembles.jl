
"""
The module DistanceEnsembles provides functions to turn a given ensembles of points in
phase space into a histogram as discret distribution over phase space, using the package
'AverageShiftedHistograms.jl', as well as functions computing based on this histrogram the
Kolmogorov distance between two such ensembles.
"""
module DistanceEnsembles
    using AverageShiftedHistograms
    export griddensity, koldis

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
        ensemble::Matrix{Float64},
        rngx::AbstractRange,
        rngy::AbstractRange
        )
        histogram = xyz( ash(ensemble[1,:], ensemble[2,:]; rngx=rngx, rngy=rngy) )[3]
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
    ensemble of trajectories by some aparatus which defins its fixed grid of bins such that each
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

    function koldis(
        trajectories1::Array{Float64,3},
        trajectories2::Array{Float64,3},
        rngx::AbstractRange, rngy::AbstractRange
        )
        if size(trajectories1) != size(trajectories2)
            error("ensembles of trajectories are of different size")
        else
            n_times = size(trajectories1)[3]
            distances = Vector{Float64}(undef, n_times)
            for τ in 1:n_times
                dis = koldis(trajectories1[:,:,τ], trajectories2[:,:,τ], rngx,rngy)
                distances[τ] = dis
            end
            return distances
        end
    end

end
