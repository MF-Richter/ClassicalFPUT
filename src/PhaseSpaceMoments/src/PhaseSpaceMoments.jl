"""
This module provides some functions to compute phase-space mean values and (co)variances for an
ensemble of points in a multipartite phase-space. The coordinates of each point should be given
in combined ordering, i.e., x = (p1, ..., pN, q1, ..., qN). If the points are sorted subsystem-
wise, please convert the ensemble first.
"""
module PhaseSpaceMoments

    using Statistics, LinearAlgebra
    export AveragedCoordinate, meanPhSp, meandisplacement, phasespaceangle,
    CovarianceMatrix, CovarianceDeterminant, CovarianceSqueezing, Covariances

    function AveragedCoordinate(
        points::Matrix{Float64},
        coordinate::Int64
        )
        n_points = size(points)[2]
        return sum(points[coordinate,:])/n_points
    end

    function AveragedCoordinate(
        points::Matrix{Float64}
        )
        dim = size(points)[1]
        averages = Vector{Float64}(undef, dim)
        for i in eachindex(averages)
            averages[i] = AveragedCoordinate(points, i)
        end
        return averages
    end 

    function AveragedCoordinate(
        trajectories::Array{Float64,3},
        coordinate::Int64
        )
        n_time = size(trajectories)[3]
        averages = Vector{Float64}(undef, n_time)
        for i in eachindex(averages)
            points = trajectories[:,:,i]
            averages[i] = AveragedCoordinate(points, coordinate)
        end
        return averages
    end

    function AveragedCoordinate(
        trajectories::Array{Float64,3}
        )
        dim = size(trajectories)[1]
        n_time = size(trajectories)[3]
        array_averages = Matrix{Float64}(undef, dim, n_time)
        for i in 1:dim
            array_averages[i,:] .= AveragedCoordinate(trajectories, i)
        end
        return array_averages
    end

    """
    Vector of phase space means
    """
    function meanPhSp(
        points::Matrix{Float64},
        site::Int64;
        )
        N = Int(0.5*size(points)[1])
        q = AveragedCoordinate(points, N+site)
        p = AveragedCoordinate(points, site)
        return [q,p]
    end

    function meanPhSp(
        trajectories::Array{Float64,3},
        site::Int64;
        )

        n_time = size(trajectories)[3]
        means = Matrix{Float64}(undef, 2,n_time)
        for i in 1:n_time
            means[:,i] .= meanPhSp(trajectories[:,:,i], site)
        end
        return means
    end

    """
    Mean displacement in optical phase space for given states
    """    
    function meandisplacement(
        points::Matrix{Float64},
        site::Int64
        )
        return norm(meanPhSp(points, site))
    end

    function meandisplacement(
        trajectories::Array{Float64,3},
        site::Int64;
        )
        n_time = size(trajectories)[3]
        μs = Vector{Float64}(undef, n_time)
        for i in eachindex(μs)
            μs[i] = meandisplacement(trajectories[:,:,i], site)
        end
        return μs
    end

    """
    Angle of mean displacement vector in optical phase space
    """
    function phasespaceangle(
        points::Matrix{Float64},
        site::Int64
        )
        qp = meanPhSp(points, site)
        return atan(qp[1],qp[2])
    end

    function phasespaceangle(
        trajectories::Array{Float64,3},
        site::Int64
        )
        n_time = size(trajectories)[3]
        φs = Vector{Float64}(undef, n_time)
        for i in eachindex(μs)
            φs[i] = phasespaceangle(trajectories[:,:,i], site)
        end
        return φs
    end


    function CovarianceMatrix(
        points::Matrix{Float64}
        )
        return 2*cov(points; dims=2, corrected=true)
    end

    function CovarianceMatrix(
        points::Matrix{Float64},
        site::Int64
        )
        N = Int(0.5*size(points)[1])
        return 2*cov(points[[N+site,site],:]; dims=2, corrected=true)
    end

    function CovarianceMatrix(
        trajectories::Array{Float64,3}
        )
        n_time = size(trajectories)[3]
        CovMs= Vector{Matrix{Float64}}(undef, n_time)
        for i in 1:n_time
            CovMs[i] = CovarianceMatrix(trajectories[:,:,i])
        end
        return CovMs
    end

    function CovarianceMatrix(
        trajectories::Array{Float64,3},
        site::Int64
        )
        n_time = size(trajectories)[3]
        CovMs= Vector{Matrix{Float64}}(undef, n_time)
        for i in 1:n_time
            CovMs[i] = CovarianceMatrix(trajectories[:,:,i], site)
        end
        return CovMs
    end

    function Covariances(
        trajectories::Array{Float64,3},
        site::Int64
        )

        Sigmas = CovarianceMatrix(trajectories, site)
        Vars_q  = []
        Vars_p  = []
        Covs_qp = []
        for i in eachindex(Sigmas)
            Σ = Sigmas[i]
            push!(Vars_q, Σ[1,1])
            push!(Vars_p, Σ[2,2])
            push!(Covs_qp, Σ[1,2])
        end
        return Vars_q, Vars_p, Covs_qp
    end

    function CovarianceDeterminant(
        points::Matrix{Float64},
        site::Int64
        )
        return det(CovarianceMatrix(points, site))  
    end

    function CovarianceDeterminant(
        trajectories::Array{Float64,3},
        site::Int64
        )
        return det.(CovarianceMatrix(trajectories, site))  
    end

    function squeezing(Σ::Matrix{Float64})
        σ1,σ2 = abs.(eigvals(Σ))
        return 1-σ1/σ2
    end

    function CovarianceSqueezing(
        points::Matrix{Float64},
        site::Int64
        )
        return squeezing(CovarianceMatrix(points, site))
    end

    function CovarianceSqueezing(
        trajectories::Array{Float64,3},
        site::Int64
        )
        return squeezing.(CovarianceMatrix(trajectories, site))
    end

end 
