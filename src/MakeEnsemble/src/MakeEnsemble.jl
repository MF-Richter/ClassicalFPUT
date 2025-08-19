
module MakeEnsemble
    using  Random, Distributions
    export directsum, SystemEnvMeanvalue, SystemEnvCovariances, Gaussian_ensemble
    
    """
    'directsum' computes the algebraic direct sum of two matrices of same type or an array of matrices of same type, i.e. a larger
    matrix of combined dimensions with the input matrices as diagonal blocks. For example:
    
    M1 = [1 2; 2 1] and M2 = [3 0 0; 0 3 0; 0 0 3]
    directsum(M1,M2) = [1 2 0 0 0; 2 1 0 0 0; 0 0 3 0 0; 0 0 0 3 0; 0 0 0 0 3]
    """
    function directsum(M1::T, M2::T) where T<:AbstractMatrix
        Tentry = typeof(M1[1,1])
        n1,m1 = size(M1)
        n2,m2 = size(M2)
        return [M1 zeros(Tentry, n1,m2); zeros(Tentry, n2,m1) M2]
    end

    function directsum(Matrices::Vector{T}) where T<:AbstractMatrix
        M = Matrices[1]
        for i in 2:length(Matrices)
            M = directsum(M, Matrices[i])
        end
        return M
    end


    """
    This function generates the vector containing the mean values (i.e. first statisticval momeents) for
    N-partite phase space system, e.g. a FPUT chain with N oscillators. The oscillator at 'siteS' is as-
    sumed to be the system of interest while the remaining oscillators are accounted as part of the envi-
    ronment. The system oscillator has mean values 'μS' while all the environmental oscillators have mean
    values 'μE'.
    """
    function SystemEnvMeanvalue(
        μS::Vector{Float64},  # mean vaule of the system oscillator
        N::Int64;             # length of chain

        siteS::Int64=1,                 # site of the system oscillator in the chain
        μE::Vector{Float64}=[0.0, 0.0]  # mean value of the bath oscillators
        )

        μ = μS

        if siteS>1
            for i in 2:siteS
                μ = [μE; μ]
            end
        end

        if siteS<N
            for i in (siteS+1):N
                μ = [μ; μE]
            end
        end

        return μ
    end

    
    """
    This function generates the covariance matrix (i.e. first statisticval momeents) for N-partite phase
    space system, e.g. a FPUT chain with N oscillators. The oscillator at 'siteS' is assumed to be the
    system of interest while the remaining oscillators are accounted as part of the environment. The sys-
    tem oscillator is cov. matrix 'ΣS' while all the environmental oscillators have meanvalues 'ΣE'.
    """
    function SystemEnvCovariances(
        ΣS::Matrix{Float64},  # cov. matrix of the system oscillator
        N::Int64;             # length of chain

        siteS::Int64=1,              # site of the system oscillator in the chain
        ΣE::Matrix{Float64}=ΣS*1e-4  # cov. matrix of the bath oscillators
        )

        Σ = ΣS

        if siteS>1
            for i in 2:siteS
                Σ = directsum(ΣE, Σ)
            end
        end

        if siteS<N
            for i in (siteS+1):N
                Σ = directsum(Σ, ΣE)
            end
        end

        return Σ
    end


    """
    A function to create an ensemble of points according to a multivariant normal distribution with mean values 'μ' and
    covariance matrix 'Σ'.
    """
    function Gaussian_ensemble(
        μ::Vector{Float64},  # mean value of ensemble
        Σ::Matrix{Float64},  # cov. matrix of ensemble
        n_points::Int64      # number of points within the ensemble
        )
        return rand(MvNormal(μ, 0.5*Σ), n_points) 
    end

end