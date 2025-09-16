"""
'EquationsOfMotionFPUT' is a module providing the EOMs of a FPUT chain as well as tools for the
dissipative and fluctuating terms of thermal baths coupled to single oscillators of the chain.
"""
module EquationsOfMotionFPUT
    using LinearAlgebra
    export HamiltonianFPUT, HamEOM_FPUT, KineticDissipation, NoiseAtSites


    """
    Hamiltonian function of a single harmonic oscillator in 1D
    """
    function HamHO(p, q, params)
        m = params[1]
        ω = params[2]
        return 0.5*(p^2/m + m*ω^2*q^2)
    end


    """
    Kinetic term of the FPUT Hamiltonian function. The number of chain elements is defined by the
    dimesnion of the momentum vector 'p'.
    """
    function KineticFPUT(
        p
        )
        N = length(p)
        T = 0.5*p[1]^2
        for i in 2:N
            T += 0.5*p[i]^2
        end
        return T
    end


    """
    Potential term of the FPUT Hamiltonian function. The number of chain elements is defined by
    the dimesnion of the position vector 'q'.
    """
    function PotentialFPUT(
        q,
        κ::Float64, α::Float64, β::Float64
        )
        N = length(q)

        potentialFPUT(r) = κ/2*r^2 + α/6*r^3 + β/24*r^4

        V = potentialFPUT(-q[1])
        for i in 1:(N-1)
            V += potentialFPUT(q[i] - q[i+1])
        end
        V += potentialFPUT(q[N])
        return V
    end

    """
    Hamiltonian function of a FPUT chain with quadratic (i.e. harmonic) factor of κ/2, cubic fac-
    tor of α/6 and tetric factor of β/24 within the interaction potential.
    """
    function HamFPUT(
        p, q,
        κ::Float64, α::Float64, β::Float64
        )
        return KineticFPUT(p) + PotentialFPUT(q, κ,α,β)
    end


    function EOMposition_FPUT(
        p
        )
        N = length(p)
        qdot = Vector{Float64}()
        for i in 1:N
            push!(qdot, p[i])
        end
        return qdot
    end

    function sglEOMposition_FPUT(
        p
        )
        return p
    end

    function EOMmomentum_FPUT(
        q,
        κ::Float64, α::Float64, β::Float64
        )
        N = length(q)

        gradFPUT(r) = κ*r + α/2*r^2 + β/6*r^3

        pdot = Vector{Float64}()
        push!(pdot, gradFPUT(-q[1]) - gradFPUT(q[1] - q[2]))
        for i in 2:(N-1)
            push!(pdot, gradFPUT(q[i-1] - q[i]) - gradFPUT(q[i] - q[i+1]))
        end
        push!(pdot, gradFPUT(q[N-1] - q[N]) - gradFPUT(q[N]))
        return pdot
    end

    function sglEOMmomentum_FPUT(
        q,
        κ::Float64, α::Float64, β::Float64
        )
        gradFPUT(r) = κ*r + α/2*r^2 + β/6*r^3
        return gradFPUT(-q) - gradFPUT(q)
    end

    """
    Generating the DiffEq of the momentum and position coordinates of a FPUT chain with coupling pot-
    ential 
    """
    function HamEOM_FPUT(
        p, q,
        κ::Float64, α::Float64, β::Float64       
        )
        if length(p)==1
            p = p[1]
            q = q[1]
            return [sglEOMmomentum_FPUT(q, κ,α,β), sglEOMposition_FPUT(p)]
        else
            return [EOMmomentum_FPUT(q, κ,α,β); EOMposition_FPUT(p)]
        end
    end


    """
    The dissipation term for the DiffEq when a damping is attached to one 'site' of the chain
    """
    function KineticDissipation(
        p, q,
        ξ::Float64,
        sites::Vector{Int64}
        )

        N = length(p)
        DissDiagonal = zeros(2*N)
        for i in 1:N
            if i in sites
                DissDiagonal[i] = -ξ
            end
        end

        Dissipator = diagm(DissDiagonal)

        return  Dissipator*[p;q]
    end

    """
    Fluctuative noise term of a bath of temperature 'kbT' attached to several 'sites' of the chain
    """
    function NoiseAtSites(
        N::Int64,
        ξ::Float64,
        kbT::Float64,
        sites::Vector{Int64}
        )
        dx = zeros(2*N)
        for i in 1:N
            if i in sites
                dx[i] = sqrt(2*ξ*kbT)
            end
        end
        return dx
    end

end