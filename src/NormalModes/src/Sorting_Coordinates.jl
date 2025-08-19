"""
Wihtin this module we need to sortings of multi-system phase-space vectors:
i)  Sort x such that for each subsystem position is followed by momentum and so on for the next subsystem, like
        x = (q1,p1, q2,p2, ... qN,pN)
    This sorting we abreviate with 'subsys'; it is used to easily read out coordinates for individual subsystems
ii) Sort x by first all momentum coordinates and than all position coordinates, like
        pq = (p1,p2,...pN, q1,q2,...qN)
    This sorting is abbreviated by 'combined' and used when passing coordinates of the combined phase-space to the ODE solver
The function 'sortPhSpVector_combined_to_subsys()' permutes a phase-space vector of combined sorting to subsys sorting and the
functon 'sortPhSpVector_subsys_to_combined()' permutes it back to combined sorting.
"""
function sortPhSpVector__combined_to_subsys(pq::Vector{Float64})
    N = Int64(length(pq)/2)
    x = Array{Float64,1}(undef,2*N)
    for k in 1:N
        x[2*k-1] = pq[k+N]
        x[2*k]   = pq[k]
    end
    return x
end

function sortPhSpVector__subsys_to_combined(x::Vector{Float64})
    N = Int64(length(x)/2)
    pq = Array{Float64,1}(undef,2*N)
    for k in 1:N
        pq[k+N] = x[2*k-1]
        pq[k]   = x[2*k]
    end
    return pq
end

function PermMatrix__combined_to_subsys(N::Int64)
    M = Array{Float64,2}(undef, 2*N,2*N)
    for i in 1:2*N
        e = zeros(2*N)
        e[i] = 1.0
        M[:,i] = sortPhSpVector__combined_to_subsys(e)
    end
    return M
end

function PermMatrix__subsys_to_combined(N::Int64)
    M = Array{Float64,2}(undef, 2*N,2*N)
    for i in 1:2*N
        e = zeros(2*N)
        e[i] = 1.0
        M[:,i] = sortPhSpVector__subsys_to_combined(e)
    end
    return M
end