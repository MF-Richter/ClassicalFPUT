module EvolveEnsemble
    using DifferentialEquations
    using Base.Threads
    export ensembleevolutionEOM, ensembleevolutionSEOM, evolutionEOM, evolutionSEOM


    include("Sorting_Coordinates.jl")
    
    """
    Reshaping the solution (i.e. the trajectory) such that it is a 2D Array of Float64 where for each subsystems
    the position coordinate is followed by its momentum coordinate, i.e. the array trajectory[:,t] corresponds ToolsBuildFPUT
    the phase-space vector

    x(t) = ( q1(t),p1(t), q2(t),p2(t), ... , qN(t),pN(t) ).
    """
    function sortPhSpTraj__combined_to_subsys(
        solArray  # Array with solution (i.e. trajectory) es returned by Hamiltonian problem ODEsolver
        )
        # read out number of subsystems 'N' and time steps 'n_times'
        N = Int64(size(solArray)[1]/2)
        n_times = size(solArray)[2]
    
        # set array for reshaped solution array
        traj = Array{Float64,2}(undef, 2*N,n_times)
    
        for t in 1:n_times
            traj[:,t] .= sortPhSpVector__combined_to_subsys(solArray[1:2*N,t])
        end
    
        return traj
    end

    """
    Converting the solution array to an 2D Array of Float64 - same as sortPhSpTraj__combined_to_subsys() - yet without reshaping
    the coordinate sequence.
    """
    function convert_solution(
        solArray  # Array with solution (i.e. trajectory) es returned by Hamiltonian problem ODEsolver
        )
        # read out number of subsystems 'N' and time steps 'n_times'
        N = Int64(size(solArray)[1]/2)
        n_times = size(solArray)[2]
    
        # set array for reshaped solution array
        traj = Array{Float64,2}(undef, 2*N,n_times)
    
        for t in 1:n_times
            traj[:,t] .= solArray[1:2*N,t]
        end
    
        return traj
    end



    """
    Evolves an ensemble of initial points 'xpoints' (in subsys sorting) according to the deterministic equations of motion
    defined by 'EOM!' in accordance to the ODE problem and solver requirements of DifferentialEquations.jl for a range 'times'
    of time steps. The result is stroed as the 3D array 'trajectories' where the first index refers to the phase-space coor-
    dinate, the second index marks the specific trajectory wihtin the ensemble and the third index is the time step; thus the
    vector trajectories[:,i,j] are all coordinates of the i-th realisation in the ensemble at the j-th time step.
    The keyword 'combined-sorting' returns the phase-space coordinates in combined sorting if set to 'true', otherwise the
    phase-space arrays are sorted subsys wise (see Sorting_Coordinates.jl for more information about the sortation). The func-
    tion returns the array 'trajectories' as well as a string marking the sorting of the phase-space coordinates.
    Please check out the links

    https://github.com/SciML/DiffEqPhysics.jl/tree/master
    https://docs.sciml.ai/DiffEqDocs/stable/

    for more information about the differential equation solver and handling of DifferentialEquations.jl.
    """
    function ensembleevolutionEOM(
        EOMs!,                      # Equations of motion as required by DifferentialEquations.jl
        times,                     # time steps at which to save the evolution
        xpoints::Matrix{Float64};  # ensemble of initial points in phase-space in subsys sorting
        alg = DP8(),
        combined_sorting::Bool=true,
        kwargs...
        )

        # number of phase-space dimensions, chain elements, points in ensemble and time steps in solution
        dim, n_traj = size(xpoints)
        N = Int64(dim/2)
        n_times = length(times)
        dt = times[2]-times[1]

        # initialise array to store trajectories
        trajectories = Array{Float64,3}(undef, dim,n_traj,n_times)

        # lock to prevent data race when compute trajectories multi-threaded
        trajlock = ReentrantLock()

        # compute trajectories in a multi-threaded loop over init. points
        @threads for i in 1:n_traj
            pq = sortPhSpVector__subsys_to_combined(xpoints[:,i])
            prob = ODEProblem(EOMs!, pq, (first(times),last(times)))
            if combined_sorting==true
                traj = convert_solution( solve(prob, alg, dt=dt, saveat=times; kwargs...) )
            else
                traj = sortPhSpTraj__combined_to_subsys( solve(prob, alg, dt=dt, saveat=times; kwargs...) )
            end
            @lock trajlock trajectories[:,i,:] = traj
        end
        
        if combined_sorting==true
            return trajectories, "combined"
        else
            return trajectories, "subsys"
        end
    end

    
    function evolutionEOM(
        EOMs!,                      # Equations of motion as required by DifferentialEquations.jl
        times,                     # time steps at which to save the evolution 
        xpoint::Vector{Float64};   # initial point in phase-space in subsys sorting
        alg=DP8(),
        combined_sorting::Bool=true,
        kwargs...
        )

        dt = times[2]-times[1]
        pq = sortPhSpVector__subsys_to_combined(xpoint)
        prob = ODEProblem(EOMs!, pq, (first(times),last(times)))
        if combined_sorting==true
            trajectory = convert_solution( solve(prob, alg, dt=dt, saveat=times; kwargs...) )
        else
            trajectory = sortPhSpTraj__combined_to_subsys( solve(prob, alg, dt=dt, saveat=times; kwargs...) )
        end
        
        if combined_sorting==true
            return trajectory, "combined"
        else
            return trajectory, "subsys"
        end
    end

    function ensembleevolutionSEOM(
        Drift!,                     # Deterministic equations of motion as required by DifferentialEquations.jl for StochasticDE
        Diffusion!,                 # Noise term as required by DifferentialEquations.jl for StochasticDE
        times,                     # time steps at which to save the evolution
        xpoints::Matrix{Float64};  # ensemble of initial points in phase-space in subsys sorting
        alg = EulerHeun(),
        combined_sorting::Bool=true,
        kwargs...
        )

        # number of phase-space dimensions, chain elements, points in ensemble and time steps in solution
        dim, n_traj = size(xpoints)
        n_times = length(times)
        dt = times[2]-times[1]

        # initialise array to store trajectories
        trajectories = Array{Float64,3}(undef, dim,n_traj,n_times)

        # lock to prevent data race
        trajlock = ReentrantLock()

        alg = EulerHeun() # gut für kurze Ketten bis N=3 und manchmal auch N=5
        # compute trajectories in multi-threaded loop over init. points
        @threads for i in 1:n_traj
            pq = sortPhSpVector__subsys_to_combined(xpoints[:,i])
            prob = SDEProblem(Drift!, Diffusion!, pq, (first(times),last(times)))
            if combined_sorting==true
                traj = convert_solution( solve(prob, alg, dt=dt, saveat=times; kwargs...) )
            else
                traj = sortPhSpTraj__combined_to_subsys( solve(prob, alg, dt=dt, saveat=times; kwargs...) )
            end
            @lock trajlock trajectories[:,i,:] = traj
        end
        
        if combined_sorting==true
            return trajectories, "combined"
        else
            return trajectories, "subsys"
        end
    end

    function evolutionSEOM(
        Drift!,                     # Deterministic equations of motion as required by DifferentialEquations.jl for StochasticDE
        Diffusion!,                 # Noise term as required by DifferentialEquations.jl for StochasticDE
        times,                     # time steps at which to save the evolution
        xpoint::Vector{Float64};   # initial point in phase-space in subsys sorting
        alg=EulerHeun(),
        combined_sorting::Bool=true,
        kwargs...
        )
        # alg = SKenCarp()
        alg = EulerHeun() # gut für kurze Ketten bis N=3 und manchmal auch N=5
        dt = 0.1+(times[2]-times[1])
        pq = sortPhSpVector__subsys_to_combined(xpoint)
        prob = SDEProblem(Drift!, Diffusion!, pq, (first(times),last(times)))
        if combined_sorting==true
            trajectory = sortPhSpTraj__combined_to_subsys( solve(prob, alg, dt=dt, saveat=times; kwargs...) )
        else
            trajectory = convert_solution( solve(prob, alg, dt=dt, saveat=times; kwargs...) )
        end 

        if combined_sorting==true
            return trajectory, "combined"
        else
            return trajectory, "subsys"
        end
    end

end
