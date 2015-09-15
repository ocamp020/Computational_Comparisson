function NGM_sol()

    #=
    Solves the Neo-Classical Growth Model for a given parameterization

    =#

    include("../../Toolbox/Julia_Tools/toolbox.jl")

    # Parameters for optimization
    const MAX_ITER = 5000
    const TOL_V = 1e-10
    
    # Production Function Parameters
    const A, ALPHA, DELTA = 10.0, 1/3.0, 0.025

    # Utility Parameters
    const BETA = 0.99
    const SIGMA = 2.0
    const NU = 0.5
    chi = 1.0

    # Productivity Parameters
    const RHO_Z = 0.9
    const SIGMA_Z = 0.05

    # Steap State Objective
    const H_SS = 1.0/3

    #Grids
    const N_K = 50
    const N_H = 100
    const N_Z = 3
    const K_MIN = 0.5
    const K_MAX = 20.0
    const K_STEP = (K_MAX - K_MIN)/(N_K - 1)
    const K_CURVE = 1
    const H_MIN = .01
    const H_MAX = .99
    const H_STEP = (H_MAX - H_MIN)/(N_H - 1)
    const H_CURVE = 1
    const LAMBDA_Z = 3

    # Initializing variables
    log_z_grid = Array(Float64, N_Z)
    k_grid  = Array(Float64, N_K)
    h_grid  = Array(Float64, N_H)
    v       = Array(Float64, N_Z, N_K)
    v_new   = Array(Float64, N_Z, N_K)
    g_k     = Array(Float64, N_Z, N_K)
    g_h     = Array(Float64, N_Z, N_K)
    g_c     = Array(Float64, N_Z, N_K)
    v_mat   = Array(Float64, N_Z, N_K)
    g_k_ind = 1
    g_h_ind = 1
    k_mat   = Array(Float64, N_Z, N_K, N_K, N_H)
    kp_mat  = Array(Float64, N_Z, N_K, N_K, N_H)
    h_mat   = Array(Float64, N_Z, N_K, N_K, N_H)
    z_mat   = Array(Float64, N_Z, N_K, N_K, N_H)
    u_h     = Array(Float64, N_Z, N_K, N_K, N_H)

    # Productivity grid and transition matrix
    p_z, log_z_grid = MarkovAR_95(N_Z,RHO_Z,SIGMA_Z)
    z_grid = exp(log_z_grid)

    # Grids for capital and labor
    k_grid = grid(K_MIN,K_MAX,N_K,K_CURVE)
    h_grid = grid(H_MIN,H_MAX,N_H,H_CURVE)

    # Auxiliary matrix for VFI on discrete state space
    for i=1:N_K
        k_mat[:,i,:,:] = k_grid[i]
        kp_mat[:,:,i,:] = k_grid[i]
    end

    for i=1:N_Z
        z_mat[i,:,:,:] = z_grid[i]
    end

    for i=1:N_H
        h_mat[:,:,:,i] = h_grid[i]
    end

    # Evaluate utility at all possible tuples (z,k,k',h)
    if SIGMA == 1.0
        u_h = (log(max(A*z_mat.*(k_mat.^ALPHA).*(h_mat.^(1-ALPHA))+(1-DELTA).*k_mat-kp_mat, 1e-20))
               - chi.*h_mat.^(1-NU)/(1-NU))
    else
        u_h = ((max( A*z_mat.*(k_mat.^ALPHA).*(h_mat.^(1-ALPHA))+(1-DELTA).*k_mat-kp_mat, 1e-20 )
               ).^(1-SIGMA) - 1)/(1-SIGMA) - chi.*h_mat.^(1-NU)/(1-NU)
    end

    # Maximize with respect to h - Obtain tridimensional array of form (z,k,k')
    u = squeeze(maximum(u_h, 4),4)

    # Auxiliary transition matrix such that p_z_k(i,j,l) = p_z(l,i) = Pr(z'=z_i | z=z_l)
    p_z_k = permutedims(repeat(p_z, outer=[1,1,N_K]), [2,3,1])

    # Initialize variables for VFI
    v_mat = zeros(N_Z,N_K)
    diff_v = 1
    iter = 0
    v = zeros(N_Z,N_K)

    # VFI Loop
    while (diff_v > TOL_V) && (iter < MAX_ITER)
        # v_mat is size z, k
        v_new = squeeze(maximum( u .+ BETA * v_mat , 3 ),3)
        diff_v = maximum( abs(v -v_new) )
        v = v_new
        v_mat = permutedims(sum(p_z_k .* v, 1), [3, 1, 2])
        iter += 1

        if iter%50==0
            println("Iteration " , iter , " diff_v " , diff_v)
        end 
    end

    if (iter<MAX_ITER) 
                println("\n DSS - Convergence was achieved at iteration: ", iter)
    else
                println("\n DSS - Convergence failed - Max_Iter achieved")
    end

    # Policy Functions
    for i=1:N_Z
    for j=1:N_K
        g_k_ind = indmax( [u .+ BETA * v_mat][i,j,:] )
        g_h_ind = indmax( [u_h][i,j,g_k_ind,:] )

        g_k[i,j] = k_grid[g_k_ind]
        g_h[i,j] = h_grid[g_h_ind]
        g_c[i,j] = A*z_grid[i]*(k_grid[j]^ALPHA)*(g_h[i,j]^(1-ALPHA))+(1-DELTA)*k_grid[j]-g_k[i,j]
    end
    end
    


    println("\n k_grid")
    for i=1:N_K
        println(i, " ", k_grid[i])
    end
    println("\n ===============================================")
    println("===============================================")

    println("\n Value Function discrete state space")
    println("Capital ", " Value by z")
    for i=1:N_K
        println(k_grid[i], " ", v[:,i])
    end

    println("\n Capital policy function discrete state space")
    println("Capital ", " K' by z")
    for i=1:N_K
        println(k_grid[i], g_k[:,i])
    end


end

@time NGM_sol()