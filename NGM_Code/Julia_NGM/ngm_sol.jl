
#=
Add comments here describing this code

=#

# Include toolbox and libraries
include("../../Toolbox/Julia_Tools/toolbox.jl")

# Global parameters and global variables
    # Parameters for optimization
    const MAX_ITER = 5000
    const TOL_V = 1e-10

    # Production Function Parameters
    const A, ALPHA, DELTA = 0.50, 1/3.0, 0.025

    # Utility Parameters
    const BETA = 0.99
    const SIGMA = 1.0
    chi = 1.0

    # Productivity Parameters
    const RHO_Z = 0.9
    const SIGMA_Z = 0.05

    # Steap State Objective
    const H_SS = 1.0/3
          k_ss = 1  
    
    #Grids
    const N_K      = 50
    const N_H      = 100
    const N_Z      = 3
    const K_MIN    = 0.001
          k_max    = 20.0
    const K_CURVE  = 4
    const H_MIN    = .01
    const H_MAX    = .99
    const H_CURVE  = 1

    # Initializing variables
    log_z_grid = Array(Float64, N_Z)
    k_grid     = Array(Float64, N_K)


# ==============================================================================
# Contains
# ==============================================================================

function NGM_dss(N_K,N_Z)
    # Initializing variables
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


    # Labor grid
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
               - chi.*log(1-h_mat))
    else
        u_h = ((max( A*z_mat.*(k_mat.^ALPHA).*(h_mat.^(1-ALPHA))+(1-DELTA).*k_mat-kp_mat, 1e-20 )
               ).^(1-SIGMA) - 1)/(1-SIGMA) - chi.*log(1-h_mat)
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

        #println("z=",z_grid[i]," k=",k_grid[j]," g_k_ind=",g_k_ind)
    end
    end

    return v, g_k, g_h, g_c
end


# =============================================================================
# FOC_SS: Computes the diffence between the steady state of hours and 
#         the objective value h_ss for a given parameter chi=x
#
# Usage:  diff = FOC_SS(x)
#
# Input:  x   , real, the value of parameter chi
#
# Output: diff  , real, h_ss - Steady_State
#
# Note:   This function uses the parameters defined for the program
#
    function  FOC_SS(x)
        psi_1 = ( (1/BETA + DELTA - 1) / (A*ALPHA) )^(1/(1-ALPHA))
        psi_2 = (1.0-ALPHA) * A * psi_1^(-ALPHA) / ( A*psi_1^(1-ALPHA)-DELTA )
        k_ss  = psi_2/(x+psi_1*psi_2)

        h_ss  = psi_1 * k_ss

        diff  = h_ss - H_SS
        return diff
    end







# ==============================================================================
# Execution
# ==============================================================================

# Steady State
    psi_1 = ( (1/BETA + DELTA - 1) / (A*ALPHA) )^(1/(1-ALPHA))
    psi_2 = (1.0-ALPHA) * A * psi_1^(-ALPHA) / ( A*psi_1^(1-ALPHA)-DELTA )
    k_ss  = psi_2/(chi+psi_1*psi_2)
    println("\nSteady State")
    println("k_ss=", k_ss)

# Capital and Productivity grid
    # Productivity grid and transition matrix
    p_z, log_z_grid = MarkovAR_95(N_Z,RHO_Z,SIGMA_Z)
    z_grid = exp(log_z_grid)

    # Grids for capital and labor
    k_max  = 2*k_ss
    k_grid = grid(K_MIN,k_max,N_K,K_CURVE)

# Solution for discrete state space
@time (v_dss,gk_dss,gh_dss,gc_dss)=NGM_dss()


# Reporting Solutions
    println("\n===============================================")
    println("===============================================")
    println("\n k_grid")
    println( [ [1:N_K]  k_grid ] )

    println("\n===============================================")
    println("===============================================")

    println("\n Value Function discrete state space")
    println("Capital ", " Value by z")
    println( [ k_grid v_dss' ] )

    println("\n Capital policy function discrete state space")
    println("Capital ", " K' by z")
    println( [ k_grid gk_dss' ] )

    