#=============================================================================
! MarkovAR_95: Approximates a contiuous AR(1) process with a discrete Markov process - Rowenhorst 95
!
! Usage: MarkovAR(n_z,rho,sigma)
!
! Input: n_z    , integer , Grid dimension for Markov process
!        rho    , real(dp), Persistence of x where x'=rho*x+sigma*e and e~N(0,1)
!        sigma  , real(dp), Standard Deviation of shocks to x where x'=rho*x+sigma*e and e~N(0,1)
!
! Output: z_grid , real(dp), dimension(n_z),     Grid of values of the Markov process. Size n_z.
!         P      , real(dp), dimension(n_z,n_z), Transition Matrix of z. Prob sum across rows.
!
! Remarks: This routine generates a grid for z (z_grid) and a transition probability matrix (P)
!          The transition probability matrix is organized (z,z') so that it sums to 1 by rows
!          P(i,j) = p_ij is the transition prob of going from z=z_i to z'=z_j
!          Note that this creates matrixes in the opposite order than MarkovAR that follows Tauchen
!          The method follows Rouwenhorst (1995) as shown in Kopecky and Suen (2010)
=#

function MarkovAR_95(n_z, rho, sigma)
    z_grid = Array(Float64, n_z)
    P_z = Array(Float64, n_z, n_z)
    P = Array(Float64n_z-1,n_z-1)

    # Compute z_grid
    psi = sqrt(n_z-1)*sigma/sqrt(1-rho^2)
    z_grid[1] = -psi
    for i=2:n_z
        step = 2*psi/(n_z-1)
        z_grid[i] = z_grid[i-1] + step
    end

    p = (1+rho)/2
    q = (1+rho)/2

    if n_z == 2
        return [p 1-q; 1-p q], z_grid
    else
        P,  = MarkovAR_95(n_z - 1, rho, sigma)
        # To create P_n you take P which is n_z -1 by n_z - 1 and create
        # 4 matrixes which have P in one corner and 0's in the row and column
        # touching the other corner.  For example,
        # [1 2; 3 4] => [1 2 0; 3 4 0; 0 0 0], [ 0 1 2; 0 3 4; 0 0 0] ...
        # plus
        P_a = zeros(n_z+1, n_z+1)
        P_a[2:n_z, 2:n_z] = P
        P_n = (p*P_a[2:n_z+1,2:n_z+1] + (1-p)*P_a[2:n_z+1,1:n_z] +
              (1-q)*P_a[1:n_z,2:n_z+1] + q*P_a[1:n_z,1:n_z])
        P_half = ones(n_z,1)
        P_half[2:n_z-1,1] = .5
        return P_n.*P_half, z_grid
    end
end

#=============================================================================
! grid: Creates grid between x_min and x_max with curvature c_grid
!
! Usage: grid(x_min,x_max,n_grid,c_grid)
!
! Input: x_min   , real(dp), the lower number to consider in the grid to search for zero
!        x_max   , real(dp), the higher number to consider in the grid to search for zero
!        n_grid  , integer , the number of grid points between x_min and x_max
!        c_grid  , real(dp), the curvature of the grid (how close are points near x_min)
!
=#

function grid(x_min,x_max,n_grid,c_grid)
    x=linspace(0.0,1,n_grid)
    return x.^c_grid*(x_max-x_min)+x_min
end



