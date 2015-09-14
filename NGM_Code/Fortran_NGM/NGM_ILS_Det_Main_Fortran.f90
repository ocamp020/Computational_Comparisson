Program NGM_ILS_Det_Main_Fortran
! This program solves the neoclassical growth model (NGM) for a given parametrization


! Load toolbox
	use Toolbox
	implicit none 

! Paramters 
	! Production function
	real(dp), parameter :: A=1/2.0d0, alpha=1.0d0/3.0d0, delta=0.0250d0
	! Utility
	real(dp), parameter :: beta=0.980d0
	! Steady state objective
	real(dp)            :: k_ss
	! Grids
	integer , parameter :: n_k=50
    real(dp), parameter :: k_curve=4.0d0
    real(dp)            :: k_min, k_max
    real(dp), dimension(n_k) :: k_grid=0
    ! Tolerance
    real(dp), parameter :: Tol_V    = 1e-10
    integer , parameter :: Max_Iter = 5000
    ! Discrete state space solution 
    real(dp), dimension(n_k)     :: V_dss=0, gk_dss=0, gc_dss=0
    ! Continuous state space solution 
    real(dp), dimension(n_k)     :: V_css=0, gk_css=0, gc_css=0
    ! Auxiliary variable for CSS solution
    real(dp), dimension(n_k)     :: V_old=0, ddV_old=0
    ! Counters and indices
    integer  :: i=1, j=1, l=1, m=1
! Auxiliary variables for Brent
	real(dp), parameter :: Tol_Brent = 1e-10
	real(dp) :: f_hat
	integer  :: iter_brent

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Capital grid	
	! Steady State
	k_ss = ( (1.0_dp/beta + delta - 1.0_dp) / (A*alpha) )**(-1.0_dp/(1.0_dp-alpha))
		write(*,*) "k_ss=", k_ss
	! Grid
	k_min  = 0.00001d0
	k_max  = 2*k_ss
	
	k_grid = grid(k_min,k_max,n_k,k_curve)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Solution for discrete state space
	print*, " "
	write(*,*) "Solution of NGM in discrete state space"
	call NGM_ILS_Det_DSS(V_dss,gk_dss,gc_dss)
	!V_dss = 0
    	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Solution for continuous state space (in capital and labor)
	print*, " "
	write(*,*) "Solution of NGM in continuous state space"
	call NGM_ILS_Det_CSS(V_dss,V_css,gk_css,gc_css)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Reporting Solutions
	print*, "k_grid"
	do i=1,n_k
		write(*,*) i, k_grid(i)
	end do

	print*, " "
	write(*,*) "Value Function discrete state space"
	do i=1,n_k
		write(*,*) k_grid(i), V_dss(i)
	end do

	print*, " "
	write(*,*) "Capital policy function discrete state space"
	do i=1,n_k
		write(*,*) k_grid(i), gk_dss(i)
	end do

	print*, " "
	write(*,*) "Value Function continuous state space"
	do i=1,n_k
		write(*,*) k_grid(i), V_css(i)
	end do

	print*, " "
	write(*,*) "Capital policy function continuous state space"
	do i=1,n_k
		write(*,*) k_grid(i), gk_css(i)
	end do

	print*, " "
	write(*,*) "Comparisson of Solutions"
	do i=1,n_k
		write(*,*) k_grid(i), V_css(i)-V_dss(i)
	end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Contains

! =============================================================================
! NGM_ILS_Det_DSS: Computes value and policy functions for the neoclassical growth model
!          with log utility in consumption and inelastica labor supply and cobb-douglas
!          production function. No productivity shocks.
!          The solution is on a discrete state space, so k' is chosen from a grid.
!
! Usage:  call NGM_DSS(V_dss,gk_dss,gc_dss)
!
! Input:  There are no inputs but the routine uses all parameters defined at the top
!
! Output: V_dss,  real(dp), dim(n_k), Value function at grid (k)
!         gk_dss, real(dp), dim(n_k), Capital policy function at grid (k)
!         gc_dss, real(dp), dim(n_k), Consumption policy function at grid (k)
!
! Note:   This function uses the parameters defined for the program
!
	Subroutine NGM_ILS_Det_DSS(V,g_k,g_c)
		real(dp), dimension(n_k), intent(out) :: V, g_k, g_c
		real(dp), dimension(n_k)     	 :: V_new=0
	    real(dp), dimension(n_k,n_k) 	 :: V_mat=0, U=0
	    integer , dimension(n_k)     	 :: g_k_ind=0
	    real(dp), dimension(n_k,n_k)     :: K_mat=0, Kp_mat=0
	    real(dp) :: diff_V=10
	    integer  :: iter=1

		! Matrices of form (z,k,k',h')
		do i=1,n_k
			K_mat(i,:)  = k_grid(i)
			Kp_mat(:,i) = k_grid(i)
		end do


		! Utility at each combination
			! c   = max( Az(k^alpha) + (1-delta)k - k' , 0 )
			U = log( max( A*(K_mat**alpha)+(1-delta)*K_mat-Kp_mat , 1e-20 ) ) 
		
		! VFI on discrete state space
			! Initialization
				V      = 0
				V_mat  = 0
				diff_V = 1
				iter   = 0
			! VFI
		    do while ((diff_V>Tol_V).and.(iter<Max_Iter))
		        V_new = maxval( U + beta*V_mat , 2 )
		        diff_V = maxval( abs(V-V_new) )
		        V = V_new
		        V_mat = transpose(spread( V , 2 , n_k))
		        
		        iter = iter+1
	    	end do

	    	if (iter<Max_Iter) then
	    		print*, " "
	    		write(*,*) "DSS - Convergence was achieved at iteration:", iter
	    	else
	    		print*, " "
	    		write(*,*) "DSS - Convergence failed - Max_Iter achieved"
	    	end if 

    	! Policy Functions
	    	! Indexes of policy function for capital and labor
		    g_k_ind = maxloc( U + beta*V_mat , 2 ) 

			! Policy Functions in levels
		    g_k = k_grid(g_k_ind)

	        do j=1,n_k
	            g_c(j) = A*(k_grid(j)**alpha)+(1-delta)*k_grid(j)-g_k(j)
	        end do

	end Subroutine NGM_ILS_Det_DSS


! =============================================================================
! NGM_ILS_Det_CSS: Computes value and policy functions for the neoclassical growth model
!          with log utility in consumption, inelastic labor supply and cobb-douglas
!          production function. No Productivity shocks.
!          The solution is on a continuous state space on capital and labor
!
! Usage:  call NGM_CSS(V_0,V_css,gk_css,gc_css)
!
! Input:  V_0, real(dp), dim(n_k), Initial value function
!
! Output: V_css,  real(dp), dim(n_k), Value function at grid (k)
!         gk_css, real(dp), dim(n_k), Capital policy function at grid (k)
!         gc_css, real(dp), dim(n_k), Consumption policy function at grid (k)
!
! Note:   This function uses the parameters defined for the program. Interpolation 
!		  uses cubic splines with secant method for endpoint derivatives.
!
	Subroutine NGM_ILS_Det_CSS(V_0,V,g_k,g_c)
		real(dp), dimension(n_k), intent(in)  :: V_0
		real(dp), dimension(n_k), intent(out) :: V, g_k, g_c
		real(dp), dimension(n_k) :: V_new=0
		real(dp) :: diff_V=10
	    integer  :: iter=0
	    ! Variables for brent
	    real(dp) :: k_min, k_a, k_b, k_c
	    ! Variables for interpolation
	    real(dp), dimension(n_k) :: ddV 

	    V     = V_0

	    do while ((diff_V>Tol_V).and.(iter<Max_Iter))
			! Input for neg_Bellman_T
	    	V_old = V
	    	
	    	! Coefficients for spline interpolation (this is used by neg_Bellman_T)
	    	! The last two elements are the derivatives at end nodes. They are needed only if method is 3
	    	ddV = cspline_coeff(n_k,k_grid,V,2,0.0_dp,0.0_dp)  
	    	ddV_old = ddV

	    	! Otpimization
	    	do j=1,n_k

	    		call Min_Bracket(k_grid(1),k_grid(n_k/2),k_grid(n_k),neg_Bellman_T,k_a,k_b,k_c)
	    		V_new(j) = brent_opt(k_a,k_b,k_c,neg_Bellman_T,Tol_Brent,k_min)
		    	
			end do 

			! The value of V_new is multiplied by -1 since brent minimizes
		    V_new = -V_new 

		    ! Difference between value functions and updating
	    	diff_V = maxval( abs(V-V_new) )
		    V      = V_new

	    	iter = iter+1

	    	write(*,*) "CSS: Iter=", iter, "Diff_V=", diff_V
	    end do

    	if (iter<Max_Iter) then
    		print*, " "
    		write(*,*) "CSS - Convergence was achieved at iteration:", iter
    	else
    		print*, " "
    		write(*,*) "CSS - Convergence failed - Max_Iter achieved"
    	end if 

    	! Policy Functions
	    	do j=1,n_k
	    		! Optimization for kp
	    		call Min_Bracket(k_grid(1),k_grid(n_k/2),k_grid(n_k),neg_Bellman_T,k_a,k_b,k_c)
	    		V_new(j) = brent_opt(k_a,k_b,k_c,neg_Bellman_T,Tol_Brent,k_min)

				! Allocation of optimal policies
				g_k(j) = k_min
			end do 


	end Subroutine NGM_ILS_Det_CSS

! =============================================================================
! Utility: Computes value of utility at given values of states and controls
!          (states are k and controls are kp)
!          Utility is sepparable in consumption and labor and is of the form:
!          U(c,h) = log(c)
!
! Usage:  U = Utility(k,kp)
!
! Input:  k , real(dp), value of current capital k
!		  kp, real(dp), value of future capital kp
!
! Output: Utility,  real(dp), Value of utility at (k,kp)
!
! Note:   If the consumption implied by (k,kp) is non-positive utilty is set
!         to a large negative number.
!		  This function uses the global parameters of the utiltiy function
!
	Function Utility(k,kp)
		real(dp), intent(in)  :: k, kp
		real(dp)              :: Utility
		real(dp)              :: c 

		!write(*,*) "kp in utility", kp

		c = A*(k**alpha)+(1.0_dp-delta)*k-kp

		if (c.ge.0.0_dp) then
			Utility = log(c)
		else
			Utility =  -1e20
		end if

	end Function Utility


! =============================================================================
! neg_Bellman_T: Computes the value of the negative of the bellman operator T
!                evaluated for a given candidate value k'.
!                The function takes as given V_old, and the current state (z,k)
!  				 Values of V_old at k' are computed with cubic spline interpolation
!                Bellman's operator T is:  
!					T[V_old](kp) = max{h}[ U(z,k,kp,h) + beta*E[ V_old(z',k') |z] ]
!                
!
! Usage:  TV = neg_Bellman_T(kp)
!
! Input:  kp, real(dp), value of next period capital kp
!		  
! Implicit Inputs: i      , integer , Index of current z in z_grid
!				   j      , integer , Index of current k in k_grid
!                  V_old  , real(dp), dimension(n_k), Current value function
!                  ddV_old, real(dp), dimension(n_k), Second derivatives of current value function
!
! Output: neg_Bellman_T,  real(dp), Negatie of the bellman operator evaluated at kp
!
! Note:  If the maximum consumption implied by (z,k,kp) is non-positive the candidate kp is 
!        discarded and the value of the operator is set to a large positive number. 
!        If this is not achieved by the specified bounds then a corner solution is assumed.
!		 This function uses the global parameters of the utiltiy function
!        This function is intended to be used with brent_opt in the VFI procedure since 
!        brent_opt is a minimization routine this function returns the negative of the operator. 
!		 Since brent_opt does not allow the function to have more than one argument the 
!        values of (z,k,V_old,ddV_old) are passed through global variables.
!
	Function neg_Bellman_T(kp)
		real(dp), intent(in)         :: kp
		real(dp)                     :: neg_Bellman_T, Bellman_T, h_opt, c_max
		real(dp)                     :: V_tilde
		real(dp), dimension(1)       :: kp_aux, V_aux, dV_aux, ddV_aux

		! Solve the problem only if positive consumption is possible
			! Evaluate consumption if working full time
			c_max = A*(k_grid(j)**alpha)+(1.0_dp-delta)*k_grid(j)-kp
		if (c_max.gt.0.0_dp) then
			! Interpolation of V(k',z') for each value of z'
			kp_aux = kp
			call cspline_IP(n_k,k_grid,V_old,ddV_old,3,kp_aux,1, V_aux, dV_aux, ddV_aux)
			V_tilde = V_aux(1)

			!write(*,*) "V(kp)=", V_tilde

			! Bellman operator on current function V
			Bellman_T = Utility(k_grid(j),kp) + beta*V_tilde
			neg_Bellman_T = -Bellman_T
		else
			neg_Bellman_T = 1e20
		end if

	end Function neg_Bellman_T


end Program NGM_ILS_Det_Main_Fortran


