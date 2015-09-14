Program NGM_Main_Fortran
! This program solves the neoclassical growth model (NGM) for a given parametrization


! Load toolbox
	use Toolbox
	implicit none 

! Paramters 
	! Production function
	real(dp), parameter :: A=1/2.0d0, alpha=1.0d0/3.0d0, delta=0.0250d0
	! Utility
	real(dp), parameter :: beta=0.990d0
	real(dp)            :: chi=1
	! Productivity process
	real(dp), parameter :: rho_z=0.90d0, sigma_z=0.050d0 
	! Steady state objective
	real(dp), parameter :: h_ss_obj = 1.0d0/3.0d0
	real(dp)            :: k_ss, psi_1, psi_2
	! Grids
	integer , parameter :: n_k=50, n_h=100, n_z=3
    real(dp), parameter :: k_curve=4.0d0
    real(dp)            :: k_min, k_max
    real(dp), parameter :: h_min=0.01, h_max=0.99, h_step=(h_max-h_min)/(n_h-1), h_curve=1.0d0
    ! Tolerance
    real(dp), parameter :: Tol_V    = 1e-10
    integer , parameter :: Max_Iter = 5000
	! Grids and such
    real(dp), dimension(n_z)         :: log_z_grid=0, z_grid=0
    real(dp), dimension(n_z,n_z)     :: P_z
	real(dp), dimension(n_k)         :: k_grid=0
    ! Discrete state space solution 
    real(dp), dimension(n_z,n_k)     :: V_dss=0, gk_dss=0, gh_dss=0, gc_dss=0, res_dss=0
    ! Continuous state space solution (Optimization method or root finding)
    real(dp), dimension(n_z,n_k)     :: V_css=0, gk_css=0, gh_css=0, gc_css=0, res_css=0
    real(dp), dimension(n_z,n_k)     :: V_css_rf=0, gk_css_rf=0, gh_css_rf=0, gc_css_rf=0, res_css_rf=0
    ! Auxiliary variable for CSS solution
    real(dp), dimension(n_z,n_k)     :: V_old=0, ddV_old=0
    ! Counters and indices
    integer  :: i=1, j=1, l=1, m=1
    ! Timers
    real(dp) :: start_time, finish_time, dss_time, css_time, css_rf_time
! Auxiliary variables for Brent
	real(dp), parameter :: Tol_Brent = 1e-10
	real(dp) :: f_hat, k_opt
	integer  :: iter_brent

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Steady State calibration
	print*, " "
	write(*,*) "Calibration of Steady State - Objective: h_ss=", h_ss_obj
	call brent(chi,f_hat,iter_brent,0.01_dp,10.0_dp,Tol_Brent,FOC_SS)
	write(*,*) "Value of chi", chi, "Error is", f_hat, "Ierations", iter_brent
	
	! Steady State
	psi_1 = ( (1.0_dp/beta + delta - 1.0_dp) / (A*alpha) )**(1.0_dp/(1.0_dp-alpha))
    psi_2 = (1.0_dp-alpha) * A * psi_1**(-alpha) / ( A*psi_1**(1.0_dp-alpha)-delta )
	k_ss = psi_2/(chi+psi_1*psi_2)
	write(*,*) "k_ss=", k_ss


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Productivity grid and transition matrix
	call MarkovAR_95(n_z,rho_z,sigma_z,log_z_grid,P_z)
	z_grid = exp(log_z_grid)
		!write(*,*) "z_grid", z_grid
		!write(*,*) "P_z", P_z(1,:)

! Capital grid	
	k_min  = 0.001d0
	k_max  = 2*k_ss
	
	k_grid = grid(k_min,k_max,n_k,k_curve)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Solution for discrete state space
	print*, " "
	write(*,*) "Solution of NGM in discrete state space"
	call cpu_time(start_time)
	call NGM_DSS(V_dss,gk_dss,gh_dss,gc_dss)
	call cpu_time(finish_time)
	dss_time = (finish_time-start_time)/60.0_dp
	!V_dss = 0
    	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Solution for continuous state space (in capital and labor)
! Optimization method
	print*, " "
	write(*,*) "Solution of NGM in continuous state space"
	call cpu_time(start_time)
	call NGM_CSS(V_dss,V_css,gk_css,gh_css,gc_css)
	call cpu_time(finish_time)
	css_time = (finish_time-start_time)/60.0_dp

! Solution for continuous state space (in capital and labor)
! Root finding method
	print*, " "
	write(*,*) "Solution of NGM in continuous state space"
	call cpu_time(start_time)
	call NGM_CSS_RF(V_dss,V_css_rf,gk_css_rf,gh_css_rf,gc_css_rf)
	call cpu_time(finish_time)
	css_rf_time = (finish_time-start_time)/60.0_dp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Reporting Solutions

	print*, "k_grid"
	do i=1,n_k
		write(*,*) i, k_grid(i)
	end do
	print*, "==============================================="
	print*, "==============================================="

	print*, " "
	write(*,*) "Value Function discrete state space"
	write(*,*) "Capital", "Value by z"
	do i=1,n_k
		write(*,*) k_grid(i), V_dss(:,i)
	end do

	print*, " "
	write(*,*) "Capital policy function discrete state space"
	write(*,*) "Capital", "K' by z"
	do i=1,n_k
		write(*,*) k_grid(i), gk_dss(:,i)
	end do

	print*, "==============================================="
	print*, "==============================================="

	print*, " "
	write(*,*) "Value Function continuous state space"
	write(*,*) "Capital", "Value by z"
	do i=1,n_k
		write(*,*) k_grid(i), V_css(:,i)
	end do

	print*, " "
	write(*,*) "Capital policy function continuous state space"
	write(*,*) "Capital", "K' by z"
	do i=1,n_k
		write(*,*) k_grid(i), gk_css(:,i)
	end do

	print*, "==============================================="
	print*, "==============================================="

	print*, " "
	write(*,*) "Value Function continuous state space"
	write(*,*) "Capital", "Value by z"
	do i=1,n_k
		write(*,*) k_grid(i), V_css_rf(:,i)
	end do

	print*, " "
	write(*,*) "Capital policy function continuous state space"
	write(*,*) "Capital", "K' by z"
	do i=1,n_k
		write(*,*) k_grid(i), gk_css_rf(:,i)
	end do
	
	print*, "==============================================="
	print*, "==============================================="

	print*, " "
	write(*,*) "Comparisson of Solutions"
	write(*,*) "Capital", "Diff in Value by z"
	do i=1,n_k
		write(*,*) k_grid(i), V_css(:,i)-V_dss(:,i)
	end do

	print*, "==============================================="
	print*, "==============================================="

	print*, " "
	write(*,*) "Comparisson of Solutions 2"
	write(*,*) "Capital", "Diff in Value by z"
	do i=1,n_k
		write(*,*) k_grid(i), V_css(:,i)-V_css_rf(:,i)
	end do

	print*, "==============================================="
	print*, "==============================================="

	! Evaluate residual of euler equations
	V_old = V_dss
	do i=1,n_z
		ddV_old(i,:) = cspline_coeff(n_k,k_grid,V_old(i,:),2,0.0_dp,0.0_dp)  
	end do
	do i=1,n_z
	do j=1,n_k
		res_dss(i,j) = FOC_Euler(gk_dss(i,j))
	end do
	end do 

	V_old = V_css
	do i=1,n_z
		ddV_old(i,:) = cspline_coeff(n_k,k_grid,V_old(i,:),2,0.0_dp,0.0_dp)  
	end do
	do i=1,n_z
	do j=1,n_k
		res_css(i,j) = FOC_Euler(gk_css(i,j))
	end do
	end do 

	V_old = V_css_rf
	do i=1,n_z
		ddV_old(i,:) = cspline_coeff(n_k,k_grid,V_old(i,:),2,0.0_dp,0.0_dp)  
	end do
	do i=1,n_z
	do j=1,n_k
		res_css_rf(i,j) = FOC_Euler(gk_css_rf(i,j))
	end do
	end do 

	print*, "==============================================="
	print*, "==============================================="
	write(*,*) "Capital", "Residual Discrete"
	do i=1,n_k
		write(*,*) k_grid(i), res_dss(:,i)
	end do

	print*, "==============================================="
	print*, "==============================================="
	write(*,*) "Capital", "Residual Continuous"
	do i=1,n_k
		write(*,*) k_grid(i), res_css(:,i)
	end do

	print*, "==============================================="
	print*, "==============================================="
	write(*,*) "Capital", "Residual Continuous RF" 
	do i=1,n_k
		write(*,*) k_grid(i), res_css_rf(:,i)
	end do

	print*, "==============================================="
	print*, "==============================================="
	write(*,*) "Times ", "DSS", dss_time, "CSS", css_time, "CSS_RF", css_rf_time 



! 		! Input for FOC_Euler
! 	    	V_old = V_css
	    	
! 	    	! Coefficients for spline interpolation (this is used by neg_Bellman_T)
! 	    	! The last two elements are the derivatives at end nodes. They are needed only if method is 3
! 	    	do i=1,n_z
! 	    		ddV_old(i,:) = cspline_coeff(n_k,k_grid,V_old(i,:),2,0.0_dp,0.0_dp)  
! 	    	end do

! 	    	! Otpimization
! 	    	print*, " "
! 	    	print*, " "
! 	    	do i=1,n_z
! 	    		print*," "
! 	    		print*,"New Z"
! 	    	do j=1,n_k
! 	    		call brent(k_opt,f_hat,iter_brent,k_grid(1),2*k_grid(n_k),Tol_Brent,FOC_Euler)
! 	    		write(*,*) "Capital Diff", k_opt-gk_css(i,j), "Residual", FOC_Euler(k_opt), FOC_Euler(gk_css(i,j)), f_hat
! 	    	end do 
! 	    	end do 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Contains

! =============================================================================
! FOC_SS: Computes the diffence between the steady state of hours and 
!         the objective value h_ss for a given parameter chi=x
!
! Usage:  diff = FOC_SS(x)
!
! Input:  x   , real(dp), the value of parameter chi
!
! Output: diff  , real(dp), h_ss - Steady_State
!
! Note:   This function uses the parameters defined for the program
!
	Function  FOC_SS(x)
	    real(dp) :: x, FOC_SS
	    real(dp) :: k_ss, h_ss, psi_1, psi_2

	    psi_1 = ( (1.0_dp/beta + delta - 1.0_dp) / (A*alpha) )**(1.0_dp/(1.0_dp-alpha))

	    psi_2 = (1.0_dp-alpha) * A * psi_1**(-alpha) / ( A*psi_1**(1.0_dp-alpha)-delta )

	    k_ss   = psi_2/(x+psi_1*psi_2)

	    h_ss   = psi_1 * k_ss

	    FOC_SS = h_ss - h_ss_obj

	end Function FOC_SS

! =============================================================================
! NGM_DSS: Computes value and policy functions for the neoclassical growth model
!          with separable CRRA utility in consumption and labor and cobb-douglas
!          production function. Productivity shocks are a discrete Markov process
!          The solution is on a discrete state space, so both k' and h are chosen
!          from a grid.
!
! Usage:  call NGM_DSS(V_dss,gk_dss,gh_dss,gc_dss)
!
! Input:  There are no inputs but the routine uses all parameters defined at the top
!
! Output: V_dss,  real(dp), dim(n_z,n_k), Value function at grid (z,k)
!         gk_dss, real(dp), dim(n_z,n_k), Capital policy function at grid (z,k)
!         gh_dss, real(dp), dim(n_z,n_k), Hours policy function at grid (z,k)
!         gc_dss, real(dp), dim(n_z,n_k), Consumption policy function at grid (z,k)
!
! Note:   This function uses the parameters defined for the program
!
	Subroutine NGM_DSS(V,g_k,g_h,g_c)
		real(dp), dimension(n_z,n_k), intent(out) :: V, g_k, g_h, g_c
		real(dp), dimension(n_h)             :: h_grid=0 
		real(dp), dimension(n_z,n_k)     	 :: V_new=0
		real(dp), dimension(n_z,n_k,n_z) 	 :: P_mat
	    real(dp), dimension(n_z,n_k,n_k) 	 :: V_mat=0, U=0
	    integer , dimension(n_z,n_k)     	 :: g_k_ind=0
	    integer , dimension(n_z,n_k,n_k) 	 :: g_h_ind=0
	    real(dp), dimension(n_z,n_k,n_k,n_h) :: K_mat=0, Kp_mat=0, H_mat=0, Z_mat=0, U_h=0
	    real(dp) :: diff_V=10
	    integer  :: iter=1

	    ! Hours grid
			h_grid = grid(h_min,h_max,n_h,h_curve)
			!write(*,*) "h_grid", h_grid

		! Auxiliary matrices
			P_mat  = spread(transpose(P_z),2,n_k)
				! P_mat is of the form P(z',k',z)=Pr(z'|z)
				!write(*,*) "P_mat", P_mat(1,:,2)

			! Matrices of form (z,k,k',h')
			do i=1,n_k
				K_mat(:,i,:,:)  = k_grid(i)
				Kp_mat(:,:,i,:) = k_grid(i)
			end do

			do i=1,n_z
				Z_mat(i,:,:,:)  = z_grid(i)
			end do

			do i=1,n_h
				H_mat(:,:,:,i)  = h_grid(i)
			end do

		! Utility at each combination
			! U_h = (c^(1-sigma)-1)/(1-sigma) - chi*h^(1-nu)/(1-nu)
			! c   = max( Az(k^alpha)(h^(1-alpha)) + (1-delta)k - k' , 0 )
				U_h = log( max( A*Z_mat*(K_mat**alpha)*(H_mat**(1-alpha))+(1-delta)*K_mat-Kp_mat , 1e-20 ) )  &
						+ chi*log(1.0_dp-H_mat)
		! Utility maximization 
			! U is the utility function evaluated at optimum h for every triple (z,k,k')
			U = maxval(U_h,4)
				!print*, "Utiltiy for (z,k)"
				!do i=1,n_k
				!	write(*,*) U(n_z/2+1,i,:)
				!end do

		! VFI on discrete state space
			! Initialization
				V      = 0
				V_mat  = 0
				diff_V = 1
				iter   = 0
			! VFI
		    do while ((diff_V>Tol_V).and.(iter<Max_Iter))
		        V_new = maxval( U + beta*V_mat , 3 )
		        diff_V = maxval( abs(V-V_new) )
		        V = V_new
		        V_mat = spread( transpose(sum(P_mat*spread( V , 3 , n_z),1)) , 2 , n_k )
		            ! Interpret V as V(z',k') so that rows are z' (realization of z tomorrow) and columns k' (capital tomorrow)
		            ! spread( V , 3 , n_z) replicates V in a third dimension n_z times, each dimension is taken as a z (current z)
		            ! P_aux*spread( V , 3 , n_z),1) multiplies each "slice" of the replicated matrix by the transition prob of z
		            ! sum(P_aux*spread( V , 3 , n_z),1) is of shape (n_k,n_z) and takes expected value across z' (the first dim)
		            ! The last result has rows for k' and columns for z (current)
		            ! The result is transposed to get (z,k') and is replicated in the second dimension (moving k' to the third)
		            ! Is repicated in this dimension since the current value of k doens't affect value fun. tomorrow
		            ! Finally V_mat is of shape (n_z,n_k,n_k) corresponding to (z,k,k')
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
		    g_k_ind = maxloc( U + beta*V_mat , 3 ) 
		    g_h_ind = maxloc(U_h,4) 

			! Policy Functions in levels
		    do i=1,n_z
		    	g_k(i,:) = k_grid(g_k_ind(i,:))

		        do j=1,n_k
		            g_h(i,j) = h_grid(g_h_ind(i,j,g_k_ind(i,j)))
		            g_c(i,j) = A*z_grid(i)*(k_grid(j)**alpha)*(g_h(i,j)**(1-alpha))+(1-delta)*k_grid(j)-g_k(i,j)
		        end do
	    	end do

	end Subroutine NGM_DSS


! =============================================================================
! NGM_CSS: Computes value and policy functions for the neoclassical growth model
!          with separable CRRA utility in consumption and labor and cobb-douglas
!          production function. Productivity shocks are a discrete Markov process
!          The solution is on a continuous state space on capital and labor
!
! Usage:  call NGM_CSS(V_0,V_css,gk_css,gh_css,gc_css)
!
! Input:  V_0, real(dp), dim(n_z,n_k), Initial value function
!
! Output: V_css,  real(dp), dim(n_z,n_k), Value function at grid (z,k)
!         gk_css, real(dp), dim(n_z,n_k), Capital policy function at grid (z,k)
!         gh_css, real(dp), dim(n_z,n_k), Hours policy function at grid (z,k)
!         gc_css, real(dp), dim(n_z,n_k), Consumption policy function at grid (z,k)
!
! Note:   This function uses the parameters defined for the program. Interpolation 
!		  uses cubic splines with secant method for endpoint derivatives.
!
	Subroutine NGM_CSS(V_0,V,g_k,g_h,g_c)
		real(dp), dimension(n_z,n_k), intent(in)  :: V_0
		real(dp), dimension(n_z,n_k), intent(out) :: V, g_k, g_h, g_c
		real(dp), dimension(n_z,n_k) :: V_new=0
		real(dp) :: diff_V=10
	    integer  :: iter=0
	    ! Variables for brent
	    real(dp) :: k_min, k_a, k_b, k_c, h_opt, h_min, h_max
	    real(dp), dimension(3) :: p
	    ! Variables for interpolation
	    real(dp), dimension(n_z,n_k) :: ddV 

	    V     = V_0

	    do while ((diff_V>Tol_V).and.(iter<Max_Iter))
			! Input for neg_Bellman_T
	    	V_old = V
	    	
	    	! Coefficients for spline interpolation (this is used by neg_Bellman_T)
	    	! The last two elements are the derivatives at end nodes. They are needed only if method is 3
	    	do i=1,n_z
	    		ddV(i,:) = cspline_coeff(n_k,k_grid,V(i,:),2,0.0_dp,0.0_dp)  
	    	end do
	    		ddV_old = ddV

	    	! Otpimization
	    	do i=1,n_z
	    	do j=1,n_k

	    		call Min_Bracket(k_grid(1),k_grid(n_k/2),k_grid(n_k),neg_Bellman_T,k_a,k_b,k_c)
	    		V_new(i,j) = brent_opt(k_a,k_b,k_c,neg_Bellman_T,Tol_Brent,k_min)
		    	
		    end do 
			end do 

			! The value of V_new is multiplied by -1 since brent minimizes
		    V_new = -V_new 

		    ! Difference between value functions and updating
	    	diff_V = maxval( abs(V-V_new) )
		    V      = V_new

	    	iter = iter+1

	    	if (modulo(iter,50).eq.0) then
	    		write(*,*) "CSS: Iter=", iter, "Diff_V=", diff_V
	    	end if
	    end do

    	if (iter<Max_Iter) then
    		print*, " "
    		write(*,*) "CSS - Convergence was achieved at iteration:", iter
    	else
    		print*, " "
    		write(*,*) "CSS - Convergence failed - Max_Iter achieved"
    	end if 

    	! Policy Functions
	    	do i=1,n_z
	    	do j=1,n_k
	    		! Optimization for kp
	    		call Min_Bracket(k_grid(1),k_grid(n_k/2),k_grid(n_k),neg_Bellman_T,k_a,k_b,k_c)
	    		V_new(i,j) = brent_opt(k_a,k_b,k_c,neg_Bellman_T,Tol_Brent,k_min)

	    		! Optimization for h
				h_min = 0.00001_dp
				h_max = 1-h_min

				p(1) = z_grid(i)
				p(2) = k_grid(j)
				p(3) = k_min
					
				call brent_p(h_opt,f_hat,iter_brent,h_min,h_max,Tol_Brent,FOC_H,p,3)

				! Allocation of optimal policies
				g_k(i,j) = k_min
				g_h(i,j) = h_opt
		    end do 
			end do 


	end Subroutine NGM_CSS

! =============================================================================
! NGM_CSS_RF: Computes value and policy functions for the neoclassical growth model
!          with separable CRRA utility in consumption and labor and cobb-douglas
!          production function. Productivity shocks are a discrete Markov process
!          The solution is on a continuous state space on capital and labor
!		   Maximization uses a root finding method on the Euler equation
!
! Usage:  call NGM_CSS(V_0,V_css,gk_css,gh_css,gc_css)
!
! Input:  V_0, real(dp), dim(n_z,n_k), Initial value function
!
! Output: V_css,  real(dp), dim(n_z,n_k), Value function at grid (z,k)
!         gk_css, real(dp), dim(n_z,n_k), Capital policy function at grid (z,k)
!         gh_css, real(dp), dim(n_z,n_k), Hours policy function at grid (z,k)
!         gc_css, real(dp), dim(n_z,n_k), Consumption policy function at grid (z,k)
!
! Note:   This function uses the parameters defined for the program. Interpolation 
!		  uses cubic splines with secant method for endpoint derivatives.
!
	Subroutine NGM_CSS_RF(V_0,V,g_k,g_h,g_c)
		real(dp), dimension(n_z,n_k), intent(in)  :: V_0
		real(dp), dimension(n_z,n_k), intent(out) :: V, g_k, g_h, g_c
		real(dp), dimension(n_z,n_k) :: V_new=0
		real(dp) :: diff_V=10
	    integer  :: iter=0
	    ! Variables for brent
	    real(dp) :: k_opt, k_a, k_b, k_c, h_opt, h_min, h_max, res_min, res_max
	    real(dp), dimension(3) :: p
	    ! Variables for interpolation
	    real(dp), dimension(n_z,n_k) :: ddV 
	    real(dp), dimension(n_z)     :: V_tilde
		real(dp), dimension(1)       :: kp_aux, V_aux, dV_aux, ddV_aux

	    V     = V_0

	    do while ((diff_V>Tol_V).and.(iter<Max_Iter))
			! Input for FOC_Euler
	    	V_old = V
	    	
	    	! Coefficients for spline interpolation (this is used by neg_Bellman_T)
	    	! The last two elements are the derivatives at end nodes. They are needed only if method is 3
	    	do i=1,n_z
	    		ddV(i,:) = cspline_coeff(n_k,k_grid,V_old(i,:),2,0.0_dp,0.0_dp)  
	    	end do
	    		ddV_old = ddV

	    	! Otpimization
	    	do i=1,n_z
	    	do j=1,n_k
	    		! Optimization for kp
	    		res_min = FOC_Euler(k_grid(1))
	    		res_max = FOC_Euler(2*k_grid(n_k))
	    		!write(*,*) "res_min", res_min, "res_max", res_max
	    		if (res_min*res_max.lt.0.0_dp) then
	    			call brent(k_opt,f_hat,iter_brent,k_grid(1),2*k_grid(n_k),Tol_Brent,FOC_Euler)
	    		elseif ((res_min.gt.0.0_dp).and.(res_max.gt.0.0_dp)) then
	    			k_opt = k_grid(n_k)
	    		else
	    			k_opt = k_grid(1)
	    		end if
	    		!write(*,*) "Capital", k_opt, "Residual", FOC_Euler(k_opt), f_hat
	    		
	    		! Optimization for h
				h_min = 0.00001_dp
				h_max = 1-h_min

				p(1) = z_grid(i)
				p(2) = k_grid(j)
				p(3) = k_opt
	    		if (FOC_H(h_min,p)*FOC_H(h_max,p).lt.0.0_dp) then
					call brent_p(h_opt,f_hat,iter_brent,h_min,h_max,Tol_Brent,FOC_H,p,3)
				else
					if (FOC_H(h_max,p).gt.0.0_dp) then
						write(*,*) "Corner solution in hours", &
						         & "z=", z_grid(i), "k=", k_grid(j), "kp=", k_opt, "FOC", FOC_H(h_max,p)
						h_opt = h_max
					else !For some reason both bounds give negative FOC_H
						print *, "problem at p !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
					end if 
				end if 

				! Interpolation of V(k',z') for each value of z'
				kp_aux = k_opt
				do l=1,n_z
					call cspline_IP(n_k,k_grid,V_old(l,:),ddV_old(l,:),3,kp_aux,1,V_aux,dV_aux,ddV_aux)
					V_tilde(l) = V_aux(1)
				end do

	    		! Value Function
	    		V_new(i,j) = Utility(z_grid(i),k_grid(j),k_opt,h_opt) + beta*sum(P_z(i,:)*V_tilde)
		    	
		    end do 
			end do

		    ! Difference between value functions and updating
	    	diff_V = maxval( abs(V-V_new) )
		    V      = V_new

	    	iter = iter+1

	    	if (modulo(iter,50).eq.0) then
	    		write(*,*) "CSS: Iter=", iter, "Diff_V=", diff_V
	    	end if
	    end do

    	if (iter<Max_Iter) then
    		print*, " "
    		write(*,*) "CSS - Convergence was achieved at iteration:", iter
    	else
    		print*, " "
    		write(*,*) "CSS - Convergence failed - Max_Iter achieved"
    	end if 

    	! Policy Functions
	    	do i=1,n_z
	    	do j=1,n_k
	    		! Optimization for kp
	    		call brent(k_opt,f_hat,iter_brent,k_grid(1),2*k_grid(n_k),Tol_Brent,FOC_Euler)
	    		
	    		! Optimization for h
				h_min = 0.00001_dp
				h_max = 1-h_min

				p(1) = z_grid(i)
				p(2) = k_grid(j)
				p(3) = k_opt
	    		call brent_p(h_opt,f_hat,iter_brent,h_min,h_max,Tol_Brent,FOC_H,p,3)

				! Allocation of optimal policies
				g_k(i,j) = k_opt
				g_h(i,j) = h_opt
		    end do 
			end do 


	end Subroutine NGM_CSS_RF


! =============================================================================
! Utility: Computes value of utility at given values of states and controls
!          (states are z and k and controls are kp and h)
!          Utility is sepparable in consumption and labor and is of the form:
!          U(c,h) = log(c) - chi*log(h)
!
! Usage:  U = Utility(z,k,kp,h)
!
! Input:  z , real(dp), value of productivity z
!		  k , real(dp), value of current capital k
!		  kp, real(dp), value of future capital kp
!		  h , real(dp), value of labor h
!
! Output: Utility,  real(dp), Value of utility at (z,k,kp,h)
!
! Note:   If the consumption implied by (z,k,kp,h) is non-positive utilty is set
!         to a large negative number.
!		  This function uses the global parameters of the utiltiy function
!
	Function Utility(z,k,kp,h)
		real(dp), intent(in)  :: z, k, kp, h
		real(dp)              :: Utility
		real(dp)              :: c 

		!write(*,*) "kp in utility", kp

		c = A*z*(k**alpha)*(h**(1.0_dp-alpha))+(1.0_dp-delta)*k-kp

		if (c.ge.0.0_dp) then
			Utility = log(c) + chi*log(1.0_dp-h)
		else
			Utility =  -1e20
		end if

	end Function Utility

! =============================================================================
! FOC_H: Computes the residual of the labor FOC of the agent's problem
!        The function takes in a value for hours (h) and a vector p.
!        p contains the values of (z,k,kp) necessary to evaluate the FOC
!        The residual is  FOC_H = w/c - chi/(1-h)
!        Where: c=Az(k^alpha)(h^(1-alpha))+(1-delta)k-k' and w=(1-alpha)Az((k/h)^alpha)
!
! Usage:  Res = FOC_H(h,p)
!
! Input:  h, real(dp)              , value of labor h
!		  p, real(dp), dimension(3), vector with (z,k,kp) - in that order
!
! Output: FOC_H,  real(dp), Residual of the FOC of labor
!
! Note:   If the consumption implied by (z,k,kp,h) is non-positive the residual is set
!         to a large positive number. Indicating that a higher h is necessary.
!		  This function uses the global parameters of the utiltiy function
!
	Function FOC_H(h,p)
		real(dp), intent(in)   :: h
		real(dp), dimension(3), intent(in) :: p
		real(dp) :: FOC_H
		real(dp) :: c, w, k, z, kp

		z  = p(1)
		k  = p(2)
		kp = p(3)

		c = A*z*(k**alpha)*(h**(1.0_dp-alpha))+(1.0_dp-delta)*k-kp
		w = (1-alpha)*A*z*(k**alpha)*(h**(-alpha))

		if (c.gt.0.0_dp) then
			FOC_H = w*(1.0_dp-h) - chi*c
		else ! If consumption is negative or zero then the rewards for working are high
			FOC_H = 1e10
		end if

	end Function FOC_H

! =============================================================================
! FOC_Euler: Computes the residual of the Euler equation of the agent's problem
!        The function takes in a value for hours (kp).
!        The residual is  FOC_H = -1/c - beta*E[V'(k',z')|z]
!        Where: c=Az(k^alpha)(h^(1-alpha))+(1-delta)k-k' 
!
! Usage:  Res = FOC_Euler(kp)
!
! Input:  kp, real(dp)              , value of next period's capital
!
! Output: FOC_Euler,  real(dp), Residual of the Euler equation
!
! Note:   If the consumption implied by (z,k,kp,h) is non-positive the residual is set
!         to a large negative number. Indicating that a lower kp is necessary.
!		  This function uses the global parameters of the utiltiy function
!
	Function FOC_Euler(kp)
		real(dp), intent(in)     :: kp
		real(dp)                 :: FOC_Euler
		real(dp)                 :: c, k, z, h_opt, c_max, h_min, h_max
		real(dp), dimension(n_z) :: dV_tilde
		real(dp), dimension(1)   :: kp_aux, V_aux, dV_aux, ddV_aux
		real(dp), dimension(3)   :: p

		! Set current state
			z = z_grid(i)
			k = k_grid(j)

		! Solve the problem only if positive consumption is possible
			! Set limits for labor choice
			h_min = 0.00001_dp
			h_max = 1-h_min
			! Evaluate consumption if working full time
			c_max = A*z_grid(i)*(k_grid(j)**alpha)*(h_max**(1-alpha))+(1.0_dp-delta)*k_grid(j)-kp
			
		if (c_max.gt.0.0_dp) then
			!write(*,*) "c_max=", c_max
			! Optimal choicee of labor
			p(1) = z
			p(2) = k
			p(3) = kp
			
			if (FOC_H(h_min,p)*FOC_H(h_max,p).lt.0.0_dp) then
				call brent_p(h_opt,f_hat,iter_brent,h_min,h_max,Tol_Brent,FOC_H,p,3)
			else
				if (FOC_H(h_max,p).gt.0.0_dp) then
					write(*,*) "Corner solution in hours", "z=", z, "k=", k, "kp=", kp, "FOC", FOC_H(h_max,p)
					h_opt = h_max
				else !For some reason both bounds give negative FOC_H
					print *, "problem at p !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
					!write(*,*) "Problem at p=", p
					!print *, A*z_grid(i)*(k_grid(j)**alpha)*(0.001**(1-alpha))+(1.0_dp-delta)*k_grid(j)-kp
					!print *, A*z_grid(i)*(k_grid(j)**alpha)*(0.999**(1-alpha))+(1.0_dp-delta)*k_grid(j)-kp
					!print *, A*z_grid(i)*(k_grid(j)**alpha)*(1.0_dp**(1-alpha))+(1.0_dp-delta)*k_grid(j)-kp
					!print *, (1-alpha)*A*p(1)*(p(2)**alpha)*(0.999**(-alpha))/&
					!         & (A*p(1)*(p(2)**alpha)*(0.999**(1-alpha))+(1.0_dp-delta)*p(2)-kp) - chi*1/(1-0.999)
				end if 
			end if 

			! Set current conumption
			c = A*z*(k**alpha)*(h_opt**(1.0_dp-alpha))+(1.0_dp-delta)*k-kp

			! Interpolation of V(k',z') for each value of z'
			kp_aux = kp
			do l=1,n_z
				if ((kp.le.k_grid(n_k)).and.(kp.ge.k_grid(1))) then
					call cspline_IP(n_k,k_grid,V_old(l,:),ddV_old(l,:),3,kp_aux,1,V_aux,dV_aux,ddV_aux)
				else
					if (kp.gt.k_grid(n_k)) then
						call cspline_IP(n_k,k_grid,V_old(l,:),ddV_old(l,:),3,k_grid(n_k),1,V_aux,dV_aux,ddV_aux)
					else
						call cspline_IP(n_k,k_grid,V_old(l,:),ddV_old(l,:),3,k_grid(1),1,V_aux,dV_aux,ddV_aux)
					end if
				end if  
				dV_tilde(l) = dV_aux(1)
			end do
			!print*,"dV"
			!print*,dV_aux

			! Set Euler equation
			FOC_Euler = -1 + c*beta*sum( P_z(i,:)*dV_tilde ) 


		else ! If consumption is negative or zero then the cost of saving is too high
			FOC_Euler = -1E20
		end if

	end Function FOC_Euler

! =============================================================================
! neg_Bellman_T: Computes the value of the negative of the bellman operator T
!                evaluated for a given candidate value k'.
!                The function takes as given V_old, and the current state (z,k)
!				 In order to evaluate T at k' the optimal labor decision is computed
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
!                  V_old  , real(dp), dimension(n_z,n_k), Current value function
!                  ddV_old, real(dp), dimension(n_z,n_k), Second derivatives of current value function
!
! Output: neg_Bellman_T,  real(dp), Negatie of the bellman operator evaluated at kp
!
! Note:  If the maximum consumption implied by (z,k,kp) is non-positive the candidate kp is 
!        discarded and the value of the operator is set to a large positive number. 
!        When optimizing for labor it is necessary to bracket a solution of FOC_H=0
!        If this is not achieved by the specified bounds then a corner solution is assumed.
!		 This function uses the global parameters of the utiltiy function
!        This function is intended to be used with brent_opt in the VFI procedure since 
!        brent_opt is a minimization routine this function returns the negative of the operator. 
!		 Since brent_opt does not allow the function to have more than one argument the 
!        values of (z,k,V_old,ddV_old) are passed through global variables.
!
	Function neg_Bellman_T(kp)
		real(dp), intent(in)         :: kp
		real(dp)                     :: neg_Bellman_T, Bellman_T, h_opt, c_max, h_min, h_max
		real(dp), dimension(n_z)     :: V_tilde
		real(dp), dimension(1)       :: kp_aux, V_aux, dV_aux, ddV_aux
		real(dp), dimension(3)       :: p

		! Solve the problem only if positive consumption is possible
			! Set limits for labor choice
			h_min = 0.00001_dp
			h_max = 1-h_min
			! Evaluate consumption if working full time
			c_max = A*z_grid(i)*(k_grid(j)**alpha)*(h_max**(1-alpha))+(1.0_dp-delta)*k_grid(j)-kp
			
		if (c_max.gt.0.0_dp) then
			!write(*,*) "c_max=", c_max
			! Optimal choicee of labor
			p(1) = z_grid(i)
			p(2) = k_grid(j)
			p(3) = kp
			
			if (FOC_H(h_min,p)*FOC_H(h_max,p).lt.0.0_dp) then
				call brent_p(h_opt,f_hat,iter_brent,h_min,h_max,Tol_Brent,FOC_H,p,3)
			elseif (FOC_H(h_max,p).gt.0.0_dp) then
				write(*,*) "Corner solution in hours choice"
				h_opt = h_max
			else !For some reason both bounds give negative FOC_H
				write(*,*) "Problem at p=", p
				print *, A*z_grid(i)*(k_grid(j)**alpha)*(0.001**(1-alpha))+(1.0_dp-delta)*k_grid(j)-kp
				print *, A*z_grid(i)*(k_grid(j)**alpha)*(0.999**(1-alpha))+(1.0_dp-delta)*k_grid(j)-kp
				print *, A*z_grid(i)*(k_grid(j)**alpha)*(1.0_dp**(1-alpha))+(1.0_dp-delta)*k_grid(j)-kp
				print *, (1-alpha)*A*p(1)*(p(2)**alpha)*(0.999**(-alpha))/&
				         & (A*p(1)*(p(2)**alpha)*(0.999**(1-alpha))+(1.0_dp-delta)*p(2)-kp) - chi*1/(1-0.999)
			end if 

			
			! Interpolation of V(k',z') for each value of z'
			kp_aux = kp
			do l=1,n_z
				call cspline_IP(n_k,k_grid,V_old(l,:),ddV_old(l,:),3,kp_aux,1,V_aux,dV_aux,ddV_aux)
				V_tilde(l) = V_aux(1)
			end do

			!write(*,*) "V(kp)=", V_tilde

			! Bellman operator on current function V
			Bellman_T = Utility(z_grid(i),k_grid(j),kp,h_opt) + beta*sum(P_z(i,:)*V_tilde)
			neg_Bellman_T = -Bellman_T
		else
			neg_Bellman_T = 1e10
		end if

	end Function neg_Bellman_T


end Program NGM_Main_Fortran


