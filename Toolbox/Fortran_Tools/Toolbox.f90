module Toolbox
	implicit none
	! Options for size and precision in the program
    integer, parameter          ::  DP      = KIND(1.0D0)
    integer, parameter          ::  QP      = selected_real_kind(32)
    integer, parameter          ::  LG      = KIND(.true.)
    real(dp), parameter         ::  eps_p   = EPSILON(1.0_dp)
    real(dp), parameter         ::  small_p = 1E-7_dp
    real(dp), parameter         ::  big_p   = HUGE(1.0_dp)

Contains
    ! diff_f:  Computes (centered) first derivative of real valued function f:R->R
    ! CDJac:   Compute Jacobian of f using central differences
    ! CDHesse: Compute approximation of the Hesse matrix of a real valued function 
    
    ! Inv:    Compute inverse of square matrix of real(8) elements
    ! Sort:   Sorts in ascending the elements of a one dimensional array of type real(8) 
    ! cumsum: Gives a matrix with the cumulative sum of an original matrix along a desired dimension
    ! Vec:    The vec-operator of matrix A(n,m)
    ! Kron:   Computes the kronecker product of matrices A (n_a by m_a) and B (n_b by m_b)
    ! StDev:  Computes the standard deviation of the elemtens in the row vector x(n)

    ! NCDF:      Evaluates the normal CDF in value x
    ! randn:     Compute pseudo random numbers from a standard normal distribution.
    ! randn_mat: Computes an n by m matrix of standard normal random numbers

    ! MarkovAR:    Approximates a contiuous AR(1) process with a discrete Markov process - Tauchen    86
    ! MarkovAR_95: Approximates a contiuous AR(1) process with a discrete Markov process - Rowenhorst 95

    ! bracket:     Brackets the zero of a continous function f(x) between x_min and x_max
    ! bisection:   Computes x_hat in [a,b] such that abs(f(x_hat))<=Tol
    ! Secant:      Computes x_hat in [a,b] such that abs(f(x_hat))<=Tol using secant method
    ! Newton:      Computes x_hat in [a,b] such that abs(f(x_hat))<=Tol using Newton's method
    ! Brent:       Computes x_hat in [a,b] such that abs(f(x_hat))<=Tol Brent's method

    ! grid:          Creates grid between x_min and x_max with curvature c_grid
    ! Linear_IP:     Conducts Linear Interpolation to compute f(x) given a grid
    ! cspline_coeff: Computes second derivates for cubic splines. Does not compute interpolation.
    ! cspline_IP:    Computes cubic spline approximation for points x_eval. Uses dff from cspline_coeff.

    ! Min_Bracket: Brackets a minimum of a continous function f(x) between a_0 and c_0
    ! GSS:         Golden section search to find the maximizer of y=f(x)
    ! brent_opt:   Computes the minimum and minimizer of a real valued function using Brent's method.
    ! dbrent_opt:  Computes the minimum and minimizer of a real valued function using D-Brent's method.
    ! Nelder_Mead: Computes the minimum and minimizer of a real valued function using Nelder-Mead method.



! =============================================================================
!
! diff_f: Computes (centered) first derivative of real valued function f:R->R
!
! Usage: df = diff_f(x,f)
!
! Input: x, real(dp), Point at which derivative is taken
!        f,           Pointer to the function f(x)
!
! Output: df, real(dp), First derivative 
!
! Remarks: Derivatives computed with central differences

    Function diff_f(x,f)
        implicit none
        real(dp) :: x, diff_f
        real(dp), external :: f
        real(dp) :: eps, eps1, h, fx1, fx2, x1, x2

        eps1 = Epsilon(eps1)
        eps  = Epsilon(eps)**(1.0_dp/3.0_dp)
          
        h    = eps*dabs(x)
        if (h .le. eps1) h=eps
        x1 = x + h
        x2 = x - h
        h  = x1 - x

        fx1 = f(x1)
        fx2 = f(x2)

        diff_f = (fx1-fx2)/(2.0_dp*h)
        return

    End Function diff_f  



! ***********************************************************************************
!
! CDJac: Compute Jacobian of f using central differences
!
! Usage: Call CDJac(m,n,x,f,df)
!
! Input: m, integer(4), the number of elements in x
!        n, integer(4), the number of elements returned by f(m,n,x,fx) in fx
!        x, real(8), m times 1 vector
!        f, pointer to the function f(m,n,x,fx)
!
! Output: df, real(8) n times n matrix
!
! Remarks: This is algorithm A.5.6.4. from Dennis and Schnabel (1983), p. 323

    Subroutine CDJac(m,n,x,f,df)
		
		Integer(4) n, m
		Real(8) x(m), df(n,m)

	    Integer(4) i, j
		Real(8) eps, eps1, h, temp, fx1(n), fx2(n), x1(m), x2(m)
		External f

	    x1=x      
	    x2=x

	    eps1 = Epsilon(eps1)
		eps  = Epsilon(eps)**(1.0/3.0)
	      
	    do j=1,m
		   temp=x1(j)
		   h=eps*dabs(temp)
		   if (h .le. eps1) h=eps
		   x1(j)=temp+h
		   x2(j)=temp-h
		   h=x1(j)-temp
		   call f(m,n,x1,fx1)

 		   ! if (FError.ne.0) return
		   call F(m,n,x2,fx2)     
 		   ! if (FError.ne.0) return    
		   do i=1,n
		      df(i,j)=(fx1(i)-fx2(i))/(2.*h)
		   end do
		   x1(j)=x(j)
		   x2(j)=x(j)
		end do

		return

	End Subroutine CDJac

!******************************************************************************************
!
! CDHesse: compute approximation of the Hesse matrix of a real valued function
!          fx=f(n,x), where x is a n vector
!
! Usage:  Call CDHesse(n,x0,f,hmat)
!
! input:  n: Integer(4), the number of elements of x0
!        x0: Real(8),    the point at which the derivatives are to be computed
!         f: pointer to the procedure that returns the value of f at x in fx
!
! output: hmat: Real(8), n by n matrix
!
! Remarks:  algorithm: based on the formulas (A.2.10) and (A.2.11) in Heer and Maussner
!
!           f must be declared as external in the calling program
!
	Subroutine CDHesse(m,x0,f,hmat)
		implicit none 
		Integer(4) n, m
		Real(8) x0(1:m), hmat(1:m,1:m)
		Integer(4) i, j
		Real(8) eps, temp, x1(1:m), x2(1:m), x3(1:m), x4(1:m), h(1:m)
		real(8) f0, f1, f2, f3, f4
		external f
	 
	 	n = 1 
		eps=Epsilon(eps)**(1.0/3.0)
	 
		Do i=1,m,1
	    
			if (x0(i).lt.0.0d0) then
				h(i)=-dmax1(dabs(x0(i)),1.d0)*eps
			else
				h(i)=dmax1(dabs(x0(i)),1.0)*eps
			endif
		
			temp=x0(i)+h(i)
			h(i)=temp-x0(i) ! the last two steps increase precision slightly, see Dennis and Schnabel (1983), p. 321 for this trick 
		End Do

	 	! print *, "Vector h for Hessian matrix"
	 	! print *, h

		Do i=1,m,1
			x1=x0
			x2=x0
			x1(i)=x1(i)+h(i)
			x2(i)=x2(i)-h(i)
				call f(m,n,x0,f0)
				call f(m,n,x1,f1)
				call f(m,n,x2,f2)
	 			! print *, "Function evaluation 1"
	 			! print *, (/f0 , f1, f2/)
			hmat(i,i)=(f1+f2-2*f0)/((h(i)**2.d0))
			Do j=i+1,m,1    
				x1=x0
				x1(i)=x1(i)+h(i)
				x1(j)=x1(j)+h(j)
	              x2=x0
				x2(i)=x2(i)-h(i)
				x2(j)=x2(j)+h(j)
	              x3=x0
				x3(i)=x3(i)+h(i)
				x3(j)=x3(j)-h(j)
	              x4=x0
				x4(i)=x4(i)-h(i)
				x4(j)=x4(j)-h(j)
					call f(m,n,x1,f1)
					call f(m,n,x2,f2)
					call f(m,n,x3,f3)
					call f(m,n,x4,f4)
	 				! print *, "Function evaluation 2"
	 				! print *, (/f0 , f1, f2, f3, f4/)
				hmat(i,j)=(f1-f2-f3+f4)/(4*h(i)*h(j))
				hmat(j,i)=hmat(i,j)           
			End Do
		End Do

		Return
	
	End Subroutine CDHesse

!******************************************************************************************
!
! Inv: Compute inverse of square matrix of real(8) elements
!
! Usage: Inv_A = Inv(A)
!
! Input: A   , real(8), dimension(n,ntor        
!
! Output: A_inv, real(8), inverse of A
! 
! Note: This function uses LAPack
!

! 	Function inv(A) result(Ainv)
! 	  real(dp), dimension(:,:), intent(in) :: A
! 	  real(dp), dimension(size(A,1),size(A,2)) :: Ainv

! 	  real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
! 	  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
! 	  integer :: n, info

! 	  ! External procedures defined in LAPACK
! 	  external DGETRF
! 	  external DGETRI

! 	  ! Store A in Ainv to prevent it from being overwritten by LAPACK
! 	  Ainv = A
! 	  n = size(A,1)

! 	  ! DGETRF computes an LU factorization of a general M-by-N matrix A
! 	  ! using partial pivoting with row interchanges.
! 	  call DGETRF(n, n, Ainv, n, ipiv, info)

! 	  if (info /= 0) then
! 	     stop 'Matrix is numerically singular!'
! 	  end if

! 	  ! DGETRI computes the inverse of a matrix using the LU factorization
! 	  ! computed by DGETRF.
! 	  call DGETRI(n, Ainv, n, ipiv, work, n, info)

! 	  if (info /= 0) then
! 	     stop 'Matrix inversion failed!'
! 	  end if
! 	End function inv

!******************************************************************************************
!
! Sort: Sorts in ascending the elements of a one dimensional array of type real(8) 
!       It also gives the original indeces of the sorted elements
!
! Usage: Call Sort(n,A,A_sort,Sort_ind)
!
! Input: n   , integer(4), the number of elements in A
!        A   , real(8), the array whose elements are to be sorted
!
! Output: A_sort, real(8), array with sorted elements (in ascending order)
!		  Sort_ind, integer(4), array with original indeces of the elements in A_sort
!

	Subroutine Sort(n,A,A_sort,Sort_ind)
		integer , intent(in) :: n    !Number of elements in A
		real(dp), intent(in) , dimension(n) :: A
		real(dp), intent(out), dimension(n) :: A_sort
		integer , intent(out), dimension(n) :: Sort_ind
		integer :: i,j

		A_sort = A
		do i=1,n
			Sort_ind(i)=i
		end do

		do i=1,(n-1)
			do j=i+1,n
				if (A_sort(i) .ge. A_sort(j)) then
					A_sort((/i,j/))   = A_sort((/j,i/))
					Sort_ind((/i,j/)) = Sort_ind((/j,i/))
				end if
			end do
		end do

		return
	End Subroutine Sort	

! =============================================================================
! cumsum: Gives a matrix with the cumulative sum of an original matrix along a desired dimension
!
! Usage: Call cumsum(A,n,m,dim,B)
!
! Input: A   , real(dp), dimension(n,m), The original matrix
!        n   , integer ,                 Number of rows
!        m   , integer ,                 Number of columns
!        dim , integer ,                 Dimension along which to sum - 1 sums rows - 2 sums cols
!
! Output: B  , real(dp), dimension(n,m), Matrix of size(A) with the cumulative sums along dim
!

	Subroutine cumsum(A,n,m,dim,B)
	    integer, intent(in)  :: n,m,dim
	    real(dp), intent(in)  :: A(n,m)
	    real(dp), intent(out) :: B(n,m)
	    integer :: i

	    if (n .ne. size(A(:,1))) write(*,*)  "Wrong input n not equal to rows of A", "Rows=", size(A(:,1))
	    if (m .ne. size(A(1,:))) write(*,*)  "Wrong input m not equal to columns of A", "Cols=", size(A(1,:))

	    if (dim.eq.1) then
	        B(1,:) = A(1,:)
	        do i=2,n
	            B(i,:) = B(i-1,:)+A(i,:)
	        end do
	    elseif (dim.eq.2) then
	        B(:,1) = A(:,1)
	        do i=2,m
	            B(:,i) = B(:,i-1)+A(:,i)
	        end do
	    else
	        write(*,*) "Wrong input dim not equat to either 1 or 2", "dim=", dim
	    end if

	end Subroutine cumsum


!************************************************************************************
!
! Vec: the vec-operator
!
! Usage Call(n,m,amat,bvec)
!
! Input: n: Integer(4) the number of rows of Amat
!        m: Integer(4) the number of columns of Amat
!     Amat: Real(8) n times m matrix
!
! Output: bvec: Real(8) n*m vector
!
	
	Function Vec(n,m,amat) result(bvec)
		implicit none
		Integer  :: n, m
		Real(dp) :: amat(n,m), bvec(n*m)
		Integer(4) i,j

		j=1
		Do i=1,m
			bvec(j:i*n)=amat(1:n,i)
			j=j+n
		End do
		Return
	End Function Vec

!******************************************************************************************
!
! Kron: Computes the kronecker product of matrices A (n_a by m_a) and B (n_b by m_b)
!
! Usage: x = Kron(n_a,m_a,A,n_b,m_b,B)
!
! input:  n_a: Integer(4), the number of rows of matrix A
!		  m_a: Integer(4), the number of columns of matrix A
! 		  A  : real(8), dimension(n_a,m_a), matrix A
! 		  n_b: Integer(4), the number of rows of matrix B
!		  m_b: Integer(4), the number of columns of matrix B
! 		  B  : real(8), dimension(n_b,m_b), matrix B
!
! output: K_AB: Real(8), n_a*n_b by m_a*m_b matrix
!
!

	Function Kron(n_a,m_a,A,n_b,m_b,B) result(K_AB)
		integer, intent(in) :: n_a,m_a,n_b,m_b
		real(8), intent(in) :: A(n_a,m_a), B(n_b,m_b)
		real(8), dimension(n_a*n_b,m_a*m_b) :: K_AB
		integer :: i, j

		do i=1,n_a
			do j=1,m_a
				K_AB( (n_b*(i-1)+1):(n_b*i) , (m_b*(j-1)+1):(m_b*j) ) = A(i,j)*B
			end do 
		end do 

		return
	end function Kron


!*********************************************************************************
! StDev: Computes the standard deviation of the elemtens in the row vector x(n)
!
! Usage: Result=StDev(n,x)
!
! Input: n, Integer(4), the number of elements in x
!        x, Real(8), n times 1 vector
!
! Output: Result, Real(8)

    Function StDev(n,x)
    	implicit none
		INTEGER  :: n, i
		REAL(dp) :: x(n), xmean, xsum, StDev
      
		xsum  = SUM(x)
		xmean = dble(Sum(x)/n)
		
		StDev = 0.0
		Do i=1,n,1
	        StDev = StDev + ((x(i)-xmean)**2)
	    End Do
		
		StDev = DSQRT(StDev/(n-1))
	    Return
	END Function StDev


! =============================================================================
! NCDF: Evaluates the normal CDF in value x
!
! Usage: phi = NCDF(x)
!
! Input: x   , real(dp), value at which to evaluate
!
! Output: NCDF , real(dp), Value of normal CDF at x
!

	real(dp) function NCDF(x)
		!     real(dp), parameter  :: pi=3.141592653589793
		real(dp), intent(in) :: x

		NCDF = ( 1+erf( x/sqrt(real(2,8)) ))/2

		return
	end function NCDF


!******************************************************************************************
!
! randn: Compute pseudo random numbers from a standard normal distribution.
!
! Usage:  x = randn(n)
!
! input:  n: Integer(4), either 1 or 2 is the number of draws to be generated
!                        The code generates 2 random numbers, either both or one is reported
!
! output: rand_n: Real(8), n array, contains the random number (or numbers) generated
!
! Remarks: The code follows Numerical Recepies (Section 7.2) but changes their uniform random 
!          number generator (rand1) for the built in routine random_number.
!          To generate the normal random numbers two uniform random numbers [-1,1] are needed.
!		   The uniform random numbers must be inside the unit circle. randn generates draws until
!		   the numbers are inside the unit circle.
! 		   Once the numbers are obtained the Box-Muller transforom is used.
!

	function randn(n) result(rand_n)
		integer, intent(in)   :: n
		real(8) :: rand_n(n), fac,rsq,v(2)

		if (n.ne.1 .and. n.ne.2) write(*,*) "Warning! n has to be either 1 or 2, currently n=", n

		rsq=0  ! This is set so that the following while statement works
		do while (rsq.ge.1..or.rsq.eq.0.)
			call random_number(v)
			v=2.*v-1
			rsq=v(1)**2+v(2)**2 
		end do
		fac=sqrt(-2.*log(rsq)/rsq) 			
		rand_n(1) = v(1)*fac
		if (n.eq.2) rand_n(2) = v(2)*fac
		return
	end function randn


!******************************************************************************************
!
! randn_mat: Computes an n by m matrix of standard normal random numbers
!
! Usage: x = randn_mat(n,m)
!
! input:  n: Integer(4), the number of rows
!		  m: Integer(4), the number of columns
!
! output: rand_mat: Real(8), n by m matrix
!
! Remarks:  This function calls randn to generate normal random numbers.
!  	    	Since normal random numbers come in pairs the matrix is filled by pairs of numbers
!

	function randn_mat(n,m) result(rand_mat)
		integer, intent(in)  :: n,m
		real(8), dimension(n,m) :: rand_mat
		real(8), dimension(1) :: temp
		integer :: i, j
		real(8) :: fac,rsq,v(2)

		do i=1,n
			! This do covers all numbers from 1 to m if m is even and all from 1 to m-1 if m is odd
			do j=1,m-1,2
				rand_mat(i,j:j+1) = randn(2)
			end do			
			if (mod(m,2).eq.1) then
				temp = randn(1)
				rand_mat(i,m) = temp(1)
			end if
		end do		
		return
	end function randn_mat



! Markov Process Subroutine


! =============================================================================
! MarkovAR: Approximates a contiuous AR(1) process with a discrete Markov process - Tauchen 86
!
! Usage: Call MarkovAR(n_z,lambda,rho,sigma,z_grid,P)
!
! Input: n_z    , integer , Grid dimension for Markov process
!        lambda , real(dp), Expansion of the grid: [-lambda*std(x) , lambda*std(x)] where x is AR(1)
!        rho    , real(dp), Persistence of x where x'=rho*x+sigma*e and e~N(0,1)
!        sigma  , real(dp), Standard Deviation of shocks to x where x'=rho*x+sigma*e and e~N(0,1)
!
! Output: z_grid , real(dp), dimension(n_z),     Grid of values of the Markov process. Size n_z.
!         P      , real(dp), dimension(n_z,n_z), Transition Matrix of z. Prob sum across cols.
!
! Remarks: This routine generates a grid for z (z_grid) and a transition probability matrix (P)
!          The transition probability matrix is organized (z',z) so that it sums to 1 by cols
! 		   P(j,i) = p_ij is the transition prob of going from z=z_i to z'=z_j
!          This method is taken from Tauchen (1986)

	subroutine MarkovAR(n_z,lambda,rho,sigma,z_grid,P)
	    integer, intent(in) :: n_z
	    real(8), intent(in) :: lambda, rho, sigma
	    real(8), intent(inout), dimension(n_z)     :: z_grid
	    real(8), intent(inout), dimension(n_z,n_z) :: P
	    integer :: i, j
	    real(8) :: step

	    ! Step of grid
	    step = 2*(lambda*sigma/sqrt(1-rho**2))/(n_z-1)
	    
	    ! Compute z_grid
	    z_grid(1) = -1*(lambda*sigma/sqrt(1-rho**2))
	    do i=2,n_z
	        z_grid(i) = z_grid(i-1) + step
	    end do
	    
	    ! Compute transition matrix    
	    do i=1,n_z
	        P(1,i) = NCDF( (z_grid(1)-rho*z_grid(i)+step/2)/sigma )
	        do j=2,n_z-1
	            P(j,i) = NCDF( (z_grid(j)-rho*z_grid(i)+step/2)/sigma ) - NCDF( (z_grid(j)-rho*z_grid(i)-step/2)/sigma )
	        end do
	        P(n_z,i)=1-sum(P(1:n_z-1,i))
	    end do

	    return
	end subroutine  MarkovAR


! =============================================================================
! MarkovAR_95: Approximates a contiuous AR(1) process with a discrete Markov process - Rowenhorst 95
!
! Usage: Call MarkovAR_95(n_z,rho,sigma,z_grid,P)
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
! 		   P(i,j) = p_ij is the transition prob of going from z=z_i to z'=z_j
!          Note that this creates matrixes in the opposite order than MarkovAR that follows Tauchen
!          The method follows Rouwenhorst (1995) as shown in Kopecky and Suen (2010)
	
	Subroutine MarkovAR_95(n_z,rho,sigma,z_grid,P_z)
	    integer, intent(in) :: n_z
	    real(dp), intent(in) :: rho, sigma
	    real(dp), intent(inout), dimension(n_z)     :: z_grid
	    real(dp), intent(inout), dimension(n_z,n_z) :: P_z
	    integer :: i, j
	    real(dp) :: step, p, q, psi, P_2(2,2)
	    real(dp) , dimension(:,:), allocatable :: zeros, P_a1, P_a2, P_a3, P_a4, P_n, P_o

	    ! Parameters p, q and psi
	        p = (1+rho)/2
	        q = (1+rho)/2
	        psi = sqrt(real(n_z)-1)*sigma/sqrt(1-rho**2)

	    ! Step of grid
	        step = 2*psi/(n_z-1)
	    
	    ! Compute z_grid
	        z_grid(1) = -psi
	        do i=2,n_z
	            z_grid(i) = z_grid(i-1) + step
	        end do
	    
	    ! Compute transition matrix for n_z=2
	        P_2 = reshape((/p,1-q,1-p,q/),(/2,2/))

	    ! Compute transition matrix for arbitrary n_z 
	        if (n_z>2) then
	            
	            allocate(P_o(2,2),P_n(1,1),P_a1(1,1),P_a2(1,1),P_a3(1,1),P_a4(1,1))
	            P_o = P_2

	            do i=3,n_z
	                deallocate(P_n,P_a1,P_a2,P_a3,P_a4)
	                allocate(P_n(i,i),P_a1(i,i),P_a2(i,i),P_a3(i,i),P_a4(i,i))
	                
	                P_a1 = 0
	                P_a1(1:i-1,1:i-1) = P_o
	                P_a2 = 0
	                P_a2(1:i-1,2:i) = P_o
	                P_a3 = 0
	                P_a3(2:i,1:i-1) = P_o
	                P_a4 = 0
	                P_a4(2:i,2:i) = P_o

	                P_n = p*P_a1+(1-p)*P_a2+(1-q)*P_a3+q*P_a4

	                ! The transition matrix is multiplied by 1/2 in all rows but the first and last
	                ! This is done so that P_z sums to 1 by rows.
	                P_a1(1,:) = 1.0d0
	                P_a1(i,:) = 1.0d0
	                P_a1(2:i-1,:) = 0.50d0
	                P_n = P_a1*P_n

	                deallocate(P_o)
	                allocate(P_o(i,i))
	                P_o = P_n

	            end do 

	            P_z = P_n
	        else
	            P_z = P_2
	        end if

	    return
	end subroutine  MarkovAR_95



! =============================================================================
! bracket: Brackets the zero of a continous function f(x) between x_min and x_max
!          bracket searchs for the points in a grid with n_grid points and curvature
!          c_grid.
!
! Usage: Call bracket(a,b,x_min,x_max,n_grid,c_grid)
!
! Input: x_min   , real(dp), the lower number to consider in the grid to search for zero
!        x_max   , real(dp), the higher number to consider in the grid to search for zero
!        n_grid  , integer , the number of grid points between x_min and x_max
!        c_grid  , real(dp), the curvature of the grid (how close are points near x_min)
!        fun     , Function
!
! Output: a , real(dp), closest grid point to the left of the zero
!         b , real(dp), closest grid point to the right of the zero
!

	Subroutine bracket(a,b,x_min,x_max,n_grid,c_grid,fun)
	    implicit none
	    integer               :: i, loc_aux(1)
	    integer , intent(in)  :: n_grid
	    real(dp), intent(in)  :: x_min, x_max, c_grid
	    real(dp), intent(out) :: a, b
	    real(dp)              :: grid_v(n_grid), fgrid(n_grid), x, y, z
	    real(dp), external    :: fun


	    ! Define the Grid
	    grid_v=0.0_dp
	    do i=1,n_grid
	        grid_v(i) = real(i-1,dp)/real(n_grid-1,dp)
	    end do

	    grid_v=grid_v**c_grid
	    grid_v=grid_v*(x_max-x_min)+x_min

	    ! Evaluate the function and save the distance between f(x) and zero
	    do i = 1,n_grid
	        fgrid(i) = abs(fun(grid_v(i))-0.0_dp)
	    end do

	    loc_aux = minloc(fgrid) ! Save the location of the element closest to zero
	    IF (loc_aux(1).eq.1) THEN
	        x = grid_v(loc_aux(1))
	        y = grid_v(loc_aux(1))
	        z = grid_v(loc_aux(1) + 1)
	    ELSEIF (loc_aux(1).eq.n_grid) THEN
	        x = grid_v(loc_aux(1) - 1)
	        y = grid_v(loc_aux(1))
	        z = grid_v(loc_aux(1))
	    ELSE
	        x = grid_v(loc_aux(1) - 1)
	        y = grid_v(loc_aux(1))
	        z = grid_v(loc_aux(1) + 1)
	    END IF

	    ! Checking if the points bracket the zero and assigning brackets
	    IF ( fun(x)*fun(z).gt.0.0_dp ) THEN
	        PRINT*, 'ssbrack did not bracket.'
	        RETURN
	    ELSE
	        IF ( fun(x)*fun(y).lt.0.0_dp ) THEN
	            a = x
	            b = y
	        ELSE
	            a = y
	            b = z
	        END IF
	    END IF

	end Subroutine bracket


! =============================================================================
! bisection: Computes x_hat in [a,b] such that abs(f(x_hat))<=Tol
!            The routine also gives f_hat and iter, the number of iterations.
!
! Usage: Call bisection(x_hat,f_hat,iter,a,b,tol,fun)
!
! Input: a   , real(dp), the lower number to consider in the to search for zero
!        b   , real(dp), the higher number to consider in theto search for zero
!        tol , real(dp), the number of grid points between x_min and x_max
!        fun , Function
!
! Output: x_hat , real(dp), number such that abs(f(x_hat))<=Tol
!         f_hat , real(dp), f(x_hat)
!         iter  , integer , number of iteration it took to find x_hat
!

! Remarks: 1) Set x= (a+b)/2.
!          2) If f(a)f(x)>0 then update a=x go to 4.
!          3) If f(a)f(x)<0 then update b=x go to 4.
!          4) If abs(f(x)-f_obj)>Tol then go back to 1.
!          5) If abs(f(x)-f_obj)<Tol then the program ends

	Subroutine bisection(x_hat,f_hat,iter,a_0,b_0,tol,fun)
	    implicit none
	    real(dp), intent(in)  :: a_0, b_0, tol
	    real(dp), intent(out) :: x_hat, f_hat
	    integer , intent(out) :: iter
	    real(dp), external    :: fun
	    real(dp)              :: fa, fb, fc, a, b, diff
	    integer               :: i, maxiter = 10000
	    
	    
	    ! Initial values
	        a = a_0
	        b = b_0

	    !Checking if the points bracket the zero
	    IF (fun(a)*fun(b).gt.0.0_dp) THEN
	        PRINT*, 'The initial guess does not bracket.'
	        RETURN
	    END IF

	    fa = fun(a)
	    fb = fun(b)

	    iter = 1
	    diff = 1000
	    do while (diff>tol)
	        ! Check maxiter
	            if ( iter.gt.maxiter ) then  
	                PRINT*, 'Bisection: Maximum number of iterations reached.'
	                RETURN
	            end if           
	        ! Compute middle point
	            x_hat  = (a+b)/2.0_dp
	            f_hat  = fun(x_hat)
	        ! Update 
	        if ( fa*f_hat.ge.0 ) then
	            a  = x_hat
	            fa = f_hat
	        else
	            b  = x_hat
	            fb = f_hat
	        end if

	        diff  = abs(f_hat)
	        iter  = iter + 1
	    end do

	end Subroutine bisection



! =============================================================================
! Secant: Computes x_hat in [a,b] such that abs(f(x_hat))<=Tol using secant method
!            The routine also gives f_hat and iter, the number of iterations.
!            This subroutine uses the secant method See NR, pp. 347 for details
!
! Usage: Call secant(x_hat,f_hat,iter,a,b,tol,fun)
!
! Input: a   , real(dp), the lower number to consider in the to search for zero
!        b   , real(dp), the higher number to consider in theto search for zero
!        tol , real(dp), the number of grid points between x_min and x_max
!        fun , Function
!
! Output: x_hat , real(dp), number such that abs(f(x_hat))<=Tol
!         f_hat , real(dp), f(x_hat)
!         iter  , integer , number of iteration it took to find x_hat
!
! Remarks: 1) Approximate f(.) between a and b with a linear function f=beta*x+m
!             beta = (fb-fa)/(b-a) and m = fb - beta*b
!          2) Set x so that the line is zero: x = -m/beta
!          3) If f(a)f(x)>0 then update a=x go to 5.
!          4) If f(a)f(x)<0 then update b=x go to 5.
!          5) If abs(f(x)-f_obj)>Tol then go back to 1.
!          6) If abs(f(x)-f_obj)<Tol then the program ends
!          This implementation is closer to the false position method. 

	Subroutine secant(x_hat,f_hat,iter,a_0,b_0,tol,fun)
	    implicit none
	    real(dp), intent(in)  :: a_0, b_0, tol
	    real(dp), intent(out) :: x_hat, f_hat
	    integer , intent(out) :: iter
	    real(dp), external    :: fun
	    real(dp)              :: fa, fb, fc, a, b, diff, beta
	    integer               :: i, maxiter = 10000
	    

	    ! Initial values
	        a = a_0
	        b = b_0

	    !Checking if the points bracket the zero
	    IF (fun(a)*fun(b).gt.0.0_dp) THEN
	        PRINT*, 'The initial guess does not bracket.'
	        RETURN
	    END IF

	    !Values of the Function and Checking Convergence
	    fa = fun(a)
	    fb = fun(b)

	    iter = 1
	    diff = 1000
	    do while (diff>tol)
	        ! Check maxiter
	        if ( iter.gt.maxiter ) then  
	            PRINT*, 'Secant: Maximum number of iterations reached.'
	            RETURN
	        end if           
	        
	        ! Compute middle point
	        beta  = (fb - fa)/(b-a)
	        x_hat = (-(fb - beta*b))/beta
	        f_hat = fun(x_hat)

	        ! Update 
	        if ( fa*f_hat.ge.0 ) then
	            a  = x_hat
	            fa = f_hat
	        else
	            b  = x_hat
	            fb = f_hat
	        end if

	        diff  = abs(f_hat)
	        iter  = iter + 1

	    end do

	end Subroutine secant


! =============================================================================
! Newton: Computes x_hat in [a,b] such that abs(f(x_hat))<=Tol using Newton's method
!         The routine also gives f_hat and iter, the number of iterations.
!         This subroutine uses the secant method. 
!         Two-Sided numerical deriatives are used.
!
! Usage: Call newton(x_hat,f_hat,iter,x0,tol,fun)
!
! Input: x0   , real(dp), Initial guess
!        tol , real(dp), the number of grid points between x_min and x_max
!        fun , Function
!
! Output: x_hat , real(dp), number such that abs(f(x_hat))<=Tol
!         f_hat , real(dp), f(x_hat)
!         iter  , integer , number of iteration it took to find x_hat
!
! Remarks: 1) Compute df(x) = (f(x+h)-f(x-h))/2h
!          2) Set x' so that the line is zero: x' = x - f(x)/df(x)
!          3) If f(x)/df(x)<=Tol then there is no change in x. Program ends.

	Subroutine newton(x_hat,f_hat,iter,x0,tol,fun)
	    implicit none
	    real(dp), intent(in)  :: x0, tol
	    real(dp), intent(out) :: x_hat, f_hat
	    integer , intent(out) :: iter
	    real(dp), external    :: fun
	    real(dp)              :: f, df, eps, eps1, x1, x2, h, x, diff
	    integer               :: i, maxiter = 10000

	    ! Set initial value
	    x = x0

	    ! Set small number for derivative
	    eps1 = Epsilon(eps1)
	    eps  = Epsilon(eps)**(1.0/3.0)

	    iter = 1
	    diff = 1000
	    do while (diff>tol)
	        ! Check maxiter
	        if ( iter.gt.maxiter ) then  
	            PRINT*, 'Newton: Maximum number of iterations reached.'
	            RETURN
	        end if 

	        f = fun(x)

	        ! Derivative
	        h    = eps*dabs(x)
	        if (h .le. eps1) h=eps
	        df = (fun(x+h)-fun(x-h))/(2.0_dp*h)
	        
	        ! Update x and f
	        x = x - f/df

	        diff  = abs(f/df)
	        iter  = iter + 1
	    end do

	    x_hat = x
	    f_hat = f

	end Subroutine newton


! =============================================================================
! Brent: Computes x_hat in [a,b] such that abs(f(x_hat))<=Tol Brent's method
!        The routine also gives f_hat and iter, the number of iterations.
!        This routine is based on the BRENT routine of: Numerical Recipes for FORTRAN 90, Ch.10      
!
! Usage: Call brent(x_hat,f_hat,iter,a,b,tol,fun)
!
! Input: a   , real(dp), the lower number to consider in the to search for zero
!        b   , real(dp), the higher number to consider in theto search for zero
!        tol , real(dp), tolerance level
!        fun , Function
!
! Output: x_hat , real(dp), number such that abs(f(x_hat))<=Tol
!         f_hat , real(dp), f(x_hat)
!         iter  , integer , number of iteration it took to find x_hat
!

	Subroutine brent(x_hat,f_hat,iter,a_0,b_0,tol,fun)
	    implicit none
	    real(dp), intent(in)  :: a_0, b_0, tol
	    real(dp), intent(out) :: x_hat, f_hat
	    integer , intent(out) :: iter
	    real(dp), external    :: fun
	    integer , parameter   :: maxiter = 10000
	    real(dp), parameter   :: CGOLD = .3819660, ZEPS = 1.0e-10
	    real(dp)              :: a,b,c,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm,fa,fb

	    ! Initial values
	        a = a_0
	        b = b_0
	        c = (a+b)/2.0_dp

	        fa = fun(a)
	        fb = fun(b)

	    !Checking if the points bracket the zero
	    IF (fa*fb.gt.0.0_dp) THEN
	        PRINT*, 'The initial guess does not bracket.'
	        write(*,*) "fun(a)=", fa, "fun(b)=", fb, "fun(a)*fun(b)", fa*fb
	        RETURN
	    END IF

	   ! Using the brecketed values for the minization (from here is BRENT function on NumRec F90)
	    v = c
	    w = v
	    x = v
	    e = 0.0_dp
	    fx = abs(fun(x))
	    fv = fx
	    fw = fx

	    DO iter = 1,maxiter
	        xm = 0.5_dp*(a+b)
	        tol1 = tol*ABS(x) + ZEPS
	        tol2 = 2.0_dp*tol1
	        IF ((ABS(x-xm)).le.(tol2- .5_dp*(b-a))) GO TO 4
	        IF (ABS(e).gt.tol1) THEN
	            r = (x-w)*(fx-fv)
	            q = (x-v)*(fx-fw)
	            p = (x-v)*q - (x-w)*r
	            q = 2.0_dp*(q-r)
	            IF (q.gt.0.0_dp) p = -q
	            q = ABS(q)
	            etemp = e
	            e = d
	            if (ABS(p).ge.ABS(.5_dp*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x)) GO TO 2
	            d = p/q
	            u = x + d
	            IF(u - a.lt.tol2 .or. b-u.lt.tol2) d = SIGN(tol1,xm-x)
	            GO TO 3
	        END IF
	2           IF (x.ge.xm) THEN
	            e = a-x
	        ELSE
	            e = b-x
	        END IF
	        d = CGOLD*e
	3           IF (ABS(d).ge.tol1) THEN
	            u = x+d
	        ELSE
	            u = x + SIGN(tol1,d)
	        END IF
	        fu = ABS(fun(u))
	        IF (fu.le.fx) THEN
	            IF (u.ge.x) THEN
	                a = x
	            ELSE
	                b = x
	            END IF
	            v = w
	            fv = fw
	            w=x
	            fw = fx
	            x = u
	            fx = fu
	        ELSE
	            IF(u.lt.x) THEN
	                a = u
	            ELSE
	                b = u
	            ENDIF
	            IF(fu.le.fw .or. ABS(w-x).lt.small_p) THEN
	                v = w
	                fv = fw
	                w = u
	                fw = fu
	            ELSE IF (fu.le.fv .or. ABS(v-x).lt.small_p .or. (v-w).lt.small_p) THEN
	                v = u
	                fv = fu
	            END IF
	        END IF
	    END DO
	    PRINT*, 'Brent: Maximum number of iterations reached.'
	4       x_hat = x
	    f_hat = fx

	end Subroutine brent

! =============================================================================
! Brent_P: Computes x_hat in [a,b] such that abs(f(x_hat,p))<=Tol Brent's method
!          The routine also gives f_hat and iter, the number of iterations.
!          The function f is allowed to depend on x and a parameter of dimension n_p
!          This routine is based on the BRENT routine of: Numerical Recipes for FORTRAN 90, Ch.10      
!
! Usage: Call brent(x_hat,f_hat,iter,a,b,tol,fun,p,n_p)
!
! Input: a   , real(dp), the lower number to consider in the to search for zero
!        b   , real(dp), the higher number to consider in theto search for zero
!        tol , real(dp), tolerance level
!        fun , Function
!        p   , real(dp), dimension(n_p), parameter vector
!        n_p , integer , dimension of parameter vector
!
! Output: x_hat , real(dp), number such that abs(f(x_hat))<=Tol
!         f_hat , real(dp), f(x_hat)
!         iter  , integer , number of iteration it took to find x_hat
!

	Subroutine brent_p(x_hat,f_hat,iter,a_0,b_0,tol,fun,par,n_p)
	    implicit none
	    real(dp), intent(in)  :: a_0, b_0, tol
	    integer , intent(in)  :: n_p
	    real(dp), dimension(n_p), intent(in) :: par
	    real(dp), intent(out) :: x_hat, f_hat
	    integer , intent(out) :: iter
	    real(dp), external    :: fun
	    integer , parameter   :: maxiter = 10000
	    real(dp), parameter   :: CGOLD = .3819660, ZEPS = 1.0e-10
	    real(dp)              :: a,b,c,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm,fa,fb

	    ! Initial values
	        a = a_0
	        b = b_0
	        c = (a+b)/2.0_dp

	        fa = fun(a,par)
	        fb = fun(b,par)

	    !Checking if the points bracket the zero
	    IF (fa*fb.gt.0.0_dp) THEN
	        PRINT*, 'The initial guess does not bracket.'
	        !write(*,*) "fun(a)=", fun(a,par), "fun(b)=", fun(b,par)
	        RETURN
	    END IF

	   ! Using the brecketed values for the minization (from here is BRENT function on NumRec F90)
	    v = c
	    w = v
	    x = v
	    e = 0.0_dp
	    fx = abs(fun(x,par))
	    fv = fx
	    fw = fx

	    DO iter = 1,maxiter
	        xm = 0.5_dp*(a+b)
	        tol1 = tol*ABS(x) + ZEPS
	        tol2 = 2.0_dp*tol1
	        IF ((ABS(x-xm)).le.(tol2- .5_dp*(b-a))) GO TO 4
	        IF (ABS(e).gt.tol1) THEN
	            r = (x-w)*(fx-fv)
	            q = (x-v)*(fx-fw)
	            p = (x-v)*q - (x-w)*r
	            q = 2.0_dp*(q-r)
	            IF (q.gt.0.0_dp) p = -q
	            q = ABS(q)
	            etemp = e
	            e = d
	            if (ABS(p).ge.ABS(.5_dp*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x)) GO TO 2
	            d = p/q
	            u = x + d
	            IF(u - a.lt.tol2 .or. b-u.lt.tol2) d = SIGN(tol1,xm-x)
	            GO TO 3
	        END IF
	2           IF (x.ge.xm) THEN
	            e = a-x
	        ELSE
	            e = b-x
	        END IF
	        d = CGOLD*e
	3           IF (ABS(d).ge.tol1) THEN
	            u = x+d
	        ELSE
	            u = x + SIGN(tol1,d)
	        END IF
	        fu = ABS(fun(u,par))
	        IF (fu.le.fx) THEN
	            IF (u.ge.x) THEN
	                a = x
	            ELSE
	                b = x
	            END IF
	            v = w
	            fv = fw
	            w=x
	            fw = fx
	            x = u
	            fx = fu
	        ELSE
	            IF(u.lt.x) THEN
	                a = u
	            ELSE
	                b = u
	            ENDIF
	            IF(fu.le.fw .or. ABS(w-x).lt.small_p) THEN
	                v = w
	                fv = fw
	                w = u
	                fw = fu
	            ELSE IF (fu.le.fv .or. ABS(v-x).lt.small_p .or. (v-w).lt.small_p) THEN
	                v = u
	                fv = fu
	            END IF
	        END IF
	    END DO
	    PRINT*, 'Brent: Maximum number of iterations reached.'
	4       x_hat = x
	    f_hat = fx

	end Subroutine brent_p

! =============================================================================
! grid: Creates grid between x_min and x_max with curvature c_grid
!
! Usage: x_grid = grid(x_min,x_max,n_grid,c_grid)
!
! Input: x_min   , real(dp), the lower number to consider in the grid to search for zero
!        x_max   , real(dp), the higher number to consider in the grid to search for zero
!        n_grid  , integer , the number of grid points between x_min and x_max
!        c_grid  , real(dp), the curvature of the grid (how close are points near x_min)
!
! Output: x_grid , real(dp), dimension(n_grid), Grid for x between x_min and x_max
!

    Function grid(x_min,x_max,n_grid,c_grid)
        implicit none
        integer               :: i
        integer , intent(in)  :: n_grid
        real(dp), intent(in)  :: x_min, x_max, c_grid
        real(dp) :: grid(n_grid)

        ! Define the Grid
        grid=0.0_dp
        do i=1,n_grid
            grid(i) = real(i-1,dp)/real(n_grid-1,dp)
        end do

        grid = grid**c_grid
        grid = grid*(x_max-x_min)+x_min

    end Function grid



! =============================================================================
!
! Linear_IP: Conducts Linear Interpolation to compute f(x) given a grid
!              Given a function y=f(x) tabulated in xvec and yvec and a point x0,
!              return the function value y0=f(x0) obtained from linear interpolations
!              between x1<x<x2.
!
! Usage: y0 = Linear_IP(n,xvec,yvec,x0)
!
! Input: n_grid, integer ,                    Size of the grid
!        x_grid, real(dp), dimension(n_grid), Values of the grid
!        f_grid, real(dp), dimension(n_grid), Values of the objective function at the grid
!        x     , real(dp),                  , Point for interpolation
!
! Output: y0, real(dp), Value of function at x.
!
! Remarks: Taken from Heer & Maussner (2nd Edition) - Original source in file Funcion.for

    Function Linear_IP(n_grid,x_grid,f_grid,x)
        implicit none
        integer  :: n_grid, j
        real(dp) :: x, x_grid(n_grid), f_grid(n_grid), Linear_IP

    
        if ((x.lt.x_grid(1)) .or. (x.gt.x_grid(n_grid))) then
            Print *, "Linear Interpolation Error: Input off of grid!"
            Linear_IP = x
            Return
        end if
    
        if (x.eq.x_grid(1)) then
            Linear_IP = f_grid(1)
        elseif (x.eq.x_grid(n_grid)) then
            Linear_IP = f_grid(n_grid)
        else
            j = count(x_grid.le.x) ! this determines the lower bracket for x
            Linear_IP=f_grid(j)+((f_grid(j+1)-f_grid(j))/(x_grid(j+1)-x_grid(j)))*(x-x_grid(j))
        end if

        Return

    End Function Linear_IP


! =============================================================================
!
! cspline_coeff: Computes second derivates for cubic splines. Does not compute interpolation.
!
! Usage: ddf = cspline_coeff(n_grid,x_grid,f_grid,method,df_1,df_n)
!
! Input: n_grid, integer ,                    Size of the grid
!        x_grid, real(dp), dimension(n_grid), Values of the grid
!        f_grid, real(dp), dimension(n_grid), Values of the objective function at the grid
!        method, integer ,                    Cubic spline method
!                                             1 - natural cubic spline (end point derivates = 0) 
!                                             2 - secant hermite spline (end point derivates = secant)
!                                             3 - User specified first derivatives (df_1 and df_n)
!        df_1   , real(dp),                   1st derivative at x(1) - Optional
!        df_n   , real(dp),                   1st derivative at x(n) - Optional
!
! Output: ddf, real(dp), dimension(n_grid), Second derivatives at nodes.
!
! Remarks: Taken from Heer & Maussner (2nd Edition) - Original source in file Funcion.for
!          They took it from Numerical Recepies Ch 3.

    Function cspline_coeff(n_grid,x_grid,f_grid,method,df_1,df_n) result(ddf)
        implicit none
        integer  :: n_grid, method, i, k
        real(dp) :: x_grid(n_grid), f_grid(n_grid), df_1, df_n, ddf(n_grid)
        real(dp) :: p, qn, sig, un, u(n_grid)

        ddf =0.0

        if ((method.ne.1) .and. (method.ne.2) .and. (method.ne.3)) then
            Print *, "Cubic Spline Error: Method has to be either 1, 2 or 3"
        end if

        Select Case(method)
            Case(1)     ! natural cubic spline
                ddf(1) = 0.0d0
                u(1)   = 0.0d0
                qn     = 0.d0
                un     = 0.d0               
            Case(2)     ! secant hermite
                !df_1   = (f_grid(2)-f_grid(1))/(x_grid(2)-x_grid(1))
                !df_n   = (f_grid(n_grid)-f_grid(n_grid-1))/(x_grid(n_grid)-x_grid(n_grid-1))
                ! When the first derivatives are replaced by the slope of the secant 
                ! the last term in u(1) and un are made zero. See case 3 for the general formula of u.
                ddf(1) = -0.5d0
                u(1)   = 0.0d0
                qn     = 0.5d0
                un     = 0.0d0
            Case(3)     ! End points supplied by user
                ddf(1) = -0.5d0
                u(1)   = (3.d0/(x_grid(2)-x_grid(1)))*((f_grid(2)-f_grid(1))/(x_grid(2)-x_grid(1)) - df_1)
                qn     = 0.5d0
                un     = (3.d0/(x_grid(n_grid)-x_grid(n_grid-1)))* &
                         (df_n - (f_grid(n_grid)-f_grid(n_grid-1))/(x_grid(n_grid)-x_grid(n_grid-1)))
        End Select          

        Do i=2,n_grid-1
            sig    = (x_grid(i)-x_grid(i-1))/(x_grid(i+1)-x_grid(i-1))
            p      = sig*ddf(i-1)+2.d0
            ddf(i) = (sig-1.d0)/p
            u(i)   = (6.d0*((f_grid(i+1)-f_grid(i))/(x_grid(i+1)-x_grid(i)) - (f_grid(i)-f_grid(i-1))/(x_grid(i)-x_grid(i-1))) &
                      /(x_grid(i+1)-x_grid(i-1))-sig*u(i-1))/p
        End Do

        ddf(n_grid) = (un-qn*u(n_grid-1))/(qn*ddf(n_grid-1)+1.d0)

        Do k=n_grid-1,1,-1
            ddf(k)=ddf(k)*ddf(k+1)+u(k)
        End Do

        return

    end Function cspline_coeff



! =============================================================================
!
! cspline_IP: Computes cubic spline approximation for points x_eval. 
!             First and second derivatives are also provided.
!             Uses dff from cspline_coeff.
!
! Usage: call cspline_IP(n,x_grid,f_grid,ddf,method,x_eval,m,f_IP,df_IP,ddf_IP)
!
! Input: n     , integer ,               Size of the grid
!        x_grid, real(dp), dimension(n), Values of the grid
!        f_grid, real(dp), dimension(n), Values of the objective function at the grid
!        ddf   , real(dp), dimension(n), Values of the second derivative of f at grid nodes
!        method, integer ,               Bracketing Method
!                                             1 - Use if x_grid nodes are equally spaced
!                                             2 - Use if log(x_grid) are equally spaced
!                                             3 - Bisection
!        x_eval, real(dp), dimension(m), Values of x at which funciton is interpolated
!        m     , integer ,               Number of points to interpolateptional
!
! Output: f_IP  , real(dp), dimension(m), Value of the function at x_eval
!         df_IP , real(dp), dimension(m), Value of the first derivative at x_eval
!         ddf_IP, real(dp), dimension(m), Value of the second derivative at x_eval
!
! Remarks: Taken from Heer & Maussner (2nd Edition) - Original source in file Funcion.for
!          They took it from Numerical Recepies Ch 3.

    Subroutine cspline_IP(n,x_grid,f_grid,ddf,method,x_eval,m,f_IP,df_IP,ddf_IP)
        implicit none
        integer  :: n, m, method
        real(dp) :: x_grid(n), f_grid(n), ddf(n), x_eval(m), f_IP(m), df_IP(m), ddf_IP(m)
        integer  :: k, klo, khi, i
        real(dp) :: A, B, C, D, Ap, Bp, Cp, Dp, h

        if ((method.ne.1) .and. (method.ne.2) .and. (method.ne.3)) then
            Print *, "Cubic Spline Error: Method has to be either 1, 2 or 3"
        end if

        do i=1,m ! For each element the interpolation is carried out

            if (x_eval(i).ge.x_grid(n)) then
                f_IP(i)=f_grid(n)
            else

            ! Bracket x_eval between two grid points of indeces klo and khi
            Select Case(method)
                Case(1)     ! equi-distant points
                      klo = floor((x_eval(i)-x_grid(1))/(x_grid(2)-x_grid(1)))+1
                      khi = 1+klo
                    if (khi.gt.n) then
                        print *, "CSpline interpolation error: x_eval outside of grid! Program stops"
                        stop
                    endif
                Case(2) ! equi-distant point in log scale
                    klo = floor( (dlog(x_eval(i))-dlog(x_grid(1)))/(dlog(x_grid(2))-dlog(x_grid(1))))+1
                    khi = 1+klo
                    if (khi.gt.n) then
                        print *, "CSpline interpolation error: x_eval outside of grid! Program stops"
                        stop
                    endif
                Case(3)     ! bisection
                    klo = 1
                    khi = n
                    do while ((khi-klo).gt.1)
                        k = (khi+klo)/2
                        if (x_grid(k).gt.x_eval(i)) then
                            khi = k
                        else 
                            klo = k
                        endif
                    Enddo
            End Select          

            ! Step between bracketing of x_eval
            h = x_grid(khi)-x_grid(klo)
            
            ! Coefficients for interpolation
            A = (x_grid(khi)-x_eval(i))/h
            B = (x_eval(i)-x_grid(klo))/h
            C = (a**3-a)*(h**2)/6.d0
            D = (b**3-b)*(h**2)/6.d0
            ! Coefficients for derivatives
            Ap = -1/h
            Bp = -ap
            Cp = (3*(A**2)-1)*Ap*(h**2)/6.d0
            Dp = (3*(B**2)-1)*Bp*(h**2)/6.d0
            ! Interpolation
            f_IP(i)   = A*f_grid(klo) + B*f_grid(khi) + C*ddf(klo) + D*ddf(khi)
            df_IP(i)  = Ap*f_grid(klo) + Bp*f_grid(khi) + Cp*ddf(klo) + Dp*ddf(khi)
            ddf_IP(i) = A*ddf(klo) + B*ddf(khi)

            endif
        end Do
        return
    end Subroutine cspline_IP    

! =============================================================================
! Min_Bracket: Brackets a minimum of a continous function f(x) between a_0 and c_0
!              The routine does not look for the smallest bracket but for the first one
!              An initial bracket is provided. If it does not work the routine uses
!              bisection to create a new candidate.
!              A bracket for a minimum is a triple (a,b,c) such that f(b)<f(a) and f(b)<f(c)
!
! Usage: Call Min_Bracket(a_0,b_0,c_0,fun,a,b,c)
!
! Input: a_0 , real(dp), the lower bound for original bracket
!        b_0 , real(dp), a number between a_0 and b_0
!        c_0 , real(dp), the upper bound for original bracket
!        fun , Function
!
! Output: a , real(dp), the lower bound for final bracket
!         b , real(dp), a number between a and b
!         a , real(dp), the upper bound for final bracket
!
! Note: It must be that a<b<c
!
	Recursive Subroutine Min_Bracket(a_0,b_0,c_0,fun,a,b,c)
		real(dp), intent(in)  :: a_0, b_0, c_0
		real(dp), intent(out) :: a, b, c
		real(dp), external    :: fun

		if (fun(a_0).lt.fun(b_0)) then
			call Min_Bracket(a_0,(a_0+b_0)/2,b_0,fun,a,b,c)
		elseif (fun(c_0).lt.fun(b_0)) then
			call Min_Bracket(b_0,(b_0+c_0)/2,c_0,fun,a,b,c)
		else
			a = a_0
			b = b_0
			c = c_0
		end if

	end Subroutine Min_Bracket


! =============================================================================
!
! GSS: Golden section search to find the maximizer of y=f(x)
!   
! Usage: x_hat = GSS(f,xl,xu)
!
! Input: f  : pointer to the function Real(dp) function y = f(x)
!        xl : Real(dp), lower bound of the interval in which the maximizer lies.
!        xu : Real(dp), upper bound of the interval in which the maximizer lies.
!
! Output: x_hat, Real(dp), the approximate maximizer
!
! Remarks: The algorithm has a default tolerance level. Approx. size of interval is 10^-8
!          This algorithm is taken from Mausner 2008 optimization.for
!          Corresponds to Algorithm 11.6.1 in Heer & Mausner (2nd Ed)
!

    Function GSS(f,xl,xu) result(x_hat)
        implicit none
        real(dp) :: xl, xu, x_hat
        real(dp) :: tol, p, q, a, b, c, d, fb, fc, test(1:2) 
        integer  :: i
        real(dp), external :: f

        ! Check bounds
        if (xl.ge.xu) print*, "Golden Section Error: Bounds for search do not satisfy xl<xu"

        ! Compute the parameters of the problem */
        tol = dsqrt(Epsilon(tol))
        p   = 0.5d0*(dsqrt(5.d0)-1.d0)
        q   = 1.d0-p

        ! Compute the initial interval [a,d] and the points b and c that divide it
        a = xl
        d = xu
        b = p*a+q*d
        c = q*a+p*d

        ! Compute the function value at b and c 
        fb = f(b)
        fc = f(c)
   
        ! Iterate so that the inverval gets small and smaller
        test(1) = 1
        test(2) = dabs(a)+dabs(c)

        ! Bisection procedure
        do while (dabs(d-a).gt.tol*MaxVal(test))
            if (fb.lt.fc) then ! choose [b,d] as next interval
                a  = b         
                b  = c
                fb = fc
                c  = p*b+q*d
                fc = f(c)
            else ! choose [a,c] as next interval
                d  = c
                c  = b
                fc = fb
                b  = p*c + q*a
                fb = f(b)
            endif
            test(2)=dabs(a)+dabs(c)
        End Do

        if (fb.gt.fc) then
            x_hat = b
        else
            x_hat = c
        endif

        Return
    End Function GSS


! =============================================================================
!
! brent_opt: Computes the minimum and minimizer of a real valued function using Brent's method.
!
! Usage: f_min = brent_opt(ax,bx,cx,func,tol,x_min) ; x_min is the minizer.
!
! Input: ax   , real(dp), Lower bound of interval bracketing a minimum
!        bx   , real(dp), Some middle point of interval bracketing a minimum
!        cx   , real(dp), Upper bound of interval bracketing a minimum
!        func ,           Pointer to a Real(dp) function to be minimized 
!        tol  , real(dp), Tolerance for stopping criteria. Don't set lower than sqrt(eps)
!        x_min, real(dp), Variable where minizer will be stored
!
! Output: f_min, real(dp), Minimum value of func between ax and cx.
!         x_min, real(dp), Minimizer of the function.
!
! Remarks:  Taken from Numerical Recepies, Section 10.2
!           Given a function func, and given a bracketing triplet of abscissas ax, bx, cx 
!           (such that bx is between ax and cx, and func(bx) is less than both func(ax) and func(cx)), 
!           this routine isolates the minimum to a fractional precision of about tol using Brent’s method. 
!           The abscissa of the minimum is returned as xmin, and the minimum function value is returned as brent_opt, 
!           the returned function value.
!           Interal Parameters: Maximum allowed number of iterations; golden ratio; 
!           and a small number that protects against trying to achieve fractional accuracy for a minimum that happens to be exactly zero.

    Function brent_opt(ax,bx,cx,func,tol,x_min)
        implicit none
        integer               :: iter
        real(dp), intent(in)  :: ax, bx, cx, tol
        real(dp), intent(out) :: x_min
        real(dp), external    :: func
        real(dp), parameter   :: CGOLD=0.3819660_dp, ZEPS=1.0e-3_dp*epsilon(ax)
        integer , parameter   :: ITMAX=100
        real(dp) :: brent_opt
        real(dp) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm,brent

        if ((func(ax).lt.func(bx)).or.(func(cx).lt.func(bx))) then
            print*, "Brent Error: A minimum is not bracketed. It must be that f(bx)<f(ax) & f(bx)<f(cx)"
        end if

        ! Initializations
        a  = min(ax,cx)
        b  = max(ax,cx)
        v  = bx
        w  = v
        x  = v
        e  = 0.0
        fx = func(x) 
        fv = fx
        fw = fx


        ! Main Loop
        do iter=1,ITMAX
            xm   = 0.5_dp*(a+b)
            tol1 = tol*abs(x)+ZEPS
            tol2 = 2.0_dp*tol1

            ! Test for stopping criterion
            if (abs(x-xm) <= (tol2-0.5_dp*(b-a))) then
                x_min     = x 
                brent     = fx 
                brent_opt = brent
                RETURN
            end if 

            ! Construct a trial parabolic fit - This gives new value for "d"
            if (abs(e) > tol1) then
                r = (x-w)*(fx-fv) 
                q = (x-v)*(fx-fw) 
                p = (x-v)*q-(x-w)*r 
                q = 2.0_dp*(q-r)
                if (q > 0.0) p=-p 
                q     = abs(q)
                etemp = e
                e     = d

                ! Check acceptability of parabolic fit. 
                ! If not fit change to golden section. Else take parabolic step.
                if (abs(p) >= abs(0.5_dp*q*etemp) .or. p <= q*(a-x) .or. p >= q*(b-x)) then
                    e=merge(a-x,b-x, x >= xm )
                    d=CGOLD*e
                else 
                    d=p/q
                    u=x+d
                    if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
                end if 

            else ! Take the golden section step into the larger of the two segments
                e=merge(a-x,b-x, x >= xm )
                d=CGOLD*e
            end if

            u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
            fu=func(u)

            ! The following lines decide which way to go for the next iteration
            if (fu <= fx) then
                if (u >= x) then
                    a=x 
                else
                    b=x 
                end if
                ! call shft(v,w,x,u)
                v = w
                w = x
                x = u
                ! call shft(fv,fw,fx,fu)
                fv = fw
                fw = fx
                fx = fu
            else
                if (u < x) then 
                    a=u
                else 
                    b=u
                end if
                if (fu <= fw .or. w == x) then
                    v=w
                    fv=fw
                    w=u
                    fw=fu
                else if (fu <= fv .or. v == x .or. v == w) then 
                    v=u
                    fv=fu 
                end if
            end if
        end do

        print *, "Brent Error: Maximum number of iterations reached"

        !Contains
        !    Subroutine shft(a,b,c,d) 
        !        REAL(dp), INTENT(OUT) :: a 
        !        REAL(dp), INTENT(INOUT) :: b,c 
        !        REAL(dp), INTENT(IN) :: d
        !        a=b
        !        b=c
        !        c=d
        !    end Subroutine shft 
    end Function brent_opt


! =============================================================================
!
! dbrent_opt: Computes the minimum and minimizer of a real valued function using D-Brent's method.
!
! Usage: f_min = dbrent_opt(ax,bx,cx,func,dfunc,tol,x_min) ; x_min is the minizer.
!
! Input: ax    , real(dp), Lower bound of interval bracketing a minimum
!        bx    , real(dp), Some middle point of interval bracketing a minimum
!        cx    , real(dp), Upper bound of interval bracketing a minimum
!        func  ,           Pointer to a Real(dp) function to be minimized 
!        dfunc ,           Pointer to a Real(dp) function that gives the first derivative of func
!        tol   , real(dp), Tolerance for stopping criteria. Don't set lower than sqrt(eps)
!        x_min , real(dp), Variable where minizer will be stored
!
! Output: f_min, real(dp), Minimum value of func between ax and cx.
!         x_min, real(dp), Minimizer of the function.
!
! Remarks:  Taken from Numerical Recepies, Section 10.3
!           Given a function func and its derivative function dfunc, and given a bracketing triplet of abscissas 
!           ax, bx, cx [such that bx is between ax and cx, and func(bx) is less than both func(ax) and func(cx)], 
!           this routine isolates the minimum to a fractional precision of about tol using a modification of 
!           Brent’s method that uses derivatives. 
!           The abscissa of the minimum is returned as xmin, and the minimum function value is returned as dbrent.
!           Inernal Parameters: Maximum allowed number of iterations, and a small number that protects against trying to 
!           achieve fractional accuracy for a minimum that happens to be exactly zero.

    Function dbrent_opt(ax,bx,cx,func,dfunc,tol,x_min)
        implicit none
        integer  :: iter
        real(dp), intent(in)  :: ax, bx, cx, tol 
        real(dp), intent(out) :: x_min
        real(dp) :: dbrent_opt
        real(dp) :: a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm,dbrent
        real(dp), external  :: func, dfunc
        real(dp), parameter :: CGOLD=0.3819660_dp, ZEPS=1.0e-3_dp*epsilon(ax)
        integer , parameter :: ITMAX=100
        logical  :: ok1,ok2 ! These are used as flags for whether proposed steps are accepted or not

        if ((func(ax).lt.func(bx)).or.(func(cx).lt.func(bx))) then
            print*, "D-Brent Error: A minimum is not bracketed. It must be that f(bx)<f(ax) & f(bx)<f(cx)"
        end if
        ! Initializations
        a  = min(ax,cx)
        b  = max(ax,cx)
        v  = bx
        w  = v
        x  = v
        e  = 0.0
        fx = func(x) 
        fv = fx
        fw = fx
        dx = dfunc(x) 
        dv = dx
        dw = dx


        ! Main Loop
        do iter=1,ITMAX
            xm   = 0.5_dp*(a+b)
            tol1 = tol*abs(x)+ZEPS
            tol2 = 2.0_dp*tol1

            ! Test for stopping criterion
            if (abs(x-xm) <= (tol2-0.5_dp*(b-a))) exit

            ! Construct a trial parabolic fit - This gives new value for "d"
            if (abs(e) > tol1) then
                ! Initialize d's to an out-of-bracket value
                d1 = 2.0_dp*(b-a)
                d2 = d1

                ! Check and apply secant method
                if (dw /= dx) d1=(w-x)*dx/(dx-dw)
                if (dv /= dx) d2=(v-x)*dx/(dx-dv)

                ! Value of d's must be within the bracket, and on the side pointed to by the derivative at x:
                u1   = x+d1
                u2   = x+d2
                ok1  = ((a-u1)*(u1-b) > 0.0) .and. (dx*d1 <= 0.0) 
                ok2  = ((a-u2)*(u2-b) > 0.0) .and. (dx*d2 <= 0.0)
                olde = e
                e    = d

                ! This "if" takes only acceptable values for d, and if both are acceptable takes the smallest.
                if (ok1 .or. ok2) then
                    if (ok1 .and. ok2) then
                        d = merge(d1,d2, abs(d1) < abs(d2))
                    else 
                        d = merge(d1,d2,ok1)
                    end if
                    
                    if (abs(d) <= abs(0.5_dp*olde)) then
                        u = x+d
                        if (u-a < tol2 .or. b-u < tol2) d = sign(tol1,xm-x)
                    else
                        e = merge(a,b, dx >= 0.0)-x
                        d = 0.5_dp*e ! Note that this is bisection, not golden section
                    end if
                else
                    e = merge(a,b, dx >= 0.0)-x
                    d = 0.5_dp*e ! Note that this is bisection, not golden section
                end if

            else

                e = merge(a,b, dx >= 0.0)-x
                d = 0.5_dp*e ! Note that this is bisection, not golden section

            end if

            if (abs(d) >= tol1) then 
                u  = x+d
                fu = func(u) 
            else
                u  = x+sign(tol1,d) 
                fu = func(u)
                if (fu > fx) exit
            end if

            ! Update Derivative
            du=dfunc(u)


            if (fu <= fx) then
                
                if (u >= x) then 
                    a=x
                else 
                    b=x
                end if
                
                !call mov3(v,fv,dv,w,fw,dw)
                v  = w
                fv = fw
                dv = dw
                !call mov3(w,fw,dw,x,fx,dx) 
                w  = x
                fw = fx
                dw = dx
                !call mov3(x,fx,dx,u,fu,du)
                x  = u
                fx = fu
                dx = du

            else

                if (u < x) then
                    a=u 
                else
                    b=u 
                end if

                if (fu <= fw .or. w == x) then
                    !call mov3(v,fv,dv,w,fw,dw)
                    v  = w
                    fv = fw
                    dv = dw
                    !call mov3(w,fw,dw,u,fu,du)
                    w  = u
                    fw = fu
                    dw = du
                else if (fu <= fv .or. v == x .or. v == w) then 
                    !call mov3(v,fv,dv,u,fu,du)
                    v  = u
                    fv = fu
                    dv = du
                end if 
            end if

        end do

        if (iter>ITMAX) print *, "D-Brent Error: Maximum number of iterations reached"

        x_min = x
        dbrent_opt = fx

    !   Contains
    !       Subroutine mov3(a,b,c,d,e,f) 
    !           REAL(SP), INTENT(IN) :: d,e,f 
    !           REAL(SP), INTENT(OUT) :: a,b,c 
    !               a=d
    !               b=e
    !               c=f
    !       END SUBROUTINE mov3 

    end Function dbrent_opt


! =============================================================================
!
! Nelder_Mead: Computes the minimum and minimizer of a real valued function using Nelder-Mead method.
!
! Usage: call Nelder_Mead(x_0,ftol,func,x_min,y_min,iter)
!
! Input: x_0   , real(dp), dimension(N), Initial point for algorithm
!        ftol  , real(dp),               Tolerance for stopping criteria
!        func  ,                         Pointer to a Real(dp) function to be minimized 
!
! Output: iter , integer ,               Number of iterations the algorithm used
!         x_min, real(dp), dimension(N), Minimizer of the function.
!         f_min, real(dp),               Minimum value of func at x_min
!
! Remarks:  Taken from Numerical Recepies, Section 10.4
!           Internal Parameters: Maximum allowed number of iterations
!                                A small number that protects against trying to 
!                                   achieve accuracy when function does not change
!                                Lambda: Spread of the initial simplex around x_0
    Subroutine Nelder_Mead(x_0,ftol,func,x_min,f_min,iter)
        implicit none
        
        ! Options for size and precision in the program
        integer , parameter ::  DP      = KIND(1.0D0)

        ! Nelder Mead Inputs
        real(dp), dimension(:)  , intent(in)    :: x_0
        real(dp),                 intent(in)    :: ftol
        real(dp), external :: func

        ! Nelder Mead Parameters
        integer , parameter :: ITMAX = 5000
        real(dp), parameter :: tiny = 1.0e-10
        real(dp), parameter :: lambda = 3
        integer             :: ndim

        ! Nelder Mead Ouputs
        integer , intent(out) :: iter
        real(dp), intent(out) :: f_min, x_min(size(x_0))

        ! Global Variables to be used by subfunctions
        !   p: A (N+1,N) matrix with the points of the simplex. It has N+1 points, each of dimension N
        !   y: A (N) vector with the value of the function at each point of the simplex
        integer  :: ihi, i
        real(dp), dimension(:) :: y(size(x_0)+1), psum(size(x_0)), p(size(x_0)+1,size(x_0))

        ! Variables for central processing
        integer  :: ilo,inhi
        real(dp) :: rtol,ysave,ytry,ytmp,dum,dum_vec(size(x_0))
        
        ! Variables for amotry
        real(dp) :: fac,fac1,fac2,ptry(size(x_0))

        ! Initial values
            if (rank(x_0).eq.1) then
                ndim = size(x_0)
            else
                print *, "Nelder Mead Error: x_0 must be a vector"
                stop
            end if 

        ! Initial simplex
            do i=1,ndim
                p(i,:) = x_0
                p(i,i) = x_0(i)+lambda
                y(i)   = func(p(i,:))
            end do
                p(ndim+1,:) = x_0
                y(ndim+1)   = func(p(ndim+1,:))

        ! Operation
        ! This function is the body of the program
        ! The function first checks the stopping criteria and then computes new points for simplex
        ! Best result is returned in the first element of the p and y objects

        ! Initialization
        iter = 0
        psum = sum(p,dim=1)

        ! Iteration Loop
        do
            ! Determine which point is the highest (worst) and lowest (best)
            ilo = minloc(y,dim=1)
            ihi = maxloc(y,dim=1)

            ! Determine which point is the second highest (inhi)
            ytmp   = y(ihi)
            y(ihi) = y(ilo)     
            inhi   = maxloc(y,dim=1)
            y(ihi) = ytmp

            ! Compute the fractional range from highest to lowest and return if satisfactory.
            !   If returning, put best point and value in slot 1
            rtol = 2.0_dp*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
            if (rtol < ftol) then
                ! call swap(y(1),y(ilo))
                    dum    = y(1)
                    y(1)   = y(ilo)
                    y(ilo) = dum 
                ! call swap(p(1,:),p(ilo,:))
                    dum_vec  = p(1,:)
                    p(1,:)   = p(ilo,:)
                    p(ilo,:) = dum_vec
                exit
            end if

            ! Check for maximum number of iterations
            if (iter >= ITMAX) then
                print *, "Nelder_Mead Error: Maximum number of iterations reached"
                stop
            end if

            ! New iteration
            !   First extrapolate by a factor ?1 through the face of the simplex across from the high point, 
            !   i.e., reflect the simplex from the high point.
            
            ! Update iteration
            iter = iter+1
            
            ! Get new try point
                ! ytry = amotry(-1.0_dp)
                fac  = -1.0_dp
                    fac1 = (1.0_dp-fac)/ndim
                    fac2 = fac1-fac 
                    ptry = psum*fac1 - p(ihi,:)*fac2 
                    ytry = func(ptry)

                    ! Evaluate the function at the trial point.
                    !   If it?s better than the highest, then replace the highest.
                    if (ytry < y(ihi)) then
                        y(ihi)   = ytry 
                        psum     = psum - p(ihi,:) + ptry 
                        p(ihi,:) = ptry
                    end if
                
            ! If new try point gives a result better than the best point then
            !   try an additional extrapolation by a factor of 2.
            ! If the reflected point is worse than the second highest, then 
            !    look for an intermediate lower point, i.e., do a one-dimensional contraction.
            if (ytry <= y(ilo)) then
                ! ytry = amotry(2.0_dp)
                fac  = 2.0_dp
                    fac1 = (1.0_dp-fac)/ndim
                    fac2 = fac1-fac 
                    ptry = psum*fac1 - p(ihi,:)*fac2 
                    ytry = func(ptry)

                    ! Evaluate the function at the trial point.
                    !   If it?s better than the highest, then replace the highest.
                    if (ytry < y(ihi)) then
                        y(ihi)   = ytry 
                        psum     = psum - p(ihi,:) + ptry 
                        p(ihi,:) = ptry
                    end if
                iter = iter+1
            else if (ytry >= y(inhi)) then
                ysave = y(ihi)
                ! ytry  = amotry(0.5_dp)
                fac  = 0.5_dp
                    fac1 = (1.0_dp-fac)/ndim
                    fac2 = fac1-fac 
                    ptry = psum*fac1 - p(ihi,:)*fac2 
                    ytry = func(ptry)

                    ! Evaluate the function at the trial point.
                    !   If it?s better than the highest, then replace the highest.
                    if (ytry < y(ihi)) then
                        y(ihi)   = ytry 
                        psum     = psum - p(ihi,:) + ptry 
                        p(ihi,:) = ptry
                    end if
                iter  = iter+1
                
                ! If can?t seem to get rid of that high point then
                ! Better contract around the lowest (best) point.
                if (ytry >= ysave) then
                    p = 0.5_dp*( p + spread( p(ilo,:) , 1 , ndim+1 ))
                    
                    do i=1,ndim+1
                        if (i /= ilo) y(i)=func(p(i,:)) 
                    end do
                    
                    iter = iter+ndim

                    psum = sum(p,dim=1)

                end if
            end if
        end do

        do i=1,ndim+1
                    print *, p(i,:)
                end do
                print *, y
        ! Final values
            x_min = p(1,:)
            f_min = y(1)

    end Subroutine Nelder_Mead 

end module Toolbox