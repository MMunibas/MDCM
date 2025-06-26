! Mike Devereux, 5.6.2025
! This routine aims to replace differential evolution with hybrid G-NES, exploiting the
! sampling capabilities of G-NES and the efficiency of analytical gradients for local
! optimization

module gnes_optimizer
  use omp_lib
  implicit none
  private
  public :: initialize_gnes_optimizer, optimize_gnes

  integer, parameter :: dp = kind(0d0)
  integer :: D, pop_size, max_iter, refine_interval
  real(dp) :: eta_nes, sigma_init, sigma_decay
  logical :: initialized = .false.

  real(dp), allocatable :: mu(:)

contains

  subroutine initialize_gnes_optimizer(initialize_mu, dim, population, iterations, &
                                         refine_every, learning_rate, sigma0, decay)
    implicit none
    integer, intent(in) :: dim, population, iterations, refine_every
    real(dp), intent(in) :: learning_rate, sigma0, decay

    ! Population initializer subroutine passed as reference
    interface
        subroutine initialize_mu(mu)
          real*8, intent(out) :: mu(:)
        end subroutine
    end interface

    D = dim ! number of parameters to fit (x,y,z,q of chgs)
    pop_size = population
    max_iter = iterations
    refine_interval = refine_every
    eta_nes = learning_rate
    sigma_init = sigma0
    sigma_decay = decay

    if(allocated(mu)) deallocate(mu)
    allocate(mu(D))
    call initialize_mu(mu)

    initialized = .true.
  end subroutine initialize_gnes_optimizer

  subroutine optimize_gnes(func, dfunc, par_range, theta_opt)
    implicit none
    integer :: i, j, iter, idx_best, no_improve_counter, max_no_improve
    real(dp) :: theta_opt(:) ! to receive optimized parameters
    real(dp) :: improve_tol, prev_best
    real(dp) :: par_range(:,:)
    real(dp) :: sigma ! represents diagonal covariance matrix
    real(dp), allocatable :: theta(:,:), epsilon(:,:), fitness(:), grad(:,:)
    real(dp), allocatable :: grad_mu(:), grad_tmp(:), x_best(:)

    ! for derivatives check
    real(dp), dimension(D) :: tx, fm, fp, gm, gp
    real(dp) :: tdx
    integer :: l

    ! for soft constraints
    real(dp) :: penalty, lambda

    ! Loss function passed to optimizer
    interface
        real*8 function func(x)
          real*8, intent(in) :: x(:)
        end function func
    end interface

    ! Gradient of loss function passed to optimizer
    interface
        function dfunc(x,n) result(y)
          integer, intent(in) :: n
          real*8, intent(in) :: x(:)
          real*8 :: y(n)
        end function dfunc
    end interface

    if (.not. initialized) then
       print *, 'Error: Optimizer not initialized.'
       stop
    end if

    allocate(theta(D, pop_size), epsilon(D, pop_size), fitness(pop_size))
    allocate(grad(D, pop_size), grad_mu(D), grad_tmp(D), x_best(D))

    sigma = sigma_init

    no_improve_counter=0
    max_no_improve=10 ! exit optimization if no improvement from last 10 iterations
    improve_tol=0.00001_dp
    prev_best=10._dp

    lambda=10._dp ! soft constraint weighting factor

    do iter = 1, max_iter
       grad_mu = 0.0_dp

#ifdef _OPENMP
       !$OMP PARALLEL DO PRIVATE(i) SHARED(mu, sigma, epsilon, theta, fitness, grad)
#endif
       do i = 1, pop_size
          ! generate new samples:
          call random_normal_vector(epsilon(:, i))
          theta(:, i) = mu + sigma * epsilon(:, i)

          fitness(i) = func(theta(:, i))
          grad(:, i) = dfunc(theta(:, i),size(theta,dim=1))
          call constraints(fitness(i),grad(:,i),theta(:,i),D,lambda,par_range)

!          ! debug: uncomment to test derivatives:
!          gp(:)=0._dp
!          gm(:)=0._dp
!          do l=1,D
!            tx(:)=theta(:, i)
!            tx(l)=tx(l)+0.0001_dp
!            fp(l)=func(tx)
!            call constraints(fp(l),gp,tx(:),D,lambda,par_range)
!          enddo
!          do l=1,D
!            tx(:)=theta(:, i)
!            tx(l)=tx(l)-0.0001_dp
!            fm(l)=func(tx)
!            call constraints(fm(l),gm,tx(:),D,lambda,par_range)
!          enddo
!          do l=1,D
!            tdx=(fp(l)-fm(l))/0.0002_dp
!            print*,l,': numerical = ',tdx,', analytical = ',grad(l,i)
!          enddo
!          stop
!          ! end test derivatives
       end do
#ifdef _OPENMP
       !$OMP END PARALLEL DO
#endif

       ! === Average Gradient ===
       do i = 1, pop_size
          grad_mu = grad_mu + grad(:, i)
       end do
       grad_mu = grad_mu / pop_size

       ! === Update Mean ===
       mu = mu - eta_nes * grad_mu

       ! === Sigma Decay ===
       sigma = max(1.e-4_dp, sigma * sigma_decay)

       ! === L-BFGS Refinement ===
       if (mod(iter, refine_interval) == 0) then
          idx_best = 1
          do i = 2, pop_size
             if (fitness(i) < fitness(idx_best)) idx_best = i
          end do
          x_best = theta(:, idx_best)

          call lbfgs_refine(func, dfunc, x_best, D, lambda, par_range)
          mu = x_best
       end if

       ! === Exit if Converged ===
       if (abs(minval(fitness)-prev_best) <= improve_tol) then
          no_improve_counter=no_improve_counter+1
       else
          no_improve_counter=0
       endif
       if (no_improve_counter >= max_no_improve) then
          print*,'   No improvement for last ',no_improve_counter,' cycles, exiting!'
          exit
       endif
       prev_best=minval(fitness)

       ! === Logging ===
       write(*, '(A,I4,A,F12.6,A,F8.4)') 'Iter:', iter, ' Best Loss:', minval(fitness), ' Sigma:', sigma
    end do

    ! copy best fit to solution array
    theta_opt(:) = x_best(:)

    deallocate(theta, epsilon, fitness, grad, grad_mu, grad_tmp, x_best)
  end subroutine optimize_gnes

  ! === Placeholder: Initialization of mu ===
!  subroutine initialize_mu(mu)
!    implicit none
!    real(dp), intent(out) :: mu(:)
!    integer :: i
!    do i = 1, size(mu)
!       mu(i) = 0.01_dp  ! Example: small initialization
!    end do
!  end subroutine initialize_mu

  ! === Placeholder: Random Normal Vector (Box-Muller) ===
  subroutine random_normal_vector(v)
    implicit none
    real(dp), intent(out) :: v(:)
    integer :: i
    real(dp) :: u1, u2
    do i = 1, size(v), 2
       call random_number(u1)
       call random_number(u2)
       v(i) = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * acos(-1.0_dp) * u2)
       if (i+1 <= size(v)) v(i+1) = sqrt(-2.0_dp * log(u1)) * sin(2.0_dp * acos(-1.0_dp) * u2)
    end do
  end subroutine random_normal_vector

  ! === Placeholder: L-BFGS Wrapper ===
  subroutine lbfgs_refine(func, dfunc, x, n, lambda, par_range)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    integer, intent(in) :: n
    real(dp), intent(inout) :: x(n)
    real(dp) :: lambda
    real(dp) :: par_range(:,:)
  
    integer :: m, iprint(2), iflag, i, icall
    integer :: k
    real(dp) :: f
    real(dp), allocatable :: g(:), diag(:), w(:)

    ! for derivatives check
    real(dp), dimension(n) :: tx, fm, fp, gm, gp
    real(dp) :: tdx

    ! Loss function passed to optimizer
    interface
        real*8 function func(x)
          real*8, intent(in) :: x(:)
        end function func
    end interface
    ! Gradient of loss function passed to optimizer
    interface
        function dfunc(x, n) result(y)
          integer, intent(in) :: n
          real*8, intent(in) :: x(:)
          real*8 :: y(n)
        end function dfunc
    end interface

    m = 5
    allocate(g(n), diag(n), w(3*n*m + 2*n))
  
    iprint(1) = 10 !-1  ! Suppress output
    iprint(2) = 0  !0

    iflag = 0
    icall = 0

    f = func(x) ! initial rmse and drmse "f" and "g"
    g = dfunc(x, n)
    call constraints(f,g,x,D,lambda,par_range)
!    ! debug: uncomment to test derivatives:
!    do i=1,n
!      tx(:)=x(:)
!      tx(i)=x(i)+0.0001_dp
!      fp(i)=func(tx)
!      call constraints(fp(i),gp,tx,D,lambda,par_range)
!    enddo
!    do i=1,n
!      tx(:)=x(:)
!      tx(i)=x(i)-0.0001_dp
!      fm(i)=func(tx)
!      call constraints(fm(i),gm,tx,D,lambda,par_range)
!    enddo
!    do i=1,n
!      tdx=(fp(i)-fm(i))/0.0002_dp
!      print*,i,': numerical = ',tdx,', analytical = ',g(i)
!    enddo
!    ! end test derivatives
  
    do
       icall = icall + 1
       call lbfgs(n, m, x, f, g, .false., diag, iprint, 1.0d-5, 1.0d-9, w, iflag)
!       call lbfgs(n, m, x, f, g, .false., diag, iprint, 1.0d-4, 1.0d-5, w, iflag)
       if (iflag == 1 .or. iflag == 2) then
         f = func(x)
         g = dfunc(x, n)
         call constraints(f,g,x,D,lambda,par_range)
       else
          exit
       end if
    end do
  
    if (allocated(diag)) deallocate(diag)
    if (allocated(w)) deallocate(w)
    if (allocated(g)) deallocate(g)
  end subroutine lbfgs_refine

  subroutine constraints(f,g,x,D,lambda,par_range)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(inout) :: f
    real(dp), intent(in) :: lambda
    real(dp), intent(inout) :: g(:)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: par_range(:,:)
    integer :: D

    real(dp) :: penalty
    integer :: j

    ! apply soft constraints:
    do j = 1, D
      if(x(j) < par_range(1,j)) then
        penalty = (x(j) - par_range(1,j))
        f = f + lambda * penalty**2
        g(j) = g(j) + 2._dp * lambda * penalty
      endif
      if(x(j) > par_range(2,j)) then
        penalty = (x(j) - par_range(2,j))
        f = f + lambda * penalty**2
        g(j) = g(j) + 2._dp * lambda * penalty
      endif
    end do
!          x(j) = max(par_range(1,j), min(par_range(2,j), x(j)))  ! for hard constraints

  end subroutine constraints

end module gnes_optimizer

