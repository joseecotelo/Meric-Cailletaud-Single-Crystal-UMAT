      program MericCailletaudTest      

      use umat_subroutines
      use math_subroutines
      use MericCailletaud_module
      implicit none      

      integer itermax
      integer i,j,k
      real*8 slip_systems(nslip,6)
      real*8 taucrss(nslip),eps_el(3,3),tau(nslip)
      real*8 deps(3,3),deps_el(3,3),gamma_n(nslip),eps_el_n(3,3)
      real*8 dt,epsdot,tol,fyield(nslip),dgamma(nslip),gamma(nslip)
      real*8 cauchy(3,3),tolNewton,gamma_accum(nslip),eps(3,3),
     &       eps_n(3,3),gamma_accum_n(nslip)

150   format(100(e20.10,1x))
      open(150,file='data.txt')

      call read_file(67,'fcc.txt',nslip,6,slip_systems)
      ! define the matprops
      eta = 1.d0
      nexp = 1.d0
      tau0 = 0.d0
      taucrss = tau0
      mathard = 0.d0
      ! store the slip normals and shear directions
      do k=1,nslip
          nschmidt(k,:) = slip_systems(k,1:3)
          sschmidt(k,:) = slip_systems(k,4:6)
      enddo
      ! initialize identity tensors
      call init_delta_f()
      ! build elasticity tensor
      cEl = (2.d0)*dev4 + (5.d0)*hydro4
      ! construc the schmidt tensor
      call init_schmidt_f()
      dt = 1.e-1
      epsdot = 1.e-2
      deps = 0.d0
      eps = 0.d0
      eps_n = 0.d0
      deps(1,2) = epsdot*dt
      deps(2,1) = epsdot*dt
      eps_el_n = 0.d0
      gamma_accum_n = 0.d0
      gamma_n = 0.d0
      tol = 1.e-12
      tolNewton = 1.e-8
      itermax = 20
      ! Solve the linear equation for the change in the increment variation
      do i=1,1
        ! update strain
        eps = deps*dt + eps_n
        ! compute trial solution
        deps_el = deps
        eps_el = eps_el_n + deps_el
        call tens4_prod_tens2_f(cauchy,cEl,eps_el)
        write(*,*) 'cauchy trial = ',cauchy
        ! check yield conditions
        do k=1,nslip
            call tens2_inner_product_f(tau(k),cauchy,schmidt(k,:,:))
        enddo
        ! evaluate yield condition
        call eval_fyield_f(fyield,tau,taucrss)
        if (any(fyield .ge. tol)) then
            write(*,*) 'solution is viscoplastic'
            call solve_nonlinear_equations_f(deps_el,dgamma,deps,
     &                            eps_el_n,gamma_n,gamma_accum_n,
     &                             dt,tolNewton,itermax)
        else
            write(*,*) 'solution is elastic'

            eps_el = eps_el_n + deps_el
            gamma = gamma_n + dgamma
            gamma_accum = gamma_accum_n + dabs(dgamma)
            ! compute updated stress
            call tens4_prod_tens2_f(cauchy,cEl,eps_el)

        end if
        ! write out data
        write(150,150) eps(1,2),cauchy(1,2)
        call flush(150)
        eps_n = eps
        gamma_n = gamma
        gamma_accum_n = gamma_accum
        eps_el_n = eps_el
        write(*,*) 'cauchy final = ',cauchy
      enddo


      end program

      subroutine read_file (UnitNum, FileName, NumRows, NumCols, Array )      

          integer, intent (in) :: UnitNum
          character (len=*), intent (in) :: FileName
          integer, intent (in) :: NumRows, NumCols
          real*8, dimension(1:NumRows, 1:NumCols),intent (out) :: Array      

          character (len=300) :: line
          integer :: i, j      

          open (unit=UnitNum,file=FileName,status='old',action='read')      

          ReadComments: do
             read (UnitNum, '(A)') line
             if (line (1:1) /= "#") exit ReadComments
          end do ReadComments      

          backspace (UnitNum)      

          do i=1,NumRows
             read (UnitNum, *) (Array (i, j), j=1,NumCols)
          end do      

          close (UnitNum)      

          return      

       end subroutine read_file