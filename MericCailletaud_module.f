      module MericCailletaud_module

      use umat_subroutines
      implicit none
      
      integer,parameter :: nslip = 12
      real*8  eta,nexp,tau0
      real*8  schmidt(nslip,3,3)
      real*8  nschmidt(nslip,3),sschmidt(nslip,3),mathard(nslip,nslip)
      real*8  cEl(3,3,3,3),delta2_schmidt(nslip,nslip)
      real*8  cEl_T_x_schmidt(nslip,3,3)
      
      contains

      subroutine init_schmidt_f()
          integer k,i,j
          ! schmidt tensors
          do k=1,nslip
              schmidt(k,:,:) = 0.d0
              do i=1,3
                  do j=1,3
                    schmidt(k,i,j) = .5d0*(nschmidt(k,i)*sschmidt(k,j)+
     &                          sschmidt(k,i)*nschmidt(k,j))
                  enddo
              enddo
          enddo
      
          ! delta2-schmidt
          delta2_schmidt = 0.d0
          do k=1,nslip
              delta2_schmidt(k,k) = 1.d0
          enddo

          do k=1,nslip
            call tens4_prod_tens2_f(cEl_T_x_schmidt(k,:,:),
     &                                  cEl,schmidt(k,:,:))
          enddo

      end subroutine
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine solve_nonlinear_equations_f(deps_el,dgamma,deps,
     &          eps_el_n,gamma_n,gamma_accum_n,dt,tol,itermax)
      
      use math_subroutines

      implicit none
      
      integer iter,itermax,k,kp,indx(nslip)
      real*8 deps_el(3,3),dgamma(nslip)
      real*8 gamma_accum_n(nslip)
      real*8 eps_el_n(3,3),gamma_n(nslip),tol,res,dt
      real*8 dgamma_i(nslip),deps_el_i(3,3),deps(3,3)
      real*8 f_epsel_i(3,3),f_gamma_i(nslip),jac(nslip,nslip)
      real*8 rhs(nslip),gamma_accum(nslip)
      real*8 f_epsel(3,3),f_gamma(nslip),gamma_accum_i(nslip)
      real*8 dfgamma_ddgamma_i,dfgamma_ddeps_el_i(3,3)
      real*8 delta_dgamma(nslip)
      real*8 deps_vp(3,3),eps_el_i(3,3),cauchy_i(3,3),tau_i(nslip)
      real*8 taucrss_i(nslip),fyield_i(nslip)
      real*8 dgammadot_dfyield_i(nslip)
      real*8 dRhard_ddgamma_i(nslip,nslip),jacInv(nslip,nslip)
      real*8 aux
      
      ! set the starting point for newton iterations
      dgamma_i = 0.d0
      deps_el_i = deps
      ! initialize variables
      res = 1.d0
      ! start the newton iteration
      iter = 0
      do while(dabs(res) .gt. dabs(tol))
          iter = iter+1
        ! if iterations are greater than max number of iterations, exit the loop
          if (iter .gt. itermax) then
              write(*,*) 'max iterations reached, exiting'
              exit
          endif
          ! update the strain-rates
          ! Update the accumulated plastic slips
          gamma_accum_i = dabs(gamma_accum_n) + dabs(dgamma_i)
          ! Evaluate the nonlinear equations
          call eval_incremental_equations_f(f_epsel_i,f_gamma_i,
     &          deps_el_i,deps,dgamma_i,gamma_accum_i,dt)
          ! Construct matrix and rhs vector
          eps_el_i = deps_el_i + eps_el_n
          call tens4_prod_tens2_f(cauchy_i,cEl,eps_el_i)
          do k=1,nslip
              call tens2_inner_product_f(tau_i(k),cauchy_i,
     &                                      schmidt(k,:,:))
          enddo
          call eval_fyield_f(fyield_i,tau_i,taucrss_i)
          call eval_dgammadot_dfyield_f(dgammadot_dfyield_i,fyield_i)
          call eval_dRhard_ddgamma_f(dRhard_ddgamma_i,dgamma_i)
      
          rhs = 0.d0
          jac = 0.d0
          aux = 0.d0
          do k=1,nslip
      
              dfgamma_ddeps_el_i = (-dt*dgammadot_dfyield_i(k)*
     &                   dsign(1.d0,tau_i(k)))*cEl_T_x_schmidt(k,:,:)
      
              do kp=1,nslip
      
                  dfgamma_ddgamma_i = delta2_schmidt(k,kp) + 
     &             dt*dgammadot_dfyield_i(k)*dRhard_ddgamma_i(k,kp)
      
                  call tens2_inner_product_f(aux,dfgamma_ddeps_el_i,
     &                                             schmidt(kp,:,:))
                  jac(k,kp) = dfgamma_ddgamma_i - aux
      
              enddo
      
              call tens2_inner_product_f(aux,dfgamma_ddeps_el_i,
     &                                                  f_epsel_i)
              rhs(k) = aux - f_gamma_i(k)
      
          enddo
          ! Solve the linear equation for the change in the increment variation
          call KLUDCMP(jac,nslip,indx)
          call KLUBKSB(jac,nslip,indx,rhs)
          do k=1,nslip
            delta_dgamma(k) = rhs(k)
            dgamma(k) = delta_dgamma(k) + dgamma_i(k)
          enddo
        !   write(*,*) 'rhs = ',rhs
          gamma_accum = gamma_accum_n + dabs(dgamma)
        !   write(*,*) 'gamma_accum = ',gamma_accum
          ! construct the increment in viscoplastic strain
          deps_vp = 0.d0
          do k=1,nslip
              deps_vp = deps_vp + dgamma(k)*schmidt(k,:,:)
          enddo
          ! Compute the new deps_el
          deps_el = deps-deps_vp
          ! Evaluate the incremental equations with the newly updated variables
          f_epsel = 0.d0
          f_gamma = 0.d0
          call eval_incremental_equations_f(f_epsel,f_gamma,
     &          deps_el,deps,dgamma,gamma_accum,dt)
          do k=1,nslip
            write(*,*) 'f_gamma new = ',f_gamma(k)
          enddo
          ! Compute the norm of the residual
          call eval_norm_residual_incremental_equations_f(res,
     &                              f_epsel,f_gamma)
          ! Update the iteration variables
          dgamma_i = dgamma
          deps_el_i = deps_el
        !   res = dabs(f_gamma(2))
          write(*,*) 'res = ',res
      enddo

      end subroutine
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !   subroutine constituve_update_f(deps,eps_el_n,gamma_n,dt,tol,
    !  &                                  itermax)
     
    !      implicit none
         
    !      integer itermax
    !      real*8 deps(3,3),eps_el_n(3,3),gamma_n(nslip),dt,tol
    !      real*8 deps_el(3,3),dgamma(nslip),gamma(nslip),eps_el(3,3)
    !      real*8 cauchy(3,3),gamma_accum_n(nslip)
         
    !      ! solve nonlinear equations
    !      call solve_nonlinear_equations_f(deps_el,dgamma,deps,eps_el_n,
    !  &   gamma_n,gamma_accum_n,dt,tol,itermax)
    !      ! update eps_el and gamma
    !      eps_el = eps_el_n + deps_el
    !      gamma = gamma_n + dgamma
    !      ! update the Cauchy stress
    !      call tens4_prod_tens2_f(cauchy,cEl,eps_el)
    !      ! set previous variables to current
    !      gamma_n = gamma
    !      eps_el_n = eps_el
         
    !   end subroutine     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
      subroutine eval_fyield_f(fyield,tau,taucrss)

          implicit none
      
          integer k
          real*8 fyield(nslip),tau(nslip),taucrss(nslip)
      
          do k=1,nslip
              fyield(k) = dabs(tau(k))-taucrss(k)
          enddo
      
      end subroutine
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine eval_dRhard_ddgamma_f(dRhard_ddgamma,dgamma)
      
          implicit none
      
          integer k,kp
          real*8 dRhard_ddgamma(nslip,nslip),dgamma(nslip)
          real*8 sgn_dgamma(nslip)
      
          do k=1,nslip
              sgn_dgamma(k) = dsign(1.d0,dgamma(k))
          enddo
      
          dRhard_ddgamma = 0.d0
          do k=1,nslip
              do kp=1,nslip
                  dRhard_ddgamma(k,kp) = mathard(k,kp)*sgn_dgamma(kp)
              enddo
          enddo
      
      end subroutine
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine eval_norm_residual_incremental_equations_f(norm,
     &                                              f_epsel,f_gamma)

          implicit none
      
          integer i,j
          real*8 f_epsel(3,3),f_gamma(nslip),norm
      
          norm = 0.d0
          do i=1,3
              do j=i,3
                  norm = norm + (f_epsel(i,j))**2.d0
              enddo
          enddo
      
          do i=1,nslip
              norm = norm + (f_gamma(i))**2.d0
          enddo
      
          norm = dsqrt(norm)
      
      end subroutine
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine eval_incremental_equations_f(f_epsel,f_gamma,
     &      deps_el,deps,dgamma,gamma_accum,dt)

      implicit none
      
      integer k
      real*8 dt,eps_el_n(3,3),gamma_accum(nslip)
      real*8 deps_el(3,3),deps_vp(3,3),gammadot(nslip)
      real*8 deps(3,3)
      real*8 dgamma(nslip)
      real*8 f_epsel(3,3),f_gamma(nslip),fyield(nslip)
      
      ! equation for plastic slip increments
      call eval_gammadot_f(gammadot,fyield,deps_el,dgamma,eps_el_n,
     &                     gamma_accum)
      do k=1,nslip
          f_gamma(k) = dgamma(k) - dt*gammadot(k)
      enddo
      ! equation for elastic strain increment
      deps_vp = 0.d0
      do k=1,nslip
        deps_vp = deps_vp + dgamma(k)*schmidt(k,:,:)
      enddo
      f_epsel = deps_el + deps_vp - deps
      
      end subroutine
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine eval_gammadot_f(gammadot,fyield,deps_el,dgamma,
     &      eps_el_n,gamma_accum)
      
          implicit none
      
          integer k
          real*8 tol
          real*8 eps_el_n(3,3),eps_el(3,3),gamma_accum(nslip)
          real*8 deps_el(3,3),deps_vp(3,3),gammadot(nslip),dgamma(nslip)
          real*8 fyield(nslip),tau(nslip),taucrss(nslip),cauchy(3,3)
          
          ! evaluate TAUCRSS
          call eval_taucrss_f(taucrss,gamma_accum)
          ! set the tolerance for evaluating the plastic strain-rates
          tol = 1.e-12
          eps_el = eps_el_n + deps_el
        !   write(*,*) 'eps_el = ',eps_el
          call tens4_prod_tens2_f(cauchy,cEl,eps_el)
          do k=1,nslip
              call tens2_inner_product_f(tau(k),cauchy,schmidt(k,:,:))
            !   write(*,*) 'tau = ',tau(k)
              fyield(k) = dabs(tau(k)) - taucrss(k)
            !   write(*,*) 'fyield = ',fyield(k)
              if (fyield(k) .le. tol) then
                  gammadot(k) = 0.d0
              else
                  gammadot(k) = (fyield(k)/eta)**(nexp)
              endif
              write(*,*) 'gammadot = ',gammadot(k)
          enddo
      
      end subroutine
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine eval_dgammadot_dfyield_f(dgammadot_dfyield,fyield)
      
          implicit none
      
          integer k
          real*8 dgammadot_dfyield(nslip),fyield(nslip),tol
      
          tol = 1.e-12
          do k=1,nslip
              if (fyield(k) .le. tol) then
                  dgammadot_dfyield(k) = 0.d0
              else
                  dgammadot_dfyield(k) = (nexp/eta)*(
     &                             (fyield(k)/eta)**(nexp-1.d0))

              endif
          enddo
      
      end subroutine
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine eval_taucrss_f(taucrss,gamma_accum)
        
          implicit none
      
          integer k,kp
          real*8 taucrss(nslip),gamma_accum(nslip),Rhard(nslip),sum
      
          Rhard = 0.d0
          do k=1,nslip
            sum = 0.d0
            do kp=1,nslip
                sum = sum + mathard(k,kp)*gamma_accum(kp)
            enddo
            Rhard(k) = sum
        enddo
      
          do k=1,nslip
            taucrss(k) = tau0 + Rhard(k)
          enddo
      
      end subroutine

      end module MericCailletaud_module
