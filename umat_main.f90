module global_variables

implicit none

integer,parameter :: nslip = 12
real*8, parameter :: pi = 3.1415926535d0
integer nexp
real*8  eta(nslip)
real*8  delta2(3,3), delta4(3,3,3,3), delta4s(3,3,3,3)
real*8  schmidt(nslip,3,3),tau0(nslip)
real*8  nschmidt(nslip,3),sschmidt(nslip,3),mathard(nslip,nslip)
real*8  cel(3,3,3,3)

contains
subroutine init_schmidt_f()
    integer k,i,j
    do k=1,nslip
        schmidt(k,:,:) = 0.
        do i=1,3
            do j=1,3
                schmidt(k,i,j) = .5*(nschmidt(k,i)*sschmidt(k,j)+sschmidt(k,i)*nschmidt(k,j))
            enddo
        enddo
    enddo
end subroutine

end module global_variables

subroutine solve_nonlinear_equations_f()


end subroutine

subroutine update_nonlinear_equations_f(fdepsel,fdgamma,deps,epsel_n,p_n)

use global_variables
implicit none

integer m,k
real*8 dt
real*8 epsel(3,3),epsel_n(3,3),depsel(3,3)
real*8 sig(3,3),depsv(3,3),deps(3,3)
real*8 tau(nslip),dgamma(nslip),dp(nslip),p(nslip),p_n(nslip),r(nslip),&
gdot(nslip),f(nslip)
real*8 fdepsel(3,3),fdgamma(nslip)

! update epsel
epsel = epsel_n + depsel
! update stress
call mat_x_vec_f(sig,cel,epsel)
! update the resolved shear stresses on each slip system
do k=1,m
    call tens2_dot_tens2_f(tau(k),sig,schmidt(k,:,:))
enddo
! update accumulated plastic strain
dp = dabs(dgamma)
p = p_n + dp
! update hardening stress
call shard_update_f(r,p)
! update fk
do k=1,m
    f(k) = dabs(tau(k))-r(k)-tau0(k)
enddo
! update gammadot
call gdot_update_f(gdot,f,tau)
! update the viscous strain
depsv = 0.
do k=1,m
    depsv = depsv + dgamma(k)*schmidt(k,:,:)
enddo
! update nonlinear equations
fdepsel = depsel + depsv - deps
do k=1,m
    fdgamma(k) = dgamma(k) - dt*gdot(k)
enddo

end subroutine

subroutine shard_update_f(r,p)

use global_variables
implicit none

real*8 r(nslip),p(nslip)

call mat_x_vec_f(r,mathard,p)

end subroutine

subroutine gdot_update_f(gdot,f,tau)

use global_variables
implicit none

integer k
real*8 gdot(nslip),f(nslip),tau(nslip)

do k=1,nslip
    if (f(k) .lt. 0.) then
        gdot(k) = 0.
    else
        gdot(k) = sign(tau(k),1.d0)*(f(k)/eta(k))**(nexp)
    endif
enddo

end subroutine