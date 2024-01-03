      module umat_subroutines

        implicit none
        
        real*8, parameter :: pi = 3.1415926535d0
        real*8  delta2(3,3), delta4(3,3,3,3), delta4s(3,3,3,3)
        real*8  hydro4(3,3,3,3), dev4(3,3,3,3)

        contains

        subroutine init_delta_f()
            implicit none
            integer i,j,k,l
            
            delta2 = 0.d0
            delta4 = 0.d0
            delta4s = 0.d0
            do i=1,3
                  delta2(i,i) = 1.d0
            enddo
            do i=1,3
                  do j=1,3
                        do k=1,3
                              do l=1,3
                              delta4(i,j,k,l) = delta2(i,k)*delta2(k,l)
                              delta4s(i,j,k,l) = 0.5d0*(
     &                               delta2(i,k)*delta2(j,l) +
     &                               delta2(i,l)*delta2(j,k))
                              hydro4(i,j,k,l) = 
     &                              (1.d0/3.d0)*delta2(i,j)*delta2(k,l)
                              enddo
                        enddo
                  enddo
            enddo
            dev4 = delta4s - hydro4

        end subroutine
C ---------------------------------------------------------------------
      subroutine pk_mooney_rivlin_compressible_f(pk,f,f_inv_T,
     &detf,matprops)

      implicit none

      real*8 matprops(3), f(3,3), pk(3,3), f_inv(3,3), f_inv_T(3,3)
      real*8 detf, c1, c2, kappa
      real*8 rcg(3,3), rcg_sqrd(3,3), iden2(3,3)
      real*8 f_T(3,3), tr_rcg, tr_rcg_sqrd, i1, i2, i1Bar, i2Bar
      real*8 lcg(3,3),lcg_x_f(3,3),cauchy(3,3),lcg_x_lcg(3,3)
      real*8 detfMinusThird,detfMinusTwoThird,detfMinusFourThird
      real*8 di1_df(3,3),di2_df(3,3),di1bar_df(3,3),di2bar_df(3,3)
      real*8 f_x_rcg(3,3),b_x_f(3,3),b(3,3),d1,bbar(3,3),detb
      real*8 detbMinusThird,i1bbar,b_sqrd(3,3),tr_b_sqrd
      integer i,j
C     Store matprops
      c1 = matprops(1)
      c2 = matprops(2)
      d1 = matprops(3)
C     Compute jacobian factors
      detfMinusThird = detf**(-1.d0/3.d0)
      detfMinusTwoThird = detfMinusThird**2.d0
      detfMinusFourThird = detfMinusTwoThird**2.d0
C     Compute auxillary tensors to compute invariants
      call tens2_transpose_f(f_T,f)
      call tens2_prod_tens2_f(rcg,f_T,f)
      call tens2_prod_tens2_f(rcg_sqrd,rcg,rcg)
      call tens2_prod_tens2_f(b,f,f_T)
      detb = det_tens2_f(b)
      detbMinusThird = detb**(-1.d0/3.d0)
      bbar = detbMinusThird*b
      i1 = 0.d0
      do i=1,3
            i1 = i1 + b(i,i)
      enddo
      i2 = 0.d0
      call tens2_prod_tens2_f(b_sqrd,b,b)
      tr_b_sqrd = 0.d0
      do i=1,3
            tr_b_sqrd = tr_b_sqrd + b_sqrd(i,i)
      enddo
      i2 = .5d0*(i1**2.d0 - tr_b_sqrd)
      i1Bar = detfMinusTwoThird*i1
      i2Bar = detfMinusFourThird*i2
      ! tr_rcg = 0.
      ! tr_rcg_sqrd = 0.
      ! do i=1,3
      !       tr_rcg = tr_rcg + rcg(i,i)
      !       tr_rcg_sqrd = tr_rcg_sqrd + rcg_sqrd(i,i)
      ! enddo
      ! i1 = tr_rcg
      ! i2 = .5d0*(i1**2.d0 - tr_rcg_sqrd)
      ! i1Bar = detfMinusTwoThird*i1
      ! i2Bar = detfMinusFourThird*i2
C     Compute b*f
      call tens2_prod_tens2_f(b_x_f,b,f)
C     Compute PK
      pk = 2.d0*(detfMinusTwoThird*(c1+i1Bar*c2)*f - 
     &      detfMinusFourThird*c2*b_x_f) + 
     &  detf*(2.d0*(detf-1.d0)/d1 - 
     & (2.d0/3.d0/detf)*(c1*i1Bar+2.d0*c2*i2Bar))*f_inv_T
     
      return
      end subroutine
C ---------------------------------------------------------------------
      subroutine pk_neo_hookean_compressible_f(pk,f,f_inv_T,
     &detf,matprops)

      implicit none
      real*8 matprops(2), f(3,3), pk(3,3), f_inv(3,3), f_inv_T(3,3)
      real*8 detf, c10, d1, pCauchy, pk_dev(3,3), pk_vol(3,3)
      real*8 pk_dev_bar(3,3),f_bar(3,3),f_T(3,3),c(3,3)
      real*8 detfMinusThird,detfMinusTwoThird,trc,i1,i1bar
      integer i,j

C     Store matprops
      c10 = matprops(1)
      d1 = matprops(2)

      ! pk = 2.d0*c10*(f-f_inv_T)+(2.d0*detf/d1)*(detf-1.d0)*f_inv_T
C     
      call tens2_transpose_f(f_T,f)
      call tens2_inv_f(f_inv_T,f_T)
      call tens2_prod_tens2_f(c,f_T,f)
      detfMinusThird = detf**(-1.d0/3.d0)
      detfMinusTwoThird = detfMinusThird**2.d0

      trc = c(1,1)+c(2,2)+c(3,3)
      i1 = trc
      i1bar = detfMinusTwoThird*i1

      pk = c10*((-2.d0/3.d0)*i1bar*f_inv_T+detfMinusTwoThird*2.d0*f) + 
     & (2.d0/d1)*detf*(detf-1.d0)*f_inv_T

      end subroutine
! C ---------------------------------------------------------------------
!       subroutine pk_compressible_f(pk,f,f_inv_T,
!      &detf,matprops)

!       

!       implicit none
!       real*8 matprops(nprops_he), f(3,3), pk(3,3), f_inv_T(3,3), detf

!       if (heOpt .eq. 1) then
!       call pk_neo_hookean_compressible_f(pk,f,f_inv_T,detf,matprops)
!       elseif (heOpt .eq. 2) then
!       call pk_mooney_rivlin_compressible_f(pk,f,f_inv_T,detf,matprops)
!       else
!       call pk_neo_hookean_compressible_f(pk,f,f_inv_T,detf,matprops)
!       endif

!       end subroutine     
C ---------------------------------------------------------------------
      subroutine dPdF_neo_hookean_compressible_f(dPdF,matprops,
     &f_inv,detf)
      

      implicit none
      real*8 matprops(2), f_inv(3,3)
      real*8 detf, c10, d1
      real*8 dPdF(3,3,3,3), finv_prod_finv(3,3,3,3)
      real*8 finv_dyad_finv(3,3,3,3)
      integer i,j,k,l

      c10 = matprops(1)
      d1 = matprops(2)
C     Construct auxillary tensors
	do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
       finv_dyad_finv(i,j,k,l) = f_inv(j,i)*f_inv(l,k)
       finv_prod_finv(i,j,k,l) = .5d0*(f_inv(l,i)*f_inv(j,k)+
     &                                 f_inv(j,k)*f_inv(l,i))
      enddo
      enddo
      enddo
      enddo
C     NOTE: The expression below has major symmetry
      dPdF = 2.d0*c10*delta4 
     &+ (2.d0*c10-(2.d0*detf/d1)*(detf-1.d0))*finv_prod_finv
     &+ (2.d0*detf/d1)*(2.d0*detf-1.d0)*finv_dyad_finv
      
      end subroutine
C ---------------------------------------------------------------------
      subroutine pk_dev_f(pk_dev,pk_dev_bar,f,f_inv_T,detF)

      implicit none
      real*8 f(3,3), pk(3,3), f_inv(3,3), f_inv_T(3,3)
      real*8 detF, detfMinusThird,pk_dev(3,3)
      real*8 scalarProd, pk_dev_bar(3,3)
      integer i,j

      detfMinusThird = detf**(-1.d0/3.d0)

      pk_dev = 0.d0
      scalarProd = 0.d0
      do i=1,3
      do j=1,3
            scalarProd = scalarProd + pk_dev_bar(i,j)*f(i,j)
      enddo
      enddo

      pk_dev = detfMinusThird*(pk_dev_bar - (scalarProd/3.)*f_inv_T)

      end subroutine
C ---------------------------------------------------------------------
      subroutine pk_vol_f(pk_vol,pCauchy,detf,f_inv_T)

      implicit none
      real*8 pk_vol(3,3),detf,f_inv_T(3,3),pCauchy

      pk_vol = pCauchy*detf*f_inv_T

      end subroutine
C ---------------------------------------------------------------------
      subroutine pk2_vol_f(pk2_vol,pCauchy,detf,c_inv)

      implicit none
      real*8 pk2_vol(3,3),detf,c_inv(3,3),pCauchy

      pk2_vol = pCauchy*detf*c_inv

      end subroutine
C ---------------------------------------------------------------------
      subroutine tens2_inner_product_f(val,tens1,tens2)
	
	integer i,j
	real*8 val,tens1(3,3),tens2(3,3)
	
	val = 0.d0
	do i=1,3
		  do j=1,3
				 val = val + tens1(i,j)*tens2(i,j)
		  enddo
	enddo
	
	end
C ---------------------------------------------------------------------
      subroutine pk2_dev_f(pk2_dev,pk2bar_dev,detF,c,c_inv)

      

      implicit none
      real*8 detF,pk2bar_dev_dot_c,detfMinusThird,detfMinusTwoThird
      real*8 pk2_dev(3,3),pk2bar_dev(3,3),c(3,3),
     &                 c_inv(3,3)     
c Construct fourth-order deviatoric projection tensor
c Compute trace of product pk2bar_dev and c
      detfMinusThird = (detf)**(-1.d0/3.d0)
      detfMinusTwoThird = detfMinusThird**2.d0
      call tens2_inner_product_f(pk2bar_dev_dot_c,pk2bar_dev,c)
      pk2_dev = pk2bar_dev - (1.d0/3.d0)*pk2bar_dev_dot_c*c_inv
      pk2_dev = pk2_dev*detfMinusTwoThird   
      
      end subroutine
C ---------------------------------------------------------------------
      subroutine tens_dyadic_product_f(tens_out,tens1,tens2)

	integer i,j,k,l
	double precision tens1(3,3),tens2(3,3),tens_out(3,3,3,3)
	
	do i=1,3
	do j=1,3
	do k=1,3
	do l=1,3
	  tens_out(i,j,k,l) = tens1(i,j)*tens2(k,l)
	enddo
	enddo
	enddo
	enddo
	
      end subroutine
C ---------------------------------------------------------------------
	subroutine tens_circ_product_f(tens_out,tens1,tens2)
	

	integer i,j,k,l
	double precision tens1(3,3),tens2(3,3),tens_out(3,3,3,3)
	
	do i=1,3
	do j=1,3
	do k=1,3
	do l=1,3
	  tens_out(i,j,k,l) = .5d0*(tens1(i,k)*tens2(j,l)+
     &                        tens1(i,l)*tens2(j,k))
	enddo
	enddo
	enddo
	enddo
	
      end subroutine
C ---------------------------------------------------------------------
      subroutine tangent_material_vol_f(tangent_vol,detf,pCauchy,
     &dpCauchydj,c_inv)    
	
	integer i,j,k,l
	
	real*8 detf,pCauchy,dpCauchydj
	real*8 tangent_vol(3,3,3,3),c_inv(3,3)
      real*8 c_inv_circ_c_inv(3,3,3,3),c_inv_dyad_c_inv(3,3,3,3)
	
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
            c_inv_circ_c_inv(i,j,k,l) = .5*(c_inv(i,k)*c_inv(j,l)
     &            + c_inv(i,l)*c_inv(j,k))
            c_inv_dyad_c_inv(i,j,k,l) = c_inv(i,j)*c_inv(k,l)
      enddo
      enddo
      enddo
      enddo

      tangent_vol = (detf*pCauchy+(detf**2.)*dpCauchydj)*
     &c_inv_dyad_c_inv -2.d0*detf*pCauchy*c_inv_circ_c_inv
	
	end subroutine 
C ---------------------------------------------------------------------
      subroutine tangent_material_dev_f(tangent_dev,cbarDev,pk2bar_dev,
     &                                    pk2_dev,c,c_inv,detF)
     
      
      
      integer i,j,k,l
      
      real*8           detf,detfMinusTwoThird,
     &                 pk2bar_dev(3,3),pk2_dev(3,3),detfMinusThird,
     &                 tangent_dev(3,3,3,3),c(3,3),c_inv(3,3),
     &                 pBar(3,3,3,3),projDevRef(3,3,3,3),
     &                 projDevRef_T(3,3,3,3),term1(3,3,3,3),
     &                 term2(3,3,3,3),term3(3,3,3,3),aux1(3,3,3,3),
     &                 cbarDev(3,3,3,3),c_inv_circ_c_inv(3,3,3,3),
     &                 c_inv_dyad_c_inv(3,3,3,3)
      real*8 detfMinusFourThird,auxScalar
     
      detfMinusThird = detf**(-1.d0/3.d0)
      detfMinusTwoThird = detfMinusThird**2.d0
      detfMinusFourThird = detfMinusTwoThird**2.d0
C     Term 1
C     Construct the dev. projection operator
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
            projDevRef(i,j,k,l) = delta4(i,j,k,l)
     &      - (1.d0/3.d0)*c_inv(i,j)*c(k,l)
            projDevRef_T(i,j,k,l) = delta4(k,l,i,j) 
     &      - (1.d0/3.d0)*c_inv(k,l)*c(i,j)
      enddo
      enddo
      enddo
      enddo
      call tens4_prod_tens4_f(aux1,cbarDev,
     &projDevRef_T)
      call tens4_prod_tens4_f(term1,detfMinusFourThird*projDevRef,
     &aux1)
C     Term 2
      call tens2_inner_product_f(auxScalar,pk2bar_dev,c)
C     Construct dyadic products/circle products
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
            c_inv_circ_c_inv(i,j,k,l) = .5*(c_inv(i,k)*c_inv(j,l)
     &            + c_inv(i,l)*c_inv(j,k))
            c_inv_dyad_c_inv(i,j,k,l) = c_inv(i,j)*c_inv(k,l)
      enddo
      enddo
      enddo
      enddo
      pBar = c_inv_circ_c_inv - (1.d0/3.d0)*c_inv_dyad_c_inv
      term2 = (2./3.)*detfMinusTwoThird*auxScalar*pBar
C     Term 3
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
            term3(i,j,k,l) = -(2./3.)*(c_inv(i,j)*pk2_dev(k,l)
     &                        +pk2_dev(i,j)*c_inv(k,l))
      enddo
      enddo
      enddo
      enddo

      tangent_dev = term1 + term2 + term3

      end subroutine
C ---------------------------------------------------------------------
      subroutine material_tangent_f(mat_tangent,mat_tangent_dev_bar,
     &   pk2_dev_bar,pk2_dev,c,c_inv,detf,pCauchy,dpCauchydj)
	
	implicit none
      real*8 mat_tangent(3,3,3,3),mat_tangent_dev(3,3,3,3),
     &mat_tangent_dev_bar(3,3,3,3),mat_tangent_vol(3,3,3,3)
      real*8 pk2_dev_bar(3,3),pk2_dev(3,3),c(3,3),c_inv(3,3)
      real*8 detf,pCauchy,dpCauchydj

C     Call deviatoric tangent subroutine
      call tangent_material_dev_f(mat_tangent_dev,mat_tangent_dev_bar,
     &                          pk2_dev_bar,pk2_dev,c,c_inv,detF)
C     Call volumetric tangent subroutine
      call tangent_material_vol_f(mat_tangent_vol,detF,pCauchy,
     &                          dpCauchydj,c_inv)
C     Compute total tangent
      mat_tangent = mat_tangent_dev + mat_tangent_vol
	
	end subroutine 
C ---------------------------------------------------------------------
      subroutine dpk_df_f(dPdF,f_inv,detF,pkDev_bar,pkDev,pkVol,pCauchy,
     &dpCauchydj,dpkDev_bar_df_bar,f_bar,f)

      

      implicit none
      real*8 matprops(2), dPdF(3,3,3,3), f_bar(3,3), f(3,3)
      real*8 f_inv(3,3), f_inv_T(3,3), dpkDev_bar_df_bar(3,3,3,3)
      real*8 f_inv_T_dyad_f_inv_T(3,3,3,3), f_inv_prod_f_inv(3,3,3,3)
      real*8 dpkDev_df(3,3,3,3), dpkVol_df(3,3,3,3), pkDev_bar(3,3)
      real*8 detF, pCauchy, dpCauchydj, pkDev(3,3), pkVol(3,3)
      integer i,j,k,l

C     Kinematics
      call tens2_transpose_f(f_inv_T,f_inv)

C     Compute dpkDev/df
      call dpkDev_df_f(dpkDev_df,dpkDev_bar_df_bar,
     &pkDev,pkDev_bar,detf,f,f_inv_T,f_bar)
C     Compute dpkVol/df
      call dpkVol_df_f(dpkVol_df,pkVol,detf,pCauchy,dpCauchydj,f_inv_T)
C     Compute dpk/df
      dPdF = dpkDev_df + dpkVol_df
C     Apply major symmetry
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
            dPdF(k,l,i,j) = dPdF(i,j,k,l)
      enddo
      enddo
      enddo
      enddo

      end subroutine
C ---------------------------------------------------------------------
      subroutine dpkVol_df_f(dpkVol_df,pkVol,detf,pCauchy,dpCauchydj,
     &f_inv_T)

      implicit none
      real*8 f_inv_T_dyad_f_inv_T(3,3,3,3),f_inv_T_prod_f_inv_T(3,3,3,3)
      real*8 dpkVol_df(3,3,3,3), detf, pCauchy, dpCauchydj, f_inv_T(3,3)
      real*8 pkVol(3,3)
      integer i,j,k,l

      dpkVol_df = 0.d0
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
            dpkVol_df(i,j,k,l) = pkVol(i,j)*f_inv_T(k,l) 
     & + (detf**2.)*dpCauchydj*f_inv_T(i,j)*f_inv_T(k,l)
     & - pkVol(k,j)*f_inv_T(i,l)
      enddo
      enddo
      enddo
      enddo

      end subroutine
C ---------------------------------------------------------------------
      subroutine dpkDev_df_f(dpkDev_df,dpkDev_bar_df_bar,
     &pkDev,pkDev_bar,detf,f,f_inv_T,f_bar)

      
      implicit none
      real*8 dpkDev_df(3,3,3,3),f(3,3),f_inv_T(3,3),detf,detfMinusThird
      real*8 dpkDev_dpkDev_bar(3,3,3,3),dpkDev_bar_df(3,3,3,3)
      real*8 dpkDev_bar_df_bar(3,3,3,3), auxScalar,detfMinusTwoThird
      real*8 f_bar(3,3),dpkDev_bar_df_bar_prod_f(3,3),pkDev(3,3)
      real*8 auxProd(3,3,3,3),dpkDev_dj(3,3),ddetf_df(3,3)
      real*8 pkDev_bar(3,3),aux1_tens2(3,3),aux2_tens2(3,3),aux3_scalar
      real*8 term3(3,3,3,3),term2(3,3,3,3),term1(3,3,3,3)
      integer i,j,k,l

      detfMinusThird = (detf)**(-1./3.)
      detfMinusTwoThird = detfMinusThird**2.
C     Term 1
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
         term1(i,j,k,l) = (-1./3.)*pkDev(i,j)*f_inv_T(k,l)
      enddo
      enddo
      enddo
      enddo
C     Term 2
      aux1_tens2 = 0.
      aux2_tens2 = 0.
      aux3_scalar = 0.
      do i=1,3
      do j=1,3
            do k=1,3
            do l=1,3
                  aux1_tens2(i,j) = aux1_tens2(i,j) + 
     &                  dpkDev_bar_df_bar(i,j,k,l)*f(k,l)

                  aux2_tens2(i,j) = aux2_tens2(i,j) + 
     &                  dpkDev_bar_df_bar(k,l,i,j)*f(k,l)

                aux3_scalar = aux3_scalar + 
     &                   dpkDev_bar_df_bar(i,j,k,l)*f(i,j)*f(k,l)  
            enddo
            enddo
      enddo
      enddo
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
      term2(i,j,k,l) = detfMinusTwoThird*(dpkDev_bar_df_bar(i,j,k,l) 
     & -(1.d0/3.d0)*aux1_tens2(i,j)*f_inv_T(k,l) 
     & -(1.d0/3.d0)*aux2_tens2(k,l)*f_inv_T(i,j)
     & +(1.d0/9.d0)*aux3_scalar*f_inv_T(i,j)*f_inv_T(k,l))
      enddo
      enddo
      enddo
      enddo
C     Term 3
      auxScalar = 0.
      do i=1,3
      do j=1,3
            auxScalar = auxScalar + pkDev_bar(i,j)*f(i,j)
      enddo
      enddo
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
          term3(i,j,k,l)=(-detfMinusThird/3.)*(
     & f_inv_T(i,j)*pkDev_bar(k,l)-auxScalar*f_inv_T(i,l)*f_inv_T(k,j))
      enddo
      enddo
      enddo
      enddo
C     Compute dpkDev/df
      dpkDev_df = term1 + term2 + term3

      end subroutine
C ---------------------------------------------------------------------
      subroutine mat_tangent_from_dPdF_f(mat_tangent,dPdF,pk2,
     & f_inv,c_inv)

      implicit none
      real*8 mat_tangent(3,3,3,3), dPdF(3,3,3,3), pk2(3,3)
      real*8 auxSum
      real*8 f_inv(3,3), c_inv(3,3)
      integer i,j,k,l,alpha,beta
C     NOTE: The expression below is FULLY symmetrized.
      mat_tangent = 0.d0
      do alpha=1,3
      do j=1,3
      do l=1,3
      do beta=1,3
            auxSum = 0.d0
            do i=1,3
            do k=1,3
                  auxSum = auxSum
     &               + f_inv(beta,i)*f_inv(alpha,k)*dPdF(i,l,k,j)
     &               + f_inv(j,i)*f_inv(l,k)*dPdF(i,alpha,k,beta)
     &               + f_inv(beta,i)*f_inv(j,k)*dPdF(i,l,k,alpha)
     &               + f_inv(alpha,i)*f_inv(l,k)*dPdF(i,j,k,beta)
     &               + f_inv(l,i)*f_inv(alpha,k)*dPdF(i,beta,k,j)
     &               + f_inv(j,i)*f_inv(beta,k)*dPdF(i,alpha,k,l)
     &               + f_inv(l,i)*f_inv(j,k)*dPdF(i,beta,k,alpha)
     &               + f_inv(alpha,i)*f_inv(beta,k)*dPdF(i,j,k,l)                     
            enddo
            enddo
            mat_tangent(alpha,j,l,beta) = (1.d0/8.d0)*(auxSum  
     &                                  - pk2(j,l)*c_inv(beta,alpha)
     &                                  - pk2(beta,alpha)*c_inv(j,l)
     &                                  - pk2(alpha,l)*c_inv(beta,j)
     &                                  - pk2(beta,j)*c_inv(alpha,l)
     &                                  - pk2(j,beta)*c_inv(l,alpha)
     &                                  - pk2(l,alpha)*c_inv(j,beta)
     &                                  - pk2(alpha,beta)*c_inv(l,j)
     &                                  - pk2(l,j)*c_inv(alpha,beta))
      enddo
      enddo
      enddo
      enddo

      end subroutine
C ---------------------------------------------------------------------
      subroutine spatial_tangent_from_dPdF_f(spatial_tangent,dPdF,tau,f)

      

      implicit none
      real*8 spatial_tangent(3,3,3,3), dPdF(3,3,3,3), tau(3,3)
      real*8 auxSum
      real*8 f(3,3)
      integer i,j,k,l,alpha,beta,n,p

      spatial_tangent = 0.d0
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
            auxSum = 0.d0
            do n=1,3
            do p=1,3
                  auxSum = auxSum + dPdF(i,n,l,p)*f(j,n)*
     &                     f(k,p)
            enddo
            enddo
            spatial_tangent(alpha,j,l,beta) = auxSum - 
     &                                    tau(k,j)*delta2(i,l)
      enddo
      enddo
      enddo
      enddo

      end subroutine
C ---------------------------------------------------------------------
      subroutine spatial_tangent_from_mat_tangent_f(spatial_tangent,
     & mat_tangent,f)

      implicit none
      real*8 spatial_tangent(3,3,3,3), mat_tangent(3,3,3,3), f(3,3)
      real*8 auxSum
      integer i,j,k,l,m,n,p,q

      spatial_tangent = 0.d0

      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
            auxSum = 0.d0
            do m=1,3
            do n=1,3
            do p=1,3
            do q=1,3
                  auxSum = auxSum + mat_tangent(m,n,p,q)*f(i,m)*f(j,n)*
     &                                     f(k,p)*f(l,q)
            enddo
            enddo
            enddo
            enddo
            spatial_tangent(i,j,k,l) = auxSum
      enddo
      enddo
      enddo
      enddo

      end subroutine
C ---------------------------------------------------------------------
      subroutine tens2_to_vec_abaqus_f(vec,tens)

      implicit none
      real*8 vec(6)
      real*8 tens(3,3)

      vec = 0.d0

      vec(1) = tens(1,1)
      vec(2) = tens(2,2)
      vec(3) = tens(3,3)
      vec(4) = tens(1,2)
      vec(5) = tens(1,3)
      vec(6) = tens(2,3)

      end subroutine
C ---------------------------------------------------------------------
      subroutine vec_to_tens2_abaqus_f(tens,vec)

      implicit none
      integer i,j
      real*8 vec(6)
      real*8 tens(3,3)

      tens = 0.d0
      
      tens(1,1) = vec(1)
      tens(2,2) = vec(2)
      tens(3,3) = vec(3)
      tens(1,2) = vec(4)
      tens(1,3) = vec(5)
      tens(2,3) = vec(6)

C     Symmetrize
      do i=1,3
            do j=i,3
                  tens(j,i) = tens(i,j)
            enddo
      enddo

      end subroutine
C ---------------------------------------------------------------------
      subroutine tens4_to_mat_abaqus_f(mat,tens)

      implicit none
      real*8 mat(6,6)
      real*8 tens(3,3,3,3)
      integer i,j,k,l

      mat = 0.d0

      mat(1,1) = tens(1,1,1,1)
      mat(2,2) = tens(2,2,2,2)
      mat(3,3) = tens(3,3,3,3)
      mat(4,4) = tens(1,2,1,2)
      mat(5,5) = tens(1,3,1,3)
      mat(6,6) = tens(2,3,2,3)

      mat(1,2) = tens(1,1,2,2)
      mat(1,3) = tens(1,1,3,3)
      mat(1,4) = tens(1,1,1,2)
      mat(1,5) = tens(1,1,1,3)
      mat(1,6) = tens(1,1,2,3)

      mat(2,3) = tens(2,2,3,3)
      mat(2,4) = tens(2,2,1,2)
      mat(2,5) = tens(2,2,1,3)
      mat(2,6) = tens(2,2,2,3)

      mat(3,4) = tens(3,3,1,2)
      mat(3,5) = tens(3,3,1,3)
      mat(3,6) = tens(3,3,2,3)

      mat(4,5) = tens(1,2,1,3)
      mat(4,6) = tens(1,2,2,3)
      
      mat(5,6) = tens(1,3,2,3)

C     Symmetrize the matrix
      do i=1,6
            do j=i,6
                  mat(j,i) = mat(i,j)
            enddo
      enddo

      end subroutine
C ---------------------------------------------------------------------
      subroutine tens2_transpose_f(TensorT,Tensor)

      implicit none
      real*8 TensorT(3,3),Tensor(3,3) 

      integer i,j       

      do i=1,3 
      do j=1,3 
      TensorT(i,j)=Tensor(j,i) 
      end do 
      end do 

      end subroutine
C ---------------------------------------------------------------------
	subroutine tens2_prod_vec_f(vecout,tens,vecin)   
	
      implicit none
	real*8 tens(3,3),vecin(3),vecout(3)
	
	integer i,j,k
	
	do i=1,3
            vecout(i) = 0.d0
            do j=1,3
                  vecout(i) = vecout(i) + tens(i,j)*vecin(j)
            end do
	end do

      END subroutine
C ---------------------------------------------------------------------
	subroutine tens2_prod_tens2_f(TensRes,Tens1,Tens2)   
	
      implicit none
	real*8 TensRes(3,3),Tens1(3,3),Tens2(3,3)
	
	integer i,j,k
	
	do i=1,3
	do j=1,3
            TensRes(i,j)=0.d0
            do k=1,3
            TensRes(i,j)=TensRes(i,j)+Tens1(i,k)*Tens2(k,j)
            end do
	end do
	end do

      END subroutine
C ---------------------------------------------------------------------
	subroutine tens2_prod_tens2_prod_tens2_f(TensRes,Tens1,
     &Tens2,Tens3)   
	
      implicit none
	real*8 TensRes1(3,3),Tens1(3,3),Tens2(3,3),Tens3(3,3),TensRes(3,3)
	
	integer i,j,k
	
	do i=1,3
	do j=1,3
            TensRes1(i,j)=0.0d0
            do k=1,3
            TensRes1(i,j)=TensRes1(i,j)+Tens1(i,k)*Tens2(k,j)
            end do
	end do
	end do

	do i=1,3
	do j=1,3
            TensRes(i,j)=0.0d0
            do k=1,3
            TensRes(i,j)=TensRes(i,j)+TensRes1(i,k)*Tens3(k,j)
            end do
	end do
	end do

	END
C ---------------------------------------------------------------------
	subroutine  tens4_prod_tens2_f(Tres,T4,T2)
	
      implicit none
	real*8 T4(3,3,3,3),T2(3,3),Tres(3,3)
	integer i,j,k,l
	
	do i=1,3
	do j=1,3
		  Tres(i,j)=0.0d0
		  do k=1,3
		  do l=1,3
				 Tres(i,j)=Tres(i,j)+T4(i,j,k,l)*T2(k,l)
		  end do
		  end do
	end do
	end do
	END
C --------------------------------------------------------------------
!       subroutine mat_inv_lapack_f(matInv,mat,dimInv)

!       implicit none
!       integer dimInv
!       integer ipiv(dimInv),info,lda,lwork
!       real*8 mat(dimInv,dimInv),matInv(dimInv,dimInv),
!      &       Amat(dimInv,dimInv),work(dimInv)
! C     Define variables associated with LAPACK subroutine
!       external dgetri, dgetrf

!       lda = dimInv
!       lwork = dimInv
! C     Compute inverse using LAPACK
!       Amat = mat
! C     Compute LU factorization
!       call dgetrf(dimInv,dimInv,Amat,lda,ipiv,info)
!       if (info .lt. 0) then
!             write(*,*) 'LU factorization has NAN'
!       elseif (info .gt. 0) then
!             write(*,*) 'LU Factorization is singular'
!       endif
! C     Compute inverse of LU factorization
!       call dgetri(dimInv,Amat,lda,ipiv,work,lwork,info)
!       if (info .ne. 0) then
!             write(*,*) 'matrix inv error'
!       endif
! C     Convert matrix back to fourth-order tensor
!       matInv = Amat

!       return
!       end subroutine
C ---------------------------------------------------------------------
      subroutine matmul_f(res,a,b,dim)

      implicit none
      integer dim
      integer i,j,k,l
      real*8 res(dim,dim),a(dim,dim),b(dim,dim)

      res = 0.d0
      do i=1,dim
      do j=1,dim
            do k=1,dim
                  res(i,j) = res(i,j) + a(i,k)*b(k,j)
            enddo
      enddo
      enddo

      end subroutine
C ---------------------------------------------------------------------
      subroutine  tens4_prod_tens4_f(Tres,T1,T2)
  	
      implicit none

      real*8 Tres(3,3,3,3),T1(3,3,3,3),T2(3,3,3,3)
      integer i,j,k,l,m,n
      
      Tres = 0.d0
      do i=1,3
      do j=1,3
            do k=1,3
            do l=1,3
            Tres(i,j,k,l) = 0.d0
                  do m=1,3
                  do n=1,3
                        Tres(i,j,k,l) = Tres(i,j,k,l) + 
     &                        T1(i,j,m,n)*T2(m,n,k,l)
                  enddo
                  enddo
            end do
            end do
      enddo
      enddo
  
      END
C ---------------------------------------------------------------------
	subroutine tens2_inv_f(TensRes,Tensr)
	
      implicit none
	real*8 TensRes(3,3),Tensr(3,3),det
	
	integer i,j
	
      TensRes = 0.d0
	TensRes(1,1)=Tensr(2,2)*Tensr(3,3)-Tensr(2,3)*Tensr(3,2)
	TensRes(1,2)=Tensr(1,3)*Tensr(3,2)-Tensr(3,3)*Tensr(1,2)
	TensRes(1,3)=Tensr(1,2)*Tensr(2,3)-Tensr(2,2)*Tensr(1,3)
	TensRes(2,1)=Tensr(2,3)*Tensr(3,1)-Tensr(3,3)*Tensr(2,1)
	TensRes(2,2)=Tensr(1,1)*Tensr(3,3)-Tensr(1,3)*Tensr(3,1)
	TensRes(2,3)=Tensr(1,3)*Tensr(2,1)-Tensr(2,3)*Tensr(1,1)
	TensRes(3,1)=Tensr(2,1)*Tensr(3,2)-Tensr(2,2)*Tensr(3,1)
	TensRes(3,2)=Tensr(1,2)*Tensr(3,1)-Tensr(1,1)*Tensr(3,2)
	TensRes(3,3)=Tensr(2,2)*Tensr(1,1)-Tensr(2,1)*Tensr(1,2)
	
	det=det_tens2_f(Tensr)
      if (det /= det) then
            write(*,*) 'det = ',det
      endif
	
      tensres = tensres/det
	
	return
	
	end
C ---------------------------------------------------------------------
      subroutine tens2_to_vec_f(vec,tens)

      implicit none

      real*8 tens(3,3), vec(9)
      integer i,j,k,l,p,q

      vec = 0.d0
      p = 1
      do i=1,3
      do j=1,3
            vec(p) = tens(i,j)
            p = p+1
      enddo
      enddo

      end subroutine
C ---------------------------------------------------------------------
      subroutine vec_to_tens2_f(tens,vec)

      implicit none

      real*8 tens(3,3), vec(9)
      integer i,j,k,l,p,q

      tens = 0.d0
      p = 1
      do i=1,3
      do j=1,3
            tens(i,j) = vec(p)
            p = p+1
      enddo
      enddo

      end subroutine
C ---------------------------------------------------------------------
      subroutine tens4_to_mat_MandelBasis_f(mat,tens)

      implicit none

      real*8 tens(3,3,3,3),mat(9,9),vec(9)
      real*8 b(9,3,3), comp
      integer i,j,k,l,p,q,opt
C     Define the basis tensors
      b = 0.d0
      b(1,1,1) = 1.d0
      b(2,2,2) = 1.d0
      b(3,3,3) = 1.d0
      b(4,2,3) = 1.d0/dsqrt(2.d0)
      b(4,3,2) = b(4,2,3)
      b(5,1,3) = 1.d0/dsqrt(2.d0)
      b(5,3,1) = b(5,1,3)
      b(6,1,2) = 1.d0/dsqrt(2.d0)
      b(6,2,1) = b(6,1,2)
      b(7,2,3) = -1.d0/dsqrt(2.d0)
      b(7,3,2) = 1.d0/dsqrt(2.d0)
      b(8,1,3) = 1.d0/dsqrt(2.d0)
      b(8,3,1) = -1.d0/dsqrt(2.d0)
      b(9,1,2) = -1.d0/dsqrt(2.d0)
      b(9,2,1) = 1.d0/dsqrt(2.d0)
C     Project the fourth-order tensor onto Mandel basis and store components
      mat = 0.d0
      do i=1,9
      do j=1,9
            comp = 0.d0
            call quadratic_product_tens4_f(comp,tens,b(i,:,:),b(j,:,:))
            mat(i,j) = comp
      enddo
      enddo
                  
      end subroutine
C ---------------------------------------------------------------------
      subroutine mat_to_tens4_MandelBasis_f(tens,mat)

      implicit none

      real*8 tens(3,3,3,3),mat(9,9),vec(9)
      real*8 b(9,3,3), b4(9,3,3,3,3), aux(3,3,3,3)
      integer i,j,k,l,p,q,opt
C     Define the basis tensors
      b = 0.d0
      b(1,1,1) = 1.d0
      b(2,2,2) = 1.d0
      b(3,3,3) = 1.d0
      b(4,2,3) = 1.d0/dsqrt(2.d0)
      b(4,3,2) = b(4,2,3)
      b(5,1,3) = 1.d0/dsqrt(2.d0)
      b(5,3,1) = b(5,1,3)
      b(6,1,2) = 1.d0/dsqrt(2.d0)
      b(6,2,1) = b(6,1,2)
      b(7,2,3) = -1.d0/dsqrt(2.d0)
      b(7,3,2) = 1.d0/dsqrt(2.d0)
      b(8,1,3) = 1.d0/dsqrt(2.d0)
      b(8,3,1) = -1.d0/dsqrt(2.d0)
      b(9,1,2) = -1.d0/dsqrt(2.d0)
      b(9,2,1) = 1.d0/dsqrt(2.d0)
C     Convert matrix back to fourth-order tensor
      tens = 0.d0
      do p=1,9
      do q=1,9
            do i=1,3
            do j=1,3
            do k=1,3
            do l=1,3
                  aux(i,j,k,l) = b(p,i,j)*b(q,k,l)
            enddo
            enddo
            enddo
            enddo
            tens = tens + mat(p,q)*aux
      enddo
      enddo
                  
      end subroutine
C ---------------------------------------------------------------------
      subroutine quadratic_product_tens4_f(val,a,b,c)

      implicit none
      real*8 val,a(3,3,3,3),b(3,3),c(3,3)
      integer i,j,k,l

      val = 0.d0
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
            val = val + c(i,j)*a(i,j,k,l)*b(k,l)
      enddo
      enddo
      enddo
      enddo

      end subroutine
C ---------------------------------------------------------------------
      subroutine mat_to_tens4_f(tens,mat)

      implicit none

      real*8 tens(3,3,3,3), mat(9,9)
      integer i,j,k,l,p,q

      tens = 0.d0
      p = 1
      q = 1
      do i=1,3
      do j=1,3
            q = 1
            do k=1,3
            do l=1,3
                  tens(i,j,k,l) = mat(p,q)
                  q = q+1
            enddo
            enddo
            p = p+1
      enddo
      enddo

      end subroutine
C ---------------------------------------------------------------------
	FUNCTION det_tens2_f(T)
      
      implicit none
	
	real*8 det_tens2_f,T(3,3)
	
	det_tens2_f=T(1,1)*T(2,2)*T(3,3)+T(2,1)*T(3,2)*T(1,3) + 
     &  T(1,2)*T(2,3)*T(3,1)-T(1,2)*T(2,1)*T(3,3) - 
     &  T(2,3)*T(3,2)*T(1,1)-T(3,1)*T(2,2)*T(1,3)
	
	RETURN
	
	END
C ---------------------------------------------------------------------
      subroutine rot_tens2_f(tensRot,tens,rot,opt)

      implicit none
      integer i,j,m,n,opt
      real*8 tensRot(3,3),tens(3,3),rot(3,3)
C     Opt = 1 --- rotate to new config
C     Opt = 2 --- un-rotate to prev. config
      tensRot = 0.d0
      if (opt .eq. 1) then
            do i=1,3
            do j=1,3
                  do m=1,3
                  do n=1,3
                        tensRot(i,j) = tensRot(i,j) + 
     &                  rot(i,m)*rot(j,n)*tens(m,n)
                  enddo
                  enddo
            enddo
            enddo
      elseif (opt .eq. 2) then
            do i=1,3
            do j=1,3
                  do m=1,3
                  do n=1,3
                        tensRot(i,j) = tensRot(i,j) + 
     &                  rot(m,i)*rot(n,j)*tens(m,n)
                  enddo
                  enddo
            enddo
            enddo
      endif

      end subroutine
C ---------------------------------------------------------------------
      subroutine vec_norm_f(xNorm,x)

      implicit none

      real*8 x(3),xNorm

      xNorm = dsqrt(x(1)**2.d0+x(2)**2.d0+x(3)**2.d0)

      end subroutine
C ---------------------------------------------------------------------
      subroutine vec_dot_vec_f(val,x,y,ndim)

      implicit none

      integer ndim,i
      real*8 x(ndim),y(ndim),val

      val = 0.
      do i=1,ndim
            val = val + x(i)*y(i)
      enddo
      
      return
      
      end subroutine
C --------------------------------------------------------------------
      function cauchy_eq_f(sig)
      
      implicit none
      real*8 cauchy_eq_f,sigDotsig
      real*8 sig(3,3),sigDev(3,3),sigHydro

      sigHydro = (sig(1,1) + sig(2,2) + sig(3,3))/3.d0
      sigDev = sig - sigHydro*delta2
      call tens2_inner_product_f(sigDotsig,sigDev,sigDev)
      cauchy_eq_f = dsqrt((3.d0/2.d0)*sigDotsig)

      return 

      end function
C --------------------------------------------------------------------
      function strain_eq_f(eps)
      
      implicit none
      real*8 strain_eq_f,epsDoteps
      real*8 eps(3,3),epsDev(3,3),epsHydro

      epsHydro = (eps(1,1) + eps(2,2) + eps(3,3))/3.d0
      epsDev = eps - epsHydro*delta2
      call tens2_inner_product_f(epsDoteps,epsDev,epsDev)
      strain_eq_f = dsqrt((2.d0/3.d0)*epsDoteps)

      return 

      end function
C --------------------------------------------------------------------
      function cauchy_hydro_f(sig)

      implicit none
      real*8 cauchy_hydro_f
      real*8 sig(3,3)

      cauchy_hydro_f = (sig(1,1) + sig(2,2) + sig(3,3))/3.d0

      return 

      end function
C --------------------------------------------------------------------
      subroutine strainEul_f(strain,f)

      implicit none
      real*8 f(3,3),f_T(3,3),b(3,3),b_inv(3,3),strain(3,3)

      call tens2_transpose_f(f_T,f)
      call tens2_prod_tens2_f(b,f_T,f)
      call tens2_inv_f(b_inv,b)

      strain = .5d0*(delta2 - b_inv)

      end subroutine

      end module umat_subroutines