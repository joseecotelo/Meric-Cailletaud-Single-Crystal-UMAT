      MODULE math_subroutines 
      implicit none
      INTEGER, PARAMETER :: PRE = KIND(1.D0)
      CONTAINS		 
C***********************************************************************
!     ORDERS THE ELEMENTS AV(I) OF THE VECTOR AV FROM MAXIMUM TO MINIMUM, I.E., AV(1)>AV(2)>...AV(N)
!       SUBROUTINE KORDERVecEL(AV,N)
!       IMPLICIT NONE
!       REAL(KIND=PRE) AV,Amax,AUX

!       INTEGER N,I,J,Jmax
!       DIMENSION AV(N)

!       DO I=1,N
!       Amax=-10.D0**(100.D0)
!       Jmax=I
!         DO J=I,N
!           IF (AV(J).GE.Amax) THEN
!             Amax=AV(J)
!             Jmax=J
!           ENDIF
!         ENDDO
!       AUX=AV(I)
!       AV(I)=Amax
!       AV(Jmax)=AUX
!       ENDDO

!       RETURN
!       END SUBROUTINE
C************************************************************************
      SUBROUTINE KINITIA(A,N)
      IMPLICIT NONE
C
      REAL(KIND=PRE) A
      INTEGER N,I
C
      DIMENSION A(N)
C
      DO I=1,N
        A(I)=0.D0
      ENDDO
C
      RETURN
      END SUBROUTINE
C***********************************************************************
      SUBROUTINE KCOPY(A,B,N)
C
      IMPLICIT NONE
C
      REAL(KIND=PRE) A,B
      INTEGER N,I
      DIMENSION A(N),B(N)
C
      DO I=1,N
      B(I)=A(I)
      ENDDO
C
      RETURN
      END SUBROUTINE
! C***********************************************************************
! !     MATRIX MULTIPLICATION C(L,N)=A(L,M)*B(M,N)
!       SUBROUTINE KMULT(A,B,C,L,M,N)
! C
!       IMPLICIT NONE
!       REAL(KIND=PRE) A,B,C,AUX
!       INTEGER L,M,N,I,J,K

!       DIMENSION A(L,M),B(M,N),C(L,N)
! C
!       DO 10 I=1,L
!       DO 10 J=1,N
!       AUX=0.D0
!       DO 20 K=1,M
!  20   AUX=AUX+A(I,K)*B(K,J)
!  10   C(I,J)=AUX
! C
!       RETURN
!       END SUBROUTINE
! C***********************************************************************
! C     DOT PRODUCT OF TWO THREE DIMENSIONAL VECTORS, C=A(1:3).B(1:3)
!       SUBROUTINE KADOTB(A,B,C)
!       IMPLICIT NONE

!       REAL(KIND=PRE) A,B,C
!       DIMENSION A(3),B(3)
! C
!       C=A(1)*B(1)+ A(2)*B(2)+ A(3)*B(3)
! C
!       RETURN
!       END SUBROUTINE
! C***********************************************************************
! !     DOT PRODUCT OF TWO VECTORIZED SYMMETRIC TENSORS, C=A(1:NTENS)*B(1:NTENS)
!       SUBROUTINE KAB(A,B,C,NTENS,NDI)

!       IMPLICIT NONE
!       REAL(KIND=PRE) A,B,C
!       INTEGER NTENS,NDI,I
!       DIMENSION A(NTENS),B(NTENS)

!       C=0.D0
!       DO I=1,NDI
!         C=C+A(I)*B(I)
!       ENDDO
!       DO I=NDI+1,NTENS
!         C=C+A(I)*B(I)*2.D0
!       ENDDO

!       RETURN
!       END SUBROUTINE
C***********************************************************************
!     CALCULATES THE INVERSE USING THE LU DECOMPOSITION
      SUBROUTINE KINV(A,AINV,N)
C
      IMPLICIT NONE
      REAL(KIND=PRE) A,AINV,Y
      INTEGER N,INDX,I,J
      DIMENSION A(N,N),AINV(N,N),Y(N,N)
      DIMENSION INDX(1000)
C
      CALL KINITIA(Y,N*N)
      DO I=1,N
      Y(I,I)=1.D0
      ENDDO
      CALL KCOPY(A,AINV,N*N)
      CALL KLUDCMP(AINV,N,INDX)
      DO J=1,N
      CALL KLUBKSB(AINV,N,INDX,Y(1,J))
      ENDDO
      CALL KCOPY(Y,AINV,N*N)
C
      RETURN
      END SUBROUTINE
C***********************************************************************
!     Replaces an NxN matrix A by its LU decomposition.
!     INDX: Output vector which records the row permutation effected by the partial pivoting.
      SUBROUTINE KLUDCMP(A,N,INDX)
      IMPLICIT NONE
      REAL(KIND=PRE) A,TINY,VV,AAMAX,SUM,DUM
      INTEGER N,INDX,NMAX,I,J,K,IMAX
      PARAMETER (NMAX=1000,TINY=1.0E-20)
      DIMENSION A(N,N),INDX(N),VV(NMAX)
C
      DO I=1,N
      AAMAX=0.D0
      DO J=1,N
      IF (DABS(A(I,J)).GT.AAMAX) AAMAX=DABS(A(I,J))
         ENDDO
       IF (AAMAX.EQ.0.D0) THEN
        WRITE(*,*) 'Singular matrix'
        STOP
         ENDIF
         VV(I)=1.D0/AAMAX
            ENDDO
          DO J=1,N
          DO I=1,J-1
	      SUM=A(I,J)
	      DO K=1,I-1
		    SUM=SUM-A(I,K)*A(K,J)
            ENDDO
              A(I,J)=SUM
         ENDDO
          AAMAX=0.D0
	      DO I=J,N
	      SUM=A(I,J)
	      DO K=1,J-1
          SUM=SUM-A(I,K)*A(K,J)
             ENDDO
              A(I,J)=SUM
	      DUM=VV(I)*DABS(SUM)
	      IF (DUM.GE.AAMAX) THEN
          IMAX=I
          AAMAX=DUM
	      ENDIF
         ENDDO
          IF (J.NE.IMAX) THEN
	      DO K=1,N
          DUM=A(IMAX,K)
           A(IMAX,K)=A(J,K)
           A(J,K)=DUM
             ENDDO
	      VV(IMAX)=VV(J)
	      ENDIF
	     INDX(J)=IMAX
	      IF (A(J,J).EQ.0.D0) A(J,J)=TINY
	      IF (J.NE.N) THEN
	      DUM=1.D0/A(J,J)
	      DO I=J+1,N
          A(I,J)=A(I,J)*DUM
             ENDDO
          ENDIF
      ENDDO
C
      RETURN
      END SUBROUTINE
C***********************************************************************
!     Solves the linear system A.x=b.
!     A & INDX are inputs, returned by SUBROUTINE LUDCMP.
!     B is input as the right-hand side vector b, and returns with the solution vector x.
      SUBROUTINE KLUBKSB(A,N,INDX,B)
      IMPLICIT NONE
      REAL(KIND=PRE) A,B,SUM
      INTEGER N,INDX,II,I,LL,J
      DIMENSION A(N,N),INDX(N),B(N)
      II=0
      DO I=1,N
      LL=INDX(I)
      SUM=B(LL)
      B(LL)=B(I)
	  IF (II.NE.0) THEN
	      DO J=II,I-1
	          SUM=SUM-A(I,J)*B(J)
              ENDDO
          ELSE IF (SUM.NE.0.D0) THEN
	      II=I
	     ENDIF
	     B(I)=SUM
      ENDDO
       DO I=N,1,-1
	     SUM=B(I)
	     IF (I.LT.N) THEN
	      DO J=I+1,N
	          SUM=SUM-A(I,J)*B(J)
              ENDDO
          ENDIF
	     B(I)=SUM/A(I,I)
      ENDDO
C
      RETURN
      END SUBROUTINE
! C***********************************************************************
! !     I HAVE SKIPPED KDetA6x6,K2DMatTrans, AS APPARENTLY IT HAS NOT BEEN USED IN THE
! !     CALCULATIONS
! C***********************************************************************
! !     COMPUTE THE PRODUCT OF TWO 4-TH ORDER TENSORS

!       SUBROUTINE TENSPROD44(A,B,C)
!       IMPLICIT NONE
!       REAL(KIND=PRE) A,B,C,SUM
!       INTEGER I,J,K,L,M,N
!       DIMENSION A(3,3,3,3),B(3,3,3,3),C(3,3,3,3)
! C
!       DO I=1,3
! 	     DO J=1,3
! 	     DO K=1,3
! 	      DO L=1,3
!            SUM=0.D0
! 		    DO M=1,3
! 		    DO N=1,3
! 		    SUM=SUM+A(I,J,M,N)*B(M,N,K,L)
!             END DO
! 		    END DO
! 		     C(I,J,K,L)=SUM
! 	      END DO
! 	     END DO
! 	     END DO
!        END DO
! C
!        RETURN
!        END SUBROUTINE
! C
! C************************************************************************
! !    COMPUTE THE PRODUCT OF A 4-TH ORDER TENSOR AND 2ND ORDER TENSOR
! C
!       SUBROUTINE TENSPROD42(A,B,C)
!          IMPLICIT NONE
! 	     REAL(KIND=PRE) A,B,C,SUM
! 	     INTEGER I,J,K,L
! 	     DIMENSION A(3,3,3,3),B(3,3),C(3,3)
! C
!           DO I=1,3
! 	      DO J=1,3
! 	      SUM=0.D0
! 	       DO K=1,3
! 	       DO L=1,3
! 	       SUM=SUM+A(I,J,K,L)*B(K,L)
! 	       END DO
! 	       END DO
! 	       C(I,J)=SUM
!            END DO
! 	      END DO
!          RETURN
! 	     END SUBROUTINE
! !
! C************************************************************************
! !   I HAVE SKIPED THE  EXPONENTIAL OF AN ANTI-SYMMETRIC TENSOR,KExpW, REQUIRED FOR ROTATION
! C***************************************************************************
! !     I HAVE SKIPPED DET3X3,JACOBI,,gauleg,KINVERSEAP,
! !     AS APPARENTLY IT HAS NOT BEEN USED IN THE CALCULATIONS
! C *****************************************************************************
! !     CALCULATES INNER PRODUCT OF TWO MATRICES REPRESENTING 4th order tensors, C=A.B
! !     A = Matrix 1, B = Matrix 2, C = Scalar product
! !     M = # of rows, N = # of columns
!       SUBROUTINE KINNERPROD(A,B,C,M,N)
!       IMPLICIT NONE
! C
!       INTEGER I,J,M,N
!       REAL(KIND=PRE) A,B,C
!       DIMENSION A(M,N),B(M,N)
! C
!       C=0.D0
! 	     DO I=1,M
! 	     DO J=1,N
! 	      C=C+A(I,J)*B(I,J)
!         END DO
!       END DO
!       RETURN
!       END SUBROUTINE
! C
! !       SUBROUTINE SCALARPROD4(A,B,C,N)
! !       IMPLICIT NONE
! ! C
! !       INTEGER I,J,M,N
! !       REAL(KIND=PRE) A,B,C
! !       DIMENSION A(N,N),B(N,N)
! ! C
! !       C=0.D0
! ! 	      DO I=1,M
! ! 	      DO J=1,N
! ! 	            C=C+A(I,J)*B(I,J)
! !             END DO
! !             END DO
! !       RETURN
! !       END SUBROUTINE
! C******************************************************************************
! !     Calculates product of a third order tensor, a matrix and a vector
!       subroutine KHCOa(S,C,a,NDIM,B)
!       implicit none
!       integer i,j,k,l,t,NDIM
!       real(kind=PRE) C(3,3),B(3,NDIM), a(3),
!      + S(3,3,NDIM),dSum
!       do 10 k=1,NDIM
!         do 10 t=1,3
!          dSum=0.d0
!           do 20 i=1,3
!             do 20 j=1,3
!   20       dSum=dSum+C(t,i)*S(i,j,k)*a(j)
!        B(t,k)= dSum
!   10  continue

!       return
!       end subroutine
! C******************************************************************************
! !     Calculates product of a matrix and a vector
!       subroutine KBa(B,a,c,M,N)
!       implicit none
!       integer i,j,M,N
!       real(kind=PRE) B(M,N), a(N),c(M),dSum
!       do  i=1,M
!       dSum=0.d0
!       do  j=1,N
!         dSum = dSum+B(i,j)*a(j)
!       end do
!       c(i)=dSum
!       end do
!       return
!       end subroutine
! C******************************************************************************
! !     Calculates product of a vector and amatrix
!       subroutine KaDB(a,B,c,M,N)
!       implicit none
!       integer i,j,M,N
!       real(kind=PRE) B(M,N), a(M),c(N),dSum
!       do  j=1,N
!       dSum=0.d0
!       do  i=1,M
!         dSum =dSum+a(i)*B(i,j)
!       end do
!       c(j)=dSum
!       end do
!       return
!       end subroutine
! C******************************************************************************
! !     Calculates dyadic of two vectors
!       subroutine KabC(a,b,C,M,N)
!       implicit none
!       integer i,j,M,N
!       real(kind=PRE) C(M,N), a(M),b(N)
!       do 10 i=1,M
!       do 10 j=1,N
!   10      C(i,j)=a(i)*b(j)
!        return
!        end subroutine
! C******************************************************************************
! !     Calculates 'O' dyadic of two Matrices of same dimension H(M,N,M,N)=A(M,N)*B(M,N)
!       subroutine KABH(A,B,H,M,N)
!       implicit none
!       integer i,j,k,l,M,N
!       real(kind=PRE) H(M,N,M,N), A(M,N),B(M,N)

!        do 10 i=1,M
!        do 10 j=1,N
!        do 10 k=1,M
!        do 10 l=1,N
!   10     H(i,j,k,l)=A(i,k)*B(j,l)
!        return
!        end subroutine
! C******************************************************************************
! !      Calculates derivative of a symmetric matrix with respect to its inverse
!        subroutine CINVCSYM(CSYM,NDIM,DCSYM)
!        IMPLICIT NONE
!        integer NDIM
!        INTEGER (KIND=PRE) i,j,k,l,r,s,t,u,v,w
!        REAL (KIND=PRE) CSYM(NDIM,NDIM),H(3,3,NDIM),
!      + DCSYM(NDIM,NDIM,NDIM,NDIM),B(3,3,6)

!        DCSYM=0.d0
!         DO i=1,5
!         DO j=1,5
!         DO k=1,5
!         DO l=1,5
!                   DCSYM(i,j,k,l)=-1.d0/2.d0*(CSYM(i,k)*CSYM(j,l)
!      +            +CSYM(i,l)*CSYM(j,k))
! C  10      continue
!                     enddo
!                 enddo
!              enddo
!          enddo

!        return
!        end subroutine
! C******************************************************************************
! ! Converts a SIXTH Order Symmetric Tensor to FOURTH Order one.
!         subroutine PDH(P,NDIM,PH)
!         IMPLICIT NONE
!         integer  NDIM
!         INTEGER (KIND=PRE) i,j,k,l,r,s,t,u,v,w
!         REAL (KIND=PRE) P(3,3,3,3,5,5),PH(NDIM,NDIM,5,5),B(3,3,6),dsum
!         B=0.d0
!         B(1,1,2)=-1.D0/dsqrt(6.d0)
!         B(2,2,2)=-1.D0/dsqrt(6.d0)
!         B(3,3,2)= 2.D0/dsqrt(6.d0)

!         B(1,1,1)=-1.D0/dsqrt(2.d0)
!         B(2,2,1)= 1.D0/dsqrt(2.d0)

!         B(2,3,3)=1.D0/dsqrt(2.d0)
!         B(3,2,3)=1.D0/dsqrt(2.d0)

!         B(1,3,4)=1.D0/dsqrt(2.d0)
!         B(3,1,4)=1.D0/dsqrt(2.d0)

!         B(1,2,5)=1.D0/dsqrt(2.d0)
!         B(2,1,5)=1.D0/dsqrt(2.d0)

!         B(1,1,6)=1.D0/dsqrt(3.d0)
!         B(2,2,6)=1.D0/dsqrt(3.d0)
!         B(3,3,6)=1.D0/dsqrt(3.d0)

!        DO 10 i=1,NDIM
!        DO 10 j=i,NDIM
!         DO 10  r=1,5
!         DO 10  s=r,5
!         dSum=0.d0
!           DO 20 t=1,3
!           DO 20 u=1,3
!           DO 20 v=1,3
!           DO 20 w=1,3
!    20     dsum=dsum+P(t,u,v,w,r,s)*B(t,u,i)*B(v,w,j)

!           PH(i,j,r,s)=dsum
!           PH(j,i,r,s)=PH(i,j,r,s)
!           PH(i,j,s,r)=PH(i,j,r,s)
!           PH(j,i,s,r)=PH(i,j,r,s)
!    10  CONTINUE

!        RETURN
!        END SUBROUTINE
! C******************************************************************************
! !      Calculates derivative of a symmetric TENSOR with respect to its inverse
!        subroutine DCM(CE4,NDIM,IOPT,B,CES,DCMS,DCSYM,SCE)
!        IMPLICIT NONE
!        INTEGER NDIM,IOPT
!        INTEGER (KIND=PRE) I,J,K,L,R,S,N,M
!        REAL (KIND=PRE) CE4(NDIM,NDIM),DCSYM(NDIM,NDIM,NDIM,NDIM),
!      + DCMS(3,3,3,3,NDIM,NDIM),B(3,3,6),dSum,CES(3,3,5),SCE(5,3,3)


!         B=0.d0                                     ! MAY BE REDUNDANT
!         B(1,1,2)=-1.d0/sqrt(6.d0)
!         B(2,2,2)=-1.d0/sqrt(6.d0)
!         B(3,3,2)= 2.D0/sqrt(6.d0)

!         B(1,1,1)=-1.d0/sqrt(2.d0)
!         B(2,2,1)= 1.d0/sqrt(2.d0)

!         B(2,3,3)=1.d0/sqrt(2.d0)
!         B(3,2,3)=1.d0/sqrt(2.d0)

!         B(1,3,4)=1.d0/sqrt(2.d0)
!         B(3,1,4)=1.d0/sqrt(2.d0)

!         B(1,2,5)=1.d0/sqrt(2.d0)
!         B(2,1,5)=1.d0/sqrt(2.d0)

!         B(1,1,6)=1.d0/sqrt(3.d0)
!         B(2,2,6)=1.d0/sqrt(3.d0)
!         B(3,3,6)=1.d0/sqrt(3.d0)


!         IF(IOPT.EQ.1) THEN
!         B=0.d0
!         B(1,1,2)=-1.d0/sqrt(6.d0)
!         B(2,2,2)=-1.d0/sqrt(6.d0)
!         B(3,3,2)= 2.D0/sqrt(6.d0)

!         B(1,1,1)=-1.d0/sqrt(2.d0)
!         B(2,2,1)= 1.d0/sqrt(2.d0)

!         B(2,3,3)=1.d0/sqrt(2.d0)
!         B(3,2,3)=1.d0/sqrt(2.d0)

!         B(1,3,4)=1.d0/sqrt(2.d0)
!         B(3,1,4)=1.d0/sqrt(2.d0)

!         B(1,2,5)=1.d0/sqrt(2.d0)
!         B(2,1,5)=1.d0/sqrt(2.d0)

!         B(1,1,6)=1.d0/sqrt(3.d0)
!         B(2,2,6)=1.d0/sqrt(3.d0)
!         B(3,3,6)=1.d0/sqrt(3.d0)
!         ENDIF


! ! *** Converts FOURTH ORDER SYMMETRIC TENSOR INTO A THIRD ORDER
!        IF(IOPT.EQ.2) THEN
!         DO 140 I=1,3
!         DO 140 J=I,3
!         DO 140 M=1,NDIM
!         CES(I,J,M)=0.d0
!         DO 150 N=1,NDIM
!         CES(I,J,M)=CES(I,J,M)+CE4(N,M)*B(I,J,N)
!   150   CONTINUE
!         CES(J,I,M)=CES(I,J,M)
!   140   CONTINUE
!       ENDIF

!        IF(IOPT.EQ.3) THEN

!         DO 10 i=1,NDIM
!         DO 10 j=i,NDIM
!         DO 10 k=1,NDIM
!         DO 10 l=k,NDIM
!          DCSYM(i,j,k,l)=-1.d0/2.d0*(CE4(i,k)*CE4(j,l)
!      +                           + CE4(i,l)* CE4(j,k))

!          DCSYM(j,i,k,l)=DCSYM(i,j,k,l)
!          DCSYM(i,j,l,k)=DCSYM(i,j,k,l)
!          DCSYM(j,i,l,k)=DCSYM(i,j,k,l)

!   10   continue
!        DCMS=0.d0
!        DO 30 R=1,NDIM
!        DO 30 S=R,NDIM
!         DO 30 I=1,3
!         DO 30 J=I,3
!         DO 30 K=1,3
!         DO 30 L=K,3
!         dSum=0.d0
!         DO 20 N=1,NDIM
!         DO 20 M=1,NDIM
!         dSum=dSum+DCSYM(N,M,R,S)*B(I,J,N)*B(K,L,M)
!  20     CONTINUE

!         DCMS(I,J,K,L,R,S)=dSum
!         DCMS(J,I,K,L,R,S)=dSum
!         DCMS(I,J,L,K,R,S)=dSum
!         DCMS(J,I,L,K,R,S)=dSum

!         DCMS(I,J,K,L,S,R)=dSum
!         DCMS(J,I,K,L,S,R)=dSum
!         DCMS(I,J,L,K,S,R)=dSum
!         DCMS(J,I,L,K,S,R)=dSum
!  30     CONTINUE

!        ENDIF

!          IF(IOPT.EQ.4) THEN
!         DO 50 i=1,NDIM
!         DO 50 j=i,NDIM
!         DO 50 k=1,NDIM
!         DO 50 l=k,NDIM
!           DCSYM(i,j,k,l)=-1.d0/2.d0*(CE4(i,k)*CE4(j,l)
!      +                           + CE4(i,l)* CE4(j,k))

!          DCSYM(j,i,k,l)=DCSYM(i,j,k,l)
!          DCSYM(i,j,l,k)=DCSYM(i,j,k,l)
!          DCSYM(j,i,l,k)=DCSYM(i,j,k,l)
!   50   CONTINUE

!        ENDIF

!         IF(IOPT.EQ.5) THEN
!         DO 60 I=1,3
!         DO 60 J=I,3
!         DO 60 M=1,NDIM
!         SCE(M,I,J)=0.d0
!          DO 70 N=1,NDIM
!          SCE(M,I,J)=SCE(M,I,J)+CE4(M,N)*B(I,J,N)
!   70     CONTINUE
!         SCE(M,J,I)=SCE(M,I,J)
!   60   CONTINUE

!        ENDIF

!        RETURN
!        END SUBROUTINE
! C******************************************************************************
!        SUBROUTINE DotFourth (BINVM,SD,NDIM,BINDotSD)
!        IMPLICIT NONE
!        INTEGER  NDIM
!        INTEGER (KIND=PRE) i,j,r,s
!        REAL (KIND=PRE) BINVM(3,3,NDIM,NDIM),SD(3),dSum,
!      +  BINDotSD(3,NDIM,NDIM)

!        BINDotSD=0.D0
!        DO 10 i=1,3
!        DO 10 r=1,NDIM
!        DO 10 s=r,NDIM
!      	   dSum=0.d0
!            DO j=1,3
!            dSum=dSum+BINVM(i,j,r,s)*SD(j)
!            END DO
!            BINDotSD(i,r,s)=dSum
!            BINDotSD(i,s,r)=dSum
!    10  CONTINUE

!        RETURN
!        END SUBROUTINE
! C********************************************************************************
! !   BASIS FOR 5x5 SYM MATRIX AND 5X5X5X5 SYM TENSOR
!        SUBROUTINE CHG_BASIS6D(CE2,C2,CE4,C4,IOPT,NDIM,KDIM)
!        IMPLICIT NONE
!        INTEGER IOPT,I,J,K,L,M,N,KDIM,NDIM
!        REAL (KIND=PRE) SQR2,RSQ2,RSQ3,RSQ6
!        PARAMETER (SQR2=1.41421356237309   )
!        PARAMETER (RSQ2=0.70710678118654744)

!       REAL (KIND=PRE) C2(NDIM,NDIM),C4(NDIM,NDIM,NDIM,NDIM),
!      + B5D(NDIM,NDIM,KDIM),CE2(KDIM),CE4(KDIM,KDIM)
! C
! C *** CALCULATES BASIS TENSORS B(N)

!         B5D=0.D0

!         B5D(1,1,1)=1.D0
!         B5D(2,2,2)=1.D0
!         B5D(3,3,3)=1.D0
!         B5D(4,4,4)=1.D0
!         B5D(5,5,5)=1.D0

!         B5D(1,2,6)=RSQ2
!         B5D(2,1,6)=RSQ2
!         B5D(1,3,7)=RSQ2
!         B5D(3,1,7)=RSQ2
!         B5D(1,4,8)=RSQ2
!         B5D(4,1,8)=RSQ2
!         B5D(1,5,9)=RSQ2
!         B5D(5,1,9)=RSQ2

!         B5D(2,3,10)=RSQ2
!         B5D(3,2,10)=RSQ2
!         B5D(2,4,11)=RSQ2
!         B5D(4,2,11)=RSQ2
!         B5D(2,5,12)=RSQ2
!         B5D(5,2,12)=RSQ2

!         B5D(3,4,13)=RSQ2
!         B5D(4,3,13)=RSQ2
!         B5D(3,5,14)=RSQ2
!         B5D(5,3,14)=RSQ2

!         B5D(4,5,15)=RSQ2
!         B5D(5,4,15)=RSQ2

!         IF (NDIM.EQ.6) THEN
!           B5D(6,6,16)=1.D0
!           B5D(1,6,17)=RSQ2
!           B5D(6,1,17)=RSQ2
!           B5D(2,6,18)=RSQ2
!           B5D(6,2,18)=RSQ2
!           B5D(3,6,19)=RSQ2
!           B5D(6,3,19)=RSQ2
!           B5D(4,6,20)=RSQ2
!           B5D(6,4,20)=RSQ2
!           B5D(5,6,21)=RSQ2
!           B5D(6,5,21)=RSQ2
!         END IF
! C
! C *** CALCULATES CARTESIAN SECOND ORDER TENSOR FROM b-COMPONENTS VECTOR.
!       IF(IOPT.EQ.1) THEN
!         DO 40 I=1,NDIM
!         DO 40 J=1,NDIM
!         C2(I,J)=0.D0
!         DO N=1,KDIM
!         C2(I,J)=C2(I,J)+CE2(N)*B5D(I,J,N)
!         END DO
!   40    CONTINUE
!       ENDIF
! C *** CALCULATES KDIMx1 b-COMPONENTS VECTOR FROM SECOND ORDER TENSOR.
!       IF(IOPT.EQ.2) THEN
!        DO N=1,KDIM
!         CE2(N)=0.D0
!         DO 50 I=1,NDIM
!         DO 50 J=1,NDIM
!         CE2(N)=CE2(N)+C2(I,J)*B5D(I,J,N)
!   50    CONTINUE
!        END DO
!       ENDIF
! C *** CALCULATES FOURTH ORDER TENSOR FROM b-COMPONENTS MATRIX.
!       IF(IOPT.EQ.3) THEN
!        DO 20 I=1,NDIM
!        DO 20 J=1,NDIM
!        DO 20 K=1,NDIM
!        DO 20 L=1,NDIM
!         C4(I,J,K,L)=0.D0
!         DO N=1,KDIM
!         DO M=1,KDIM
!         C4(I,J,K,L)=C4(I,J,K,L)+CE4(N,M)*B5D(I,J,N)*B5D(K,L,M)
!         END DO
!         END DO
!   20   CONTINUE
!       ENDIF
! C *** CALCULATES KDIMxKDIM b-COMPONENTS MATRIX FROM FOURTH ORDER TENSOR.
!       IF(IOPT.EQ.4) THEN
!        DO N=1,KDIM
!        DO M=1,KDIM
!         CE4(N,M)=0.D0
!         DO 30 I=1,NDIM
!         DO 30 J=1,NDIM
!         DO 30 K=1,NDIM
!         DO 30 L=1,NDIM
!         CE4(N,M)=CE4(N,M)+C4(I,J,K,L)*B5D(I,J,N)*B5D(K,L,M)
!   30    CONTINUE
!        END DO
!        END DO
!       ENDIF

!       RETURN
!       END SUBROUTINE
! C******************************************************************************
! ! This subroutine expands tensors in terms of the Lequeu basis.
! !     SUBROUTINE CHG_BASIS    --->   VERSION 19/JUL/01
! C
! C     (modif. 06/FEB/98 - same convention as SELFPOLY - C.N.T.)
! C     (modif. 16/JUN/99 - same convention as Maudlin  - C.N.T.)
! C     (modif. 10/MAY/01 - KDIM version - R.L.)
! C
! C     KDIM=5 or 6, FOR DEVIATORIC or DEV+HYDROST TENSORS, RESPECTIVELY.
! C     IOPT=0: DEFINES A BASIS OF 6 SECOND ORDER TENSORS B(N).
! C     IOPT=1: CALCULATES SECOND ORDER TENSOR 'C2' AS AN EXPANSION IN TERMS
! C             OF VECTOR COMPONENTS CE2(KDIM) AND THE BASIS TENSORS B(KDIM).
! C     IOPT=2: CALCULATES COMPONENTS OF C2 AS A VECTOR CE2(KDIM).
! C     IOPT=3: CALCULATES FOURTH ORDER TENSOR 'C4' AS AN EXPANSION IN TERMS
! C             OF MATRIX COMPONENTS CE4(K,K) AND THE BASIS TENSORS B(KDIM).
! C     IOPT=4: CALCULATES MATRIX COMPONENTS CE4(K,K) OF TENSOR 'C4'.
!       SUBROUTINE CHG_BASIS(CE2,C2,CE4,C4,IOPT,KDIM)
! 	      IMPLICIT NONE
! 	     INTEGER KDIM,IOPT,I,J,K,L,M,N

!       REAL (KIND=PRE) CE2(KDIM),C2(3,3),CE4(KDIM,KDIM),C4(3,3,3,3),
!      + B(3,3,6)

! C     DIMENSION B(3,3,6)
! C     DATA B /RSQ6,0,   0,   0,   RSQ6,0,   0,   0,  -2*RSQ6,
! C    #        RSQ2,0,   0,   0,  -RSQ2,0,   0,   0,   0,
! C    #        0,   0,   0,   0,   0,   RSQ2,0,   RSQ2,0,
! C    #        0,   0,   RSQ2,0,   0,   0,   RSQ2,0,   0,
! C    #        0,   RSQ2,0,   RSQ2,0,   0,   0,   0,   0,
! C    #        RSQ3,0,   0,   0,   RSQ3,0,   0,   0,   RSQ3/
! C
! C
!         IF(IOPT.EQ.0) THEN
! ! *** CALCULATES BASIS TENSORS B(N)

!         B=0.d0

!         B(1,1,2)=-1/sqrt(6.d0)
!         B(2,2,2)=-1/sqrt(6.d0)
!         B(3,3,2)= 2.D0/sqrt(6.d0)

!         B(1,1,1)=-1/sqrt(2.d0)
!         B(2,2,1)= 1/sqrt(2.d0)

!         B(2,3,3)=1/sqrt(2.d0)
!         B(3,2,3)=1/sqrt(2.d0)

!         B(1,3,4)=1/sqrt(2.d0)
!         B(3,1,4)=1/sqrt(2.d0)

!         B(1,2,5)=1/sqrt(2.d0)
!         B(2,1,5)=1/sqrt(2.d0)

!         B(1,1,6)=1/sqrt(3.d0)
!         B(2,2,6)=1/sqrt(3.d0)
!         B(3,3,6)=1/sqrt(3.d0)

!       ENDIF

! ! *** CALCULATES CARTESIAN SECOND ORDER TENSOR FROM b-COMPONENTS VECTOR.
!       IF(IOPT.EQ.1) THEN
!         B=0.d0
!         B(1,1,2)=-1/sqrt(6.d0)
!         B(2,2,2)=-1/sqrt(6.d0)
!         B(3,3,2)= 2.D0/sqrt(6.d0)

!         B(1,1,1)=-1/sqrt(2.d0)
!         B(2,2,1)= 1/sqrt(2.d0)

!         B(2,3,3)=1/sqrt(2.d0)
!         B(3,2,3)=1/sqrt(2.d0)

!         B(1,3,4)=1/sqrt(2.d0)
!         B(3,1,4)=1/sqrt(2.d0)

!         B(1,2,5)=1/sqrt(2.d0)
!         B(2,1,5)=1/sqrt(2.d0)

!         B(1,1,6)=1/sqrt(3.d0)
!         B(2,2,6)=1/sqrt(3.d0)
!         B(3,3,6)=1/sqrt(3.d0)
!         DO 40 I=1,3
!         DO 40 J=1,3
!         C2(I,J)=0.0
!         DO 40 N=1,KDIM
!    40   C2(I,J)=C2(I,J)+CE2(N)*B(I,J,N)
!       ENDIF

! ! *** CALCULATES KDIMx1 b-COMPONENTS VECTOR FROM SECOND ORDER TENSOR.
!       IF(IOPT.EQ.2) THEN
!         B=0.d0
!         B(1,1,2)=-1/sqrt(6.d0)
!         B(2,2,2)=-1/sqrt(6.d0)
!         B(3,3,2)= 2.D0/sqrt(6.d0)

!         B(1,1,1)=-1/sqrt(2.d0)
!         B(2,2,1)= 1/sqrt(2.d0)

!         B(2,3,3)=1/sqrt(2.d0)
!         B(3,2,3)=1/sqrt(2.d0)

!         B(1,3,4)=1/sqrt(2.d0)
!         B(3,1,4)=1/sqrt(2.d0)

!         B(1,2,5)=1/sqrt(2.d0)
!         B(2,1,5)=1/sqrt(2.d0)

!         B(1,1,6)=1/sqrt(3.d0)
!         B(2,2,6)=1/sqrt(3.d0)
!         B(3,3,6)=1/sqrt(3.d0)
!         DO 50 N=1,KDIM
!         CE2(N)=0.0
!         DO 50 I=1,3
!         DO 50 J=1,3
!    50   CE2(N)=CE2(N)+C2(I,J)*B(I,J,N)
!       ENDIF

! ! *** CALCULATES FOURTH ORDER TENSOR FROM b-COMPONENTS MATRIX.
!       IF(IOPT.EQ.3) THEN
!          B=0.d0
!         B(1,1,2)=-1/sqrt(6.d0)
!         B(2,2,2)=-1/sqrt(6.d0)
!         B(3,3,2)= 2.D0/sqrt(6.d0)

!         B(1,1,1)=-1/sqrt(2.d0)
!         B(2,2,1)= 1/sqrt(2.d0)

!         B(2,3,3)=1/sqrt(2.d0)
!         B(3,2,3)=1/sqrt(2.d0)

!         B(1,3,4)=1/sqrt(2.d0)
!         B(3,1,4)=1/sqrt(2.d0)

!         B(1,2,5)=1/sqrt(2.d0)
!         B(2,1,5)=1/sqrt(2.d0)

!         B(1,1,6)=1/sqrt(3.d0)
!         B(2,2,6)=1/sqrt(3.d0)
!         B(3,3,6)=1/sqrt(3.d0)
!         DO 20 I=1,3
!         DO 20 J=1,3
!         DO 20 K=1,3
!         DO 20 L=1,3
!         C4(I,J,K,L)=0.0
!         DO 20 N=1,KDIM
!         DO 20 M=1,KDIM
!    20   C4(I,J,K,L)=C4(I,J,K,L)+CE4(N,M)*B(I,J,N)*B(K,L,M)
!       ENDIF

! ! *** CALCULATES KDIMxKDIM b-COMPONENTS MATRIX FROM FOURTH ORDER TENSOR.
!       IF(IOPT.EQ.4) THEN
!          B=0.d0
!         B(1,1,2)=-1/sqrt(6.d0)
!         B(2,2,2)=-1/sqrt(6.d0)
!         B(3,3,2)= 2.D0/sqrt(6.d0)

!         B(1,1,1)=-1/sqrt(2.d0)
!         B(2,2,1)= 1/sqrt(2.d0)

!         B(2,3,3)=1/sqrt(2.d0)
!         B(3,2,3)=1/sqrt(2.d0)

!         B(1,3,4)=1/sqrt(2.d0)
!         B(3,1,4)=1/sqrt(2.d0)

!         B(1,2,5)=1/sqrt(2.d0)
!         B(2,1,5)=1/sqrt(2.d0)

!         B(1,1,6)=1/sqrt(3.d0)
!         B(2,2,6)=1/sqrt(3.d0)
!         B(3,3,6)=1/sqrt(3.d0)
!         DO 30 N=1,KDIM
!         DO 30 M=1,KDIM
!         CE4(N,M)=0.0
!         DO 30 I=1,3
!         DO 30 J=1,3
!         DO 30 K=1,3
!         DO 30 L=1,3
!    30   CE4(N,M)=CE4(N,M)+C4(I,J,K,L)*B(I,J,N)*B(K,L,M)
!       ENDIF

!       RETURN
!       END SUBROUTINE
! C*********************************************************************************
!       SUBROUTINE KBH(B,H,A,N)
! C
!       IMPLICIT NONE
!       INTEGER I,J,K,L,M,N
!       REAL (KIND=PRE) B(N,N),H(N,N,N,N),A(N,N)
! C
!        A=0.D0

!        DO I=1,N
!         DO J=1,N
!         A(I,J)=0.D0
!         DO K=1,N
!         DO L=1,N
! 	     A(I,J)=A(I,J)+B(K,L)*H(K,L,I,J)
!          END DO
!         END DO
!         END DO
!        END DO
! C
!       RETURN
!       END SUBROUTINE
! C
! C*******************************************************************************
!       SUBROUTINE KCOHOB(C,H,B,D,N)
!       IMPLICIT NONE

!       INTEGER N,I,J,K,L,M1,M2
!       REAL(KIND=PRE) H(N,N,N,N),B(N,N),C(N,N),D(N,N,N,N)
! C
!       D=0.D0
!       DO I=1,N
!         DO J=1,N
!          DO K=1,N
! 	     DO L=1,N
! 	        D(I,J,K,L)=0.D0
!            DO M1=1,N
! 	          DO M2=1,N
! 			   D(I,J,K,L)=D(I,J,K,L)+C(I,M1)*H(M1,M2,K,L)*B(M2,J)
!                 END DO
!               END DO
!            END DO
!          END DO
!         END DO
!       END DO
! C
!       RETURN
!       END SUBROUTINE
! C
! C******************************************************************************
! C   GAUSSIAN QUADRATURE CALCULTION
! C******************************************************************************
!       SUBROUTINE GAUSS(ngp,a,b,x,w)

!       implicit none
!       integer, parameter :: PRE2 = kind(1.d0)
!       INTEGER ngp
!       REAL(KIND=PRE2) a, b, x(ngp),w(ngp)
!       CALL gauleg ( ngp, x, w )

!       CALL rescale ( ngp, a, b, x, w )           !     Rescale the data.
!       END SUBROUTINE

! C****************************************************************************************
!       SUBROUTINE  gauleg(ngp, xabsc, weig)

!       IMPLICIT NONE
!       integer, parameter :: PRE2 = kind(1.d0)
!       INTEGER  i, j, m
!       REAL(KIND=PRE2)  p1, p2, p3, pp, z, z1,API,EPS
!       INTEGER ngp           !# of Gauss Points
!       REAL(KIND=PRE2) xabsc(ngp), weig(ngp)
!       PARAMETER (EPS=3.0d-15)
! C
!       API=4.D0*DATAN(1.D0)
! C
! 	   m = (ngp + 1) / 2
! !* Roots are symmetric in the interval - so only need to find half of them  */
! 	     do i = 1, m				! Loop over the desired roots */
!      		z = Dcos( API * (i-0.25d0) / (ngp+0.5d0) )
! !*   Starting with the above approximation to the ith root,
! !*          we enter the main loop of refinement by NEWTON'S method   */
!   100       p1 = 1.0d0
!         	p2 = 0.0d0
! !*  Loop up the recurrence relation to get the Legendre
! !*  polynomial evaluated at z                 */

!         	do j = 1, ngp
!            	p3 = p2
!            	p2 = p1
!            	p1 = ((2.0d0*j-1.0d0) * z * p2 - (j-1.0d0)*p3) / j
!         	enddo

! !* p1 is now the desired Legendre polynomial. We next compute pp,
! !* its derivative, by a standard relation involving also p2, the
! !* polynomial of one lower order.      */
!         	pp = ngp*(z*p1-p2)/(z*z-1.0d0)
!         	z1 = z
!         	z = z1 - p1/pp             ! Newton's Method  */

!         	if (dabs(z-z1) .gt. EPS) GOTO  100

!       	xabsc(i) =  - z                    	! Roots will be bewteen -1.0 & 1.0 */
!       	xabsc(ngp+1-i) =  + z                	! and symmetric about the origin  */
!       	weig(i) = 2.0d0/((1.0d0-z*z)*pp*pp) ! Compute the weight and its       */
!       	weig(ngp+1-i) = weig(i)               ! symmetric counterpart         */

!       end do     ! i loop
! C
!       RETURN
!       END SUBROUTINE gauleg

! C***********************************************************************************
!       SUBROUTINE rescale ( ngp, a, b, x, w )
!       implicit none
!       integer, parameter :: PRE2 = kind(1.d0)
!       INTEGER i,ngp
!       double precision a
!       double precision b
!       double precision w(ngp)
!       double precision x(ngp)

!       do i = 1, ngp
!         x(i) = ( ( a + b ) + ( b - a ) * x(i) ) / 2.0
!       end do
!       do i = 1, ngp
!         w(i) = ( b - a ) * w(i) / 2.0
!       end do
!       RETURN
!       End subroutine
! C***********************************************************************************
!       SUBROUTINE ADOTB (A,B,DOTSUM,N)
!       IMPLICIT NONE

!       INTEGER N,I,J
!       REAL(KIND=PRE) B(N),A(N),DOTSUM

!       DOTSUM=0.d0
!       DO I=1,N
!        DOTSUM=DOTSUM+A(I)*B(I)
!       END DO
!       RETURN
!       END  SUBROUTINE
! C***********************************************************************************
! !===========================================================
! ! Evaluate eigenvalues and eigenvectors
! ! of a real symmetric matrix a(n,n): a*x = lambda*x
! ! method: Jacoby method for symmetric matrices
! ! Alex G. (December 2009)
! !-----------------------------------------------------------
! ! input ...
! ! a(n,n) - array of coefficients for matrix A
! ! n      - number of equations
! ! abserr - abs tolerance [sum of (off-diagonal elements)^2]
! ! output ...
! ! a(i,i) - eigenvalues
! ! x(i,j) - eigenvectors
! ! comments ...
! !===========================================================
!        SUBROUTINE JACOBI(a,x,abserr,n)
!        implicit none
!        integer i, j, k, n
!        double precision a(n,n),x(n,n)
!        double precision abserr, b2, bar
!        double precision beta, coeff, c, s, cs, sc

! ! initialize x(i,j)=0, x(i,i)=1
! ! *** the array operation x=0.0 is specific for Fortran 90/95
!        x = 0.0
!        do i=1,n
!        x(i,i) = 1.0
!        end do

! ! find the sum of all off-diagonal elements (squared)
!        b2 = 0.0
!        do i=1,n
!        do j=1,n
!        if (i.ne.j) b2 = b2 + a(i,j)**2
!        end do
!        end do

!       if (b2 <= abserr) return
! ! average for off-diagonal elements /2
!       bar = 0.5*b2/float(n*n)

!       do while (b2.gt.abserr)
!        do i=1,n-1
!        do j=i+1,n
!        if (a(j,i)**2 <= bar) cycle  ! do not touch small elements
!        b2 = b2 - 2.0*a(j,i)**2
!        bar = 0.5*b2/float(n*n)
! ! calculate coefficient c and s for Givens matrix
!        beta = (a(j,j)-a(i,i))/(2.0*a(j,i))
!        coeff = 0.5*beta/sqrt(1.0+beta**2)
!        s = sqrt(max(0.5+coeff,0.0))
!        c = sqrt(max(0.5-coeff,0.0))
! ! recalculate rows i and j
!           do k=1,n
!          cs =  c*a(i,k)+s*a(j,k)
!          sc = -s*a(i,k)+c*a(j,k)
!          a(i,k) = cs
!          a(j,k) = sc
!            end do
! ! new matrix a_{k+1} from a_{k}, and eigenvectors
!            do k=1,n
!         cs =  c*a(k,i)+s*a(k,j)
!         sc = -s*a(k,i)+c*a(k,j)
!         a(k,i) = cs
!         a(k,j) = sc
!         cs =  c*x(k,i)+s*x(k,j)
!         sc = -s*x(k,i)+c*x(k,j)
!         x(k,i) = cs
!         x(k,j) = sc
!            end do
!           end do
!          end do
!         end do

!       RETURN
!       END SUBROUTINE
! C***********************************************************************************
!       SUBROUTINE KExpW(WMX,EXPWMX)
! C
!       IMPLICIT DOUBLE PRECISION (A-H, O-Z)
!       DOUBLE PRECISION WMX(3,3),EXPWMX(3,3),DELTA(3,3)
!       DOUBLE PRECISION  AUX
!       INTEGER I,J
! C
!       CALL KINITIA(DELTA,3*3)
!       DELTA(1,1)=1.D0
!       DELTA(2,2)=1.D0
!       DELTA(3,3)=1.D0
! C
!       AUX=0.D0
!       DO I=1,3
!       DO J=1,3
!       AUX=AUX+WMX(I,J)*WMX(I,J)
!       ENDDO
!       ENDDO
!       AUX=DSQRT(AUX/2.D0)
! C
!        IF (WMX(1,2).EQ.0.D0.AND.WMX(1,3).EQ.0.D0.AND.WMX(2,3).EQ.0.D0)
!      + THEN
!        EXPWMX=DELTA
!        ELSE
!        EXPWMX=DELTA+(DSIN(AUX)/AUX)*WMX+
!      +	((1.D0-DCOS(AUX))/(AUX**2.D0))*MATMUL(WMX,WMX)
!        END IF
! C
!       RETURN
!       END SUBROUTINE
! C***************************************************************************
! !     COMPUTES THE "O" PRODUCT OF A FOURTH ORDER TENSOR AND SECOND ORDER TENSOR     
! C     
!       SUBROUTINE KHOB_SYM(H,B,C,N)
!       IMPLICIT NONE
!       INTEGER N,I,J,K,L,M
!       REAL(KIND=PRE) H(N,N,N,N),B(N,N),C(N,N,N,N)
! C
!       C=0.D0

!       DO I=1,N
! 	  DO J=1,N
!          DO K=1,N
! 	     DO L=K,N
! 	        C(I,J,K,L)=0.D0
! 			DO M=1,N
! 			  C(I,J,K,L)=C(I,J,K,L)+H(I,M,K,L)*B(M,J) 
!               END DO
!               C(I,J,L,K)=C(I,J,K,L)
!            END DO
!          END DO
!         END DO
!       END DO
! C
!       RETURN
!       END SUBROUTINE
! C
! C******************************************************************************
! !     COMPUTES THE "O" PRODUCT OF A SECOND ORDER TENSOR AND FOURTH ORDER TENSOR     
!       SUBROUTINE KBOH_SYM(B,H,C,N)
! C
!       IMPLICIT NONE
!       INTEGER N,I,J,K,L,M
!       REAL(KIND=PRE) H(N,N,N,N),B(N,N),C(N,N,N,N)
! C
!       C=0.D0

!       DO I=1,N
! 	  DO J=1,N
!          DO K=1,N
! 	     DO L=K,N
! 	        C(I,J,K,L)=0.D0
! 			DO M=1,N
! 			  C(I,J,K,L)=C(I,J,K,L)+H(M,J,K,L)*B(I,M) 
!               END DO
!               C(I,J,L,K)=C(I,J,K,L)
!            END DO
!          END DO
!         END DO
!       END DO
! C
!       RETURN
!       END SUBROUTINE
! C
! C******************************************************************************
!       SUBROUTINE KCOHOB_SYM(C,H,B,D,N)
! C
!       IMPLICIT NONE
!       INTEGER N,I,J,K,L,M1,M2
!       REAL(KIND=PRE) H(N,N,N,N),B(N,N),C(N,N),D(N,N,N,N)
! C
!       D=0.D0

!       DO I=1,N
! 	  DO J=I,N
!          DO K=1,N
! 	     DO L=K,N
! 	        D(I,J,K,L)=0.D0
! 			DO M1=1,N
! 	        DO M2=1,N
! 			   D(I,J,K,L)=D(I,J,K,L)+C(I,M1)*H(M1,M2,K,L)*B(M2,J)
!             END DO
!             END DO
! 		    D(J,I,K,L)=D(I,J,K,L)
! 	        D(I,J,L,K)=D(I,J,K,L)
! 	        D(J,I,L,K)=D(I,J,L,K)
!          END DO
!          END DO
!       END DO
!       END DO
! C
!       RETURN
!       END SUBROUTINE
! C
! C******************************************************************************
! !     COMPUTE THE "O" PRODUCT OF A VECTOR AND A FOURTH ORDER TENSOR     
!       SUBROUTINE KaOH_SYM(A,H,S,N)
! C
!       IMPLICIT NONE     
! 	  INTEGER J,K,L,N,M
!       REAL(KIND=PRE) A(N),H(N,N,N,N),S(N,N,N)
! C
!       S=0.D0

!       DO J=1,N
! 	  DO K=1,N
! 	    DO L=K,N
! 	      S(J,K,L)=0.D0
! 	      DO M=1,N
! 	        S(J,K,L)=S(J,K,L)+A(M)*H(M,J,K,L)
!           END DO
! 	      S(J,L,K)=S(J,K,L)
!           END DO
!         END DO
!       END DO
! C
!       RETURN
!       END SUBROUTINE
! C
! C******************************************************************************
! !  COMPUTE THE "O" PRODUCT OF A FOURTH ORDER TENSOR AND A VECTOR     
!       SUBROUTINE KHOa(H,A,S,N)
! C
!       IMPLICIT NONE    
! 	  INTEGER J,K,L,N,M
!       REAL(KIND=PRE) A(N),H(N,N,N,N),S(N,N,N)
! C
!       S=0.D0

!       DO J=1,N
! 	  DO K=1,N
! 	    DO L=1,N
! 	      S(J,K,L)=0.D0
! 	      DO M=1,N
! 	        S(J,K,L)=S(J,K,L)+A(M)*H(J,M,K,L)
!           END DO
!           END DO
!         END DO
!       END DO
! C
!       RETURN
!       END SUBROUTINE
! C
! C*******************************************************************************
! C
!       SUBROUTINE KaOHOb_SYM(A,H,B,C,N)
! C
!       IMPLICIT NONE      
! 	  INTEGER I,J,K,L,N,M
!       REAL(KIND=PRE) A(N),H(N,N,N,N),B(N),C(N,N)
! C
!       C=0.D0
! 	  DO K=1,N
! 	  DO L=K,N
! 	    C(K,L)=0.D0
!           DO I=1,N
! 		  DO M=1,N
! 		    C(K,L)=C(K,L)+A(I)*H(I,M,K,L)*B(M)
! 		  END DO
! 		END DO
! 	    C(L,K)=C(K,L)
! 	  END DO
! 	  END DO
! C	
!       RETURN
! 	  END SUBROUTINE	
! C	
! C*******************************************************************************
!       SUBROUTINE KaOS(A,S,C,N)
! C
!       IMPLICIT NONE
! 	  INTEGER I,J,K,L,N,M
!       REAL(KIND=PRE) A(N),S(N,N,N),C(N,N)
! C
!       C=0.D0

! 	  DO K=1,N
! 	  DO L=1,N
! 	    C(K,L)=0.D0
! 	    DO M=1,N
! 	    C(K,L)=C(K,L)+A(M)*S(M,K,L)
! 	    END DO
!         END DO
!       END DO
! C
! 	  RETURN
! 	  END SUBROUTINE    
! C*******************************************************************************C    
!       SUBROUTINE KSOa_SYM(S,A,C,N)
!       IMPLICIT NONE
! 	  INTEGER I,J,K,L,N,M

! 	  REAL(KIND=PRE) A(N),S(N,N,N),C(N,N)
! C
!       C=0.D0

! 	  DO K=1,N
! 	  DO L=K,N
! 	    C(K,L)=0.D0
! 	    DO M=1,N
! 	    C(K,L)=C(K,L)+A(M)*S(M,K,L)
! 	    END DO
! 	    C(L,K)=C(K,L)
!         END DO
!       END DO
! C
! 	  RETURN
! 	  END SUBROUTINE    
! C
! C*******************************************************************************C
!       SUBROUTINE KSOB_SYM(S,B,C,N)
!       IMPLICIT NONE
!       INTEGER I,J,K,L,N,M

!       REAL(KIND=PRE) S(N,N,N),B(N,N),C(N,N,N)
! C
!       C=0.D0
      
! 	  DO I=1,N
! 	  DO K=1,N
! 	    DO L=K,N
! 	      C(I,K,L)=0.D0
! 	      DO J=1,N
! 	        C(I,K,L)=C(I,K,L)+S(J,K,L)*B(J,I)
!           END DO
! 	      C(I,L,K)=C(I,K,L)
!           END DO
!         END DO
!        END DO
! C
!       RETURN
!       END SUBROUTINE
! C
! C*******************************************************************************
!       SUBROUTINE KBH_SYM(B,H,A,M,N)
!       IMPLICIT NONE     
! 	  INTEGER I,J,K,L,M,N
! 	  REAL (KIND=PRE) B(M,M),H(M,M,N,N),A(N,N)
! C
!       A=0.D0

!       DO I=1,N
! 	  DO J=I,N

! 	   A(I,J)=0.D0
! 	   DO K=1,M
! 	   DO L=1,M
! 	     A(I,J)=A(I,J)+B(K,L)*H(K,L,I,J)
!          END DO
! 	   END DO
!          A(J,I)=A(I,J)

!         END DO
!        END DO
! C
!        RETURN
! 	   END SUBROUTINE
! C***********************************************************************
! !     COMPUTE THE PRODUCT OF TWO 4-TH ORDER TENSORS IN
!       SUBROUTINE TENSPROD_SYM(A,B,C,P)
!       IMPLICIT NONE
! 	  REAL(KIND=PRE) A,B,C,DSUM
! 	  INTEGER I,J,K,L,M,N,P
! C
!       DIMENSION A(P,P,P,P),B(P,P,P,P),C(P,P,P,P)
! C
!        DO I=1,P
! 	   DO J=I,P
! 	   DO K=1,P
! 	   DO L=K,P
!            DSUM=0.D0
! 		 DO M=1,P
! 		 DO N=1,P
! 		   DSUM=DSUM+A(I,J,M,N)*B(M,N,K,L)
! 		 END DO
! 		 END DO
! 		 C(I,J,K,L)=DSUM
! 	     C(J,I,K,L)=DSUM
! 	     C(I,J,L,K)=DSUM
! 	     C(J,I,L,K)=DSUM
! 	   END DO
! 	   END DO
! 	   END DO
! 	   END DO
! C
! 	   RETURN
! 	   END SUBROUTINE
! C********************************************************************************
! !     COMPUTES THE DERIVATIVE OF A  SYMMETRIC TENSOR WITH RESPECT TO ITS INVERSE
!        SUBROUTINE DLDMX(LMX,NDIM,DLXDMX,DLDM)
!        INTEGER NDIM,I,J,K,L,R,S,N,M
!        REAL(KIND=PRE) LMX(NDIM,NDIM),DLXDMX(NDIM,NDIM,NDIM,NDIM),
!      + DLDM(3,3,3,3,NDIM,NDIM),B(3,3,6),dSum

!         B=0.d0
!         B(1,1,2)=-1.D0/dsqrt(6.d0)
!         B(2,2,2)=-1.D0/dsqrt(6.d0)
!         B(3,3,2)= 2.D0/dsqrt(6.d0)

!         B(1,1,1)=-1.D0/dsqrt(2.d0)
!         B(2,2,1)= 1.D0/dsqrt(2.d0)

!         B(2,3,3)=1.D0/dsqrt(2.d0)
!         B(3,2,3)=1.D0/dsqrt(2.d0)

!         B(1,3,4)=1.D0/dsqrt(2.d0)
!         B(3,1,4)=1.D0/dsqrt(2.d0)

!         B(1,2,5)=1.D0/dsqrt(2.d0)
!         B(2,1,5)=1.D0/dsqrt(2.d0)

!         B(1,1,6)=1.D0/dsqrt(3.d0)
!         B(2,2,6)=1.D0/dsqrt(3.d0)
!         B(3,3,6)=1.D0/dsqrt(3.d0)


!          DO 30 I=1,NDIM
!          DO 30 J=I,NDIM
!          DO 30  K=1,NDIM
!          DO 30  L=K,NDIM
!          DLXDMX(I,J,K,L)=-1.D0/2.D0*(LMX(I,K)*LMX(J,L)     ! DLXDMX
!      +                             + LMX(I,L)*LMX(J,K))
!          DLXDMX(J,I,K,L)=DLXDMX(I,J,K,L)
!          DLXDMX(I,J,L,K)=DLXDMX(I,J,K,L)
!          DLXDMX(J,I,L,K)=DLXDMX(I,J,K,L)

!   30     CONTINUE

!         DLDM=0.d0
!         DO 40 R=1,NDIM
!         DO 40 S=R,NDIM
!         DO 50 I=1,3
!         DO 50 J=I,3
!         DO 50 K=1,3
!         DO 50 L=K,3
!          dSum=0.D0
!          DO 60 N=1,NDIM
!          DO 60 M=1,NDIM
!          dSum=dSum+DLXDMX(N,M,R,S)*B(I,J,N)*B(K,L,M)
!  60      CONTINUE
!          DLDM(I,J,K,L,R,S)=dSum

!          DLDM(J,I,K,L,R,S)=DLDM(I,J,K,L,R,S)            ! SYMMETRICALLY ASSIGNING THE COMPONENTS
!          DLDM(I,J,L,K,R,S)=DLDM(I,J,K,L,R,S)
!          DLDM(J,I,L,K,R,S)=DLDM(I,J,K,L,R,S)
		 
!     	 DLDM(I,J,K,L,S,R)=DLDM(I,J,K,L,R,S)
!          DLDM(J,I,K,L,S,R)=DLDM(J,I,K,L,R,S)
!          DLDM(I,J,L,K,S,R)=DLDM(I,J,L,K,R,S)
!          DLDM(J,I,L,K,S,R)=DLDM(J,I,L,K,R,S)

!  50    CONTINUE	 

!  40    CONTINUE

!       RETURN
!       END SUBROUTINE
C********************************************************************************
      END MODULE math_subroutines
		 