!****************************************************************************
! *   FILE         = DGGEV.f90                                        *
! *   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 *
! *   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        *
!****************************************************************************
! *  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.          *
!****************************************************************************
      PROGRAM DGGEV
      IMPLICIT NONE
      INTEGER, PARAMETER :: N=5
      DOUBLE PRECISION :: A(N,N), B(N,N) 
      DOUBLE PRECISION :: ER(N), EI(N), EVR(N,N), EVL(N,N)
      INTEGER ::  I, J
      
       

      A(1,:) = (/1.0D0,5.0D0,1.0D0,2.0D0,3.0D0/)
      A(2,:) = (/5.0D0,2.0D0,13.0D0,5.0D0,6.0D0/)
      A(3,:) = (/1.0D0,13.0D0,4.0D0,1.0D0,5.0D0/)
      A(4,:) = (/2.0D0,5.0D0,1.0D0,1.0D0,7.0D0/)
      A(5,:) = (/3.0D0,6.0D0,5.0D0,7.0D0,3.0D0/)
      
      DO I = 1, N
        DO J = 1, N
        IF (I .EQ. J) THEN
           B(I,J) = 1.0D0
        ELSE
           B(I,J) = 0.0D0
        END IF
        END DO
      END DO
           
      
      
!      CALL LPCK (A, N, ER, EI, EVR, EVL)
      CALL LPCKB(A, B, N, ER, EI, EVR, EVL)
      
      PRINT*, 'Eigenvalues (REAL, IMAGINARY)'
      DO I = 1, N
          PRINT*, ER(I), EI(I)
      END DO
      
      PRINT*, 'Right eigenvectors (arranged as columns)'
      
      DO I = 1, N
         PRINT*, (EVR(I,J), J = 1, N)
      END DO
     
      PRINT*, 'Left eigenvectors (arranged as columns)'
      
      DO I = 1, N
         PRINT*, (EVL(I,J), J = 1, N)
      END DO
       
      
      END PROGRAM
      
      
      SUBROUTINE LPCK(a,n, WR, WI, VR, VL)
!====================================================================
!      LAPACK SUBROUTINE TO SOLVE EIGENVALUE PROBLEM- DGEEV
!====================================================================
     IMPLICIT NONE


     INTEGER :: n, INFO, I, LWORK
     
     DOUBLE PRECISION :: a(n, n), WR(n), WI(n), VL(n,n), VR(n,n), WORK(4*n)
     CHARACTER*1 :: V
     
     LWORK = 4*n

     CALL DGEEV('V', 'V', n, a, n, WR, WI, VL, n, VR, n, WORK, LWORK, INFO)

     IF (INFO .EQ. 0) THEN
       PRINT*, 'Successful exit'
     ELSE
       PRINT*, 'Error'
     END IF
     RETURN

     END SUBROUTINE
     
     SUBROUTINE LPCKB(a, b, n, WR, WI, VR, VL)
!====================================================================
!      LAPACK SUBROUTINE TO SOLVE EIGENVALUE PROBLEM - DGGEV
!====================================================================
     IMPLICIT NONE

     INTEGER :: n, INFO, I, LWORK
     DOUBLE PRECISION :: a(n,n), b(n,n), ALPHAR(n), ALPHAI(n), BETA(n) 
     DOUBLE PRECISION :: VL(n,n), VR(n,n), WORK(8*n), WR(n), WI(n)
     CHARACTER*1 :: V
     
     LWORK = 8*n

     CALL DGGEV('V', 'V', n, a, n, b, n, ALPHAR, ALPHAI, BETA, VL, n, VR, n, WORK, LWORK, INFO)

     IF (INFO .EQ. 0) THEN
       PRINT*, 'Successful exit'
     ELSE
       PRINT*, 'Error'
     END IF
      
      DO I = 1, N
        WR(I) = ALPHAR(I)/BETA(I)
        WI(I) = ALPHAI(I)/BETA(I)
      END DO
      
     RETURN

     END SUBROUTINE
