!****************************************************************************
! *   FILE         = ZGGEV.f90                                              *
! *   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 *
! *   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        *
!****************************************************************************
! *  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.          *
!****************************************************************************    
      PROGRAM ZGGEV
      IMPLICIT NONE
      INTEGER, PARAMETER :: N=5
      COMPLEX*16 :: A(N,N), B(N,N)
      COMPLEX*16 :: ER(N), EI(N), EVR(N,N), EVL(N,N)
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
           

      CALL LPCKB(A, B, N, ER, EVR, EVL)
      
      PRINT*, 'Eigenvalues '
      DO I = 1, N
          PRINT*, ER(I)
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
      
     
     SUBROUTINE LPCKB(a, b, n, WR, VR, VL)
!====================================================================
!      LAPACK SUBROUTINE TO SOLVE EIGENVALUE PROBLEM - zGGEV
!====================================================================
     IMPLICIT NONE

     INTEGER :: n, INFO, I, LWORK
     COMPLEX*16 :: a(n,n), b(n,n), ALPHA(n), BETA(n) 
     COMPLEX*16 :: VL(n,n), VR(n,n), WR(n), WORK(8*n)
     DOUBLE PRECISION :: RWORK(8*n)
     CHARACTER*1 :: V
     
     LWORK = 8*n

     CALL ZGGEV('V', 'V', n, a, n, b, n, ALPHA, BETA, VL, n, VR, n, WORK, LWORK, RWORK, INFO)

     IF (INFO .EQ. 0) THEN
       PRINT*, 'Successful exit'
     ELSE
       PRINT*, 'Error'
     END IF
      
      DO I = 1, N
        WR(I) = ALPHA(I)/BETA(I)
      END DO
      
     RETURN

     END SUBROUTINE
