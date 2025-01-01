!****************************************************************************
! *   FILE         = qr_hessenbergmatrix.f90                                *
! *   AUTHOR       = CHANDAN KUMAR (CHANDANKR@IITKGP.AC.IN)                 *
! *   INSTITUTE    = INDIAN INSTITUTE OF TECHNOLOGY (IIT), KHARAGPUR        *
!****************************************************************************
! *  THIS PROGRAM IS DISTRIBUTED IN A HOPE THAT IT WILL BE USEFUL.          *
!****************************************************************************
      PROGRAM QR
      IMPLICIT NONE
      INTEGER, PARAMETER :: N=5
      DOUBLE PRECISION :: A(N,N), D(N), V(N,N)
      DOUBLE PRECISION :: WR(N), WI(N)
      INTEGER :: NROT, I, J
      
       

      A(1,:) = (/1.0D0,5.0D0,1.0D0,2.0D0,3.0D0/)
      A(2,:) = (/5.0D0,2.0D0,13.0D0,5.0D0,6.0D0/)
      A(3,:) = (/0.0D0,13.0D0,4.0D0,1.0D0,5.0D0/)
      A(4,:) = (/0.0D0,0.0D0,1.0D0,1.0D0,7.0D0/)
      A(5,:) = (/0.0D0,0.0D0,0.0D0,7.0D0,3.0D0/)
      
      CALL HQR(A, N, WR, WI)
      
      DO I =1, N
          PRINT*, WI(I), WR(I)
      END DO
      
      END PROGRAM
      
      SUBROUTINE HQR(A,N,WR, WI)
      IMPLICIT NONE
      INTEGER :: N
      DOUBLE PRECISION :: A(N,N), WI(N), WR(N)
      INTEGER :: I, ITS, J,K , L, M, NN
      DOUBLE PRECISION :: ANORM, P, Q, R, S, T, U, V, W, X, Y, Z
      ANORM = 0.0D0
      DO I = 1, N
        DO J = MAX(I-1,1), N
            ANORM = ANORM+DABS(A(I,J))
        END DO
     END DO
     
     NN = N
     T = 0.0D0
1     IF (NN .GE. 1) THEN
        ITS = 0
2        DO L = NN , 2, -1
            S = DABS(A(L-1,L-1))+DABS(A(L,L))
            IF (S .EQ. 0.0D0) S = ANORM
            IF (DABS(A(L,L-1)) + S .EQ. S) GOTO 3
         END DO
         L =1
3        X = A(NN,NN)
        IF (L .EQ. NN) THEN
            WR(NN) = X+T
            WI(NN) = 0.0D0
            NN = NN-1
         ELSE
            Y = A(NN-1, NN-1)
            W = A(NN, NN-1)*A(NN-1,NN)
            IF (L .EQ. NN-1) THEN
                P = 0.5D0*(Y-X)
                Q = P**2+W
                Z = DSQRT(DABS(Q))
                X = X+T
                IF (Q  .GE. 0.0D0) THEN
                   Z = P +DSIGN(Z,P)
                   WR(NN) =X+Z
                   WR(NN-1) = WR(NN)
                   IF (Z .NE. 0.0D0)WR(NN) = X -W/Z
                   WI(NN) = 0.0D0
                   WI(NN-1) = 0.0D0
                ELSE
                    WR(NN) = X+P
                    WR(NN-1) = WR(NN)
                    WI(NN) = Z
                    WI(NN-1) = -Z
                END IF
                NN = NN-2
             ELSE
                IF (ITS .EQ. 30) THEN
                   STOP
                  PRINT*, 'TOO MANY ITERATION'
                END IF
                IF (ITS .EQ. 10 .OR. ITS .EQ. 20) THEN
                    T = T+X
                     DO I =1, NN
                         A(I,I) = A(I,I)-X
                     END DO
                     S = DABS(A(NN, NN-1))+DABS(A(NN-1, NN-2))
                     X = 0.75D0*S
                     Y = X
                     W = -0.4375*S**2
                 END IF
                 ITS = ITS+1
                 
                 DO M =NN-2, 1, -1
                    Z = A(M,M)
                    R = X-Z
                    S = Y-Z
                    P =(R*S-W)/A(M+1,M)+A(M,M+1)
                    Q = A(M+1, M+1)-Z-R-S
                    R = A(M+2, M+1)
                    S = DABS(P) + DABS(Q) + DABS(R)
                    P = P/S
                    Q = Q/S
                    R = R/S
                    IF (M.EQ.1) GOTO 4
                    U = DABS(A(M,M-1))*(DABS(Q)+DABS(R))
                    V = DABS(P)*(DABS(A(M-1,M-1))+DABS(Z) + DABS(A(M+1, M+1)))
                    IF (U+V .EQ. V) GOTO 4
4                 END DO
                 DO I =M+2, NN
                   A(I,I-2)  = 0.0D0
                   IF (I .NE. M+2) A(I,I-3) = 0.0D0
                 END DO
               
                 DO K =M, NN-1
                  
                    IF (K .NE. M) THEN
                       P = A(K,K-1)
                       Q = A(K+1, K-1)
                       R = 0.0D0
                       IF (K .NE. NN-1) R= A(K+2, K-1)
                       X =DABS(P) + DABS(Q) +DABS(R)
                       IF (X .NE. 0.0D0) THEN
                           P = P/X
                           Q = Q/X 
                           R = R/X
                       END IF
                    END IF
                    S = DSIGN(DSQRT(P**2+ Q**2+R**2),P)
                    IF (S .NE. 0.0D0) THEN
                       IF (K .EQ. M) THEN
                         IF (L .NE. M) A(K,K-1) = -A(K,K-1)
                       ELSE
                         A(K,K-1)=-S*X
                       END IF
                       P =P+S
                       X=P/S
                       Y=Q/S
                       Z=R/S
                       Q=Q/P
                       R=R/P
                       DO J = K, NN
                          P = A(K,J)+Q*A(K+1, J)
                          IF (K .NE. NN-1) THEN
                             P = P+R*A(K+2,J)
                             A(K+2, J) = A(K+2,J)-P*Z
                          END IF
                          A(K+1, J) = A(K+1,J)-P*Y
                          A(K,J) = A(K,J)-P*X
                       END DO
                       DO I =1, MIN(NN, K+3)
                          P = X*A(I,K)+Y*A(I,K+1)
                          IF ( K .NE. NN-1) THEN
                              P = P+Z*A(I,K+2)
                              A(I,K+2) = A(I,K+2)-P*R
                          END IF
                          A(I,K+1) = A(I,K+1)-P*Q
                          A(I,K) = A(I,K)-P
                       END DO
                       
                     END IF
                  END DO
                  GOTO 2
                END IF
              END IF
            GOTO 1
            END IF
            RETURN
            END SUBROUTINE
                 
                       
                         
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
