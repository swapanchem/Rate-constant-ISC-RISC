       PROGRAM RISC_DSO_T1_S1
       IMPLICIT NONE
       INTEGER, PARAMETER :: N = 159  ! N IS THE NUMBER OF NORMAL MODES 
       INTEGER, PARAMETER :: N1  = 2*N
 
       REAL*8, DIMENSION (N,N) :: RJ, RJT, OMEGAS, OMEGAT, &
                                      BS, AS, AT1,P,P2
       REAL*8, DIMENSION (N,1) :: WS, WT, D, D1
       REAL*8, DIMENSION (1,N) :: DT
       COMPLEX*16, DIMENSION (N,1) :: D_COMPLEX, V1
       COMPLEX*16, DIMENSION (1,N) :: DT_COMPLEX, V3
       COMPLEX*16, DIMENSION (N,N) :: A, B, AT, BT, BB, &
                                          B_INV, AA, B2, AA_INV, E2, &
                                          E1, AB_INV, RJ_COMPLEX, &
                                          RJT_COMPLEX, X1, X2, X3, &
                                          X4, X5,Y1,BS_COMPLEX, &
                                          OMEGAT_COMPLEX,AS_COMPLEX, &
                                          E, BT2, AT2, P2_COMPLEX, BBB

       COMPLEX*16, DIMENSION (N1,N1) :: AK, AK_INV
       COMPLEX*16, DIMENSION (N1,1) :: F
       COMPLEX*16, DIMENSION (1,N1) :: FT, V2
       REAL*8 :: TMIN, TMAX, HH, CM_AU, S_AU, T, SOC, AMU_AU, & 
                 CONST, AU_HZ, BETA, KB, KBT, TEMP, EV_AU, PF, &
                 AMP, DET_I, DET_R, DELE, NUM1R, NUM2R, NUM1I, &
                 NUM2I, NUMR, NUMI, K_RISC_DSO1_REAL, ETA, &
                 K_RISC_DSO1, THETA, T1, SCAL, DET_P2
       COMPLEX*16 :: DET, RISC_DSO, NUM
       COMPLEX*16, DIMENSION (1,1) :: NUM1, NUM2, NUM3, NUM4
       INTEGER :: I, J, K, M, K1, M1, INFO, M2
       INTEGER, DIMENSION (N1) :: IPIV
       COMPLEX*16, DIMENSION (N1) :: WORK
      
!...............................................................................................................................................\\

! REQUIRED INPUT FILES

       !OPEN (30, FILE = 'DUSCHINSKY_MATRIX.DAT')      ! FILE FOR DUSCHINSKY MATRIX BETWEEN INITIAL AND FINAL ELECTRONIC STATE 
       !OPEN (3, FILE = 'WS1.DAT')             ! FILE FOR FINAL ELECTRONIC STATE FREQUENCY
       !OPEN (2, FILE = 'WT1.DAT')             ! FILE FOR INITIAL ELECTRONIC STATE FREQUENCY 
       !OPEN (10, FILE = 'DISPLACEMENT_VECTOR.DAT')    ! FILE FOR DISPLACEMENT VECTOR BETWEEN INITIAL AND FINAL ELECTRONIC STATES

! OUTPUT FILE FOR STORING THE TIME DEPENDENT CORRELATION FUNCTION

       OPEN (34, FILE = 'RISC_DSO_T1_S1.TXT')              

!...................................................................................................................................................\\

! INPUT PARAMETERS

       WRITE (*,*) "ENTER THE NUMBER OF TIME INTERVAL"
       READ (*,*) M1

       WRITE (*,*) "ENTER THE LOWER LIMIT OF TIME IN SEC"
       READ (*,*) TMIN

       WRITE (*,*) "ENTER THE UPPER LIMIT OF TIME IN SEC"
       READ (*,*) TMAX

       WRITE (*,*) "ENTER THE TEMPERATURE IN KELVIN &
                   (NOT LESS THAN 20K)"
       READ (*,*) TEMP

       WRITE (*,*) "ENTER THE DAMPING PARAMETER IN CM^(-1)"
       READ (*,*) ETA

       WRITE (*,*) "ENTER THE ENERGY GAP BETWEEN S1 AND T1 IN eV"
       READ (*,*) DELE

       WRITE (*,*) "ENTER THE SPIN-ORBIT COUPLING BETWEEN S1 AND T1 IN &
                    CM^(-1)"
       READ (*,*) SOC

!........................................................................................................................................................\\

       !READ THE INPUT FILES

       READ (30,*) ((RJ(I, J), J = 1, N), I = 1, N)       !RJ = DUSCHINSKY ROTATION MATRIX     
       READ (3,*) (WS(I, 1), I = 1, N)                !WS = FREQUENCY VECTOR OF SINGLET STATE
       READ (2,*) (WT(I, 1), I = 1, N)                !WT = FREQUENCY VECTOR OF TRIPLET STATE
       READ (10,*) (D(I, 1), I = 1, N)                !D = DISPLACEMENT VECTOR 


!..........................................................................................................................................................\\

       !CONVERSION FACTOR TO ATOMIC UNIT

       AU_HZ = 6.579D0 * ( 10.0D0 ** 15.0D0 )
       CM_AU = 4.55634D0 * ( 10.0D0 ** (-6.0D0 ))
       S_AU = 0.4134D0 * ( 10.0D0 ** 17.0D0 )
       AMU_AU = ( 1.82289D0 * ( 10.0D0 ** 3.0D0 ))
       AMU_AU = ( ( AMU_AU ) ** 0.5D0 )
       EV_AU = 0.036749844D0
       KB = 1.3806452D0 * ( 10.0D0 ** (-23.0D0))
       KBT = KB * TEMP
       KBT = KBT * (6.242D0 * (10.0D0 ** (18.0D0)))    !JOULE To EV
       KBT = KBT * EV_AU                           !EV TO AU
       BETA = (1.0D0 / KBT) 

!DELE = ENERGY GAP BETWEEN S1 T1

       DELE = DELE * EV_AU

! SOC = SPIN-ORBIT COUPLING BETWEEN S1 AND T1

       SOC = SOC * CM_AU


       DO I = 1, N
          WS(I, 1) = WS(I, 1) * CM_AU     
          WT(I, 1) = WT(I, 1) * CM_AU
       ENDDO 

       TMAX = TMAX * S_AU
       TMIN = TMIN * S_AU
       ETA = ETA * CM_AU


!...........................................................................................................................................\\

!GENERATE THE DIAGONAL FREQUENCY MATRIX

       DO I = 1, N
         DO J = 1, N
           IF (I .EQ. J) THEN
           OMEGAS(I, J) = WS(I, 1)
           OMEGAT(I, J) = WT(I, 1)
           OMEGAT_COMPLEX(I, J) = COMPLEX(WT(I,1), 0.0D0)
           ELSE
           OMEGAS(I,J) = 0.0D0
           OMEGAT(I,J) = 0.0D0
           OMEGAT_COMPLEX(I,J) = (0.0D0, 0.0D0)
           ENDIF
         ENDDO
       ENDDO


!GENERATE THE TRANSPOSE OF RJ AND D

       DO I = 1, N
          DO J = 1, N
          RJT(J, I) = RJ(I, J)                                      !RJT = TRANSPOSE OF DUSCHINSKY ROTATION MATRIX RJ
          RJ_COMPLEX(I, J) = COMPLEX(RJ(I, J), 0.0D0)
          RJT_COMPLEX(J, I) = COMPLEX(RJT(J, I), 0.0D0)
          ENDDO
       ENDDO   


!GENERATE THE TRANSPOSE OF THE DISPLACEMENT VECTOR

       DO I = 1, N
         DT(1, I) = D(I, 1)                                          !DT=TRANSPOSE OF THE DISPLACEMENT VECTOR
         D_COMPLEX(I, 1) = COMPLEX(D(I, 1), 0.0D0)
         DT_COMPLEX(1, I) = COMPLEX(DT(1, I), 0.0D0)
       ENDDO
      
         DO I = 1, N
          DO J = 1, N
          IF (I .EQ. J) THEN
                  P(I, J) = 2.0D0 * SINH(OMEGAT(I, J) * BETA / 2.0D0)
          ELSE
                  P(I, J) = 0.0D0
          ENDIF
          ENDDO
       ENDDO

       CALL DGEMM ('N', 'N', N, N ,N ,1.0D0, P, N, P, N, &
                0.0D0, P2, N)

        DET_P2 = 1.0D0
        DO I = 1, N
          DET_P2 = DET_P2 * P2(I, I)
        ENDDO

     
       DO I = 1, N
          DO J = 1, N
          P2_COMPLEX(I, J) = COMPLEX(P2(I, J), 0.0D0)
          ENDDO
       ENDDO



       CONST = (SOC * SOC)                  !FOR DSO
       
!........................................................................................................................................................\\

       !GENERATE THE REQUIRED MATRIX AND MATRIX MULTIPLICATION WITHIN
       !THE TIME LOOP FOR DETERMINANT CALCULATION IN THE FC REGION

        K_RISC_DSO1 = (0.0D0, 0.0D0)  
        HH = ( TMAX - TMIN ) / FLOAT(M1)  !HH=GAP BETWEEN THE TWO SUCCESSIVE TIME INTERVAL
        DO K1 = 1, M1 + 1                 ! STARTING THE TIME LOOP 
        M2 = (M1 / 2) + 1                 !M2 IS THE POINT AT T=0 IN BETWEEN POSITIVE AND NEGATIVE TIME INTERVAL
        IF (K1 .NE. M2) THEN
        T = TMIN + (HH * (K1 - 1))
        ELSE 
        T = 1.0D0*( 10.0D0 ** (- 40.0D0))
        T = T * S_AU
        ENDIF

          DO I = 1, N
           DO J = 1, N
           IF (I .EQ. J) THEN
                  
!BS=OMEGAS*TAN(OMEGAS*T/2)

          BS(I, J) = OMEGAS(I, J) / TAN(T * OMEGAS(I, J))
          BS_COMPLEX(I, J) = COMPLEX(BS(I, J), 0.0D0)
   
          
!AS=OMEGAS*SIN(OMEGAS*T/2)

          AS(I, J) = OMEGAS(I, J) / SIN(T * OMEGAS(I, J))
          AS_COMPLEX(I, J) = COMPLEX(AS(I, J), 0.0D0)

!FORMATION OF AT (EXCLUDING THE WT , WHICH CAN BE MULTIPLIED LATER)

          AT1(I, J) = (COSH( - 2.0D0 * BETA * OMEGAT(I, J))) - (COS(- & 
          2.0D0 * T * OMEGAT(I, J)))
          AT(I, J) = COMPLEX((( 2.0D0 * SIN( -OMEGAT( I, J) * T) * &
           COSH(-BETA * &
          OMEGAT(I, J))) / AT1(I, J)), (( -2.0D0 * COS(- T * &
           OMEGAT(I, J))* &
          SINH( - BETA * OMEGAT(I, J))) / AT1(I, J)))

!FORMATION OF BT (EXCLUDING THE WT , WHICH CAN BE MULTIPLIED LATER)

          BT (I, J) = COMPLEX((( SIN( - 2.0D0 * OMEGAT(I, J) * T )) / &
          AT1(I, J)), &
          (( - SINH( - 2.0D0 * BETA * OMEGAT(I, J))) / AT1(I, J)))
          ELSE
          AS(I, J) = 0.0D0
          AS_COMPLEX(I, J) = (0.0D0, 0.0D0)
          BS(I, J) = 0.0D0
          BS_COMPLEX(I, J) = (0.0D0, 0.0D0)
          AT(I, J) = (0.0D0, 0.0D0)
          AT1(I, J) = 0.0D0
          BT(I, J) = (0.0D0, 0.0D0)
          ENDIF
          ENDDO
       ENDDO



!FORMATION OF K MATRIX AND ITS INVERSE. HERE, WE DENOTE IT AS AK
      
      !FORMATION OF A AND B MATRIX

      CALL ZGEMM ( 'N', 'N', N, N, N, (1.0D0, 0.0D0), OMEGAT_COMPLEX, &
                  N, AT, N,(0.0D0, 0.0D0), AT2, N )
      CALL ZGEMM ('N', 'N', N, N, N, (1.0D0, 0.0D0), OMEGAT_COMPLEX, &
                 N, BT, N, (0.0D0, 0.0D0), BT2, N )
      CALL ZGEMM ('N', 'N', N, N, N, (1.0D0, 0.0D0), RJT_COMPLEX, &
                 N, AT2, N, (0.0D0, 0.0D0), X1, N )

      CALL ZGEMM ('N', 'N', N, N, N, (1.0D0, 0.0D0), X1, N, &
                 RJ_COMPLEX, N, (0.0D0, 0.0D0), X2, N )

      CALL ZGEMM ('N', 'N', N, N, N, (1.0D0, 0.0D0), RJT_COMPLEX, N, &
                BT2, N, (0.0D0, 0.0D0), X3, N )

      CALL ZGEMM ('N', 'N', N, N, N, (1.0D0, 0.0D0), X3, N, &
                 RJ_COMPLEX, N, (0.0D0, 0.0D0), X4, N )

      DO I = 1, N
         DO J = 1, N
         A(I,J) = (AS_COMPLEX(I, J) + X2(I, J))
         B(I,J) = (BS_COMPLEX(I, J) + X4(I, J))
         ENDDO
      ENDDO
      
      !INVERSE OF B 

      DO I = 1, N
        DO J = 1, N
        B_INV(I, J) = B(I, J)
        ENDDO
      ENDDO

      CALL ZGETRF (N, N, B_INV, N, IPIV, INFO)
      CALL ZGETRI (N, B_INV, N, IPIV, WORK, N, INFO)

!FORMATION OF THE DENOMINATOR MATRIX WITHIN THE SQUARE ROOT OF THE
!DETERMINANT WHICH IS DENOTED HERE AS AA

      CALL ZGEMM ('N', 'N', N, N, N, (1.0D0, 0.0D0), A, N, B_INV, N, &
              (0.0D0, 0.0D0), AB_INV, N)
      CALL ZGEMM ('N', 'N', N, N, N, (1.0D0, 0.0D0), AB_INV, N, A, N, &
              (0.0D0, 0.0D0), E1, N)
      CALL ZGEMM ('N', 'N', N, N, N, (1.0D0, 0.0D0), B, N, E1, N, &
              (0.0D0, 0.0D0), E2, N)
      CALL ZGEMM ('N', 'N', N, N, N, (1.0D0, 0.0D0), B, N, B, N, &
              (0.0D0, 0.0D0), B2, N)

      DO I = 1, N
        DO J = 1, N
        AA(I, J) = (B2(I, J) - E2(I, J))
        ENDDO
      ENDDO

      !INVERSE OF AA MATRIX (INVERSE OF THE DENOMINATOR WITHIN THE SQUARE
      !ROOT OF THE DETERMINANT)

      DO I = 1, N
         DO J = 1, N
         AA_INV(I, J) = AA(I, J)
         ENDDO
      ENDDO

      CALL ZGETRF (N, N, AA_INV, N, IPIV, INFO)
      CALL ZGETRI (N, AA_INV, N, IPIV, WORK, N, INFO)

     

      !FORMATION OF K MATRIX

      AK(1 : N, 1 : N) = B
      AK(1 : N, N + 1 : N1) = - A
      AK(N + 1 : N1, 1: N) = - A
      AK(N + 1 : N1, N + 1 : N1) = B

      !INVERSE OF AK 
      !STORED THE AK MATRIX AS AK_INV FOR INVERSION CALCULATION USING BLAS/LAPACK

      DO I = 1, N1
         DO J = 1, N1
         AK_INV(I, J) = AK(I, J)
         ENDDO
      ENDDO
      CALL ZGETRF (N1, N1, AK_INV, N1, IPIV, INFO)
      CALL ZGETRI (N1, AK_INV, N1, IPIV, WORK, N1, INFO)


!XX=NUMERATOR MATRIX WITHIN THE SQUARE ROOT OF THE DETERMINANT

      CALL ZGEMM ('N', 'N', N, N, N, (1.0D0, 0.0D0), AS_COMPLEX, N, &
                AT2, N, (0.0D0, 0.0D0), X5, N)
   

      !FORMATION OF THE COMPLETE DETERMINANT DENOTED BY BB

      CALL ZGEMM ('N', 'N', N, N, N, (1.0D0, 0.0D0), X5, N, AA_INV, N, &
                 (0.0D0, 0.0D0), BB, N)

      CALL ZGEMM ('N', 'N', N, N, N, (1.0D0, 0.0D0), BB, N, &
                P2_COMPLEX, N, (0.0D0, 0.0D0), BBB, N)

!DETERMINANAT OF BB WHICH IS N BY N MATRIX
     
      CALL DETERMINANT (N, BBB, DET)

!AMPLITUDE AND PHASE OF THE DETERMINANT 

      AMP = ABS(DET)
      AMP = SQRT(AMP)
      DET_I = AIMAG(DET)     !IMAGINARY PART OF DET
      DET_R = REAL(DET)      !REAL PART OF DET
      THETA = DET_I / DET_R    !THETA=PHASE FACTOR
      THETA = THETA / 2.0D0
      WRITE (* , *) T, DET

!........................................................................................................................................................\\

!........................................................................................................................................................\\

!SIMPLIFICATION OF THE EXPONENTIAL PART 

      !FORMATION OF E MATRIX OF N BY N 

      DO I = 1, N
        DO J = 1 , N
        E(I, J) = (BT2(I, J) - AT2(I, J))
        ENDDO
      ENDDO

      !FORMATION OF F VECTOR AND ITS TRANSPOSE FT CONTAINING 2N ELEMENTS

      CALL ZGEMM ('N', 'N', N, N, N, (1.0D0, 0.0D0), RJT_COMPLEX, N, &
                 E, N, (0.0D0, 0.0D0), Y1, N)

      CALL ZGEMM ('N', 'N', N, 1, N, (1.0D0, 0.0D0), Y1, N, &
                 D_COMPLEX, N, (0.0D0, 0.0D0), V1, N)

       
      F(1 : N, 1 : 1) = V1
      F(N + 1 : N1,1 : 1) = V1 

      DO I = 1, N1
        FT(1, I) = F(I, 1)
      ENDDO
     
!MATRIX MULTIPLICATION FOR THE EXPONENTIAL PART
      
      CALL ZGEMM ('N1', 'N1', 1, N1, N1, (1.0D0, 0.0D0), FT, 1, &
                 AK_INV, N1, (0.0D0, 0.0D0), V2, 1)

      CALL ZGEMM ('N1', 'N1', 1, 1, N1, (1.0D0, 0.0D0), V2, 1, F, N1, &
                 (0.0D0, 0.0D0), NUM1, 1)

      
      CALL ZGEMM ('N', 'N', 1, N, N, (1.0D0, 0.0D0), DT_COMPLEX, 1, &
                 E, N, (0.0D0, 0.0D0), V3, 1)

      CALL ZGEMM ('N', 'N', 1, 1, N, (1.0D0, 0.0D0), V3, 1, &
                 D_COMPLEX, N, (0.0D0, 0.0D0), NUM2, 1)

      NUM1R = REAL(NUM1(1, 1)) / 2.0D0
      NUM1I = AIMAG(NUM1(1, 1)) / 2.0D0
      NUM2R = REAL(NUM2(1, 1))
      NUM2I = AIMAG(NUM2(1, 1))
      NUMI = - NUM1R + NUM2R + THETA - (T * DELE)
      NUMR = NUM1I - NUM2I

      !TOTAL PART WITHIN THE EXPONENTIAL 

      NUM = COMPLEX( NUMR, NUMI)

!...........................................................................................................................................................\\      

!FINAL EXPRESSION OF REVERSE INTERSYSTEM CROSSING RATE FOR THE DSO PART

       T1 = ABS(T)          !ABSLOUTE VALUE OF TIME T

       RISC_DSO = AMP * EXP(NUM) * EXP( - T1 * ETA)

       WRITE (34, * ) T, REAL(RISC_DSO)

       !INTEGRATION USING SIMPSON'S ONE-THIRD RULE WITHIN TIME DOMAIN

      IF ((K1 .EQ. 1) .OR. (K1 .EQ. (M1 + 1))) THEN
          K_RISC_DSO1 = K_RISC_DSO1 + RISC_DSO
          ELSE
              IF (MOD(K1, 2) .EQ. 0.0D0) THEN
              K_RISC_DSO1 = K_RISC_DSO1 + (RISC_DSO * 4.0D0)
              ELSE
              K_RISC_DSO1 = K_RISC_DSO1 + (2.0D0 * RISC_DSO)
              ENDIF
        ENDIF

        ENDDO                          !END LOOP FOR T

       K_RISC_DSO1_REAL = REAL(K_RISC_DSO1)

       K_RISC_DSO1_REAL = (K_RISC_DSO1_REAL * HH * CONST) / 3.0D0

       K_RISC_DSO1_REAL = K_RISC_DSO1_REAL * AU_HZ

       WRITE (*, *) 'K_RISC_DSO1_REAL IS', K_RISC_DSO1_REAL, CONST

       END PROGRAM RISC_DSO_T1_S1


!.................................................................................................................................................................\\

!SUBROUTINE FOR DETERMINANT CALCULATION USING LAPACK-BLAS LIBRARY

      SUBROUTINE DETERMINANT (N, A, DET)
      COMPLEX*16 :: A(N, N)
      COMPLEX*16 :: DET, SGN
      INTEGER :: I, INFO, N
      INTEGER :: IPIV(N)

      CALL ZGETRF (N, N, A, N, IPIV, INFO)
      DET = (1.0D0, 0.0D0)
       DO I = 1, N
          DET = DET*A(I, I)
       ENDDO
      SGN =(1.0D0, 0.0D0)
      DO I = 1, N
      IF (IPIV(I) /= I) THEN
        SGN = - SGN
      ENDIF
      ENDDO
      DET = SGN * DET
      RETURN
      END



