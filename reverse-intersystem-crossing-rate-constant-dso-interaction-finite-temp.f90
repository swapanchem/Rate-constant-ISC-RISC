       PROGRAM RISC_DSO_FINITE_TEMP

               USE DeterminantModule_DSO_RISC  ! IMPORT THE MODULE FOR THE CALCULATION OF DETERMINANT

       IMPLICIT NONE

       REAL*8,DIMENSION(:,:), ALLOCATABLE :: RJ, RJT, OMEGAS, OMEGAT, &
                                             BS, AS, AT1, P, P2
       REAL*8,DIMENSION(:,:), ALLOCATABLE :: WS, WT, D, D1

       REAL*8,DIMENSION(:,:), ALLOCATABLE :: DT

       COMPLEX*16,DIMENSION(:,:), ALLOCATABLE :: D_COMPLEX, V1

       COMPLEX*16,DIMENSION(:,:), ALLOCATABLE :: DT_COMPLEX, V3

       COMPLEX*16,DIMENSION(:,:),ALLOCATABLE :: A, B, AT, BT, BB,  &
                                                B_INV, AA, B2, & 
       AA_INV, E2, E1, AB_INV, RJ_COMPLEX, BBB_COPY, &
       RJT_COMPLEX, X1, X2, X3, X4, X5, Y1, BS_COMPLEX, &
       OMEGAT_COMPLEX, AS_COMPLEX, E, BT2, AT2, P2_COMPLEX, BBB

       COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: AK, AK_INV

       COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: F

       COMPLEX*16, DIMENSION(:,:), ALLOCATABLE  :: FT, V2

       REAL*8 :: TMIN, TMAX, HH, CM_AU, S_AU, T, AMU_AU, & 
                 CONST, AU_HZ, BETA, KB, KBT, TEMP, EV_AU, &
                 PF, AMP, DET_I, DET_R, NUM1R, NUM2R, NUM1I, &
                 NUM2I, NUMR, NUMI, K_RISC_DSO_REAL, ETA, &
                 THETA, T1, DET_P2, DELE_S1_TL, &
                 SOC_S1_TL

       COMPLEX*16 :: DET, RISC_DSO, NUM, K_RISC_DSO

       COMPLEX*16,DIMENSION(:,:), ALLOCATABLE :: NUM1, NUM2

       INTEGER :: I, J, K, M, N, N1, K1, M1, INFO, M2, NTRUNC

       INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV

       COMPLEX*16, DIMENSION(:), ALLOCATABLE :: WORK

       !Character for reading input file data 

       CHARACTER(LEN=100) ::NAM1,NAM2,NAM3,NAM4,NAM5,NAM6,NAM7,NAM8, &
                            NAM0,NAM


       OPEN(30,FILE='input_dso_risc.txt')


       READ (30, *) NAM

       READ (30, *) NAM0
       READ (30, *) N

       READ (30, *) NAM1
       READ (30, *) M1

       READ (30, *) NAM2
       READ (30, *) TMAX

       READ (30, *) NAM3
       READ (30, *) TMIN

       READ (30, *) NAM4
       READ (30, *) TEMP

       READ (30, *) NAM5
       READ (30, *) ETA

       READ (30, *) NAM6
       READ (30, *) DELE_S1_TL


       READ(30,*)NAM7
       READ(30,*)SOC_S1_TL


       READ(30,*)NAM8

       N1 = 2 * N


       OPEN(66,FILE='J.TXT')       ! J is the Duschinsky rotation matrix

       OPEN(77,FILE='WS1.TXT')     ! frequencies of the S1 state (final state)

       OPEN(88,FILE='WTL.TXT')     ! frequencies of the TL state (initial state)

       OPEN(99,FILE='D.TXT')       ! displacement vector
         
       OPEN(119,FILE='RISC_DSO_CORR.TXT')  ! real part of the correlation function accounting DSO interaction


      ! N       NORMAL MODES OF THE INVESTIGATED MOLECULE
      ! M1      M1=NUMBER OF INTERVAL FOR INTEGRATION
      ! ETA    ETA=DAMPING PARAMETER IN CM^(-1)
      ! TEMP    TEMP=TEMPERATURE IN K
      ! TMIN    TMIN=LOWER LIMIT OF TIME IN SEC
      ! TMAX    TMAX=UPPER LIMIT OF TIME IN SEC
       
       

       !INPUT THE MATRIX AND READ THE MATRIX

       ALLOCATE ( RJ(N,N), WS(N,1), WT(N,1), D(N,1) )


       READ(66,*) ((RJ(I,J), J = 1,N),I = 1,N)       !RJ = DUSCHINSKY ROTATION MATRIX     

       READ(77,*) (WS(I,1), I = 1,N)                !WS = FREQUENCY VECTOR OF SINGLET STATE

       READ(88,*) (WT(I,1), I = 1,N)                !WT = FREQUENCY VECTOR OF TRIPLET STATE

       READ (99, *) (D(I, 1), I = 1,N)                !D = DISPLACEMENT VECTOR 


       !CONVERSION FACTOR TO  ATOMIC UNIT

       AU_HZ = 6.579D0*(10.0D0**15.0D0)

       CM_AU = 4.55634D0*(10.0D0**(-6.0D0))

       S_AU = 0.4134D0*(10.0D0**17.0D0)

       AMU_AU = (1.82289D0*(10.0D0**3.0D0))

       AMU_AU = ((AMU_AU)**0.5D0)

       EV_AU = 0.036749844D0

       KB = 1.3806452D0*(10.0D0**(-23.0D0))

       KBT = KB*TEMP

       KBT = KBT*(6.242D0*(10.0D0**(18.0D0)))    !JOULE To EV

       KBT = KBT*EV_AU                           !EV TO AU

       BETA = (1.0D0/KBT) 

!ENERGY GAP BETWEEN THE S1 AND THE LOWEST EXCITED TRIPLET (TL)

       DELE_S1_TL = DELE_S1_TL * EV_AU

!SOC GAP BETWEEN THE S1 AND THE LOWEST EXCITED TRIPLET (TL)

       SOC_S1_TL = SOC_S1_TL * CM_AU

    

       DO I = 1,N

          WS(I,1) = WS(I,1)*CM_AU     

          WT(I,1) = WT(I,1)*CM_AU

       ENDDO 

       TMAX = TMAX*S_AU
       TMIN = TMIN*S_AU
       ETA  = ETA*CM_AU


!GENERATE THE DIAGONAL FREQUENCY MATRIX

       ALLOCATE ( OMEGAS(N,N), OMEGAT(N,N), OMEGAT_COMPLEX(N,N) )


       DO I = 1, N

         DO J = 1, N

         IF (I.EQ.J) THEN

         OMEGAS(I,J) = WS(I,1)

         OMEGAT(I,J) = WT(I,1)

         OMEGAT_COMPLEX(I,J) = COMPLEX(WT(I,1),0.0D0)

         ELSE

         OMEGAS(I,J) = 0.0D0

         OMEGAT(I,J) = 0.0D0

         OMEGAT_COMPLEX(I,J) = (0.0D0,0.0D0)

         ENDIF

         ENDDO

       ENDDO


!GENERATE THE TRANSPOSE OF RJ AND D

       ALLOCATE (RJT(N,N))

       ALLOCATE (RJ_COMPLEX(N,N))

       ALLOCATE (RJT_COMPLEX(N,N))

       DO I = 1,N

          DO J = 1,N

          RJT(J,I) = RJ(I,J)                                      !RJT = TRANSPOSE OF DUSCHINSKY ROTATION MATRIX RJ

          RJ_COMPLEX(I,J) = COMPLEX(RJ(I,J),0.0D0)

          RJT_COMPLEX(J,I) = COMPLEX(RJT(J,I),0.0D0)

          ENDDO

       ENDDO   


!GENERATE THE TRANSPOSE OF THE DISPLACEMENT VECTOR

       ALLOCATE ( DT(1,N), D_COMPLEX(N,1), DT_COMPLEX(1,N) )


       DO I = 1, N

         DT(1, I) = D(I, 1)                                  !DT=TRANSPOSE OF DISPLACEMENT VECTOR

         D_COMPLEX(I, 1) = COMPLEX(D(I, 1), 0.0D0)

         DT_COMPLEX(1, I) = COMPLEX(DT(1, I), 0.0D0)

       ENDDO
      
         ALLOCATE ( P(N,N), P2(N,N), P2_COMPLEX(N,N) )

         DO I = 1, N

          DO J = 1, N

          IF (I.EQ.J) THEN

                  P(I,J) = 2.0D0 * SINH( OMEGAT(I, J) * BETA / 2.0D0 )
          ELSE 

                  P(I,J) = 0.0D0

          ENDIF

          ENDDO

       ENDDO


       CALL DGEMM( 'N','N',N, N, N, 1.0D0, P, N, P, N, &
                0.0D0, P2, N )


        DET_P2 = 1.0D0

        DO I = 1, N

         DET_P2 = DET_P2 * P2(I, I)

         ENDDO



       DO I = 1,N

         DO J = 1,N

         P2_COMPLEX(I, J) = COMPLEX(P2(I, J), 0.0D0)

         ENDDO

       ENDDO


       CONST = ( SOC_S1_TL * SOC_S1_TL )                    !SQUARE OF THE SOC


!........................................................................................................................................................\\
       !GENERATE THE REQUIRED MATRIX AND MATRIX MULTIPLICATION WITHIN
       !THE TIME LOOP FOR DETERMINANT CALCULATION IN THE FC REGION

        K_RISC_DSO = ( 0.0D0, 0.0D0 )  


        HH = ( TMAX - TMIN ) / FLOAT( M1 )  !HH=GAP BETWEEN THE TWO SUCCESSIVE POINTS

        DO K1 = 1,M1 + 1

          M2 = (M1 / 2) + 1               !M2 IS THE POINT AT T=0 IN BETWEEN POSITIVE AND NEGATIVE TIME INTERVAL

        IF (K1.NE.M2) THEN

        T = TMIN + (HH * (K1 - 1))

        ELSE

        T = 1.0D0 * ( 10.0D0 ** ( - 24.0D0 ))

        T = T * S_AU

        ENDIF

          ALLOCATE( BS(N,N), BS_COMPLEX(N,N), AS(N,N), &
                    AS_COMPLEX(N,N), AT1(N,N), AT(N,N), BT(N,N)) 

          DO I = 1,N

          DO J = 1,N

          IF (I.EQ.J) THEN
                  
!BS=OMEGAS*TAN(OMEGAS*T/2)


          BS(I, J) = OMEGAS(I, J) / TAN(T * OMEGAS( I, J ))

          BS_COMPLEX(I, J) = COMPLEX(BS(I, J), 0.0D0)
   
          
!AS=OMEGAS*SIN(OMEGAS*T/2)


          AS(I, J) = OMEGAS(I, J) / SIN(T * OMEGAS(I, J))

          AS_COMPLEX(I, J) = COMPLEX(AS(I, J), 0.0D0)


!FORMATION OF AT (EXCLUDING THE WT , WHICH CAN BE MULTIPLIED LATER)



          AT1(I, J) = (COSH( - 2.0D0 * BETA * OMEGAT(I, J))) - &
          (COS( - 2.0D0 * T * OMEGAT(I, J)))

          AT(I, J) = COMPLEX(((2.0D0 * SIN( - OMEGAT(I, J) * T) * &
          COSH( - BETA * OMEGAT(I, J))) / AT1(I, J)),(( - 2.0D0 * &
          COS( -T * OMEGAT(I, J)) * SINH( - BETA * OMEGAT(I, J))) &
          / AT1(I, J)))

!FORMATION OF BT (EXCLUDING THE WT , WHICH CAN BE MULTIPLIED LATER)

          BT(I, J) = COMPLEX(((SIN( - 2.0D0 * OMEGAT(I, J) * T)) / &
          AT1(I, J)),(( - SINH( - 2.0D0 * BETA * OMEGAT(I, J))) / &
          AT1(I, J)))

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

!       WRITE(*,*)'CONST2 IS', CONST1

!FORMATION OF K MATRIX AND ITS INVERSE. HERE, WE DENOTE IT AS AK
      
      !FORMATION OF A AND B MATRIX

      ALLOCATE(AT2(N,N), BT2(N,N), X1(N,N), X2(N,N), X3(N,N), X4(N,N))


      CALL ZGEMM('N','N',N, N, N, (1.0D0, 0.0D0), OMEGAT_COMPLEX, N, &
                 AT, N, (0.0D0,0.0D0), AT2, N)


      CALL ZGEMM('N','N',N, N, N, (1.0D0, 0.0D0), OMEGAT_COMPLEX, N, &
                 BT, N, (0.0D0, 0.0D0), BT2, N)


      CALL ZGEMM('N','N',N, N, N, (1.0D0, 0.0D0), RJT_COMPLEX, N,  &
                AT2, N, (0.0D0, 0.0D0), X1, N)


      CALL ZGEMM('N','N',N, N, N, (1.0D0, 0.0D0), X1, N, RJ_COMPLEX, &
                 N, (0.0D0, 0.0D0), X2, N)


      CALL ZGEMM('N','N',N, N, N, (1.0D0, 0.0D0), RJT_COMPLEX, N, &
                BT2, N, (0.0D0, 0.0D0), X3, N)


      CALL ZGEMM('N','N',N, N, N, (1.0D0, 0.0D0), X3, N, RJ_COMPLEX, &
                 N, (0.0D0, 0.0D0), X4, N)

   
      
     ALLOCATE(A(N,N), B(N,N), B_INV(N,N),IPIV(N1), WORK(N1), &
             AB_INV(N,N) )  


      DO I = 1, N

       DO J = 1, N

        A(I, J) = (AS_COMPLEX(I, J) + X2(I, J))

        B(I, J) = (BS_COMPLEX(I, J) + X4(I, J))

        ENDDO

      ENDDO
      
      !INVERSE OF B 
 


      DO I = 1, N

       DO J = 1, N

       B_INV(I, J) = B(I, J)

       ENDDO

      ENDDO


      CALL ZGETRF(N, N, B_INV, N, IPIV, INFO)

      CALL ZGETRI(N, B_INV, N, IPIV, WORK, N, INFO)


!FORMATION OF THE DENOMINATOR MATRIX WITHIN THE SQUARE ROOT OF THE
!DETERMINANT WHICH IS DENOTED HERE AS AA


      ALLOCATE( E1(N,N), E2(N,N), B2(N,N), AA(N,N), AA_INV(N,N) )



      CALL ZGEMM('N','N',N, N, N, (1.0D0, 0.0D0), A, N, B_INV, &
                N, (0.0D0, 0.0D0), AB_INV, N)


      CALL ZGEMM('N','N',N, N, N, (1.0D0, 0.0D0), AB_INV, N, A, &
                N, (0.0D0, 0.0D0), E1, N)


      CALL ZGEMM('N','N',N, N, N, (1.0D0, 0.0D0), B, N, E1, N, &
              (0.0D0, 0.0D0), E2, N)


      CALL ZGEMM('N','N',N, N, N, (1.0D0, 0.0D0), B, N, B, N, &
              (0.0D0,0.0D0), B2, N)



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



      CALL ZGETRF(N, N, AA_INV, N, IPIV, INFO)


      CALL ZGETRI(N, AA_INV, N, IPIV, WORK, N, INFO)


      !FORMATION OF K MATRIX
   
      ALLOCATE( AK(N1,N1), AK_INV(N1,N1) )

      AK(1 : N, 1 : N) = B

      AK( 1 : N, N + 1 : N1) = - A
      AK( N + 1 : N1, 1 : N) = - A
      AK (N + 1 : N1, N + 1 : N1) = B

      !INVERSE OF AK 
      !STORED THE AK MATRIX AS AK_INV FOR INVERSION CALCULATION USING BLAS/LAPACK


      DO I = 1, N1

        DO J = 1, N1

        AK_INV(I, J) = AK(I, J)

        ENDDO
       
      ENDDO



      CALL ZGETRF(N1, N1, AK_INV, N1, IPIV, INFO)

      CALL ZGETRI(N1, AK_INV, N1, IPIV, WORK, N1, INFO)



!XX=NUMERATOR MATRIX WITHIN THE SQUARE ROOT OF THE DETERMINANT

      ALLOCATE (X5(N,N), BB(N,N), BBB(N,N), BBB_COPY(N,N) )

      CALL ZGEMM('N','N',N, N, N, (1.0D0, 0.0D0), AS_COMPLEX, N, &
                AT2, N, (0.0D0, 0.0D0), X5, N)
   

      !FORMATION OF THE COMPLETE DETERMINANT DENOTED BY BB

      CALL ZGEMM('N','N',N, N, N, (1.0D0, 0.0D0), X5, N, AA_INV, N, &
              (0.0D0, 0.0D0), BB, N)


      CALL ZGEMM('N','N',N, N, N, (1.0D0, 0.0D0), BB, N, P2_COMPLEX, &
                 N, (0.0D0, 0.0D0), BBB, N)


       DO I = 1, N

         DO J = 1, N

         BBB_COPY(I, J) = BBB(I, J)

         ENDDO

       ENDDO
!DETERMINANAT OF BB WHICH IS N BY N MATRIX
     
      CALL DETERMINANT(N, BBB_COPY, DET)


!AMPLITUDE AND PHASE OF THE DETERMINANT 

      AMP = ABS( DET ) 

      AMP = SQRT( AMP )

      DET_I = AIMAG( DET )     !IMAGINARY PART OF DET

      DET_R = REAL( DET )      !REAL PART OF DET

      THETA = DET_I / DET_R    !THETA=PHASE FACTOR

      THETA = THETA / 2.0D0


!........................................................................................................................................................

!........................................................................................................................................................

!SIMPLIFICATION OF THE EXPONENTIAL PART 

      !FORMATION OF E MATRIX OF N BY N 

      ALLOCATE (E(N,N), Y1(N,N), V1(N,1), F(N1,1), FT(1,N1) )

      DO I = 1, N

       DO J = 1, N

       E(I,J) = (BT2(I, J) - AT2(I, J))

       ENDDO

      ENDDO

      !FORMATION OF F VECTOR AND ITS TRANSPOSE FT CONTAINING 2N ELEMENTS



      CALL ZGEMM('N','N',N, N, N, (1.0D0, 0.0D0), RJT_COMPLEX, N, E, &
                 N, (0.0D0, 0.0D0), Y1, N)


      CALL ZGEMM('N','N',N, 1, N, (1.0D0, 0.0D0), Y1, N, D_COMPLEX, &
                N, (0.0D0, 0.0D0), V1, N)

       
      
      F(1 : N, 1) = V1(:,1)
      F(N + 1 : N1, 1) = V1(:,1)

      DO I = 1, N1

       FT(1, I) = F(I, 1)

      ENDDO
     
!MATRIX MULTIPLICATION FOR THE EXPONENTIAL PART
      

      ALLOCATE (V2(1,N1), NUM1(1,1), V3(1,N), NUM2(1,1) )


      CALL ZGEMM('N','N',1, N1, N1, (1.0D0,0.0D0), FT, 1, AK_INV, N1, &
              (0.0D0, 0.0D0), V2, 1)


      CALL ZGEMM('N','N',1, 1, N1, (1.0D0, 0.0D0), V2, 1, F, N1, &
              (0.0D0, 0.0D0), NUM1, 1)

      
      CALL ZGEMM('N','N',1, N, N, (1.0D0, 0.0D0), DT_COMPLEX, 1, E, N, &
                (0.0D0, 0.0D0), V3, 1)


      CALL ZGEMM('N','N',1, 1, N, (1.0D0, 0.0D0), V3, 1, D_COMPLEX, N, &
            (0.0D0, 0.0D0), NUM2, 1)


      NUM1R = REAL( NUM1(1, 1) ) / 2.0D0

      NUM1I = AIMAG( NUM1(1, 1) ) / 2.0D0

      NUM2R = REAL(NUM2(1, 1))

      NUM2I = AIMAG(NUM2(1, 1))

      NUMI = -NUM1R + NUM2R + THETA - ( T * DELE_S1_TL ) 

      NUMR = NUM1I - NUM2I

      !TOTAL PART WITHIN THE EXPONENTIAL 

      NUM = COMPLEX (NUMR, NUMI)


!FINAL EXPRESSION OF REVERSE INTERSYSTEM CROSSING RATE FOR THE DSO PART

       T1 = ABS( T )          !ABSLOUTE VALUE OF TIME T

       RISC_DSO = AMP * EXP(NUM) * EXP( -T1 * ETA)

       RISC_DSO = RISC_DSO * CONST 

       RISC_DSO = RISC_DSO * (1.0D0 / 3.0D0 )  !AVERAGING DUE TO THREE TRIPLET SUBLEVELS

       WRITE (119, *) T, REAL (RISC_DSO)

       !INTEGRATION USING SIMPSON'S ONE-THIRD RULE WITHIN TIME DOMAIN

      IF ((K1 .EQ.1 ).OR.( K1.EQ.( M1 + 1 ))) THEN

          K_RISC_DSO = K_RISC_DSO + RISC_DSO

          ELSE

              IF (MOD(K1, 2) .EQ. 0.0D0) THEN

              K_RISC_DSO = K_RISC_DSO + (RISC_DSO * 4.0D0)

              ELSE

              K_RISC_DSO = K_RISC_DSO + (2.0D0 * RISC_DSO)

              ENDIF

        ENDIF

       DEALLOCATE(BS, BS_COMPLEX, AS, AS_COMPLEX, AT1, AT, BT, AT2, &
                 BT2, X1, X2, X3, X4, A, B, B_INV, AB_INV, E1, E2,  &
                 B2, AA, AA_INV,AK, AK_INV, X5, BB, BBB, BBB_COPY, &
                 E, Y1, F, FT, V2, V1, NUM1, V3, NUM2, IPIV, WORK ) 
 

        ENDDO                          !END LOOP FOR T

       DEALLOCATE (RJ, WS, WT, D, OMEGAS, OMEGAT, OMEGAT_COMPLEX, &
                   RJT, RJ_COMPLEX, DT, D_COMPLEX, DT_COMPLEX, P, &
                   P2, P2_COMPLEX)

       K_RISC_DSO_REAL = REAL (K_RISC_DSO)

       K_RISC_DSO_REAL = (K_RISC_DSO_REAL * HH ) / 3.0D0

       K_RISC_DSO_REAL =  ( K_RISC_DSO_REAL * AU_HZ ) 

       WRITE (*,*) 'K_RISC_DSO1_REAL IS', K_RISC_DSO_REAL

       END PROGRAM RISC_DSO_FINITE_TEMP



