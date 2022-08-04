       PROGRAM RISC_RATE_CONSTANT
!THIS CODE COMPUTES THE RATE CONSTANT OF RISC CONSIDERING BOTH     
! DIRECT  SOC AND SPIN-VIBRONIC COUPLING FOR ALL NORMAL MODES
       USE MPI
       IMPLICIT NONE
       REAL*8,DIMENSION(159,159)::RJ,RJT,OMEGAS,OMEGAT,BS,AS,AT1,TKK
       REAL*8,DIMENSION(159,1)::WS,WT,D,NAC
       REAL*8,DIMENSION(1,159)::DT,NAC_T
       COMPLEX*16,DIMENSION(159,1)::H1,H2,D_COMPLEX,V1
       COMPLEX*16,DIMENSION(1,159)::DT_COMPLEX,V3
       COMPLEX*16,DIMENSION(159,159)::A,B,AT,BT,BB,B_INV,AA,B2, &
               AA_INV,E2,E1,AB_INV,G11,G12,G21,G22,RJ_COMPLEX, &
               RJT_COMPLEX,X1,X2,X3,X4,X5,Y1,BS_COMPLEX, &
               OMEGAT_COMPLEX,AS_COMPLEX, &
               E,TKK_COMPLEX,BT2,AT2
       COMPLEX*16,DIMENSION(318,318)::AK,AK_INV,G,GK_INV
       COMPLEX*16,DIMENSION(318,1)::F,H,AKINV_F
       COMPLEX*16,DIMENSION(1,318)::FT,AKINV_F_TRANS,AKINV_T_G, &
               H_TRANS,V2
       REAL*8::TMIN,TMAX,HH,CM_AU,S_AU,T,SOC1,AMU_AU,CONST1,AU_HZ, &
               BETA,KB,KBT,TEMP,EV_AU,PF,AMP,DET_I,DET_R, &
               DELE1,NUM1R,NUM2R,NUM1I,NUM2I,NUMR,NUMI, &
               TRACE_GK_INV_REAL,TRACE_GK_INV_IMAG,K_RISC_DSO_REAL,
               CONST2,SOC2,DELE2,THETA,T1
       COMPLEX*16::DET,TRACE_GK_INV,SV,RISC_DSO,RISC_SV,SV2, &
                   NUM,TRACE,K_RISC_DSO,K_RISC_SV,SV1
       COMPLEX*16,DIMENSION(120)::SV1
       COMPLEX*16,DIMENSION(1,1)::NUM1,NUM2,NUM3,NUM4
       INTEGER::I,J,K,M,K1,N,N1,M1,INFO,M2,NP,IERROR,MYID,NP1
       INTEGER,DIMENSION(318)::IPIV
       COMPLEX*16,DIMENSION(318)::WORK

!MPI INITILIALIZATION 

       CALL MPI_INIT(IERROR)

!INITIALIZE THE NUMBER OF PROCESSOR (NP) 
!MYID =  PROCESSORS ID FOR EACH JOB 

       CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NP,IERROR)
       CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERROR)

       OPEN(30,FILE='JMAT-S1-T1-TXO-PhCz-TRANS-TOLUENE-M06.TXT')
       OPEN(3,FILE='WS1-TXO-PhCz-TOLUENE-M06.TXT')
       OPEN(2,FILE='WT1-TXO-PhCz-TOLUENE-M06.TXT')
       OPEN(10,FILE='SHIFT-VECTOR-S1-T1-TXO-PhCz-TOLUENE-M06.TXT')
       OPEN(50,FILE='NAC-T2-T1-NORM-TOLUENE-TXO-PhCz-MWT.TXT')
       OPEN(34,FILE='RISC_DSO.TXT')
       OPEN(41,FILE='RISC_SV.TXT')


       N=159                               !N=NORMAL MODES OF THE INVESTIGATED MOLECULE
       M1=20000                            !M1=NUMBER OF INTERVAL FOR INTEGRATION
       ETA1=60.0D0                         !ETA=DAMPING PARAMETER IN CM^(-1)
       TEMP=300.0D0                        !TEMP=TEMPERATURE IN K
       TMIN=-10.0D0*(10.0D0**(-12.0D0))    !TMIN=LOWER LIMIT OF TIME IN SEC
       TMAX=10.0D0*(10.0D0**(-12.0D0))     !TMAX=UPPER LIMIT OF TIME IN SEC
       
       


       !INPUT THE MATRIX AND READ THE MATRIX

       READ(30,*)((RJ(I,J),J=1,N),I=1,N)       !RJ=DUSCHINSKY ROTATION MATRIX     
       READ(3,*)(WS(I,1),I=1,N)                !WS=FREQUENCY VECTOR OF SINGLET STATE
       READ(2,*)(WT(I,1),I=1,N)                !WT=FREQUENCY VECTOR OF TRIPLET STATE
       READ(10,*)(D(I,1),I=1,N)                !D=DISPLACEMENT VECTOR 
       READ(50,*)(NAC(I,1),I=1,N)              !NAC=NONADIABATIC COUPLING VECTOR 


       !CONVERSION FACTOR TO  ATOMIC UNIT

       AU_HZ=6.579D0*(10.0D0**15.0D0)
       CM_AU=4.55634D0*(10.0D0**(-6.0D0))
       S_AU=0.4134D0*(10.0D0**17.0D0)
       AMU_AU=(1.82289D0*(10.0D0**3.0D0))
       AMU_AU=((AMU_AU)**0.5D0)
       EV_AU=0.036749844D0
       KB=1.3806452D0*(10.0D0**(-23.0D0))
       KBT=KB*TEMP
       KBT=KBT*(6.242D0*(10.0D0**(18.0D0)))    !JOULE To EV
       KBT=KBT*EV_AU                           !EV TO AU
       BETA=(1.0D0/KBT) 

       DELE1=0.21D0         !ENERGY GAP BETWEEN S1-T1
       DELE1=DELE1*EV_AU
       DELE2=0.69D0         !ENERGY GAP BETWEEN T2-T1
       DELE2=DELE2*EV_AU
       SOC1=0.153D0         !SOC FOR S1-T1
       SOC1=SOC1*CM_AU
       SOC2=0.815D0         !SOC FOR S1-T2
       SOC2=SOC2*CM_AU    


       DO I=1,N
          WS(I,1)=WS(I,1)*CM_AU     
          WT(I,1)=WT(I,1)*CM_AU
       ENDDO
       DO I=1,N
          NAC_T(1,I)=NAC(I,1)     !ROW VECTOR OF NAC
      ENDDO

       
       TMAX=TMAX*S_AU
       TMIN=TMIN*S_AU
       ETA=ETA*CM_AU


!GENERATE THE DIAGONAL FREQUENCY MATRIX

       DO I=1,N
         DO J=1,N
         IF(I.EQ.J)THEN
         OMEGAS(I,J)=WS(I,1)
         OMEGAT(I,J)=WT(I,1)
         OMEGAT_COMPLEX(I,J)=COMPLEX(WT(I,1),0.0D0)
         !RJ(I,J)=1.0D0
         ELSE
         OMEGAS(I,J)=0.0D0
         OMEGAT(I,J)=0.0D0
         OMEGAT_COMPLEX(I,J)=(0.0D0,0.0D0)
         !RJ(I,J)=0.0D0
         ENDIF
         ENDDO
       ENDDO

!FORMATION OF THE NONADIABATIC COUPLING MATRIX

       CALL DGEMM('N','N',N,N,1,1.0D0,NAC,N,NAC_T,1, &
                0.0D0,TKK,N)                                  !TKK=NONADIABATIC COUPLING MATRIX 

!GENERATE THE TRANSPOSE OF RJ AND D

       DO I=1,N
          DO J=1,N
          RJT(J,I)=RJ(I,J)                                      !RJT = TRANSPOSE OF DUSCHINSKY ROTATION MATRIX RJ
          RJ_COMPLEX(I,J)=COMPLEX(RJ(I,J),0.0D0)
          RJT_COMPLEX(J,I)=COMPLEX(RJT(J,I),0.0D0)

          TKK(I,J)=(TKK(I,J)*(SOC2*SOC2))/(DELE2*DELE2)          
          TKK_COMPLEX(I,J)=COMPLEX(TKK(I,J),0.0D0)
          ENDDO
       ENDDO   


!GENERATE THE TRANSPOSE OF THE DISPLACEMENT VECTOR

       DO I=1,N
         DT(1,I)=D(I,1)                                                 !DT=TRANSPOSE OF DISPLACEMENT VECTOR
         D_COMPLEX(I,1)=COMPLEX(D(I,1),0.0D0)
         DT_COMPLEX(1,I)=COMPLEX(DT(1,I),0.0D0)
       ENDDO

       PF=0.0D0                                               !PF = PARTITION FUNCTION
       DO I=1,N
         DO J=1,N
         IF(I.EQ.J)THEN
         PF=PF+(EXP(-OMEGAT(I,J)*BETA))                    !PF = SUMMATION OF DIAGONAL ELEMENT OF WT1 ELEMENT
         ELSE
         ENDIF
         ENDDO
       ENDDO

       CONST1=(1.0D0/PF)*(SOC1*SOC1)                  !FOR DSO
       CONST2=(1.0D0/PF)                              !FOR SPIN-VIBRONIC (SV) 
       

!........................................................................................................................................................\\
       !GENERATE THE REQUIRED MATRIX AND MATRIX MULTIPLICATION WITHIN
       !THE TIME LOOP FOR DETERMINANT CALCULATION IN THE FC REGION

        K_RISC_DSO=(0.0D0,0.0D0)  
        K_RISC_SV=(0.0D0,0.0D0)
        HH=(TMAX-TMIN)/FLOAT(M1)  !HH=GAP BETWEEN THE TWO SUCCESSIVE POINTS
        DO K1=1,M1+1
        M2=(M1/2)+1               !M2 IS THE POINT AT T=0 IN BETWEEN POSITIVE AND NEGATIVE TIME INTERVAL
        IF(K1.NE.M2)THEN
        T=TMIN+(HH*(K1-1))
        ELSE 
        T=1.0D0*(10.0D0**(-24.0D0))
        T=T*S_AU
        ENDIF

          DO I=1,N
          DO J=1,N
          IF(I.EQ.J)THEN
                  
!BS=OMEGAS*TAN(OMEGAS*T/2)

          BS(I,J)=OMEGAS(I,J)/TAN(T*OMEGAS(I,J))
          BS_COMPLEX(I,J)=COMPLEX(BS(I,J),0.0D0)
   
          
!AS=OMEGAS*SIN(OMEGAS*T/2)

          AS(I,J)=OMEGAS(I,J)/SIN(T*OMEGAS(I,J))
          AS_COMPLEX(I,J)=COMPLEX(AS(I,J),0.0D0)

!FORMATION OF AT (EXCLUDING THE WT , WHICH CAN BE MULTIPLIED LATER)

          AT1(I,J)=(COSH(-2.0D0*BETA*OMEGAT(I,J)))-(COS(-2.0D0*T &
          *OMEGAT(I,J)))
          AT(I,J)=COMPLEX(((2.0D0*SIN(-OMEGAT(I,J)*T)*COSH(-BETA* &
          OMEGAT(I,J)))/AT1(I,J)),((-2.0D0*COS(-T*OMEGAT(I,J))* &
          SINH(-BETA*OMEGAT(I,J)))/AT1(I,J)))

!FORMATION OF BT (EXCLUDING THE WT , WHICH CAN BE MULTIPLIED LATER)

          BT(I,J)=COMPLEX(((SIN(-2.0D0*OMEGAT(I,J)*T))/AT1(I,J)), &
          ((-SINH(-2.0D0*BETA*OMEGAT(I,J)))/AT1(I,J)))
          ELSE
          AS(I,J)=0.0D0
          AS_COMPLEX(I,J)=(0.0D0,0.0D0)
          BS(I,J)=0.0D0
          BS_COMPLEX(I,J)=(0.0D0,0.0D0)
          AT(I,J)=(0.0D0,0.0D0)
          AT1(I,J)=0.0D0
          BT(I,J)=(0.0D0,0.0D0)
          ENDIF
          ENDDO
       ENDDO



!FORMATION OF K MATRIX AND ITS INVERSE. HERE, WE DENOTE IT AS AK
      
      !FORMATION OF A AND B MATRIX

      CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),OMEGAT_COMPLEX,N,AT,N, &
                (0.0D0,0.0D0),AT2,N)
      CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),OMEGAT_COMPLEX,N,BT,N, &
                (0.0D0,0.0D0),BT2,N)
      CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),RJT_COMPLEX,N,AT2,N, &
                (0.0D0,0.0D0),X1,N)

      CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),X1,N,RJ_COMPLEX,N, &
                (0.0D0,0.0D0),X2,N)

      CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),RJT_COMPLEX,N,BT2,N, &
             (0.0D0,0.0D0),X3,N)

      CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),X3,N,RJ_COMPLEX,N, &
                (0.0D0,0.0D0),X4,N)

      DO I=1,N
       DO J=1,N
        A(I,J)=(AS_COMPLEX(I,J)+X2(I,J))
        B(I,J)=(BS_COMPLEX(I,J)+X4(I,J))
        ENDDO
      ENDDO
      
      !INVERSE OF B 

      DO I=1,N
       DO J=1,N
       B_INV(I,J)=B(I,J)
       ENDDO
      ENDDO

      CALL ZGETRF(N,N,B_INV,N,IPIV,INFO)
      CALL ZGETRI(N,B_INV,N,IPIV,WORK,N,INFO)

!FORMATION OF THE DENOMINATOR MATRIX WITHIN THE SQUARE ROOT OF THE
!DETERMINANT WHICH IS DENOTED HERE AS AA

      CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),A,N,B_INV,N, &
              (0.0D0,0.0D0),AB_INV,N)
      CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),AB_INV,N,A,N, &
              (0.0D0,0.0D0),E1,N)
      CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),B,N,E1,N, &
              (0.0D0,0.0D0),E2,N)
      CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),B,N,B,N, &
              (0.0D0,0.0D0),B2,N)
      DO I=1,N
       DO J=1,N
       AA(I,J)=(B2(I,J)-E2(I,J))
       ENDDO
      ENDDO

      !INVERSE OF AA MATRIX (INVERSE OF THE DENOMINATOR WITHIN THE SQUARE
      !ROOT OF THE DETERMINANT)

      DO I=1,N
       DO J=1,N
       AA_INV(I,J)=AA(I,J)
       ENDDO
      ENDDO
      CALL ZGETRF(N,N,AA_INV,N,IPIV,INFO)
      CALL ZGETRI(N,AA_INV,N,IPIV,WORK,N,INFO)

      N1=2*N

      !FORMATION OF K MATRIX

      AK(1:N,1:N)=B
      AK(1:N,N+1:N1)=-A
      AK(N+1:N1,1:N)=-A
      AK(N+1:N1,N+1:N1)=B

      !INVERSE OF AK 
      !STORED THE AK MATRIX AS AK_INV FOR INVERSION CALCULATION USING BLAS/LAPACK

      DO I=1,N1
       DO J=1,N1
        AK_INV(I,J)=AK(I,J)
       ENDDO
      ENDDO
      CALL ZGETRF(N1,N1,AK_INV,N1,IPIV,INFO)
      CALL ZGETRI(N1,AK_INV,N1,IPIV,WORK,N1,INFO)


!XX=NUMERATOR MATRIX WITHIN THE SQUARE ROOT OF THE DETERMINANT

      CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),AS_COMPLEX,N,AT2,N, &
                (0.0D0,0.0D0),X5,N)
   

      !FORMATION OF THE COMPLETE DETERMINANT DENOTED BY BB

      CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),X5,N,AA_INV,N, &
              (0.0D0,0.0D0),BB,N)

!DETERMINANAT OF BB WHICH IS N BY N MATRIX
     
      CALL DETERMINANT(N,BB,DET)

!AMPLITUDE AND PHASE OF THE DETERMINANT 

      AMP=ABS(DET)
      AMP=SQRT(AMP)
      DET_I=AIMAG(DET)     !IMAGINARY PART OF DET
      DET_R=REAL(DET)      !REAL PART OF DET
      THETA=DET_I/DET_R    !THETA=PHASE FACTOR
      THETA=THETA/2.0D0
      

!........................................................................................................................................................

!........................................................................................................................................................

!SIMPLIFICATION OF THE EXPONENTIAL PART 

      !FORMATION OF E MATRIX OF N BY N 

      DO I=1,N
       DO J=1,N
       E(I,J)=(BT2(I,J)-AT2(I,J))
       ENDDO
      ENDDO

      !FORMATION OF F VECTOR AND ITS TRANSPOSE FT CONTAINING 2N ELEMENTS

      CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),RJT_COMPLEX,N,E,N, &
                (0.0D0,0.0D0),Y1,N)
      CALL ZGEMM('N','N',N,1,N,(1.0D0,0.0D0),Y1,N,D_COMPLEX,N, &
                (0.0D0,0.0D0),V1,N)

       
      F(1:N,1:N)=V1
      F(N+1:N1,1:N)=V1 

      DO I=1,N1
       FT(1,I)=F(I,1)
      ENDDO
     
!MATRIX MULTIPLICATION FOR THE EXPONENTIAL PART
      
      CALL ZGEMM('N1','N1',1,N1,N1,(1.0D0,0.0D0),FT,1,AK_INV,N1, &
              (0.0D0,0.0D0),V2,1)
      CALL ZGEMM('N1','N1',1,1,N1,(1.0D0,0.0D0),V2,1,F,N1, &
              (0.0D0,0.0D0),NUM1,1)

      
      CALL ZGEMM('N','N',1,N,N,(1.0D0,0.0D0),DT_COMPLEX,1,E,N, &
                (0.0D0,0.0D0),V3,1)
      CALL ZGEMM('N','N',1,1,N,(1.0D0,0.0D0),V3,1,D_COMPLEX,N, &
            (0.0D0,0.0D0),NUM2,1)

      NUM1R=REAL(NUM1(1,1))/2.0D0
      NUM1I=AIMAG(NUM1(1,1))/2.0D0
      NUM2R=REAL(NUM2(1,1))
      NUM2I=AIMAG(NUM2(1,1))
      NUMI=-NUM1R+NUM2R+THETA   !-T*DELE1
      NUMR=NUM1I-NUM2I

      !TOTAL PART WITHIN THE EXPONENTIAL 

      NUM=COMPLEX(NUMR,NUMI)


!FINAL EXPRESSION OF REVERSE INTERSYSTEM CROSSING RATE FOR THE DSO PART

       T1=ABS(T)          !ABSLOUTE VALUE OF TIME T
       RISC_DSO=AMP*EXP(NUM)*EXP(-T1*ETA)
       WRITE(34,*)T,REAL(RISC_DSO1,AIMAG(RISC_DSO)

       !INTEGRATION USING SIMPSON'S ONE-THIRD RULE WITHIN TIME DOMAIN

!      IF((K1.EQ.1).OR.(K1.EQ.(M1+1)))THEN
!          K_RISC_DSO1=K_RISC_DSO1+RISC_DSO1
!          ELSE
!              IF(MOD(K1,2).EQ.0.0D0)THEN
!              K_RISC_DSO1=K_RISC_DSO1+(RISC_DSO1*4.0D0)
!              ELSE
!              K_RISC_DSO1=K_RISC_DSO1+(2.0D0*RISC_DSO1)
!              ENDIF
!        ENDIF



!SPIN-VIBRONIC(SV) PART FOR THE RISC CALCULATION      

       SV=(0.0D0,0.0D0)           !SV= TEMPORARY SPIN-VIBRONIC PART FOR EACH TIME FOR EACH PROCESSOR ID
       SV2=(0.0D0,0.0D0)          !SV2 = TOTAL SV BY COLLECTING THE SV FROM ALL PROCESSOR
  
!FOR EACH TIME SV2 WOULD BE OBTAINED BY COLLECTING ALL SV2 FROM ALL 
!PROCESSORS  


       !FORMATION OF 2N BY 2N DIMENSIONAL G MATRIX AND THE VECTOR QUANTITY H OF 2N BY 1 DIMENSION

!PARALLELIZATION OF THE MOST TIME CONSUMING LOOP BY USING ALL PROCESSORS
!THROUGH  DIVIDING THE LOOP ITERATIONS AMONG ALL PROCESSSORS

       NP1=MYID+1                    !MYID = PROCESSOR ID, NP1=TEMPORARY PROCESSOR ID
       DO K=NP1,N,NP
       DO M=1,N
             DO I=1,N
               DO J=1,N
                IF(I.EQ.K)THEN
                G11(I,J)=-BS_COMPLEX(K,K)*X2(M,J)
                G12(I,J)=BS_COMPLEX(K,K)*X4(M,J)
                G21(I,J)=AS_COMPLEX(K,K)*X2(M,J)
                G22(I,J)=-AS_COMPLEX(K,K)*X4(M,J)
                ELSE
                G11(I,J)=(0.0D0,0.0D0)
                G12(I,J)=(0.0D0,0.0D0)
                G21(I,J)=(0.0D0,0.0D0)
                G22(I,J)=(0.0D0,0.0D0)
                ENDIF
               ENDDO
               IF(I.EQ.K)THEN
               H1(I,1)=BS_COMPLEX(K,K)*V1(M,1)
               H2(I,1)=-AS_COMPLEX(K,K)*V1(M,1)
               ELSE
               H1(I,1)=(0.0D0,0.0D0)
               H2(I,1)=(0.0D0,0.0D0)
               ENDIF
              ENDDO
         G(1:N,1:N)=G11
         G(1:N,N+1:N1)=G12
         G(N+1:N1,1:N)=G21
         G(N+1:N1,N+1:N1)=G22
         H(1:N,1:N)=H1
         H(N+1:N1,1:N)=H2

!TRANSPOSE OF THE H VECTOR

         DO I=1,N1
           H_TRANS(1,I)=H(I,1)
        ENDDO

!MATRIX MULTIPLICATION FOR THE SPIN-VIBRONIC PART

        !THIS IS FOR THE Tr(GK^(-1))+ [(K^(-1)F)^(T)G(K^(-1)F)]-H^(T)K^(-1)F


     CALL ZGEMM('N1','N1',N1,N1,N1,(1.0D0,0.0D0),G,N1,AK_INV,N1, &
              (0.0D0,0.0D0),GK_INV,N1)
      CALL ZGEMM('N1','N1',N1,1,N1,(1.0D0,0.0D0),AK_INV,N1,F,N1, &
              (0.0D0,0.0D0),AKINV_F,N1)


      DO I=1,N1
       AKINV_F_TRANS(1,I)=AKINV_F(I,1)   !TRANSPOSE OF AKINV_F i.e. K^(-1)F
      ENDDO


      CALL ZGEMM('N1','N1',1,N1,N1,(1.0D0,0.0D0),AKINV_F_TRANS,1,G,N1, &
              (0.0D0,0.0D0),AKINV_T_G,1)
      CALL ZGEMM('N1','N1',1,1,N1,(1.0D0,0.0D0),AKINV_T_G,1,AKINV_F, &
              N1,(0.0D0,0.0D0),NUM3,1)
      CALL ZGEMM('N1','N1',1,1,N1,(1.0D0,0.0D0),H_TRANS,1,AKINV_F,N1, &
              (0.0D0,0.0D0),NUM4,1)



      TRACE_GK_INV=(0.0D0,0.0D0)


      DO I=1,N1
        DO J=1,N1
        IF(I.EQ.J)THEN
        TRACE_GK_INV=TRACE_GK_INV+(GK_INV(I,J))
        ELSE
        TRACE_GK_INV=TRACE_GK_INV+(0.0D0,0.0D0)
        ENDIF
        ENDDO
      ENDDO

      TRACE_GK_INV_REAL=REAL(TRACE_GK_INV) 
      TRACE_GK_INV_IMAG=AIMAG(TRACE_GK_INV) 

      TRACE=COMPLEX(-TRACE_GK_INV_IMAG,TRACE_GK_INV_REAL)

      !TOTAL SPIN-VIBRONIC PART 

      SV=SV+((TRACE+NUM3(1,1)-NUM4(1,1))*TKK_COMPLEX(K,M))
    

      ENDDO   !END LOOP FOR K 
      ENDDO   !END LOOP FOR M


!      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
      
!CALLING ALL LOOP ITERATIONS RESULTS INTO THE MASTER ID I.E. NP=0 BY
!CALLING THE MPI_GATHER FUNCTION

       CALL MPI_GATHER(SV,1,MPI_DOUBLE_COMPLEX, SV1 ,1 , &
       MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD, IERROR)
     
!SV1= TEMPORARY SV FOR EACH PROCESSOR 

!REVERSE INTERSYSTEM CROSSING RATE FOR THE SPIN-VIBRONIC PART 

!CALCULATION OF TOTAL SV FOR EACH TIME

      IF(MYID==0)THEN
      DO I=1,NP
      SV2=SV2+SV1(I)
      ENDDO

      RISC_SV=SV2*AMP*EXP(NUM)*EXP(-T1*ETA1)
     
      WRITE(41,*)T,REAL(RISC_SV),AIMAG(RISC_SV)

!INTEGRATION USING SIMPSON'S ONE-THIRD RULE WITHIN TIME DOMAIN
!      IF((K1.EQ.1).OR.(K1.EQ.(M1+1)))THEN
!          K_RISC_SV1=K_RISC_SV1+RISC_SV1
!          ELSE
!              IF(MOD(K1,2).EQ.0.0D0)THEN
!              K_RISC_SV1=K_RISC_SV1+(RISC_SV1*4.0D0)
!              ELSE
!              K_RISC_SV1=K_RISC_SV1+(2.0D0*RISC_SV1)
!              ENDIF
!      ENDIF
!      ENDIF

        ENDDO                          !END LOOP FOR T
       K_RISC_DSO1_REAL=ABS(K_RISC_DSO1)
       K_RISC_SV1_REAL=ABS(K_RISC_SV1)
       K_RISC_DSO1_REAL=(K_RISC_DSO1_REAL*HH*CONST1)/3.0D0
       K_RISC_SV1_REAL=(K_RISC_SV1_REAL*HH*CONST2)/3.0D0
       K_RISC_DSO1_REAL=K_RISC_DSO1_REAL*AU_HZ
       K_RISC_SV1_REAL=K_RISC_SV1_REAL*AU_HZ

       IF(MYID==0)THEN
       WRITE(*,*)'K_RISC_DSO1_REAL IS',K_RISC_DSO1_REAL
       WRITE(*,*)'K_RISC_SV1_REAL IS',K_RISC_SV1_REAL
       ENDIF

!DESTRUCTION OF THE MPI 

       CALL MPI_FINALIZE(IERROR)
       END PROGRAM RISC_RATE_CONSTANT


!.................................................................................................................................................................\\
!SUBROUTINE FOR DETERMINANT CALCULATION USING LAPACK-BLAS LIBRARY

      SUBROUTINE DETERMINANT(N,A,DET)
      COMPLEX*16::A(N,N)
      COMPLEX*16::DET,SGN
      INTEGER::I,INFO,N
      INTEGER::IPIV(N)

      CALL ZGETRF(N, N, A, N, IPIV,INFO)
      DET =(1.0D0,0.0D0)
       DO I = 1,N
        DET = DET*A(I,I)
       ENDDO
      SGN =(1.0D0,0.0D0)
      DO I = 1, N
      IF(IPIV(I) /= I) THEN
        SGN = -SGN
      END IF
      ENDDO
      DET= SGN*DET
      RETURN
      END



