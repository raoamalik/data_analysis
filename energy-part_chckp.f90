!     energy partitioning for A(multi-atom-system), B(atom): malik
      PROGRAM ENERGY_PARTITION
      IMPLICIT NONE

      REAL :: EREL, ERELSQ, ETOT, VZERO, EVTOT, C1, PSQ, BDC
      REAL :: DUM1, DUM2, DUM3, DIS, KINETIC, POTENTIAL, EVIB
      REAL :: VCMA(3), QCMA(3), SWA, AMA(4), EROTA, ECMRVA, EVIBA
      REAL :: VCMB(3), QCMB(3), SWB, AMB(4), EROTB, ECMRVB, EVIBB
      REAL :: VSQA, VSQB, ECMTA, ECMTB
      REAL :: RCM, VREL, VRELSQ, RDMASS, QR(3), VR(3)
      REAL, ALLOCATABLE :: P(:), Q(:), W(:), QA(:), PA(:), WA(:)
      REAL, ALLOCATABLE :: PB(:), QB(:), WB(:), PPB(:), QQB(:)      
      REAL, ALLOCATABLE :: PPA(:), QQA(:), INDA(:), INDB(:)
      REAL, ALLOCATABLE :: CMRVA(:), CMRVB(:)
      INTEGER :: N, DOF, NP, I, CC, DOFA, DOFB, NA, NB
      INTEGER :: J1, J2, J3, PATH, X, Y
      CHARACTER :: ATOM*2, L1, L2
     
      INTEGER:: COUNTA,COUNTB
!
!     Set up the parameters, need to make read statement at last
!     NA is the number of atoms in molecule A, NB = atoms in projectile
!     BDC is the bond distance criteria used in the auto group
!
      BDC = 15.0D0
      C1 = 0.04184D0
      open(5, file='energy-part.dat') 
      READ(5,*) N, L1, L2, KINETIC, L1, L2, POTENTIAL
      READ(5,*) NA, L1, L2, ETOT
      DOF = 3 * N
      NB = N - NA
 
      ALLOCATE (Q(DOF), P(DOF), W(N), INDA(N), INDB(N))
!
!     Read in the coordinates
!
      READ(5,*) L1, L2
      DO I=1,N
       J3 = I * 3
       J2 = J3 - 1
       J1 = J2 - 1
       READ(5,*) ATOM,Q(J1),Q(J2),Q(J3),P(J1),P(J2),P(J3)
       SELECT CASE (ATOM)
        CASE ('H')
         W(I) = 1.00794D0
        CASE ('C')
         W(I) = 12.01100D0
        CASE ('O')
         W(I) = 15.999D0
        CASE ('N')
         W(I) = 14.007D0
        CASE ('S')
        W(I) = 32.06D0
        CASE ('Zn')
        W(I) = 65.38D0
        CASE ('Ar')
        W(I) = 39.948D0
        CASE DEFAULT
         WRITE(6,*)'THIS ATOM HAS NOT BEEN DEFINED.PLZ DEFINE IT'
       END SELECT
      ENDDO
!
!     Auto Group     
!
      DO I=1, N
         INDA(I) = 0
         INDB(I) = 0
      ENDDO

      INDA(1) = 1 
      INDB(1) = 0
      DO I=2, N
         J3 = I * 3
         J2 = J3 - 1
         J1 = J2 - 1
         DUM1 = Q(1) - Q(J1)
         DUM2 = Q(2) - Q(J2)
         DUM3 = Q(3) - Q(J3)
!-------------CHECK--------------
         WRITE(6,*)"Q1,Q2,Q3", Q(1), Q(2), Q(3)
         WRITE(6,*)"QJ1,QJ2,QJ3", Q(J1), Q(J2), Q(J3)
         WRITE(6,*)"J1,J2,J3", J1, J2, J3
!------------CHECK_END----------
         DIS = SQRT(DUM1**2 + DUM2**2 + DUM3**2)
         IF (DIS .LE. BDC) THEN
            INDA(I) = I
         ELSE
           INDB(I) = I
         ENDIF
      ENDDO

       DO I=1, N
          WRITE(6,*) "INDA, INDB", INDA(I), INDB(I)
       ENDDO

!---  check whether the auto group matches the input classification
!---  of molecule A and B
       COUNTA =0
       DO I = 1, N
          IF (INDA(I).GT.0.5) THEN
          COUNTA = COUNTA + 1
          ENDIF
       ENDDO
!       IF (COUNTA.NE.NA) THEN
!          WRITE(*,*) "Wrong group for molecule A, COUNTA :", COUNTA
!          GO TO 100
!       ENDIF

       COUNTB =0
       DO I = 1, N
          IF (INDB(I).GT.0.5) THEN
          COUNTB = COUNTB + 1
          ENDIF
       ENDDO
       IF (COUNTB.NE.NB) THEN
          WRITE(*,*) "Wrong group for molecule B, COUNTB :", COUNTB
          GO TO 100
       ENDIF




!
!     Restore new Ps and Qs into A and B
!
      ALLOCATE (PA(DOF), QA(DOF), PPA(DOF), QQA(DOF), WA(N))
      ALLOCATE (PB(DOF), QB(DOF), PPB(DOF), QQB(DOF), WB(N)) 
!     Set initial value to zero
      DO I=1, N
         J3 = I * 3
         J2 = J3 - 1
         J1 = J2 - 1
            PA(J3) = 0
            PA(J2) = 0
            PA(J1) = 0
            QA(J3) = 0
            QA(J2) = 0
            QA(J1) = 0
            WA(I) = 0
            PB(J3) = 0
            PB(J2) = 0
            PB(J1) = 0
            QB(J3) = 0
            QB(J2) = 0
            QB(J1) = 0
            WB(I) = 0
       ENDDO

      DO I=1, N     
         J3 = I * 3
         J2 = J3 - 1
         J1 = J2 - 1
         IF (INDA(I) .GE. 0.5) THEN
            PA(J3) = P(J3)
            PA(J2) = P(J2)
            PA(J1) = P(J1)
            QA(J3) = Q(J3)
            QA(J2) = Q(J2)
            QA(J1) = Q(J1)
            WA(I) = W(I)
        write (6,*) wa(i)
         ELSE
            PB(J3) = P(J3)
            PB(J2) = P(J2)
            PB(J1) = P(J1)
            QB(J3) = Q(J3)
            QB(J2) = Q(J2)
            QB(J1) = Q(J1)
            WB(I) = W(I)
         ENDIF
       ENDDO
       SWA = 0
       SWB = 0
       DO I=1, N
          SWA = SWA + WA(I)
          SWB = SWB + WB(I)
       ENDDO
        WRITE (6,*) "SWA", SWA, "SWB", SWB
        
!      ================ Change by Malik ==============
         EROTB = 0
         DO I = 1, N
           IF (NB.LT.1.5d0.AND.INDB(I) .GE. 0.5) THEN
              J3 = INDB(I)*3
              J2 = J3 - 1
              J1 = J2 -1

              QCMB(1) = QB(J1)
              QCMB(2) = QB(J2)
              QCMB(3) = QB(J3)
               WRITE(6,*) "QCMB_MAIN", QCMB ! CHECK

              VCMB(1) = PB(J1)/SWB
              VCMB(2) = PB(J2)/SWB
              VCMB(3) = PB(J3)/SWB
           ENDIF
         ENDDO
!      ===============================================


333   IF (NA .GE. 1.5) THEN
         CALL CENMAS(DOF, SWA, QCMA, QCMB, VCMA, N, PA, QA, PPA, QQA, WA)
         CALL ROTN(DOF, PPA, QQA, N, EROTA, AMA, WA)
      ELSE
         EROTA = 0
         DO I = 1, N
           IF (NA.LT.1.5d0.AND. INDA(I) .GE. 0.5) THEN
              J3 = INDA(I)*3
              J2 = J3 - 1
              J1 = J2 -1

              QCMA(1) = QA(J1)
              QCMA(2) = QA(J2)
              QCMA(3) = QA(J3)

              VCMA(1) = PA(J1)/SWA
              VCMA(2) = PA(J2)/SWA
              VCMA(3) = PA(J3)/SWA
           ENDIF
         ENDDO
      ENDIF

!      IF (NB .GE. 1.5) THEN
!         CALL CENMAS(DOF, SWB, QCMB, VCMB, N, PB, QB, PPB, QQB, WB)
!         CALL ROTN(DOF, PPB, QQB, N, EROTB, AMB, WB)
!      ELSE
!         EROTB = 0
!         DO I = 1, N
!           IF (NB.LT.1.5d0.AND.INDB(I) .GE. 0.5) THEN
!              J3 = INDB(I)*3
!              J2 = J3 - 1
!              J1 = J2 -1
!
!              QCMB(1) = QB(J1)
!              QCMB(2) = QB(J2)
!              QCMB(3) = QB(J3)
!              WRITE(6,*) "QCMB_MAIN", QCMB ! CHECK
!
!              VCMB(1) = PB(J1)/SWB
!              VCMB(2) = PB(J2)/SWB
!              VCMB(3) = PB(J3)/SWB
!           ENDIF
!         ENDDO
!      ENDIF

      VSQA = 0.0D0
      VSQB = 0.0D0 
      ECMTA = 0.0D0
      ECMTB = 0.0D0

      DO I=1, 3
        VSQA = VSQA + VCMA(I)**2
        VSQB = VSQB + VCMB(I)**2
      ENDDO
        ECMTA = SWA*VSQA/2.0D0/C1
        ECMTB = SWB*VSQB/2.0D0/C1

!     Calculate the relative translational energy of the system
!     The center of mass information is obtained from the soubroutine
!      cenmas and stored in qcm and vcm
!
      RDMASS = SWA * SWB / (SWA + SWB)
      RCM = 0.0D0
      DO I=1, 3
         QR(I) = QCMA(I) - QCMB(I)
         VR(I) = VCMA(I) - VCMB(I)
         RCM = RCM + QR(I)**2
!check        write(*,*) "VR(I)",VR(I)
      ENDDO
      RCM = SQRT(RCM)
      VREL = 0.0D0
      VRELSQ = 0.0D0
      DO I=1, 3
         VREL = VREL + VR(I) * QR(I)
         VRELSQ = VRELSQ + VR(I)**2
      ENDDO
!        write(*,*) VRELSQ
!      write(6,*) vrelsq, rdmass, rcm
      VREL = VREL / RCM
      EREL = RDMASS * VREL**2 / 2.0D0 / C1
      ERELSQ = RDMASS * VRELSQ / 2.0D0 / C1

      ALLOCATE (CMRVA(N), CMRVB(N))
      PSQ = 0.0D0
      ECMRVA = 0.00D0
      ECMRVB = 0.00D0
      DO I=1, N
         J3 = I * 3
         J2 = J3 - 1
         J1 = J2 - 1
         IF (INDA(I) .GE. 0.5) THEN
            PSQ = PPA(J1)**2 + PPA(J2)**2 + PPA(J3)**2
            CMRVA(I) = PSQ/2.0D0/C1/W(I)
         ELSE
            CMRVA(I) = 0
         ENDIF
         ECMRVA = ECMRVA + CMRVA(I)
      ENDDO
      DO I=1, N
         J3 = I * 3
         J2 = J3 - 1
         J1 = J2 - 1
         IF (INDB(I) .GE. 0.5) THEN
            PSQ = PPB(J1)**2 + PPB(J2)**2 + PPB(J3)**2
            CMRVB(I) = PSQ/2.0D0/C1/W(I)
         ELSE
            CMRVB(I) = 0
         ENDIF
         ECMRVB = ECMRVB + CMRVB(I)
      ENDDO
!      WRITE(6,*) ECMRVA, ECMRVB
      EVIBA = ECMRVA - EROTA
      EVIBB = ECMRVB - EROTB
      WRITE(6,*) "    E_rel        Rot_A       Vib_k_A"
!      &           Rot_B       Vib_k_B"
      WRITE(6,"(6F12.4)") ERELSQ,EROTA,EVIBA ! ,EROTB,EVIBB

100   CONTINUE
      END PROGRAM ENERGY_PARTITION

!============================================================
!
!     Calculate the center of mass momenta and coordinates
!     Those are stored in PPA and QQA
!
      SUBROUTINE CENMAS(DOF, SW, QCM, QCMB, VCM, N, P, Q, PP, QQ, W)
      INTEGER :: DOF
      REAL :: PP(DOF), QQ(DOF), P(DOF), Q(DOF)
      REAL :: QCM(3), VCM(3), W(N), CM_DIS, QCMB(3)

      DO I=1,3
       VCM(I) = 0.0D0
       QCM(I) = 0.0D0
      ENDDO

      DO I=1,N
       J3 = I * 3
       J2 = J3 - 1
       J1 = J2 - 1
       VCM(1) = VCM(1) + P(J1)
       VCM(2) = VCM(2) + P(J2)
       VCM(3) = VCM(3) + P(J3)
       QCM(1) = QCM(1) + W(I) * Q(J1)
       QCM(2) = QCM(2) + W(I) * Q(J2)
       QCM(3) = QCM(3) + W(I) * Q(J3)
      ENDDO
      
      DO I=1,3
       VCM(I) = VCM(I) / SW
       QCM(I) = QCM(I) / SW
      ENDDO
!      WRITE(6,*) "QCM, ", QCM ! CHECK
      DO I=1,N
       J3 = I * 3
       J2 = J3 - 1
       J1 = J2 - 1
       PP(J1) = P(J1) - W(I) * VCM(1)
       PP(J2) = P(J2) - W(I) * VCM(2)
       PP(J3) = P(J3) - W(I) * VCM(3)
       IF (W(I) .GE. 0.5) THEN 
          QQ(J1) = Q(J1) - QCM(1)
          QQ(J2) = Q(J2) - QCM(2)
          QQ(J3) = Q(J3) - QCM(3)
       ELSE
          QQ(J1) = 0
          QQ(J2) = 0
          QQ(J3) = 0
       ENDIF
      ENDDO

!        WRITE(6,*) "QQA", QQ(1), QQ(2), QQ(3)    ! CHECK MALIK
!     ================================================
!               CHANGES BY MALIK 08/22/2019
!     ================================================
!       CALCULATE COM DISTANCE BETWEEN THE SYSTEM
!                    AND PROJECTILE
!       
       CM_DIS = 0.0D0
       CM_DIS = SQRT( (QCM(1) - QCMB(1) )**2 + (QCM(2) - QCMB(2) )**2&
     & +  ( QCM(3) - QCMB(3) )**2)
!        WRITE(6,*) "QCM", QCM(1), QCM(2), QCM(3)
!        WRITE(6,*) "QQ", QQ  
!        WRITE(6,*) "QCMB", QCMB(1), QCMB(2), QCMB(3) 
       WRITE(6,*) "COM_DIS", CM_DIS

      RETURN
!
      END
!========================================================
!========================================================
!     Calculate angular momentum, moment of inertia 
!      tensor, angular velocity and rotational energy 
!
      SUBROUTINE ROTN(DOF, PP, QQ, N, EROT, AM, W)

      INTEGER :: DOF, N
      REAL :: AIXX, AIYY, AIZZ, AIXY, AIXZ, AIYZ
      REAL :: DET, UXX, UXY, UXZ, UYY, UYZ, UZZ
      REAL :: WX, WY, WZ
      REAL :: PP(DOF), QQ(DOF), AM(4), W(N)

      C1 = 0.04184D0
      AM(1) = 0.0D0
      AM(2) = 0.0D0
      AM(3) = 0.0D0
      AM(4) = 0.0D0
      EROT = 0.0D0

!
!     Calculate angular momentum. The center of mass 
!      coordinate qq and momenta pp come from subroutine
!      cenmas in the code
!
      IF (N .EQ. 1) THEN
         WRITE(6,*) 'THIS IS AN ATOM, NO ROT'
         RETURN
      ENDIF

      IF (N .EQ. 2) THEN
         WRITE(6,*) 'DIATOMIC, NOT SET UP IN THE CODE'
         RETURN
      ENDIF

      DO I=1, N
         J3 = I * 3
         J2 = J3 - 1
         J1 = J2 - 1
         AM(1) = AM(1) + (QQ(J2) * PP(J3) - QQ(J3) * PP(J2))
         AM(2) = AM(2) + (QQ(J3) * PP(J1) - QQ(J1) * PP(J3))
         AM(3) = AM(3) + (QQ(J1) * PP(J2) - QQ(J2) * PP(J1))
      ENDDO
      AM(4) = SQRT(AM(1)**2 + AM(2)**2 + AM(3)**2)
!
!     Calculate the moment of inertia 
!
      AIXX = 0.0D0
      AIYY = 0.0D0
      AIZZ = 0.0D0
      AIXY = 0.0D0
      AIXZ = 0.0D0
      AIYZ = 0.0D0
      DO I = 1, N
         J3 = 3 * I
         J2 = J3 - 1
         J1 = J2 - 1
!         write(6,*) w(i), qq(j1), qq(j2), qq(j3)
         AIXX = AIXX + W(I) * (QQ(J2)**2 + QQ(J3)**2)
         AIYY = AIYY + W(I) * (QQ(J1)**2 + QQ(J3)**2)
         AIZZ = AIZZ + W(I) * (QQ(J1)**2 + QQ(J2)**2)
         AIXY = AIXY + W(I) * QQ(J1) * QQ(J2)
         AIXZ = AIXZ + W(I) * QQ(J1) * QQ(J3)
         AIYZ = AIYZ + W(I) * QQ(J2) * QQ(J3)
      ENDDO
      DET=AIXX*(AIYY*AIZZ-AIYZ*AIYZ)-AIXY*(AIXY*AIZZ+AIYZ*AIXZ)-AIXZ &
          *(AIXY*AIYZ+AIYY*AIXZ)
!
!     Calculate inverse of the inertia tensor
!
      IF (ABS(DET).GE.0.01D0) THEN
         UXX=(AIYY*AIZZ-AIYZ*AIYZ)/DET
         UXY=(AIXY*AIZZ+AIXZ*AIYZ)/DET
         UXZ=(AIXY*AIYZ+AIXZ*AIYY)/DET
         UYY=(AIXX*AIZZ-AIXZ*AIXZ)/DET
         UYZ=(AIXX*AIYZ+AIXZ*AIXY)/DET
         UZZ=(AIXX*AIYY-AIXY*AIXY)/DET
!
!         CALCULATE ANGULAR VELOCITIES
!
         WX=UXX*AM(1)+UXY*AM(2)+UXZ*AM(3)
         WY=UXY*AM(1)+UYY*AM(2)+UYZ*AM(3)
         WZ=UXZ*AM(1)+UYZ*AM(2)+UZZ*AM(3)
      ELSE
!
!         CALCULATE ROTATIONAL ENERGY
!
         AIXX=0.0D0
         DO I=1,N
            J=3*I+1
            SR=0.0D0
            DO K=1,3
               SR=SR+QQ(J-K)**2
            ENDDO
            AIXX=AIXX+SR*W(I)
         ENDDO
         EROT=AM(4)**2/AIXX/2.0D0/C1
         RETURN
      ENDIF
      EROT=(WX*AM(1)+WY*AM(2)+WZ*AM(3))/2.0D0/C1
      RETURN
!
      END
!=====================================================

