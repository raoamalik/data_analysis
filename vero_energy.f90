!April2018 Veronica

     PROGRAM ENERGY_TRANSFER
      IMPLICIT NONE

      REAL :: EREL, ERELSQ, ETOTB, ETOTB_F, EVTOT, C1, PSQ, BDC,
E_ZERO_B
      REAL :: DUM1, DUM2, DUM3, DIS, KINETIC, POTENTIAL, EVIB,
EN_TRANSFER
      REAL :: VCMA(3), QCMA(3), SWA, AMA(4), EROTA, ECMRVA, EVIBA
      REAL :: VCMB(3), QCMB(3), SWB, AMB(4), EROTB, ECMRVB, EVIBB
      REAL :: VSQA, VSQB, ECMTA, ECMTB, V, V_KCAL_B
      REAL :: RCM, VREL, VRELSQ, RDMASS, QR(3), VR(3)
      REAL, ALLOCATABLE :: P(:), Q(:), W(:), QA(:), PA(:), WA(:)
      REAL, ALLOCATABLE :: PB(:), QB(:), WB(:), PPB(:), QQB(:)      
      REAL, ALLOCATABLE :: PPA(:), QQA(:), INDA(:), INDB(:)
      REAL, ALLOCATABLE :: CMRVA(:), CMRVB(:)
      INTEGER, ALLOCATABLE :: SIM_R(:), TRAJ_R(:)
      INTEGER :: N, DOF, NP, I, CC, NA, NB, J, NATOM, IO=0, O=0, TJ,
SIM, NRT, M
      INTEGER :: J1, J2, J3, PATH, X, Y, STEP=0
      CHARACTER :: L1, L2, L3, ATOM
      CHARACTER, ALLOCATABLE :: ARRAY_ATOM(:)
      INTEGER:: COUNTA,COUNTB
      LOGICAL:: last, reactive
!

!     I read from 'coordinates' the number of reactive trajectories, the
!     number of atoms (N),
!     the number of A atoms (NA), the directory(SIM), and the
!     coordinates and momenta of
!     each step of each trajectory.
!     The classification of trajectories is done reading the file
!     'react_traj.dat' (generated!
!     by the script from react.dat). There will be an  energy transfer
!     output file for reactive trajectories and
!     another one for non-reactive (non-reactive + isomers)


      Open(unit=5, action='read', file= 'coordinates')
      Read(5,*) NRT
      !NRT is the number of reactive trajectories
      Read(5,*) E_ZERO_B
      Open(unit=3, action='read', file= 'react_traj.dat')
      Allocate (SIM_R(NRT))
     
      Allocate (TRAJ_R(NRT))
      Do i=1, NRT
       Read(3,*) L3, SIM_R(i), TRAJ_R(i)
      End do
     ! I red which are the reactive trajectories

      Read(5,*) SIM, N, NA

     !Part of the set up             
      M=N+5
      NB= N - NA
      DOF = 3 * N
      ALLOCATE (Q(DOF), P(DOF), W(N), INDA(N), INDB(N))
      ALLOCATE (PA(DOF), QA(DOF), PPA(DOF), QQA(DOF), WA(N))
      ALLOCATE (PB(DOF), QB(DOF), PPB(DOF), QQB(DOF), WB(N)) 
      ALLOCATE (CMRVA(N), CMRVB(N))
     
     !BDC = 2.99D0

      C1 = 0.04184D0
     !C1 is a parameter to convert venus units to kcal/mol
      
      Do i=1, NA
       INDA(I)=I
       INDB(I)=0
      End Do
      
      j=NA+1
      
      Do i=j, N
        INDB(i)=i
        INDA(i)=0
      End Do

!     I read in the file xyz the atom type and masses are assigned
      Open(unit=1, action='read', file= 'file_xyz')
      Read(1,*) NATOM
      If (NATOM .NE. N) THEN
       Write(*,*) 'NATOM in file_xyz is different from N in the script'
      ELSE
       continue
      end if 
   
      Read(1,*)
      Allocate (ARRAY_ATOM(NATOM))
      Do i=1, NATOM
         Read(1,*) ARRAY_ATOM(i)
         ATOM=ARRAY_ATOM(i)
         SELECT CASE (ATOM)                                             
          CASE ('H')
           W(I) = 1.00794D0
          CASE ('C')
           W(I) = 12.00000D0
          CASE ('O')
           W(I) = 15.9994D0
          CASE ('N')
          W(I) = 14.007D0
          CASE ('S')
          W(I) = 32.00000D0
          CASE DEFAULT
          WRITE(6,*)'THIS ATOM HAS NOT BEEN DEFINED.PLZ DEFINE IT'
         END SELECT
      End Do
      Close(1)

!Opening writing files before the cycle over the steps      
      Open(unit=8, file='energy_transfer_reactive.dat', action='write')
      Open(unit=9, file='short_energy_transfer_reactive.dat',
action='write')
      Open(unit=10, file='energy_transfer_non_reactive.dat',
action='write')
      Open(unit=11, file='short_energy_transfer_non_reactive.dat',
action='write')
!input_mopac contains cartesian coordinates and momenta

     
!**********************QUI INIZIA IL CICLO SUGLI STEPS DELLA
!TRAIETTORIA!!!********************* 
    1 FORMAT(45X,I3)
    4 FORMAT(2(A4,I3,5X),A15,X,F15.6)
    5 FORMAT(F15.6) 

      Read(5,1) TJ
    DO while(IO==0)
      O=TJ       
!Now I'll see if TJ is reactive or not
       reactive=.false.
       Do i=1, NRT
        If (TJ==TRAJ_R(i) .and. SIM==SIM_R(i)) then
         reactive=.true.
        end if
       End do

 
!Now i know TJ is reactive.
      
!     Read in the file coordinates
       Do I=1,5
        Read(5,*) 
       End Do
       STEP=STEP+1
       If (STEP==1) then
        Print*, 'Trajectory ', TJ
       end if
       Do I=1,N
        J3 = I * 3
        J2 = J3 - 1
        J1 = J2 - 1
        Read(5,*) Q(J1),Q(J2),Q(J3),P(J1),P(J2),P(J3)
       End Do
       last=.false.
       Read(5,1, IOSTAT=IO) TJ
       If (TJ.gt.O) then
          STEP=0
       End if
       If ((TJ.gt.O) .OR. (IO.ne.0)) then
          last=.true.
       endif
       
   If ((STEP==1) .OR. (last.eqv..true.))  then
       !voglio calcolare l'energia nel caso lo step sia il primo, la
       !traiettoria sia cambiata (quindi ho letto le coordinate
       !dell'ultimo step), oppure se il file sia finito (come prima).
          
!     Restore new Ps and Qs into A and B (Assegna Ps, Qs and Ws alle
!     molecole di appartenenza)
!
!     Set initial value to be zero (settandoli a 0 PA... and PB hanno
!     stesso rango (DOF))
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
       END DO

       DO I=1, NA     
         J3 = I * 3
         J2 = J3 - 1
         J1 = J2 - 1
            PA(J3) = P(J3)
            PA(J2) = P(J2)
            PA(J1) = P(J1)
            QA(J3) = Q(J3)
            QA(J2) = Q(J2)
            QA(J1) = Q(J1)
            WA(I) = W(I)
       End Do
      
       j = NA + 1
       Do i=j, N 
        J3 = I * 3
        J2 = J3 - 1
        J1 = J2 - 1    
            PB(J3) = P(J3)
            PB(J2) = P(J2)
            PB(J1) = P(J1)
            QB(J3) = Q(J3)
            QB(J2) = Q(J2)
            QB(J1) = Q(J1)
            WB(I) = W(I)
       END DO




!Assegna le masse alle molecole

       
       SWA = 0
       SWB = 0
       DO I=1, N
          SWA = SWA + WA(I)
          SWB = SWB + WB(I)
       ENDDO


       IF (NA .GE. 1.5) THEN
        CALL CENMAS(DOF, SWA, QCMA, VCMA, N, PA, QA, PPA, QQA, WA)
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

       IF (NB .GE. 1.5) THEN
         CALL CENMAS(DOF, SWB, QCMB, VCMB, N, PB, QB, PPB, QQB, WB)
         CALL EPOT(DOF, QQB, N, ARRAY_ATOM, V_KCAL_B)
         CALL ROTN(DOF, PPB, QQB, N, EROTB, AMB, WB)
       ELSE
         EROTB = 0
         DO I = 1, N
           IF (NB.LT.1.5d0.AND.INDB(I) .GE. 0.5) THEN
              J3 = INDB(I)*3
              J2 = J3 - 1
              J1 = J2 -1

              QCMB(1) = QB(J1)
              QCMB(2) = QB(J2)
              QCMB(3) = QB(J3)
!              WRITE(6,*) QCMB

              VCMB(1) = PB(J1)/SWB
              VCMB(2) = PB(J2)/SWB
              VCMB(3) = PB(J3)/SWB
           ENDIF
         ENDDO
       ENDIF

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

!      WRITE(6,*) 'THE ROTATIONAL ENERGY FOR MOLECULAR A IS'
!      WRITE(6,*) EROTA, 'KCAL/MOL'
!      WRITE(6,*) 'THE ROTATIONAL ENERGY FOR MOLECULAR B IS'
!      WRITE(6,*) EROTB, 'KCAL/MOL'
!
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
        ENDDO

        RCM = SQRT(RCM)
        VREL = 0.0D0
        VRELSQ = 0.0D0
        DO I=1, 3
         VREL = VREL + VR(I) * QR(I)
         VRELSQ = VRELSQ + VR(I)**2
        ENDDO
!      write(6,*) vrelsq, rdmass, rcm
        VREL = VREL / RCM
        EREL = RDMASS * VREL**2 / 2.0D0 / C1
        ERELSQ = RDMASS * VRELSQ / 2.0D0 / C1

!      WRITE(6,*) 'THE RELATIVE TRANSLATIONAL ENERGY IS'
!      WRITE(6,*) ERELSQ
!      WRITE(6,*) ERELSQ, 'KCAL/MOL'

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
        V_KCAL_B= V_KCAL_B - E_ZERO_B
        EVIBB = ECMRVB - EROTB
        EVIBB = EVIBB + V_KCAL_B
        ETOTB= EVIBB + EROTB + ECMTB


        if (STEP==1) then
           ETOTB_F=ETOTB
        else if (last.eqv..true.) then
           EN_TRANSFER=ETOTB_F-ETOTB
        end if

        if ((reactive.eqv..true.) .and.(last.eqv..true.)) then
           
          WRITE(8, 4)'SIM=',SIM, 'TRJ=',TJ, 'ENERGY TRANSFER',
EN_TRANSFER
          WRITE(9, 5) EN_TRANSFER

        elseif ((reactive.eqv..false.) .and.(last.eqv..true.)) then 

          WRITE(10,4) 'SIM=', SIM, 'TRJ=', TJ, 'ENERGY TRANSFER',
EN_TRANSFER
          WRITE(11,5) EN_TRANSFER
        end if

    end if
   End Do
!Questo End Do chiude il do while.

100 CONTINUE
    Close(8)
    Close(3)
    Close(10)
    Close(5)
Stop
END PROGRAM ENERGY_TRANSFER

!============================================================
!
!     Calculate the center of mass momenta and coordinates
!     Those are stored in PPA and QQA
!
      SUBROUTINE CENMAS(DOF, SW, QCM, VCM, N, P, Q, PP, QQ, W)
      INTEGER :: DOF
      REAL :: PP(DOF), QQ(DOF), P(DOF), Q(DOF)
      REAL :: QCM(3), VCM(3), W(N)

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
!      WRITE(6,*) 'hehehe'
      RETURN
!
      END 
!========================================================
!========================================================
!Calculating potential energy
!To use the subroutine you need to change I=1
SUBROUTINE EPOT(DOF, QQ, N, ARRAY_ATOM, V_KCAL_B)

INTEGER :: N, DOF, INDEX_B, JF1, JF2, JF3
CHARACTER :: ARRAY_ATOM(2), L4, L5, L6, L7, L8
REAL :: QQ(DOF), V_EV_B, V_KCAL_B, DIST_N2

Open (9, file='input_mopac', action='write')

 Write(9,*) 'RM1 1SCF UHF PRECISE CHARGE=0 GEO-OK'
 Write(9,*) 'single point'
 Write(9,*) 'calculation'





!I AM WRITING THE MOPAC INPUT IN Z-MATRIX COORDINATES
 I=N-1
 J3=I*3
 J2=J3-1
 J1=J2-1
 !Write (9, '(A2,60X )' ) ARRAY_ATOM(i), '0.000000 0   0.000000  0
 !0.000000  0   0   0   0'

 Write (9, *) ARRAY_ATOM(i), ' 0.000000 0   0.000000  0   0.000000  0
0   0   0'
 I=N
 JF3=I*3
 JF2=JF3-1
 JF1=JF2-1

 DIST_N2=(QQ(JF1)-QQ(J1))**2 + (QQ(JF2)-QQ(J2))**2 + (QQ(JF3)-QQ(J3))**2
 DIST_N2=SQRT(DIST_N2)

 Write (9, *) ARRAY_ATOM(i), DIST_N2,' 1   0.000000  0   0.000000  1   0
0   0'
! Write (9, '(A2, 3(F9.4,60X)') ARRAY_ATOM(i), DIST_N2,' 1   0.000000  0
! 0.000000  1   0   0   0'

!Call system('mopac5021mn.exe <input_mopac> output_mopac')

Call system('mopac5021mn.exe <input_mopac> output_mopac')

Call system('grep TOTAL --binary-files=text output_mopac | grep ENERGY
C--binary-files=text >temp')


Open(13, file='temp', action='read')
Read(13, *) L4, L5, L6, V_EV_B, L7, L8, V_KCAL_B


!Write(*,*) V_KCAL_B




Close (13)


Close (9)
END 
















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
!     Calculate inverse of the inertia tensor (IN PEZZI!!)
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
!=====================================================


























