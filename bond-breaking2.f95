! This program calculates the bond lengths between atoms and reports
! time when
! any bond reaches 6A.to get the coordinates in a file: grep -A104 "Q
! !P"
! output.dat | grep -v -- "^--$" |grep -vE "Q
! P"> file
! cut -c3-35 file > qcoord.dat

        implicit real*8(a-h, o-z)
        character s*2
        integer, parameter :: rows=5366,columns=104
        dimension x(rows,columns),y(rows,columns),z(rows,columns)
        dimension s(rows,columns)
        dimension r12(rows),r36(rows),r59(rows),r1213(rows),r9102(rows)
        dimension r11102(rows), r1317(rows),r1315(rows),r3034(rows)
        dimension r2930(rows),r4950(rows),r5263(rows),r6364(rows)
        dimension r6670(rows),r72102(rows),r7071(rows),r7579(rows)
        dimension r79102(rows)
        dimension r7380(rows),r8081(rows)
        dimension r8185(rows),r8183(rows),r8396(rows)
        dimension r9697(rows), r13(rows),r2455(rows)

!        write(6,*) "5=C-C backbone,  2=C-N backbone,
!     x  2=N-C backbone,
!     x  8=side chain C-C,
!       2=C-N(NH3+)"
        open(7, file='qcoord.dat')
        open(17, file='distance.dat')
!       read(5,*)s(i,j), x(i,j),y(i,j),z(i,j)
!        read(7, FMT="(3(F10.7,3x))") x(i,j),y(i,j),z(i,j)
        natom = 104
       do i=1,rows
!        read(5,*) nst
        do j=1,columns
!        open(7, file='qcoord.dat')
       read(7,*) x(i,j),y(i,j),z(i,j)
!        read(7, FMT="(3(F10.7,3x))") x(i,j),y(i,j),z(i,j)
        enddo
!       C-C bonds (backbone)
        r12(i)= sqrt ((x(i,1)-x(i,2))**2+(y(i,1)-y(i,2))**2 +(z(i,1)&
        -z(i,2))**2)
        r13(i)= sqrt ((x(i,1)-x(i,3))**2+(y(i,1)-y(i,3))**2 +(z(i,1)&
        -z(i,3))**2)
        r36(i)= sqrt ((x(i,3)-x(i,6))**2+(y(i,3)-y(i,6))**2&
        +(z(i,3)-z(i,6))**2)
        r59(i)= sqrt ((x(i,5)-x(i,9))**2+(y(i,5)-y(i,9))**2&
        +(z(i,5)-z(i,9))**2)
        r1213(i)= sqrt ((x(i,13)-x(i,12))**2+(y(i,13)-y(i,12))**2&
        +(z(i,13)-z(i,12))**2)
!       C-N    (backbone)
        r9102(i)= sqrt ((x(i,9)-x(i,102))**2+(y(i,9)-y(i,102))**2&
        +(z(i,9)-z(i,102))**2)
        r11102(i)= sqrt ((x(i,11)-x(i,102))**2+(y(i,11)-y(i,102))**2&
        +(z(i,11)-z(i,102))**2)
!       N-C     (backbone)
        r1317(i)= sqrt ((x(i,17)-x(i,13))**2+(y(i,17)-y(i,13))**2&
        +(z(i,17)-z(i,13))**2)
        r1315(i)= sqrt ((x(i,13)-x(i,15))**2+(y(i,13)-y(i,15))**2&
        +(z(i,13)-z(i,15))**2)
!       sidechain
        r3034=sqrt ((x(i,30)-x(i,34))**2+(y(i,30)-y(i,34))**2&
        +(z(i,30)-z(i,34))**2)
        r2930=sqrt ((x(i,29)-x(i,30))**2+(y(i,29)-y(i,30))**2&
        +(z(i,29)-z(i,30))**2)
        r4950=sqrt ((x(i,49)-x(i,50))**2+(y(i,49)-y(i,50))**2&
        +(z(i,49)-z(i,50))**2)
        r5263=sqrt ((x(i,52)-x(i,63))**2+(y(i,52)-y(i,63))**2&
        +(z(i,52)-z(i,63))**2)
        r6364(i)= sqrt ((x(i,63)-x(i,64))**2+(y(i,63)-y(i,64))**2&
        +(z(i,63)-z(i,64))**2)
        r6670(i)= sqrt ((x(i,66)-x(i,70))**2+(y(i,66)-y(i,70))**2&
        +(z(i,66)-z(i,70))**2)
        r72102(i)= sqrt ((x(i,72)-x(i,102))**2+(y(i,72)-y(i,102))**2&
        +(z(i,72)-z(i,102))**2)
        r7071(i)= sqrt ((x(i,70)-x(i,71))**2+(y(i,70)-y(i,71))**2&
        +(z(i,70)-z(i,71))**2)
        r7579(i)= sqrt ((x(i,75)-x(i,79))**2+(y(i,75)-y(i,79))**2&
        +(z(i,75)-z(i,79))**2)

!       C-N(NH3+)
        r79102(i)= sqrt ((x(i,79)-x(i,102))**2+(y(i,79)-y(i,102))**2&
        +(z(i,70)-z(i,102))**2)
        r7380(i)= sqrt ((x(i,73)-x(i,80))**2+(y(i,73)-y(i,80))**2&
        +(z(i,73)-z(i,80))**2)
!       Proton transfer
!       O-H bonds (C-OO-H)
!        r528(i)=sqrt ((x(i,5)-x(i,28))**2+(y(i,5)-y(i,28))**2&
!        +(z(i,5)-z(i,28))**2)
!        r2559(i)=sqrt ((x(i,25)-x(i,59))**2+(y(i,25)-y(i,59))**2&
!        +(z(i,25)-z(i,59))**2)
!       N-H bonds (NH3+)
        r8081(i)=sqrt ((x(i,80)-x(i,81))**2+(y(i,80)-y(i,81))**2&
        +(z(i,80)-z(i,81))**2)
        r8185(i)=sqrt ((x(i,81)-x(i,85))**2+(y(i,81)-y(i,85))**2&
         +(z(i,81)-z(i,85))**2)
        r8183(i)=sqrt ((x(i,81)-x(i,83))**2+(y(i,81)-y(i,83))**2&
         +(z(i,81)-z(i,83))**2)

        r8396(i)=sqrt ((x(i,83)-x(i,96))**2+(y(i,83)-y(i,96))**2&
         +(z(i,83)-z(i,96))**2)
        r9697(i)=sqrt ((x(i,96)-x(i,97))**2+(y(i,96)-y(i,97))**2&
         + (z(i,96)-z(i,97))**2)
!        write(6,91) r12(i), r23(i), r36(i), r89(i), r1112(i), r67(i),&
!        r910(i), r78(i), r1011(i), r815(i),r1516(i), r1518(i),&
!        r1617(i), r1120(i),r2021(i),r2122(i),r2223(i),r34(i),r2324(i),&
!        r528(i), r2559(i), r456(i), r457(i), r458(i), r2453(i),&
!        r2454(i) ,r2455(i)
!        write(6,*) "C-N, N-C backbone"
!        write(6,92) r67(i), r910(i), r78(i), r1011(i)
!        write (6,*) "side chain C-C"
!        write(6,93) r815(i),r1120(i), r1516(i),r1518(i), r1617(i),&
!        r2021(i), r2122(i), r2223(i)
!        write (6,*) "C-N(NH3+)"
!        write (6,*) r34(i),r2324(i)
91      format (27f9.3)
      if (r12(i) .ge. 6.0) then
        write (17,*) 'r12, ', r12(i),  i*0.130E-14
      elseif (r36(i).ge. 6.0) then
        write(17,*)'r36, ',r36(i),i*0.130E-14
        elseif (r59(i).ge. 6.0) then
        write(17,*) 'r59, ', r59(i),  i*0.130E-14
      elseif (r1213(i).ge. 6.0) then
        write(17,*)'r1213, ',r1213(i),i*0.130E-14
      elseif (r9102(i).ge. 6.0) then
        write(17,*)'r9102, ',r12(i),i*0.130E-14
      elseif (r11102(i).ge. 6.0) then
        write(17,*)'r11102, ',r11102(i),i*0.130E-14
      elseif (r1317(i).ge. 6.0)then 
        write(17,*)'r12, ',r1317(i),i*0.130E-14
      elseif (r1315(i).ge. 6.0)then 
        write(17,*)'r1315, ',r1315(i),i*0.130E-14
      elseif (r3034(i).ge. 6.0)then
        write(17,*)'r3034, ',r3034(i),i*0.130E-14
      elseif (r2930(i).ge. 6.0) then
        write(17,*)'r2930, ',r2930(i),i*0.130E-14
      elseif (r4950(i).ge. 6.0) then    
        write(17,*)'r4950, ',r4950(i),i*0.130E-14
      elseif (r5263(i).ge. 6.0) then
        write(17,*)'r5263, ',r5263(i),i*0.130E-14
      elseif (r6364(i).ge. 6.0) then
        write(17,*)'r6364, ',r6364(i),i*0.130E-14
      elseif (r6670(i).ge. 6.0) then
        write(17,*)'r6670, ',r6670(i),i*0.130E-14
      elseif (r72102(i).ge.6.0.AND.r72102(i).le.8) then
        write(17,*)'r72102, ',r72102(i),i*0.130E-14
      elseif (r7071(i).ge. 6.0) then
        write(17,*)'r7071, ',r7071(i),i*0.130E-14
      elseif (r7579(i).ge. 6.0) then
        write(17,*)'r7579, ',r7579(i),i*0.130E-14
      elseif (r79102(i).ge. 6.0) then
        write(17,*)'r79102, ',r79102(i),i*0.130E-14
      elseif (r7380(i).ge. 6.0) then
        write(17,*)'r7380, ',r7380(i),i*0.130E-14
      elseif (r8081(i).ge. 6.0) then
        write(17,*)'r8081, ',r8081(i),i*0.130E-14
      elseif (r8185(i).ge. 6.0) then
        write(17,*)'r8185, ',r8185(i),i*0.130E-14
      elseif (r8183(i).ge. 6.0) then
        write(17,*)'r8183, ',r8183(i),i*0.130E-14
      elseif (r8396(i).ge. 6.0) then
        write(17,*)'r8396, ',r8396(i),i*0.130E-14
      elseif (r9697(i).ge. 6.0) then
        write(17,*)'r9697, ',r9697(i),i*0.130E-14
      endif


!92      format (4f10.5)
!93      format (8f10.5)
!94      format (2f10.5)
!        endif
!        enddo
!
!        enddo
        enddo
!
        stop
        end
