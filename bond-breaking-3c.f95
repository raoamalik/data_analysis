! This program calculates the bond lengths between atoms and reports
! time when
! any bond reaches 6A.to get the coordinates in a file: grep -A104 "Q
! !P"
! output.dat | grep -v -- "^--$" |grep -vE "Q
! P"> file
! cut -c3-35 file > qcoord.dat

        implicit real*8(a-h, o-z)
        integer, parameter :: rows=16000,columns=104
        dimension x(rows,columns),y(rows,columns),z(rows,columns)
        dimension r12(rows),r13(rows),r26(rows),r67(rows),r73(rows)
        dimension r713(rows),r712(rows),r39(rows),r915(rows),r1517(rows)
        dimension r1719(rows),r1923(rows),r1720(rows),r3638(rows)
        dimension r4044(rows),r63103(rows),r4153(rows),r5557(rows)
        dimension r5761(rows),r91104(rows),r3435(rows),r5558(rows)
        dimension r6465(rows),r6567(rows),r6671(rows),r7172(rows)
        dimension r5863(rows), r8182(rows),r80103(rows),r7276(rows)
        dimension r7680(rows),r8286(rows),r7274(rows),r3840(rows)
        dimension r7481(rows),r60102(rows),r5355(rows),r6771(rows)
        dimension r9897(rows),r3841(rows),r2036(rows)
        dimension r164(rows),r61103(rows)
        
        open(7, file='qcoord.dat')
        open(17, file='distance')
       do i=1,rows
        do j=1,columns
       read(7,*) x(i,j),y(i,j),z(i,j)
        enddo
         r12(i)= sqrt ((x(i,1)-x(i,2))**2+(y(i,1)-y(i,2))**2+(z(i,1)&
        -z(i,2))**2)
        r13(i)= sqrt ((x(i,1)-x(i,3))**2+(y(i,1)-y(i,3))**2 +(z(i,1)&
        -z(i,3))**2)
        r26(i)= sqrt ((x(i,2)-x(i,6))**2+(y(i,2)-y(i,6))**2&
        +(z(i,2)-z(i,6))**2)
        r67(i)= sqrt ((x(i,6)-x(i,7))**2+(y(i,6)-y(i,7))**2&
        +(z(i,6)-z(i,7))**2)
        r73(i)= sqrt ((x(i,3)-x(i,7))**2+(y(i,3)-y(i,7))**2&
        +(z(i,3)-z(i,7))**2)
        r713(i)= sqrt ((x(i,13)-x(i,7))**2+(y(i,13)-y(i,7))**2&
        +(z(i,13)-z(i,7))**2)
        r712(i)= sqrt ((x(i,12)-x(i,7))**2+(y(i,12)-y(i,7))**2&
        +(z(i,12)-z(i,7))**2)
        r39(i)= sqrt ((x(i,3)-x(i,9))**2+(y(i,3)-y(i,9))**2&
        +(z(i,3)-z(i,9))**2)
        r164(i)= sqrt ((x(i,1)-x(i,64))**2+(y(i,1)-y(i,64))**2+(z(i,1)&
        -z(i,64))**2)
        r915(i)= sqrt ((x(i,9)-x(i,15))**2+(y(i,9)-y(i,15))**2&
        +(z(i,9)-z(i,15))**2)
        r1517(i)= sqrt ((x(i,15)-x(i,17))**2+(y(i,15)-y(i,17))**2&
        +(z(i,15)-z(i,17))**2)
        r1719(i)= sqrt ((x(i,19)-x(i,17))**2+(y(i,19)-y(i,17))**2&
        +(z(i,19)-z(i,17))**2)
        r1923(i)= sqrt ((x(i,19)-x(i,23))**2+(y(i,19)-y(i,23))**2&
        +(z(i,19)-z(i,23))**2)
        r1720(i)= sqrt ((x(i,17)-x(i,20))**2+(y(i,17)-y(i,20))**2&
        +(z(i,17)-z(i,20))**2)
        r2036(i)= sqrt ((x(i,36)-x(i,20))**2+(y(i,36)-y(i,20))**2&
        +(z(i,36)-z(i,20))**2)
        r3638(i)= sqrt ((x(i,36)-x(i,38))**2+(y(i,36)-y(i,38))**2&
        +(z(i,36)-z(i,38))**2)
        r3840(i)= sqrt ((x(i,38)-x(i,40))**2+(y(i,38)-y(i,40))**2&
        +(z(i,38)-z(i,40))**2)
        r4044(i)= sqrt ((x(i,40)-x(i,44))**2+(y(i,40)-y(i,44))**2&
        +(z(i,40)-z(i,44))**2)
       r63103(i)= sqrt ((x(i,63)-x(i,103))**2+(y(i,63)-y(i,103))**2&
       +(z(i,63)-z(i,103))**2)
        r61103(i)= sqrt ((x(i,61)-x(i,103))**2+(y(i,61)-y(i,103))**2&
        +(z(i,61)-z(i,103))**2)
        r4153(i)= sqrt ((x(i,41)-x(i,53))**2+(y(i,41)-y(i,53))**2&
        +(z(i,41)-z(i,53))**2)
        r5557(i)=sqrt ((x(i,55)-x(i,57))**2+(y(i,55)-y(i,57))**2&
         + (z(i,55)-z(i,57))**2)
        r5355(i)=sqrt ((x(i,55)-x(i,53))**2+(y(i,55)-y(i,53))**2&
         + (z(i,55)-z(i,53))**2)        
       r91104(i)= sqrt ((x(i,91)-x(i,104))**2+(y(i,91)-y(i,104))**2&
       +(z(i,91)-z(i,104))**2)
        r3435(i)= sqrt ((x(i,34)-x(i,35))**2+(y(i,34)-y(i,35))**2&
        +(z(i,34)-z(i,35))**2)
        r3841(i)= sqrt ((x(i,38)-x(i,41))**2+(y(i,38)-y(i,41))**2&
        +(z(i,38)-z(i,41))**2)
        r5558(i)= sqrt ((x(i,55)-x(i,58))**2+(y(i,55)-y(i,58))**2&
        +(z(i,55)-z(i,58))**2)
       r5761(i)= sqrt ((x(i,57)-x(i,61))**2+(y(i,57)-y(i,61))**2&
         + (z(i,57)-z(i,61))**2)
       r5863(i)= sqrt ((x(i,58)-x(i,63))**2+(y(i,58)-y(i,63))**2&
         + (z(i,58)-z(i,63))**2)
       r6465(i)= sqrt ((x(i,64)-x(i,65))**2+(y(i,64)-y(i,65))**2&
         + (z(i,64)-z(i,65))**2)
        r80103(i)= sqrt ((x(i,80)-x(i,103))**2+(y(i,80)-y(i,103))**2&
        +(z(i,80)-z(i,103))**2)
       r6567(i)= sqrt ((x(i,65)-x(i,67))**2+(y(i,65)-y(i,67))**2&
        +(z(i,65)-z(i,67))**2)
        r6671(i)= sqrt ((x(i,66)-x(i,71))**2+(y(i,66)-y(i,71))**2&
        +(z(i,66)-z(i,71))**2)
        r6771(i)= sqrt ((x(i,67)-x(i,71))**2+(y(i,67)-y(i,71))**2&
        +(z(i,67)-z(i,71))**2)
        r7172(i)= sqrt ((x(i,72)-x(i,71))**2+(y(i,72)-y(i,71))**2&
        +(z(i,72)-z(i,71))**2)

        r7274(i)= sqrt ((x(i,72)-x(i,74))**2+(y(i,72)-y(i,74))**2&
        +(z(i,72)-z(i,74))**2)
        r7276(i)= sqrt ((x(i,72)-x(i,76))**2+(y(i,72)-y(i,76))**2&
        +(z(i,72)-z(i,76))**2)
       r7680(i)= sqrt ((x(i,76)-x(i,80))**2+(y(i,76)-y(i,80))**2&
        +(z(i,76)-z(i,80))**2)
       r8286(i)= sqrt ((x(i,82)-x(i,86))**2+(y(i,82)-y(i,86))**2&
        +(z(i,82)-z(i,86))**2)
       r8182(i)= sqrt((x(i,81)-x(i,82))**2+(y(i,81)-y(i,82))**2&
        +(z(i,81)-z(i,82))**2)
        r7481(i)= sqrt ((x(i,81)-x(i,74))**2+(y(i,81)-y(i,74))**2&
         +(z(i,81)-z(i,74))**2)

        r9897(i)= sqrt ((x(i,98)-x(i,97))**2+(y(i,98)-y(i,97))**2&
         + (z(i,98)-z(i,97))**2)



      if (r12(i)  .ge. 6.0) then
        write (17,*) 'r12, ', r12(i),  i*0.130E-14
        endif
      if (r13(i) .ge. 6.0) then
        write(17,*)'r13, ',r13(i),i*0.130E-14
        endif
      if (r26(i) .ge. 6.0) then
        write(17,*) 'r26, ', r26(i),  i*0.130E-14
        endif
      if (r67(i)  .ge. 6.0) then
        write (17,*) 'r67, ', r67(i),  i*0.130E-14
        endif
      if (r73(i) .ge. 6.0) then
        write(17,*)'r73, ',r73(i),i*0.130E-14
        endif
      if (r713(i) .ge. 6.0) then
        write(17,*) 'r713, ', r713(i),  i*0.130E-14
        endif
      if (r712(i) .ge. 6.0) then
        write(17,*) 'r712, ', r712(i),  i*0.130E-14
        endif
      if (r39(i) .ge. 6.0) then
        write(17,*)'r39, ',r39(i),i*0.130E-14
        endif
        if (r164(i)  .ge. 6.0) then
        write (17,*) 'r164, ', r164(i),  i*0.130E-14
        endif
      if (r915(i) .ge. 6.0) then
        write(17,*)'r915, ',r915(i),i*0.130E-14
        endif
      if (r1517(i) .ge. 6.0) then
        write(17,*)'r1517, ',r1517(i),i*0.130E-14
        endif
      if (r1719(i) .ge. 6.0) then
        write(17,*)'r1719, ',r1719(i),i*0.130E-14
        endif
      if (r1923(i) .ge. 6.0) then
        write(17,*)'r1923, ',r1923(i),i*0.130E-14
        endif
      if (r1720(i) .ge. 6.0) then
        write(17,*)'r1720, ',r1720(i),i*0.130E-14
        endif
      if (r2036(i) .ge. 6.0) then
        write(17,*)'r2036, ',r2036(i),i*0.250E-14
        endif
      if (r3638(i) .ge. 6.0)then
        write(17,*)'r3638, ',r3638(i),i*0.130E-14
        endif
      if (r3840(i) .ge. 6.0)then
        write(17,*)'r3840, ',r3840(i),i*0.130E-14
        endif
      if (r4044(i) .ge. 6.0)then
        write(17,*)'r4044, ',r4044(i),i*0.130E-14
        endif
      if (r5557(i) .ge. 6.0) then
        write(17,*)'r5557, ',r5557(i),i*0.130E-14
      endif
      if (r63103(i) .ge. 6.0) then
        write(17,*)'r63103, ',r63103(i),i*0.130E-14
      endif
      if (r61103(i) .ge. 6.0) then
        write(17,*)'r61103, ',r61103(i),i*0.130E-14
      endif
      if (r91104(i) .ge. 6.0) then
        write(17,*)'r91104, ',r91104(i),i*0.130E-14
      endif
      if (r3435(i) .ge. 6.0) then
        write(17,*)'r3435, ',r3435(i),i*0.130E-14
        endif
      if (r3841(i) .ge. 6.0) then
        write(17,*)'r3841, ',r3841(i),i*0.130E-14
        endif
      if (r4153(i) .ge. 6.0) then
        write(17,*)'r4153, ',r4153(i),i*0.130E-14
        endif

      if (r5558(i) .ge. 6.0) then
        write(17,*)'r5558, ',r5558(i),i*0.130E-14
        endif
      if (r5355(i) .ge. 6.0) then
        write(17,*)'r5355, ',r5355(i),i*0.130E-14
        endif

      if (r5761(i) .ge. 6.0) then
        write(17,*)'r5761, ',r5761(i),i*0.130E-14
        endif
     if (r5863(i) .ge. 6.0) then
        write(17,*)'r5863, ',r5863(i),i*0.130E-14
        endif
      if (r6465(i) .ge. 6.0) then
        write(17,*)'r6465, ',r6465(i),i*0.130E-14
        endif
      if (r6567(i) .ge. 6.0) then
        write(17,*)'r6567, ',r6567(i),i*0.130E-14
        endif
      if (r6671(i) .ge. 6.0) then
        write(17,*)'r6671, ',r6671(i),i*0.130E-14
        endif
      if (r6771(i) .ge. 6.0) then
        write(17,*)'r6771, ',r6771(i),i*0.130E-14
        endif
      if (r80103(i) .ge. 6.0) then
        write(17,*)'r80103, ',r80103(i),i*0.130E-14
        endif
      if (r7172(i) .ge. 6.0) then
        write(17,*)'r7172, ',r7172(i),i*0.130E-14
        endif
      if (r7274(i) .ge. 6.0) then
        write(17,*)'r7274, ',r7274(i),i*0.130E-14
        endif
      if (r7276(i) .ge. 6.0) then
        write(17,*)'r7276, ',r7276(i),i*0.130E-14
        endif
      if (r7680(i) .ge. 6.0) then
        write(17,*)'r7680, ',r7680(i),i*0.130E-14
        endif
      if (r8182(i) .ge. 6.0) then
        write(17,*)'r8182, ',r8182(i),i*0.130E-14
        endif
        if (r8286(i) .ge. 6.0) then
        write(17,*)'r8286, ',r8286(i),i*0.130E-14
        endif
      if (r7481(i) .ge. 6.0) then
        write(17,*)'r7481, ',r7481(i),i*0.130E-14
        endif

      if (r9897(i) .ge. 6.0) then
        write(17,*)'r9897, ',r9897(i),i*0.130E-14
        endif  
        enddo
        stop
        end
