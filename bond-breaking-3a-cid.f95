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
        dimension r4044(rows),r49102(rows),r4152(rows),r5456(rows)
        dimension r90102(rows),r93104(rows),r3435(rows),r5457(rows)
        dimension r6466(rows),r6364(rows),r6670(rows),r7175(rows)
        dimension r5762(rows), r62103(rows),r8185(rows),r79102(rows)
        dimension r7380(rows),r8081(rows),r7071(rows),r3840(rows)
        dimension r8183(rows),r8396(rows),r60102(rows),r5254(rows)
        dimension r9697(rows),r7579(rows),r3841(rows),r2036(rows)
        dimension r163(rows),r8589(rows),r5660(rows)
!Proton Transfer
        dimension r38(rows),r5253(rows),r5455(rows),r1516(rows)
        dimension r8588(rows), r3637(rows),r8587(rows),r6669(rows)
        dimension r8386(rows),r7376(rows)

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
        r163(i)= sqrt ((x(i,1)-x(i,63))**2+(y(i,1)-y(i,63))**2+(z(i,1)&
        -z(i,63))**2)
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
        r49102(i)= sqrt ((x(i,49)-x(i,102))**2+(y(i,49)-y(i,102))**2&
        +(z(i,49)-z(i,102))**2)
        r4152(i)= sqrt ((x(i,41)-x(i,52))**2+(y(i,41)-y(i,52))**2&
        +(z(i,41)-z(i,52))**2)
        r5456(i)= sqrt ((x(i,54)-x(i,56))**2+(y(i,54)-y(i,56))**2&
         + (z(i,54)-z(i,56))**2)
        r5254(i)= sqrt ((x(i,54)-x(i,52))**2+(y(i,54)-y(i,52))**2&
         + (z(i,54)-z(i,52))**2)
        r60102(i)= sqrt ((x(i,60)-x(i,102))**2+(y(i,60)-y(i,102))**2&
        +(z(i,60)-z(i,102))**2)
        r90102(i)= sqrt ((x(i,90)-x(i,102))**2+(y(i,90)-y(i,102))**2&
        +(z(i,90)-z(i,102))**2)
        r93104(i)= sqrt ((x(i,93)-x(i,104))**2+(y(i,93)-y(i,104))**2&
        +(z(i,93)-z(i,104))**2)
        r3435(i)= sqrt ((x(i,34)-x(i,35))**2+(y(i,34)-y(i,35))**2&
        +(z(i,34)-z(i,35))**2)
        r3841(i)= sqrt ((x(i,38)-x(i,41))**2+(y(i,38)-y(i,41))**2&
        +(z(i,38)-z(i,41))**2)
        r5457(i)= sqrt ((x(i,54)-x(i,57))**2+(y(i,54)-y(i,57))**2&
        +(z(i,54)-z(i,57))**2)
       r5762(i)= sqrt ((x(i,57)-x(i,62))**2+(y(i,57)-y(i,62))**2&
         + (z(i,57)-z(i,62))**2)
       r5660(i)= sqrt ((x(i,56)-x(i,60))**2+(y(i,56)-y(i,60))**2&
         + (z(i,56)-z(i,60))**2)
       r6466(i)= sqrt ((x(i,64)-x(i,66))**2+(y(i,64)-y(i,66))**2&
         + (z(i,64)-z(i,66))**2)
        r79102(i)= sqrt ((x(i,79)-x(i,102))**2+(y(i,79)-y(i,102))**2&
        +(z(i,79)-z(i,102))**2)
       r6364(i)= sqrt ((x(i,63)-x(i,64))**2+(y(i,63)-y(i,64))**2&
        +(z(i,63)-z(i,64))**2)
        r6670(i)= sqrt ((x(i,66)-x(i,70))**2+(y(i,66)-y(i,70))**2&
        +(z(i,66)-z(i,70))**2)
        r7071(i)= sqrt ((x(i,70)-x(i,71))**2+(y(i,70)-y(i,71))**2&
        +(z(i,70)-z(i,71))**2)
        r7579(i)= sqrt ((x(i,75)-x(i,79))**2+(y(i,75)-y(i,79))**2&
        +(z(i,75)-z(i,79))**2)
        r7175(i)= sqrt ((x(i,75)-x(i,71))**2+(y(i,75)-y(i,71))**2&
        +(z(i,75)-z(i,71))**2)
        r79102(i)= sqrt ((x(i,79)-x(i,102))**2+(y(i,79)-y(i,102))**2&
        +(z(i,70)-z(i,102))**2)
        r62103(i)= sqrt ((x(i,62)-x(i,103))**2+(y(i,62)-y(i,103))**2&
        +(z(i,62)-z(i,103))**2)
       r7380(i)= sqrt ((x(i,73)-x(i,80))**2+(y(i,73)-y(i,80))**2&
        +(z(i,73)-z(i,80))**2)
       r8081(i)= sqrt ((x(i,80)-x(i,81))**2+(y(i,80)-y(i,81))**2&
        +(z(i,80)-z(i,81))**2)
       r8185(i)= sqrt((x(i,81)-x(i,85))**2+(y(i,81)-y(i,85))**2&
        +(z(i,81)-z(i,85))**2)
       r8587(i)= sqrt((x(i,85)-x(i,87))**2+(y(i,85)-y(i,87))**2&
        +(z(i,85)-z(i,87))**2)

        r8183(i)= sqrt ((x(i,81)-x(i,83))**2+(y(i,81)-y(i,83))**2&
         +(z(i,81)-z(i,83))**2)
        r8396(i)= sqrt ((x(i,83)-x(i,96))**2+(y(i,83)-y(i,96))**2&
         +(z(i,83)-z(i,96))**2)
        r9697(i)= sqrt ((x(i,96)-x(i,97))**2+(y(i,96)-y(i,97))**2&
         + (z(i,96)-z(i,97))**2)
       r7380(i)= sqrt ((x(i,80)-x(i,73))**2+(y(i,80)-y(i,73))**2&
        +(z(i,80)-z(i,73))**2)
        r8589(i)= sqrt ((x(i,85)-x(i,89))**2+(y(i,85)-y(i,89))**2&
        +(z(i,85)-z(i,89))**2)

!Proton Transfer

        r38(i)= sqrt ((x(i,3)-x(i,8))**2+(y(i,3)-y(i,8))**2 +(z(i,8)&
        -z(i,3))**2)

        r5253(i)= sqrt ((x(i,52)-x(i,53))**2+(y(i,52)-y(i,53))**2&
        +(z(i,52)-z(i,53))**2)

        r5455(i)= sqrt ((x(i,54)-x(i,55))**2+(y(i,54)-y(i,55))**2&
        +(z(i,54)-z(i,55))**2)

        r1516(i)= sqrt ((x(i,15)-x(i,16))**2+(y(i,15)-y(i,16))**2&
        +(z(i,15)-z(i,16))**2)

        r8588(i)= sqrt ((x(i,85)-x(i,88))**2+(y(i,85)-y(i,88))**2&
        +(z(i,34)-z(i,35))**2)

        r3637(i)= sqrt ((x(i,36)-x(i,37))**2+(y(i,36)-y(i,37))**2&
         + (z(i,36)-z(i,37))**2)

        r6669(i)= sqrt ((x(i,66)-x(i,69))**2+(y(i,66)-y(i,69))**2&
        +(z(i,66)-z(i,69))**2)

        r7376(i)= sqrt ((x(i,73)-x(i,76))**2+(y(i,73)-y(i,76))**2&
        +(z(i,73)-z(i,76))**2)

        r8386(i)= sqrt ((x(i,83)-x(i,86))**2+(y(i,83)-y(i,86))**2&
        +(z(i,83)-z(i,86))**2)




      if (r12(i)  .ge. 6.0) then
        write (17,*) 'r12, ', r12(i),  i*0.130E-14
        write(17,*) 'sidechain'
        endif
      if (r13(i) .ge. 6.0) then
        write(17,*)'r13, ',r13(i),i*0.130E-14
        endif
      if (r26(i) .ge. 6.0) then
        write(17,*) 'r26, ', r26(i),  i*0.130E-14
        write(17,*) 'sidechain'
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
        if (r163(i)  .ge. 6.0) then
        write (17,*) 'r163, ', r163(i),  i*0.130E-14
        endif
      if (r915(i) .ge. 6.0) then
        write(17,*)'r915, ',r915(i),i*0.130E-14
        endif
      if (r1517(i) .ge. 6.0) then
        write(17,*)'r1517, ',r1517(i),i*0.130E-14
        endif
      if (r1719(i) .ge. 6.0) then
        write(17,*)'r1719, ',r1719(i),i*0.130E-14
        write(17,*) 'sidechain'
        endif
      if (r1923(i) .ge. 6.0) then
        write(17,*)'r1923, ',r1923(i),i*0.130E-14
        write(17,*) 'sidechain'
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
        write(17,*) 'sidechain'
        endif
      if (r4044(i) .ge. 6.0)then
        write(17,*)'r4044, ',r4044(i),i*0.130E-14
        write(17,*) 'sidechain'
        endif
      if (r5456(i) .ge. 6.0) then
        write(17,*)'r5456, ',r5456(i),i*0.130E-14
        write(17,*) 'sidechain'
      endif
      if (r5660(i) .ge. 6.0) then
        write(17,*)'r5660, ',r5660(i),i*0.130E-14
        write(17,*) 'sidechain'
        endif
      if (r60102(i) .ge. 6.0) then
        write(17,*)'r60102, ',r60102(i),i*0.130E-14
      endif
      if (r93104(i) .ge. 6.0) then
        write(17,*)'r93104, ',r93104(i),i*0.130E-14
      endif
      if (r3435(i) .ge. 6.0) then
        write(17,*)'r3435, ',r3435(i),i*0.130E-14
        endif
      if (r3841(i) .ge. 6.0) then
        write(17,*)'r3841, ',r3841(i),i*0.130E-14
        endif
      if (r5457(i) .ge. 6.0) then
        write(17,*)'r5457, ',r5457(i),i*0.130E-14
        endif
      if (r5254(i) .ge. 6.0) then
        write(17,*)'r5254, ',r5254(i),i*0.130E-14
        endif
      if (r5762(i) .ge. 6.0) then
        write(17,*)'r5762, ',r5762(i),i*0.130E-14
        write(17,*) 'sidechain'
        endif
        if (r6466(i) .ge. 6.0) then
        write(17,*)'r6466, ',r6466(i),i*0.130E-14
        endif
      if (r6364(i) .ge. 6.0) then
        write(17,*)'r6364, ',r6364(i),i*0.130E-14
        endif
      if (r6670(i) .ge. 6.0) then
        write(17,*)'r6670, ',r6670(i),i*0.130E-14
        endif
      if (r79102(i) .ge. 6.0) then
        write(17,*)'r79102, ',r79102(i),i*0.130E-14
        endif
      if (r7071(i) .ge. 6.0) then
        write(17,*)'r7071, ',r7071(i),i*0.130E-14
        endif
      if (r7175(i) .ge. 6.0) then
        write(17,*)'r7175, ',r7175(i),i*0.130E-14
        write(17,*) 'sidechain'
        endif
      if (r7579(i) .ge. 6.0) then
        write(17,*)'r7579, ',r7579(i),i*0.130E-14
        write(17,*) 'sidechain'
        endif
      if (r62103(i) .ge. 6.0) then
        write(17,*)'r62103, ',r62103(i),i*0.130E-14
        endif
      if (r7380(i) .ge. 6.0) then
        write(17,*)'r7380, ',r7380(i),i*0.130E-14
        endif
      if (r8185(i) .ge. 6.0) then
        write(17,*)'r8185, ',r8185(i),i*0.130E-14
        write(17,*) 'sidechain'
        endif
        if (r8081(i) .ge. 6.0) then
        write(17,*)'r8081, ',r8081(i),i*0.130E-14
        endif
      if (r8183(i) .ge. 6.0) then
        write(17,*)'r8183, ',r8183(i),i*0.130E-14
        endif
      if (r8587(i) .ge. 6.0) then
        write(17,*)'r8587, ',r8587(i),i*0.130E-14
        endif

        if (r8589(i) .ge. 6.0) then
        write(17,*)'r8589, ',r8589(i),i*0.130E-14
        endif
      if (r8396(i) .ge. 6.0) then
        write(17,*)'r8396, ',r8396(i),i*0.130E-14
        endif
      if (r9697(i) .ge. 6.0) then
        write(17,*)'r9697, ',r9697(i),i*0.130E-14
        endif
!Proton Transfer
      if (r38(i)  .ge. 6.0) then
        write (17,*) 'r38, ', r38(i),  i*0.130E-14
        endif

      if (r5253(i) .ge. 6.0) then
        write(17,*)'r5253, ',r5253(i),i*0.130E-14
        endif

      if (r5455(i) .ge. 6.0) then
        write(17,*) 'r5455, ', r5455(i),  i*0.130E-14
        endif

      if (r1516(i) .ge. 6.0) then
        write(17,*) 'r1516, ', r1516(i),  i*0.130E-14
        endif
 
      if (r8588(i) .ge. 6.0) then
        write(17,*)'r8588, ',r8588(i),i*0.130E-14
        endif

      if (r3637(i) .ge. 6.0) then
        write(17,*)'r3637, ',r3637(i),i*0.130E-14
        endif
      if (r6669(i) .ge. 6.0)then
        write(17,*)'r6669, ',r6669(i),i*0.130E-14
        endif

      if (r7376(i) .ge. 6.0)then
        write(17,*)'r7376, ',r7376(i),i*0.130E-14
        endif

      if (r8386(i) .ge. 6.0) then
        write(17,*)'r8386, ',r8386(i),i*0.130E-14
        endif
  
        enddo
        stop
        end
