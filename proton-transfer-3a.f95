! This program calculates the bond lengths between atoms and reports
! time when
! any bond reaches 6A.to get the coordinates in a file: grep -A104 "Q
! !P"
! output.dat | grep -v -- "^--$" |grep -vE "Q
! P"> file
! cut -c3-35 file > qcoord.dat

        implicit real*8(a-h, o-z)
        integer, parameter :: rows=15000,columns=104
        dimension x(rows,columns),y(rows,columns),z(rows,columns)
        dimension s(rows,columns)
        dimension r38(rows),r5455(rows),r57(rows),r2224(rows)
        dimension r3637(rows),r3435(rows),r5253(rows),r7376(rows)
        dimension r8386(rows),r62103(rows),r93104(rows),r1516(rows)
        dimension r6669(rows),r72102(rows),r7071(rows),r7579(rows)
        dimension r79102(rows)
        dimension r7380(rows),r8081(rows)
        dimension r8185(rows),r8183(rows),r8396(rows)
        dimension r9697(rows), r13(rows)

        open(7, file='qcoord.dat')
        open(17, file='distance')
       do i=1,rows
        do j=1,columns
       read(7,*) x(i,j),y(i,j),z(i,j)
        enddo

        r38(i)= sqrt ((x(i,3)-x(i,8))**2+(y(i,3)-y(i,8))**2 +(z(i,8)&
        -z(i,3))**2)

        r5253(i)= sqrt ((x(i,52)-x(i,53))**2+(y(i,52)-y(i,53))**2&
        +(z(i,52)-z(i,53))**2)

        r5455(i)= sqrt ((x(i,54)-x(i,55))**2+(y(i,54)-y(i,55))**2&
        +(z(i,54)-z(i,55))**2)

        r1516(i)= sqrt ((x(i,15)-x(i,16))**2+(y(i,15)-y(i,16))**2&
        +(z(i,15)-z(i,16))**2)

        r3435(i)= sqrt ((x(i,34)-x(i,35))**2+(y(i,34)-y(i,35))**2&
         + (z(i,34)-z(i,35))**2)

        r3637(i)= sqrt ((x(i,36)-x(i,37))**2+(y(i,36)-y(i,37))**2&
         + (z(i,36)-z(i,37))**2)        

        r6669(i)= sqrt ((x(i,66)-x(i,69))**2+(y(i,66)-y(i,69))**2&
        +(z(i,66)-z(i,69))**2)

        r7376(i)= sqrt ((x(i,73)-x(i,76))**2+(y(i,73)-y(i,76))**2&
        +(z(i,73)-z(i,76))**2)

        r8386(i)= sqrt ((x(i,83)-x(i,86))**2+(y(i,83)-y(i,86))**2&
        +(z(i,83)-z(i,86))**2)

        r62103(i)= sqrt ((x(i,62)-x(i,103))**2+(y(i,62)-y(i,103))**2&
        +(z(i,62)-z(i,103))**2)

        r93104(i)= sqrt ((x(i,93)-x(i,104))**2+(y(i,93)-y(i,104))**2&
        +(z(i,93)-z(i,104))**2)

       
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
      if (r93104(i) .ge. 6.0) then
        write(17,*)'r93104, ',r93104(i),i*0.130E-14
        endif
 
      if (r62103(i) .ge. 6.0) then
        write(17,*)'r62103, ',r62103(i),i*0.130E-14
        endif
      if (r3435(i) .ge. 6.0) then
        write(17,*)'r3435, ',r3435(i),i*0.130E-14
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
