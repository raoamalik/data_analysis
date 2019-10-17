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
        dimension r12(rows),r58(rows),r57(rows),r2224(rows),r1518(rows)
        dimension r48104(rows),r3235(rows),r6669(rows),r7376(rows)
        dimension r8386(rows),r4950(rows),r5263(rows),r6364(rows)
        dimension r6670(rows),r72102(rows),r7071(rows),r7579(rows)
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

        r12(i)= sqrt ((x(i,1)-x(i,2))**2+(y(i,1)-y(i,2))**2 +(z(i,1)&
        -z(i,2))**2)

        r57(i)= sqrt ((x(i,5)-x(i,7))**2+(y(i,5)-y(i,7))**2&
        +(z(i,5)-z(i,7))**2)

        r58(i)= sqrt ((x(i,5)-x(i,8))**2+(y(i,5)-y(i,8))**2&
        +(z(i,5)-z(i,8))**2)

        r2224(i)= sqrt ((x(i,22)-x(i,24))**2+(y(i,22)-y(i,24))**2&
        +(z(i,22)-z(i,24))**2)

        r1518(i)= sqrt ((x(i,15)-x(i,18))**2+(y(i,15)-y(i,18))**2&
        +(z(i,15)-z(i,18))**2)

        r48104(i)= sqrt ((x(i,48)-x(i,104))**2+(y(i,48)-y(i,104))**2&
        +(z(i,48)-z(i,104))**2)

        r3235(i)= sqrt ((x(i,32)-x(i,35))**2+(y(i,32)-y(i,35))**2&
         + (z(i,32)-z(i,35))**2)        

        r6669(i)= sqrt ((x(i,66)-x(i,69))**2+(y(i,66)-y(i,69))**2&
        +(z(i,66)-z(i,69))**2)

        r7376(i)= sqrt ((x(i,73)-x(i,76))**2+(y(i,73)-y(i,76))**2&
        +(z(i,73)-z(i,76))**2)

        r8386(i)= sqrt ((x(i,83)-x(i,86))**2+(y(i,83)-y(i,86))**2&
        +(z(i,83)-z(i,86))**2)
       
      if (r12(i)  .ge. 6.0) then
        write (17,*) 'r12, ', r12(i),  i*0.130E-14
        endif

      if (r57(i) .ge. 6.0) then
        write(17,*)'r57, ',r57(i),i*0.130E-14
        endif

      if (r58(i) .ge. 6.0) then
        write(17,*) 'r58, ', r58(i),  i*0.130E-14
        endif

      if (r1518(i) .ge. 6.0) then
        write(17,*)'r1518, ',r1518(i),i*0.130E-14
        endif

      if (r2224(i) .ge. 6.0) then
        write(17,*)'r2224, ',r2224(i),i*0.130E-14
        endif

      if (r48104(i) .ge. 6.0) then
        write(17,*)'r48104, ',r48104(i),i*0.130E-14
        endif

      if (r3235(i) .ge. 6.0) then
        write(17,*)'r3235, ',r3235(i),i*0.130E-14
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
