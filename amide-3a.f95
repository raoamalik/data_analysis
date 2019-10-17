! This program calculates the bond lengths between atoms and reports
! time when
! any bond reaches 6A.to get the coordinates in a file: grep -A104 "Q
! !P"
! output.dat | grep -v -- "^--$" |grep -vE "Q
! P"> file
! cut -c3-35 file > qcoord.dat

        implicit real*8(a-h, o-z)
        integer, parameter :: rows=40000,columns=104
        dimension x(rows,columns),y(rows,columns),z(rows,columns)
        dimension r915(rows)
        dimension r4152(rows)
        dimension r6670(rows)
        dimension r7380(rows)
        dimension r8396(rows)
        dimension r2036(rows)
        dimension r163(rows)

        open(7, file='qcoord.dat')
        open(17, file='distance')
       do i=1,rows
        do j=1,columns
       read(7,*) x(i,j),y(i,j),z(i,j)
        enddo
        r163(i)= sqrt ((x(i,1)-x(i,63))**2+(y(i,1)-y(i,63))**2+(z(i,1)&
        -z(i,63))**2)
        r915(i)= sqrt ((x(i,9)-x(i,15))**2+(y(i,9)-y(i,15))**2&
        +(z(i,9)-z(i,15))**2)
        r2036(i)= sqrt ((x(i,36)-x(i,20))**2+(y(i,36)-y(i,20))**2&
        +(z(i,36)-z(i,20))**2)
        r4152(i)= sqrt ((x(i,41)-x(i,52))**2+(y(i,41)-y(i,52))**2&
        +(z(i,41)-z(i,52))**2)
        r6670(i)= sqrt ((x(i,66)-x(i,70))**2+(y(i,66)-y(i,70))**2&
        +(z(i,66)-z(i,70))**2)
        r7380(i)= sqrt ((x(i,73)-x(i,80))**2+(y(i,73)-y(i,80))**2&
        +(z(i,73)-z(i,80))**2)
        r8396(i)= sqrt ((x(i,83)-x(i,96))**2+(y(i,83)-y(i,96))**2&
        +(z(i,83)-z(i,96))**2)

        if (r163(i)  .ge. 6.0) then
        write (17,*) 'r163, ', r163(i),  i*0.130E-14
        endif
      if (r915(i) .ge. 6.0) then
        write(17,*)'r915, ',r915(i),i*0.130E-14
        endif
      if (r2036(i) .ge. 6.0) then
        write(17,*)'r2036, ',r2036(i),i*0.130E-14
        endif
      if (r6670(i) .ge. 6.0) then
        write(17,*)'r6670, ',r6670(i),i*0.130E-14
        endif
      if (r7380(i) .ge. 6.0) then
        write(17,*)'r7380, ',r7380(i),i*0.130E-14
        endif
        if (r8396(i) .ge. 6.0) then
        write(17,*)'r8396, ',r8396(i),i*0.130E-14
        endif
      if (r4152(i) .ge. 6.0) then
        write(17,*)'r4152, ',r4152(i),i*0.130E-14
        endif  
        enddo
        stop
        end
