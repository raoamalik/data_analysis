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
        dimension r90102(rows)

        open(7, file='qcoord.dat')
        open(17, file='r90102')
       do i=1,rows
        do j=1,columns
       read(7,*) x(i,j),y(i,j),z(i,j)
        enddo
        r90102(i)= sqrt ((x(i,90)-x(i,102))**2+(y(i,90)-y(i,102))**2&
        +(z(i,90)-z(i,102))**2)
 
      if (r90102(i) .ge. 6.0) then
        write(17,*)'r90102, ',r90102(i),i*0.130E-14
      endif
        enddo
        stop
        end
