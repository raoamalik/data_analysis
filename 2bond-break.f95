program bondbreaking
!This program gives the time after a bond breaks in a trajectory (bond length>6A) 
implicit none
real :: x,y,z,r
integer :: natom, i, j
natoms=104

open(unit=11, file='output.dat')
do 10 i = 1, natoms
         read (11, *) x
do j = 1, 2
         r(i)(j) = sqrt((x(i) - x(j))**2 + (y(i) - y(j))**2 + (z(i) - z(j))**2)
 
         write(*,*) 'i =', i
         write(*,*) 'sum =', sum
        
  10 enddo
close(11)

end program bondbreaking

