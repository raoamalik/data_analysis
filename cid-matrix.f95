! This program calculates the bond lengths between atoms and reports
! time when
! any bond reaches 6A.to get the coordinates in a file: grep -A104 "Q
! !P"
! output.dat | grep -v -- "^--$" |grep -vE "Q
! P"> file
! cut -c3-35 file > qcoord.dat

        implicit real*8(a-h, o-z)
        integer, parameter :: rows=15000,columns=105
        dimension x(rows,columns),y(rows,columns),z(rows,columns)
        dimension s(rows,columns)
        dimension r1051(rows),r1052(rows),r1053(rows),r1054(rows)
        dimension r1055(rows),r1056(rows),r1057(rows),r1058(rows)
        dimension r1059(rows),r10510(rows),r10511(rows),r10512(rows)
        dimension r10513(rows),r10514(rows),r10515(rows),r10516(rows)
        dimension r10517(rows),r10518(rows),r10519(rows),r10520(rows)
        dimension r10521(rows),r10522(rows),r10523(rows),r10524(rows)
        dimension r10525(rows),r10526(rows),r10527(rows),r10528(rows)
        dimension r10529(rows),r10530(rows),r10531(rows),r10532(rows)
        dimension r10533(rows),r10534(rows),r10535(rows),r10536(rows)
        dimension r10537(rows),r10538(rows),r10539(rows),r10540(rows)
        dimension r10541(rows),r10542(rows),r10543(rows),r10544(rows)
        dimension r10545(rows),r10546(rows),r10547(rows),r10548(rows)
        dimension r10549(rows),r10550(rows),r10551(rows),r10552(rows)
        dimension r10553(rows),r10554(rows),r10555(rows),r10556(rows)
        dimension r10557(rows),r10558(rows),r10559(rows),r10560(rows)
        dimension r10561(rows),r10562(rows),r10563(rows),r10564(rows)
        dimension r10565(rows),r10566(rows),r10567(rows),r10568(rows)
        dimension r10569(rows),r10570(rows),r10571(rows),r10572(rows)
        dimension r10573(rows),r10574(rows),r10575(rows),r10576(rows)
        dimension r10577(rows),r10578(rows),r10579(rows),r10580(rows)
        dimension r10581(rows),r10582(rows),r10583(rows),r10584(rows)
        dimension r10585(rows),r10586(rows),r10587(rows),r10588(rows)
        dimension r10589(rows),r10590(rows),r10591(rows),r10592(rows)
        dimension r10593(rows),r10594(rows),r10595(rows),r10596(rows)
        dimension r10597(rows),r10598(rows),r10599(rows),r105100(rows)
        dimension r105101(rows),r105102(rows),r105103(rows)
        dimension r105104(rows)


        open(7, file='qcoord.dat')
        open(17, file='matrix')
!       read(5,*)s(i,j), x(i,j),y(i,j),z(i,j)
!        read(7, FMT="(3(F10.7,3x))") x(i,j),y(i,j),z(i,j)

       do i=1,rows
!        read(5,*) nst
        do j=1,columns
!        open(7, file='qcoord.dat')
       read(7,*) x(i,j),y(i,j),z(i,j)
!        read(7, FMT="(3(F10.7,3x))") x(i,j),y(i,j),z(i,j)
        enddo
        r1051(i)= sqrt ((x(i,105)-x(i,1))**2+(y(i,105)-y(i,1))**2&
         +(z(i,105)-z(i,1))**2)
        if (r1051(i)  .le. 3.0) then
                write (17,*) 'r1051, ', r1051(i),  i*0.130E-14
                endif
        r1052(i)= sqrt ((x(i,105)-x(i,2))**2+(y(i,105)-y(i,2))**2&
         +(z(i,105)-z(i,2))**2)
        if (r1052(i)  .le. 3.0) then
                write (17,*) 'r1052, ', r1052(i),  i*0.130E-14
                endif
        r1053(i)= sqrt ((x(i,105)-x(i,3))**2+(y(i,105)-y(i,3))**2&
         +(z(i,105)-z(i,3))**2)
        if (r1053(i)  .le. 3.0) then
                write (17,*) 'r1053, ', r1053(i),  i*0.130E-14
                endif
        r1054(i)= sqrt ((x(i,105)-x(i,4))**2+(y(i,105)-y(i,4))**2&
         +(z(i,105)-z(i,4))**2)
        if (r1054(i)  .le. 3.0) then
                write (17,*) 'r1054, ', r1054(i),  i*0.130E-14
                endif
        r1055(i)= sqrt ((x(i,105)-x(i,5))**2+(y(i,105)-y(i,5))**2&
         +(z(i,105)-z(i,5))**2)
        if (r1055(i)  .le. 3.0) then
                write (17,*) 'r1055, ', r1055(i),  i*0.130E-14
                endif
        r1056(i)= sqrt ((x(i,105)-x(i,6))**2+(y(i,105)-y(i,6))**2&
         +(z(i,105)-z(i,6))**2)
        if (r1056(i)  .le. 3.0) then
                write (17,*) 'r1056, ', r1056(i),  i*0.130E-14
                endif
        r1057(i)= sqrt ((x(i,105)-x(i,7))**2+(y(i,105)-y(i,7))**2&
         +(z(i,105)-z(i,7))**2)
        if (r1057(i)  .le. 3.0) then
                write (17,*) 'r1057, ', r1057(i),  i*0.130E-14
                endif
        r1058(i)= sqrt ((x(i,105)-x(i,8))**2+(y(i,105)-y(i,8))**2&
         +(z(i,105)-z(i,8))**2)
        if (r1058(i)  .le. 3.0) then
                write (17,*) 'r1058, ', r1058(i),  i*0.130E-14
                endif
        r1059(i)= sqrt ((x(i,105)-x(i,9))**2+(y(i,105)-y(i,9))**2&
         +(z(i,105)-z(i,9))**2)
        if (r1059(i)  .le. 3.0) then
                write (17,*) 'r1059, ', r1059(i),  i*0.130E-14
                endif
        r10510(i)= sqrt ((x(i,105)-x(i,10))**2+(y(i,105)-y(i,10))**2&
         +(z(i,105)-z(i,10))**2)
        if (r10510(i)  .le. 3.0) then
                write (17,*) 'r10510, ', r10510(i),  i*0.130E-14
                endif
        r10511(i)= sqrt ((x(i,105)-x(i,11))**2+(y(i,105)-y(i,11))**2&
         +(z(i,105)-z(i,11))**2)
        if (r10511(i)  .le. 3.0) then
                write (17,*) 'r10511, ', r10511(i),  i*0.130E-14
                endif
        r10512(i)= sqrt ((x(i,105)-x(i,12))**2+(y(i,105)-y(i,12))**2&
         +(z(i,105)-z(i,12))**2)
        if (r10512(i)  .le. 3.0) then
                write (17,*) 'r10512, ', r10512(i),  i*0.130E-14
                endif
        r10513(i)= sqrt ((x(i,105)-x(i,13))**2+(y(i,105)-y(i,13))**2&
         +(z(i,105)-z(i,13))**2)
        if (r10513(i)  .le. 3.0) then
                write (17,*) 'r10513, ', r10513(i),  i*0.130E-14
                endif
        r10514(i)= sqrt ((x(i,105)-x(i,14))**2+(y(i,105)-y(i,14))**2&
         +(z(i,105)-z(i,14))**2)
        if (r10514(i)  .le. 3.0) then
                write (17,*) 'r10514, ', r10514(i),  i*0.130E-14
                endif
        r10515(i)= sqrt ((x(i,105)-x(i,15))**2+(y(i,105)-y(i,15))**2&
         +(z(i,105)-z(i,15))**2)
        if (r10515(i)  .le. 3.0) then
                write (17,*) 'r10515, ', r10515(i),  i*0.130E-14
                endif
        r10516(i)= sqrt ((x(i,105)-x(i,16))**2+(y(i,105)-y(i,16))**2&
         +(z(i,105)-z(i,16))**2)
        if (r10516(i)  .le. 3.0) then
                write (17,*) 'r10516, ', r10516(i),  i*0.130E-14
                endif
        r10517(i)= sqrt ((x(i,105)-x(i,17))**2+(y(i,105)-y(i,17))**2&
         +(z(i,105)-z(i,17))**2)
        if (r10517(i)  .le. 3.0) then
                write (17,*) 'r10517, ', r10517(i),  i*0.130E-14
                endif
        r10518(i)= sqrt ((x(i,105)-x(i,18))**2+(y(i,105)-y(i,18))**2&
         +(z(i,105)-z(i,18))**2)
        if (r10518(i)  .le. 3.0) then
                write (17,*) 'r10518, ', r10518(i),  i*0.130E-14
                endif
        r10519(i)= sqrt ((x(i,105)-x(i,19))**2+(y(i,105)-y(i,19))**2&
         +(z(i,105)-z(i,19))**2)
        if (r10519(i)  .le. 3.0) then
                write (17,*) 'r10519, ', r10519(i),  i*0.130E-14
                endif
        r10520(i)= sqrt ((x(i,105)-x(i,20))**2+(y(i,105)-y(i,20))**2&
         +(z(i,105)-z(i,20))**2)
        if (r10520(i)  .le. 3.0) then
                write (17,*) 'r10520, ', r10520(i),  i*0.130E-14
                endif
        r10521(i)= sqrt ((x(i,105)-x(i,21))**2+(y(i,105)-y(i,21))**2&
         +(z(i,105)-z(i,21))**2)
        if (r10521(i)  .le. 3.0) then
                write (17,*) 'r10521, ', r10521(i),  i*0.130E-14
                endif
        r10522(i)= sqrt ((x(i,105)-x(i,22))**2+(y(i,105)-y(i,22))**2&
         +(z(i,105)-z(i,22))**2)
        if (r10522(i)  .le. 3.0) then
                write (17,*) 'r10522, ', r10522(i),  i*0.130E-14
                endif
        r10523(i)= sqrt ((x(i,105)-x(i,23))**2+(y(i,105)-y(i,23))**2&
         +(z(i,105)-z(i,23))**2)
        if (r10523(i)  .le. 3.0) then
                write (17,*) 'r10523, ', r10523(i),  i*0.130E-14
                endif
        r10524(i)= sqrt ((x(i,105)-x(i,24))**2+(y(i,105)-y(i,24))**2&
         +(z(i,105)-z(i,24))**2)
        if (r10524(i)  .le. 3.0) then
                write (17,*) 'r10524, ', r10524(i),  i*0.130E-14
                endif
        r10525(i)= sqrt ((x(i,105)-x(i,25))**2+(y(i,105)-y(i,25))**2&
         +(z(i,105)-z(i,25))**2)
        if (r10525(i)  .le. 3.0) then
                write (17,*) 'r10525, ', r10525(i),  i*0.130E-14
                endif
        r10526(i)= sqrt ((x(i,105)-x(i,26))**2+(y(i,105)-y(i,26))**2&
         +(z(i,105)-z(i,26))**2)
        if (r10526(i)  .le. 3.0) then
                write (17,*) 'r10526, ', r10526(i),  i*0.130E-14
                endif
        r10527(i)= sqrt ((x(i,105)-x(i,27))**2+(y(i,105)-y(i,27))**2&
         +(z(i,105)-z(i,27))**2)
        if (r10527(i)  .le. 3.0) then
                write (17,*) 'r10527, ', r10527(i),  i*0.130E-14
                endif
        r10528(i)= sqrt ((x(i,105)-x(i,28))**2+(y(i,105)-y(i,28))**2&
         +(z(i,105)-z(i,28))**2)
        if (r10528(i)  .le. 3.0) then
                write (17,*) 'r10528, ', r10528(i),  i*0.130E-14
                endif
        r10529(i)= sqrt ((x(i,105)-x(i,29))**2+(y(i,105)-y(i,29))**2&
         +(z(i,105)-z(i,29))**2)
        if (r10529(i)  .le. 3.0) then
                write (17,*) 'r10529, ', r10529(i),  i*0.130E-14
                endif
        r10530(i)= sqrt ((x(i,105)-x(i,30))**2+(y(i,105)-y(i,30))**2&
         +(z(i,105)-z(i,30))**2)
        if (r10530(i)  .le. 3.0) then
                write (17,*) 'r10530, ', r10530(i),  i*0.130E-14
                endif
        r10531(i)= sqrt ((x(i,105)-x(i,31))**2+(y(i,105)-y(i,31))**2&
         +(z(i,105)-z(i,31))**2)
        if (r10531(i)  .le. 3.0) then
                write (17,*) 'r10531, ', r10531(i),  i*0.130E-14
                endif
        r10532(i)= sqrt ((x(i,105)-x(i,32))**2+(y(i,105)-y(i,32))**2&
         +(z(i,105)-z(i,32))**2)
        if (r10532(i)  .le. 3.0) then
                write (17,*) 'r10532, ', r10532(i),  i*0.130E-14
                endif
        r10533(i)= sqrt ((x(i,105)-x(i,33))**2+(y(i,105)-y(i,33))**2&
         +(z(i,105)-z(i,33))**2)
        if (r10533(i)  .le. 3.0) then
                write (17,*) 'r10533, ', r10533(i),  i*0.130E-14
                endif
        r10534(i)= sqrt ((x(i,105)-x(i,34))**2+(y(i,105)-y(i,34))**2&
         +(z(i,105)-z(i,34))**2)
        if (r10534(i)  .le. 3.0) then
                write (17,*) 'r10534, ', r10534(i),  i*0.130E-14
                endif
        r10535(i)= sqrt ((x(i,105)-x(i,35))**2+(y(i,105)-y(i,35))**2&
         +(z(i,105)-z(i,35))**2)
        if (r10535(i)  .le. 3.0) then
                write (17,*) 'r10535, ', r10535(i),  i*0.130E-14
                endif
        r10536(i)= sqrt ((x(i,105)-x(i,36))**2+(y(i,105)-y(i,36))**2&
         +(z(i,105)-z(i,36))**2)
        if (r10536(i)  .le. 3.0) then
                write (17,*) 'r10536, ', r10536(i),  i*0.130E-14
                endif
        r10537(i)= sqrt ((x(i,105)-x(i,37))**2+(y(i,105)-y(i,37))**2&
         +(z(i,105)-z(i,37))**2)
        if (r10537(i)  .le. 3.0) then
                write (17,*) 'r10537, ', r10537(i),  i*0.130E-14
                endif
        r10538(i)= sqrt ((x(i,105)-x(i,38))**2+(y(i,105)-y(i,38))**2&
         +(z(i,105)-z(i,38))**2)
        if (r10538(i)  .le. 3.0) then
                write (17,*) 'r10538, ', r10538(i),  i*0.130E-14
                endif
        r10539(i)= sqrt ((x(i,105)-x(i,39))**2+(y(i,105)-y(i,39))**2&
         +(z(i,105)-z(i,39))**2)
        if (r10539(i)  .le. 3.0) then
                write (17,*) 'r10539, ', r10539(i),  i*0.130E-14
                endif
        r10540(i)= sqrt ((x(i,105)-x(i,40))**2+(y(i,105)-y(i,40))**2&
         +(z(i,105)-z(i,40))**2)
        if (r10540(i)  .le. 3.0) then
                write (17,*) 'r10540, ', r10540(i),  i*0.130E-14
                endif
        r10541(i)= sqrt ((x(i,105)-x(i,41))**2+(y(i,105)-y(i,41))**2&
         +(z(i,105)-z(i,41))**2)
        if (r10541(i)  .le. 3.0) then
                write (17,*) 'r10541, ', r10541(i),  i*0.130E-14
                endif
        r10542(i)= sqrt ((x(i,105)-x(i,42))**2+(y(i,105)-y(i,42))**2&
         +(z(i,105)-z(i,42))**2)
        if (r10542(i)  .le. 3.0) then
                write (17,*) 'r10542, ', r10542(i),  i*0.130E-14
                endif
        r10543(i)= sqrt ((x(i,105)-x(i,43))**2+(y(i,105)-y(i,43))**2&
         +(z(i,105)-z(i,43))**2)
        if (r10543(i)  .le. 3.0) then
                write (17,*) 'r10543, ', r10543(i),  i*0.130E-14
                endif
        r10544(i)= sqrt ((x(i,105)-x(i,44))**2+(y(i,105)-y(i,44))**2&
         +(z(i,105)-z(i,44))**2)
        if (r10544(i)  .le. 3.0) then
                write (17,*) 'r10544, ', r10544(i),  i*0.130E-14
                endif
        r10545(i)= sqrt ((x(i,105)-x(i,45))**2+(y(i,105)-y(i,45))**2&
         +(z(i,105)-z(i,45))**2)
        if (r10545(i)  .le. 3.0) then
                write (17,*) 'r10545, ', r10545(i),  i*0.130E-14
                endif
        r10546(i)= sqrt ((x(i,105)-x(i,46))**2+(y(i,105)-y(i,46))**2&
         +(z(i,105)-z(i,46))**2)
        if (r10546(i)  .le. 3.0) then
                write (17,*) 'r10546, ', r10546(i),  i*0.130E-14
                endif
        r10547(i)= sqrt ((x(i,105)-x(i,47))**2+(y(i,105)-y(i,47))**2&
         +(z(i,105)-z(i,47))**2)
        if (r10547(i)  .le. 3.0) then
                write (17,*) 'r10547, ', r10547(i),  i*0.130E-14
                endif
        r10548(i)= sqrt ((x(i,105)-x(i,48))**2+(y(i,105)-y(i,48))**2&
         +(z(i,105)-z(i,48))**2)
        if (r10548(i)  .le. 3.0) then
                write (17,*) 'r10548, ', r10548(i),  i*0.130E-14
                endif
        r10549(i)= sqrt ((x(i,105)-x(i,49))**2+(y(i,105)-y(i,49))**2&
         +(z(i,105)-z(i,49))**2)
        if (r10549(i)  .le. 3.0) then
                write (17,*) 'r10549, ', r10549(i),  i*0.130E-14
                endif
        r10550(i)= sqrt ((x(i,105)-x(i,50))**2+(y(i,105)-y(i,50))**2&
         +(z(i,105)-z(i,50))**2)
        if (r10550(i)  .le. 3.0) then
                write (17,*) 'r10550, ', r10550(i),  i*0.130E-14
                endif
        r10551(i)= sqrt ((x(i,105)-x(i,51))**2+(y(i,105)-y(i,51))**2&
         +(z(i,105)-z(i,51))**2)
        if (r10551(i)  .le. 3.0) then
                write (17,*) 'r10551, ', r10551(i),  i*0.130E-14
                endif
        r10552(i)= sqrt ((x(i,105)-x(i,52))**2+(y(i,105)-y(i,52))**2&
         +(z(i,105)-z(i,52))**2)
        if (r10552(i)  .le. 3.0) then
                write (17,*) 'r10552, ', r10552(i),  i*0.130E-14
                endif
        r10553(i)= sqrt ((x(i,105)-x(i,53))**2+(y(i,105)-y(i,53))**2&
         +(z(i,105)-z(i,53))**2)
        if (r10553(i)  .le. 3.0) then
                write (17,*) 'r10553, ', r10553(i),  i*0.130E-14
                endif
        r10554(i)= sqrt ((x(i,105)-x(i,54))**2+(y(i,105)-y(i,54))**2&
         +(z(i,105)-z(i,54))**2)
        if (r10554(i)  .le. 3.0) then
                write (17,*) 'r10554, ', r10554(i),  i*0.130E-14
                endif
        r10555(i)= sqrt ((x(i,105)-x(i,55))**2+(y(i,105)-y(i,55))**2&
         +(z(i,105)-z(i,55))**2)
        if (r10555(i)  .le. 3.0) then
                write (17,*) 'r10555, ', r10555(i),  i*0.130E-14
                endif
        r10556(i)= sqrt ((x(i,105)-x(i,56))**2+(y(i,105)-y(i,56))**2&
         +(z(i,105)-z(i,56))**2)
        if (r10556(i)  .le. 3.0) then
                write (17,*) 'r10556, ', r10556(i),  i*0.130E-14
                endif
        r10557(i)= sqrt ((x(i,105)-x(i,57))**2+(y(i,105)-y(i,57))**2&
         +(z(i,105)-z(i,57))**2)
        if (r10557(i)  .le. 3.0) then
                write (17,*) 'r10557, ', r10557(i),  i*0.130E-14
                endif
        r10558(i)= sqrt ((x(i,105)-x(i,58))**2+(y(i,105)-y(i,58))**2&
         +(z(i,105)-z(i,58))**2)
        if (r10558(i)  .le. 3.0) then
                write (17,*) 'r10558, ', r10558(i),  i*0.130E-14
                endif
        r10559(i)= sqrt ((x(i,105)-x(i,59))**2+(y(i,105)-y(i,59))**2&
         +(z(i,105)-z(i,59))**2)
        if (r10559(i)  .le. 3.0) then
                write (17,*) 'r10559, ', r10559(i),  i*0.130E-14
                endif
        r10560(i)= sqrt ((x(i,105)-x(i,60))**2+(y(i,105)-y(i,60))**2&
         +(z(i,105)-z(i,60))**2)
        if (r10560(i)  .le. 3.0) then
                write (17,*) 'r10560, ', r10560(i),  i*0.130E-14
                endif
        r10561(i)= sqrt ((x(i,105)-x(i,61))**2+(y(i,105)-y(i,61))**2&
         +(z(i,105)-z(i,61))**2)
        if (r10561(i)  .le. 3.0) then
                write (17,*) 'r10561, ', r10561(i),  i*0.130E-14
                endif
        r10562(i)= sqrt ((x(i,105)-x(i,62))**2+(y(i,105)-y(i,62))**2&
         +(z(i,105)-z(i,62))**2)
        if (r10562(i)  .le. 3.0) then
                write (17,*) 'r10562, ', r10562(i),  i*0.130E-14
                endif
        r10563(i)= sqrt ((x(i,105)-x(i,63))**2+(y(i,105)-y(i,63))**2&
         +(z(i,105)-z(i,63))**2)
        if (r10563(i)  .le. 3.0) then
                write (17,*) 'r10563, ', r10563(i),  i*0.130E-14
                endif
        r10564(i)= sqrt ((x(i,105)-x(i,64))**2+(y(i,105)-y(i,64))**2&
         +(z(i,105)-z(i,64))**2)
        if (r10564(i)  .le. 3.0) then
                write (17,*) 'r10564, ', r10564(i),  i*0.130E-14
                endif
        r10565(i)= sqrt ((x(i,105)-x(i,65))**2+(y(i,105)-y(i,65))**2&
         +(z(i,105)-z(i,65))**2)
        if (r10565(i)  .le. 3.0) then
                write (17,*) 'r10565, ', r10565(i),  i*0.130E-14
                endif
        r10566(i)= sqrt ((x(i,105)-x(i,66))**2+(y(i,105)-y(i,66))**2&
         +(z(i,105)-z(i,66))**2)
        if (r10566(i)  .le. 3.0) then
                write (17,*) 'r10566, ', r10566(i),  i*0.130E-14
                endif
        r10567(i)= sqrt ((x(i,105)-x(i,67))**2+(y(i,105)-y(i,67))**2&
         +(z(i,105)-z(i,67))**2)
        if (r10567(i)  .le. 3.0) then
                write (17,*) 'r10567, ', r10567(i),  i*0.130E-14
                endif
        r10568(i)= sqrt ((x(i,105)-x(i,68))**2+(y(i,105)-y(i,68))**2&
         +(z(i,105)-z(i,68))**2)
        if (r10568(i)  .le. 3.0) then
                write (17,*) 'r10568, ', r10568(i),  i*0.130E-14
                endif
        r10569(i)= sqrt ((x(i,105)-x(i,69))**2+(y(i,105)-y(i,69))**2&
         +(z(i,105)-z(i,69))**2)
        if (r10569(i)  .le. 3.0) then
                write (17,*) 'r10569, ', r10569(i),  i*0.130E-14
                endif
        r10570(i)= sqrt ((x(i,105)-x(i,70))**2+(y(i,105)-y(i,70))**2&
         +(z(i,105)-z(i,70))**2)
        if (r10570(i)  .le. 3.0) then
                write (17,*) 'r10570, ', r10570(i),  i*0.130E-14
                endif
        r10571(i)= sqrt ((x(i,105)-x(i,71))**2+(y(i,105)-y(i,71))**2&
         +(z(i,105)-z(i,71))**2)
        if (r10571(i)  .le. 3.0) then
                write (17,*) 'r10571, ', r10571(i),  i*0.130E-14
                endif
        r10572(i)= sqrt ((x(i,105)-x(i,72))**2+(y(i,105)-y(i,72))**2&
         +(z(i,105)-z(i,72))**2)
        if (r10572(i)  .le. 3.0) then
                write (17,*) 'r10572, ', r10572(i),  i*0.130E-14
                endif
        r10573(i)= sqrt ((x(i,105)-x(i,73))**2+(y(i,105)-y(i,73))**2&
         +(z(i,105)-z(i,73))**2)
        if (r10573(i)  .le. 3.0) then
                write (17,*) 'r10573, ', r10573(i),  i*0.130E-14
                endif
        r10574(i)= sqrt ((x(i,105)-x(i,74))**2+(y(i,105)-y(i,74))**2&
         +(z(i,105)-z(i,74))**2)
        if (r10574(i)  .le. 3.0) then
                write (17,*) 'r10574, ', r10574(i),  i*0.130E-14
                endif
        r10575(i)= sqrt ((x(i,105)-x(i,75))**2+(y(i,105)-y(i,75))**2&
         +(z(i,105)-z(i,75))**2)
        if (r10575(i)  .le. 3.0) then
                write (17,*) 'r10575, ', r10575(i),  i*0.130E-14
                endif
        r10576(i)= sqrt ((x(i,105)-x(i,76))**2+(y(i,105)-y(i,76))**2&
         +(z(i,105)-z(i,76))**2)
        if (r10576(i)  .le. 3.0) then
                write (17,*) 'r10576, ', r10576(i),  i*0.130E-14
                endif
        r10577(i)= sqrt ((x(i,105)-x(i,77))**2+(y(i,105)-y(i,77))**2&
         +(z(i,105)-z(i,77))**2)
        if (r10577(i)  .le. 3.0) then
                write (17,*) 'r10577, ', r10577(i),  i*0.130E-14
                endif
        r10578(i)= sqrt ((x(i,105)-x(i,78))**2+(y(i,105)-y(i,78))**2&
         +(z(i,105)-z(i,78))**2)
        if (r10578(i)  .le. 3.0) then
                write (17,*) 'r10578, ', r10578(i),  i*0.130E-14
                endif
        r10579(i)= sqrt ((x(i,105)-x(i,79))**2+(y(i,105)-y(i,79))**2&
         +(z(i,105)-z(i,79))**2)
        if (r10579(i)  .le. 3.0) then
                write (17,*) 'r10579, ', r10579(i),  i*0.130E-14
                endif
        r10580(i)= sqrt ((x(i,105)-x(i,80))**2+(y(i,105)-y(i,80))**2&
         +(z(i,105)-z(i,80))**2)
        if (r10580(i)  .le. 3.0) then
                write (17,*) 'r10580, ', r10580(i),  i*0.130E-14
                endif
        r10581(i)= sqrt ((x(i,105)-x(i,81))**2+(y(i,105)-y(i,81))**2&
         +(z(i,105)-z(i,81))**2)
        if (r10581(i)  .le. 3.0) then
                write (17,*) 'r10581, ', r10581(i),  i*0.130E-14
                endif
        r10582(i)= sqrt ((x(i,105)-x(i,82))**2+(y(i,105)-y(i,82))**2&
         +(z(i,105)-z(i,82))**2)
        if (r10582(i)  .le. 3.0) then
                write (17,*) 'r10582, ', r10582(i),  i*0.130E-14
                endif
        r10583(i)= sqrt ((x(i,105)-x(i,83))**2+(y(i,105)-y(i,83))**2&
         +(z(i,105)-z(i,83))**2)
        if (r10583(i)  .le. 3.0) then
                write (17,*) 'r10583, ', r10583(i),  i*0.130E-14
                endif
        r10584(i)= sqrt ((x(i,105)-x(i,84))**2+(y(i,105)-y(i,84))**2&
         +(z(i,105)-z(i,84))**2)
        if (r10584(i)  .le. 3.0) then
                write (17,*) 'r10584, ', r10584(i),  i*0.130E-14
                endif
        r10585(i)= sqrt ((x(i,105)-x(i,85))**2+(y(i,105)-y(i,85))**2&
         +(z(i,105)-z(i,85))**2)
        if (r10585(i)  .le. 3.0) then
                write (17,*) 'r10585, ', r10585(i),  i*0.130E-14
                endif
        r10586(i)= sqrt ((x(i,105)-x(i,86))**2+(y(i,105)-y(i,86))**2&
         +(z(i,105)-z(i,86))**2)
        if (r10586(i)  .le. 3.0) then
                write (17,*) 'r10586, ', r10586(i),  i*0.130E-14
                endif
        r10587(i)= sqrt ((x(i,105)-x(i,87))**2+(y(i,105)-y(i,87))**2&
         +(z(i,105)-z(i,87))**2)
        if (r10587(i)  .le. 3.0) then
                write (17,*) 'r10587, ', r10587(i),  i*0.130E-14
                endif
        r10588(i)= sqrt ((x(i,105)-x(i,88))**2+(y(i,105)-y(i,88))**2&
         +(z(i,105)-z(i,88))**2)
        if (r10588(i)  .le. 3.0) then
                write (17,*) 'r10588, ', r10588(i),  i*0.130E-14
                endif
        r10589(i)= sqrt ((x(i,105)-x(i,89))**2+(y(i,105)-y(i,89))**2&
         +(z(i,105)-z(i,89))**2)
        if (r10589(i)  .le. 3.0) then
                write (17,*) 'r10589, ', r10589(i),  i*0.130E-14
                endif
        r10590(i)= sqrt ((x(i,105)-x(i,90))**2+(y(i,105)-y(i,90))**2&
         +(z(i,105)-z(i,90))**2)
        if (r10590(i)  .le. 3.0) then
                write (17,*) 'r10590, ', r10590(i),  i*0.130E-14
                endif
        r10591(i)= sqrt ((x(i,105)-x(i,91))**2+(y(i,105)-y(i,91))**2&
         +(z(i,105)-z(i,91))**2)
        if (r10591(i)  .le. 3.0) then
                write (17,*) 'r10591, ', r10591(i),  i*0.130E-14
                endif
        r10592(i)= sqrt ((x(i,105)-x(i,92))**2+(y(i,105)-y(i,92))**2&
         +(z(i,105)-z(i,92))**2)
        if (r10592(i)  .le. 3.0) then
                write (17,*) 'r10592, ', r10592(i),  i*0.130E-14
                endif
        r10593(i)= sqrt ((x(i,105)-x(i,93))**2+(y(i,105)-y(i,93))**2&
         +(z(i,105)-z(i,93))**2)
        if (r10593(i)  .le. 3.0) then
                write (17,*) 'r10593, ', r10593(i),  i*0.130E-14
                endif
        r10594(i)= sqrt ((x(i,105)-x(i,94))**2+(y(i,105)-y(i,94))**2&
         +(z(i,105)-z(i,94))**2)
        if (r10594(i)  .le. 3.0) then
                write (17,*) 'r10594, ', r10594(i),  i*0.130E-14
                endif
        r10595(i)= sqrt ((x(i,105)-x(i,95))**2+(y(i,105)-y(i,95))**2&
         +(z(i,105)-z(i,95))**2)
        if (r10595(i)  .le. 3.0) then
                write (17,*) 'r10595, ', r10595(i),  i*0.130E-14
                endif
        r10596(i)= sqrt ((x(i,105)-x(i,96))**2+(y(i,105)-y(i,96))**2&
         +(z(i,105)-z(i,96))**2)
        if (r10596(i)  .le. 3.0) then
                write (17,*) 'r10596, ', r10596(i),  i*0.130E-14
                endif
        r10597(i)= sqrt ((x(i,105)-x(i,97))**2+(y(i,105)-y(i,97))**2&
         +(z(i,105)-z(i,97))**2)
        if (r10597(i)  .le. 3.0) then
                write (17,*) 'r10597, ', r10597(i),  i*0.130E-14
                endif
        r10598(i)= sqrt ((x(i,105)-x(i,98))**2+(y(i,105)-y(i,98))**2&
         +(z(i,105)-z(i,98))**2)
        if (r10598(i)  .le. 3.0) then
                write (17,*) 'r10598, ', r10598(i),  i*0.130E-14
                endif
        r10599(i)= sqrt ((x(i,105)-x(i,99))**2+(y(i,105)-y(i,99))**2&
         +(z(i,105)-z(i,99))**2)
        if (r10599(i)  .le. 3.0) then
                write (17,*) 'r10599, ', r10599(i),  i*0.130E-14
                endif
        r105100(i)= sqrt ((x(i,105)-x(i,100))**2+(y(i,105)-y(i,100))**2&
         +(z(i,105)-z(i,100))**2)
        if (r105100(i)  .le. 3.0) then
                write (17,*) 'r105100, ', r105100(i),  i*0.130E-14
                endif
        r105101(i)= sqrt ((x(i,105)-x(i,101))**2+(y(i,105)-y(i,101))**2&
         +(z(i,105)-z(i,101))**2)
        if (r105101(i)  .le. 3.0) then
                write (17,*) 'r105101, ', r105101(i),  i*0.130E-14
                endif
        r105102(i)= sqrt ((x(i,105)-x(i,102))**2+(y(i,105)-y(i,102))**2&
         +(z(i,105)-z(i,102))**2)
        if (r105102(i)  .le. 3.0) then
                write (17,*) 'r105102, ', r105102(i),  i*0.130E-14
                endif
        r105103(i)= sqrt ((x(i,105)-x(i,103))**2+(y(i,105)-y(i,103))**2&
         +(z(i,105)-z(i,103))**2)
        if (r105103(i)  .le. 3.0) then
                write (17,*) 'r105103, ', r105103(i),  i*0.130E-14
                endif
        r105104(i)= sqrt ((x(i,105)-x(i,104))**2+(y(i,105)-y(i,104))**2&
         +(z(i,105)-z(i,104))**2)
        if (r105104(i)  .le. 3.0) then
                write (17,*) 'r105104, ', r105104(i),  i*0.130E-14
                endif
     
        enddo
        stop
        end
