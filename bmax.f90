        PROGRAM SORT_BMAX
        IMPLICIT NONE

        REAL :: BMAX
        REAL, ALLOCATABLE :: COM_DIS(:), EVIB(:), EREL(:)
        INTEGER :: COUNTA, I

        COUNTA = 150
        BMAX = 9.6D0
        OPEN(5,FILE='bmax.dat') 
        OPEN(444,FILE='temp2')

        ALLOCATE (COM_DIS(COUNTA), EVIB(COUNTA), EREL(COUNTA))

        DO I=1,COUNTA
        READ(5,*) COM_DIS(I), EVIB(I), EREL(I)
         IF (COM_DIS(I).LE.BMAX) THEN
         WRITE(444,*) COM_DIS(I), EVIB(I), EREL(I), BMAX
         ENDIF
        ENDDO

        END PROGRAM SORT_BMAX
