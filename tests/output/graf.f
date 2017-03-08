

      implicit double precision (a-h,o-z)

      dimension e(300), ared(300)
      character lixo

      open(10,file='1.dat')
      do i = 1, 140
      read(10,*) lixo, lixo, lixo, lixo, lixo, lixo, lixo, lixo, lixo,
     *          lixo, ared(i)
      end do
      close(10)
      open(10,file='2.dat')
      do i = 1, 140
      read(10,*) lixo, lixo, lixo, lixo, lixo, lixo, lixo, e(i)
      end do

      do i = 1, 140
        if(ared(i).gt.0.d0) icol = 4
        if(ared(i).le.0.d0) icol = 2
        write(*,100) i, e(i) - e(140) + 1.d-8, icol
      end do
100   format(i6,tr2,e17.10,tr2,i6)

      end
