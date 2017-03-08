c
c
c  ##########################################################
c  # 
c  #  Program 
c  #
c  #
c  #
c  #
c  #  Author: Leandro Martínez
c  #
c  ##########################################################
c
c   Subroutine fixpoint: solves the scf problem with the classic
c   fix-point optimization scheme.
c   
c

      subroutine fixpoint(c)

      implicit double precision (a-h,o-z)      

      include 'sizes.i'
      include 'arrays.i'
      include 'options.i'
      include 'files.i'
      include 'scf.i'

      integer nc(maxd)
      double precision ss(maxd,maxd)
      double precision u(maxd,maxd), ut(maxd,maxd)
      double precision temp(maxd,maxd), ftemp(maxd,maxd),
     *                 fl(maxd,maxd)

c Checking if the total number of electron is even

      if(mod(nel,2).ne.0) then
        write(*,*) ' ERROR: Odd number of electrons '
        stop
      end if

c Obtaining U and Ut such that Ut S U = I

      do i = 1, kdim
        do j = 1, kdim
          ss(i,j) = s(i,j)
        end do
      end do

      call diag(ss,temp,kdim)

      do i = 1, kdim
        ss(i,i) = 1.d0 / dsqrt(ss(i,i))
      end do

      call pmatrix(kdim, kdim, kdim, temp, ss, u)

      do i = 1, kdim
        do j = 1, kdim
          ut(i,j) = u(j,i)
        end do
      end do

c Getting initial point

c      call guess(ndim,kdim,c,iguess)

c
c Starting fix-point scf calculations
c

      nit = 1
      dehf = 1.d0
      ehf = 1.d0

c Computing the first fock matrix and energy

      call computefe(c,f,ehf)

	write(*,10) ehf
   10 format(t9,'0 Energia de Hartree-Fock = ', f10.4)  

c Starting iterative procedure

      do while (dehf.gt.ecrit.and.nit.le.maxitfix)

c Computing F' = Ut F U

        k = kdim
        call pmatrix(k, k, k, f, u, g)
        call pmatrix(k, k, k, ut, g, fl)         

c Diagonalize F'

        call diag(fl,c,kdim) 

c Ordering eigenvectors and eigenvalues from lower to higher

        do i = 1, kdim
          nc(i) = 1
          do j = 1, kdim
            if(fl(i,i).gt.fl(j,j)) nc(i) = nc(i) + 1
            if(fl(i,i).eq.fl(j,j).and.i.gt.j) nc(i) = nc(i) + 1
          end do
        end do   

        do i = 1, kdim
          ftemp(i,i) = fl(i,i)
          do j = 1, kdim
            temp(i,j) = c(i,j)
          end do
        end do
      
        do i = 1, kdim
          fl(nc(i),nc(i)) = ftemp(i,i)
          do j = 1, kdim
            c(j,nc(i)) = temp(j,i)
          end do
        end do  

c Computing the matrix C from the eigenvectors for F'

        k = kdim
        call pmatrix(k, k, k, u, c, g) 

        do i = 1, kdim
          do j = 1, kdim
            c(i,j) = g(i,j)
          end do
        end do 

c Computing the new fock and density matrix

        dehf = ehf
        call computefe(c,f,ehf)

        dehf = abs( dehf - ehf ) 
        write(*,20) nit, ehf, dehf
20      format(i9,' Energia de Hartree-Fock = ', f10.4, ' ±', e10.4 )  
      
        nit = nit + 1
	end do  

      return
      end


c Subroutine computefe: calculates the fock matrix f and the energy
c from the coeficient matrix
c

      subroutine computefe(c,f,ehf)

      implicit double precision (a-h,o-z)

      logical wtodo
      include 'sizes.i'
      include 'scf.i'
      include 'arrays.i'

c Oppening the integral's file

      open(20,file="intvalues.dat",access="direct",
     *        form="unformatted",recl=17) 

c Computing the density matrix

     	do i = 1, kdim
        do j = 1, kdim
          p(i,j) = 0.d0
          do in = 1, nel / 2
            p(i,j) = p(i,j) + 2.d0 * c(i,in) * c(j,in)
          end do
        end do
	end do
 
c Computing the G matrix ( F = Hcore + G )  

	do i1 = 1, kdim
        do i2 = 1, kdim
          g(i1,i2) = 0.d0

          do i3 = 1, kdim
            do i4 = 1, kdim

              CALL JUMP(I1,I2,I4,I3,N1,N2,N3,N4,WTODO)
              index = maxd*maxd*maxd*(n4-1) + 
     *                maxd*maxd*(n3-1) +
     *                maxd*(n2-1) + n1
              read(20,rec=index) a1e

              CALL JUMP(I1,I3,I4,I2,N1,N2,N3,N4,WTODO)
              index = maxd*maxd*maxd*(n4-1) + 
     *                maxd*maxd*(n3-1) +
     *                maxd*(n2-1) + n1
              read(20,rec=index) a2e

              g(i1,i2) = g(i1,i2) + p(i3,i4) * ( a1e - 0.5d0 * a2e )

            end do
          end do	
        end do
	end do       

c Computing the Fock matrix

      do i = 1, kdim
        do j = 1, kdim 
          f(i,j) = hcore(i,j) + g(i,j)
        end do
	end do

      close(20) 

c Computing the energy

      ehf = 0.d0      
      do i = 1, natoms - 1
        do j = i + 1, natoms
          d = (x(i) - x(j))**2 + 
     *        (y(i) - y(j))**2 +
     *        (z(i) - z(j))**2
          d = dsqrt(d)
          ehf = ehf + zatom(i)*zatom(j) / d
        end do
      end do

     	do i = 1, kdim
        do j = 1, kdim
          ehf = ehf + 0.5 * p(j,i) * (hcore(i,j) + f(i,j))
        end do
	end do

      return
      end
 
