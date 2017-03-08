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

      subroutine fixpacc(c)

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
      dimension rvecs(maxd*maxd,nvecs),fvecs(maxd,maxd,nvecs)       

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
   10 format(t9,'0 E = ', e17.10)  

c Starting iterative procedure

      do while (dehf.gt.ecrit.and.nit.le.maxitfix)

c Performing a DIIS acceleration

        call diis(kdim,nel/2,nit,ss,rvecs,fvecs,f,c)

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
20      format(i9,' E = ', e17.10, ' D = ', e17.10 )  
      
        nit = nit + 1
	end do  

      return
      end


