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
c   Subroutine readxyz: reads xyz cartesian coordinate file
c   
c

      subroutine readxyz

      include 'sizes.i'         
      include 'kill.i'
      include 'files.i'
      include 'arrays.i'
       
      integer i
      character*80 title 
     
c Opening cartesian coordinate file

      open(10, file = xyzfile, err=30, status = 'old')

c Read the number of atoms in the system

      read(10,*) natoms

      if(natoms.gt.maxatom) then
        write(*,15)
   15   format(' ERROR: Number of atoms greater than maxatom.')
        kill = .true.
        goto 50
      end if

c Read title of the system, if any

      read(10,20) title
   20 format(a80)
     
c Read elements and cartesian coordinates

      do i = 1, natoms
        read(10,*,err=30,end=30) ele(i), x(i), y(i), z(i)
      end do

      close(10)
      return

   30 write(*,40)
   40 format(' ERROR: Could not find or read cartesian',
     &                               ' coordinate file')
      kill = .true.

   50 close(10)

      return
      end

















