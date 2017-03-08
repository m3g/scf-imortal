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
c   Program main: manages all the subroutines
c
c

      program main

      include 'sizes.i'
      include 'files.i'
      include 'kill.i'
      include 'arrays.i'
      include 'scf.i'

c Printing title

      write(*,10)
   10 format(/,t2,' Program SCF',/)
      
c Initializing arrays and arguments
      call initial

c Reading input file           
      write(*,*) ' Reading input file...'
      call getinp
      if(killnow) stop

c Reading cartesian coordinate file
      write(*,*) ' Reading cartesian coordinate file...'
      call readxyz

c Reading base parameter file
      write(*,*) ' Reading basis file...'
      call readbase

c Stopping if errors were found
      if(kill) stop

c Writing problem data
      call writedata

c Calculating electron integrals
      write(*,*) ' Calculating molecular integrals...'
      call integral

c Solving with the fmm algorithm
c      write(*,*) ' Starting the FMM algorithm ... '
c      call fmm
 
c Performing a fixpoint calculation
c      write(*,*) ' Starting fixed point SCF...'
c      call guess(nel/2,kdim,c,2)
c      call fixpoint(c)

      write(*,*) ' Starting accelerated fixed point SCF...'
      call guess(nel/2,kdim,c,0)
      call fixpacc(c)
 
      end



