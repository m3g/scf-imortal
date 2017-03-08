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
c   Subroutine initial: initializes the variables and set 
c   default parameters.
c
c
      subroutine initial

      include 'kill.i'
      include 'sizes.i'
      include 'files.i'
      include 'options.i'
      include 'arrays.i'

c Setting the dimension of the space to 0
      kdim = 0

c Setting kill signal to false
      kill = .false.

c Setting printlevel default
      printlevel = 1

c Setting the use of basis as uncontracted
      contract = .false.

c Setting the presence of the coordinate file to false
      xyzfile(1:1) = '#'

c Setting the presence of the base data file to false
      basefile(1:1) = '#'

c Setting the guess parameter to nule density matrix
      iguess = 0

c Setting the default criteria for energy convergence
      ecrit = 1.d-10

c Setting the default maximum number of iterations for fixpoint method
      maxitfix = 50

c Setting the default charge to 0.d0
      icharge = 0

      return
      end
















