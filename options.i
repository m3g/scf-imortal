c
c contract: use, or not, contracted gaussian basis sets
c printlevel: well, printlevel, what else?
c iguess: type of guess of the density matrix.
c icharge: charge of the molecule, used to calculate nel
c maxitfix: maximum number of fix-point iterations
c ecrit: criteria for energy convergence
c
      logical contract
      integer printlevel
      integer iguess
      integer icharge
      integer maxitfix
      double precision ecrit
        
      common / options / ecrit, 
     *                   maxitfix, icharge, iguess, printlevel, 
     *                   contract
