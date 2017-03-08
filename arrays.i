c
c   natoms: number of atoms of the system
c   kdim:   dimension of the problem (total number of basis functions)
c   kcont:  temporary integer to store kdim for contracted basis sets
c   ntotf:  total number of basis functions not counting 3*p etc.
c   nel:    number of electrons in the system
c   x,y,z:  cartesian coordinates of the atoms.
c   zatom:  charge of atoms' nuclei
c   s:      overlap matrix
c   hcore:  Core interactions matrix (Hcore)
c   abcd:   Two electron integral array (use jump to set indices)
c   ele:    character of the element
c



      integer natoms
      integer kdim
      integer kcont
      integer ntotf
      integer nel
      double precision x(maxatom), y(maxatom), z(maxatom)
      double precision zatom(maxatom)    
      double precision s(maxd,maxd), hcore(maxd,maxd)
      character*2 ele(maxatom)

      common / arrays / s, hcore, 
     *                  x, y, z, zatom,
     *                  natoms, kcont, kdim, ntotf, nel, 
     *                  ele

