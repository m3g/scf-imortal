c
c  sizes.i: define the maximum size of each array in the program
c
c   maxline: maximum number of lines of the input file
c   maxatom: maximum number of atoms in the system
c   maxbase: maximum number of basis functions for each atom
c   maxbt:   maximum number of basis functions of each type (s,p,d,f)
c            for each atom
c   maxl:    maximum angular moment of the basis functions
c   nvecs:   number of residual vectors used in DIIS acceleration
c
c  dependent dimensions
c   
c   maxd: maximum number of basis functions in the system as a hole.
c
c

      parameter(maxline = 100)
      parameter(maxatom = 3)
      parameter(maxbase = 24)
      parameter(maxbt = 28)
      parameter(maxl = 6)
      parameter(nvecs = 10)

c
      parameter(maxd = maxbase*maxatom)




