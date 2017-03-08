c
c xyzfile:   cartesian coordinate file
c basefile:  file where the basis is writen
c densfile:  initial guess for the density matrix 
c

        character*80 xyzfile
        character*80 basefile
        character*80 densfile

        common /files/ xyzfile, basefile, densfile
