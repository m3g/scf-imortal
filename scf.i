c
c c: the matrix of the coeficients
c p: density matrix
c f: fock matrix
c g: f = hcore + g
c
      double precision c(maxd,maxd), p(maxd, maxd),
     *                 f(maxd,maxd), g(maxd, maxd)

      common / scf / p, g

