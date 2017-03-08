c
c  nang:  angular moment of the function of maximum ang. moment of each
c         atom.
c  ncgtf: number of functions of each angular moment for each atom.
c  noneint:   I don't know, ask L. N. Vidal
c  ntwopoint: I don't know, ask L. N. Vidal
c  alpha: exponents of the gaussian basis functions.
c  an:    norm of the functions
c  cont: contraction coefficients for contracted basis sets
c  ncont: number of primitive gaussians for each group of primitives
c
      integer nang(maxatom)
      integer ncgtf(maxatom,maxl)
      integer ncont(maxatom,maxl,maxbt)
      integer ngcf(maxatom,maxl)
      integer noneint, ntwoint
      double precision alpha(maxatom,maxl,maxbt),
     *                 cont(maxatom,maxl,maxbase,maxbase)

      common / integrals / cont, alpha, 
     *                     ncgtf, ngcf, nang, ncont,
     *                     noneint, ntwoint


