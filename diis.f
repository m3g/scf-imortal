c
c Subroutine DIIS: Performs a DIIs aceleration as 
c                  described by Pulay.
c

      subroutine diis(kk,ndm,nit,ssrm,rvecs,fvecs,f,xx)

      implicit double precision (a-h,o-z)
      include 'sizes.i'
      include 'arrays.i'
      include 'scf.i'

      dimension rvecs(maxd*maxd,nvecs),fvecs(maxd,maxd,nvecs)
      dimension res(maxd,maxd)
      dimension temp(maxd,maxd)
      dimension weight(nvecs)
      dimension ssrm(maxd,maxd), ssrp(maxd,maxd), ss(maxd,maxd)     
      dimension d(maxd,maxd),xx(maxd,maxd)

      if(kk.ne.kdim) then
        write(*,*) ' kk.ne.kdim '
        stop
      end if

c Computing the residual vector R = F(D)DS -SDF(D) 

      do i = 1, kdim
        do j = 1, kdim
          d(i,j) = 0.d0
          do in = 1, ndm
            d(i,j) = d(i,j) + xx(i,in) * xx(j,in)
          end do
        end do
      end do

c Computing ssrp = Qs Ds^(1/2) QsT and ssrm = Qs Ds^(-1/2) QsT

      do i = 1, kdim
        do j = 1, kdim
          ssrp(i,j) = 0.d0
          ssrm(i,j) = 0.d0
          ss(i,j) = s(i,j)
        end do
      end do

      call diag(ss,temp,kdim)

      do i = 1, kdim
        sqrss = dsqrt(ss(i,i))
        ssrp(i,i) = sqrss
        ssrm(i,i) = 1.d0 / sqrss
      end do
 
      do i = 1, kdim
        do j = 1, kdim
          res(j,i) = temp(i,j)
        end do
      end do
      
      call pmatrix(kdim, kdim, kdim, ssrp, res, ssrp)
      call pmatrix(kdim, kdim, kdim, temp, ssrp, ssrp)

      call pmatrix(kdim, kdim, kdim, ssrm, res, ssrm)
      call pmatrix(kdim, kdim, kdim, temp, ssrm, ssrm)  

c Computing the residual vector
c R = S^(-1/2) F(X) D S^(1/2) - S^(1/2) D F(X) S^{-1/2)

      call pmatrix(kdim,kdim,kdim,d,ssrp,res)
      call pmatrix(kdim,kdim,kdim,f,res,res)
      call pmatrix(kdim,kdim,kdim,ssrm,res,res)

      call pmatrix(kdim,kdim,kdim,f,ssrm,temp)
      call pmatrix(kdim,kdim,kdim,d,temp,temp)
      call pmatrix(kdim,kdim,kdim,ssrp,temp,temp)

      do i = 1, kdim
        do j = 1, kdim
          res(i,j) = res(i,j) - temp(i,j)
        end do
      end do

c      call printm(res,kdim,kdim,'s')
c      stop
      
c Moving places at fvecs and at the residual vector array rvecs

      do ivec = nvecs, 2, -1
        do i = 1, kdim
          do j = 1, kdim
            ipos = j + (i-1)*kdim
            rvecs(ipos,ivec) = rvecs(ipos,ivec-1)
            fvecs(i,j,ivec) = fvecs(i,j,ivec-1)
          end do
        end do
      end do

c Updating the arrays fvecs and rvecs at the last position

      do i = 1, kdim
        do j = 1, kdim
          rvecs(j+(i-1)*kdim,1) = res(i,j)
          fvecs(i,j,1) = f(i,j)
        end do
      end do

c Computing the weights of each matrix to form the accelerated one

      if(nit.gt.nvecs) then
        
c Checking if DIIS will be used

        write(*,*) ' Accelerating... '
        call pesos(kdim*kdim,nvecs,rvecs,weight)

c Computing the accelerated fock matrix

        do ivec = 1, nvecs
          do i = 1, kdim
            do j = 1, kdim
              f(i,j) = weight(ivec)*fvecs(i,j,ivec)
            end do
          end do
        end do

      end if

      return
      end


c
c  Computes de weights of each fock matrix for DIIS extrapolation
c

      subroutine pesos (m, n, E, w)
c     *                    ata, save, eb, eu, b, v, diag, aux)
      implicit double precision (a-h, o-z)
c   Given the n-vectors that are columns of the matrix E = (e1,...,en),
c   where n > 1,
c   this subroutine computes the weights w = (w1,...,wn) such that
c   the squared norm of [w1 e1 + ... + wn en]
c   is minimal, subject to  w1 + ... + wn = 1.
c   mlin is the number of rows given to E in the calling program
c   eb is an auxiliary array of at least m*n positions
c   ata is an auxiliary array of at least n*n positions
c   save is an auxiliary matrix of at least n*n positions
c   eu is an auxiliary vector of at least m positions
c   diag is an auxiliary vector of at least n positions
c   aux is an auxiliary vector  of at least n positions
c   b is an auxiliary vector   of at least n positions
c   v is an auxiliary vector  of at least n positions

      include 'sizes.i'
      parameter(mlin=maxd*maxd)
      parameter(mm=maxd*maxd)
      parameter(nn=nvecs)

      dimension E(mlin, nn), w(nn), eb(mm, nn), eu(mm), b(nn)
      dimension ata(nn, nn), v(nn), save(nn, nn), diag(nn), aux(nn)

      if(n.lt.2) then
        write(*, *) ' Subroutine Pesos requires n > 1, m > 0'
        stop
      endif
      iur = 0
      n1 = n - 1
      uene = 1.d0/dfloat(n)
c
c   Compute the matrix EB
      do j = 1, n1
        do i = 1, m
          eb(i, j) = E(i, j) - E(i, j+1)
        end do
      end do
c
c   Compute the vector eu
      do i = 1, m
        eu(i) = 0.d0
        do j = 1, n
          eu(i) = eu(i) + e(i, j)
        end do
        eu(i) = eu(i)/dfloat(n)
      end do
c
c  Compute b = (EB)^T eu
      do i = 1, n1
        b(i) = 0.d0
        do j = 1, m
          b(i) = b(i) - eb(j, i) * eu(j)
        end do
      end do

c
c
      do i = 1, n1
        do j = 1, n1
        ata(i, j) = 0.d0
          do k = 1, m
            ata(i, j) = ata(i, j) + eb(k, i) * eb(k, j)
          end do
        end do
      end do

10    do j = 1, n1
        do i = 1, n1
          save(i, j) = ata(i, j)
        end do
      end do

      call chole (n1, save, ier, diag, n)


      if(ier.ne.0) then
        if(iur.eq.0) then
          z = 1.d-8
        else
          z = 10.d0 * z
        endif
        do i = 1, n-1
          z = dmax1(z, dabs(ata(i, i)))
        end do
        do i = 1, n-1
          ata(i, i) = ata(i, i) + 1.d-6 * z
        end do
        iur = 1
        go to 10
      endif

      call sicho(n1, save, v, b, aux, n)

      w(1) = uene + v(1)
      w(n) = uene - v(n1)
      do i = 2, n-1
        w(i) = uene - v(i-1) + v(i)
      end do

      sw = 0.d0
      do j = 1, n
        sw = sw + w(j)
      end do
      do j = 1, n
        w(j) = w(j)/sw
      end do 

      return
      end

	subroutine chole (n, a, ier, diag, nlin)
	implicit double precision (a-h,o-z)
	dimension a(nlin, n), diag(n)	
	ier = 0
c Save the diagonal of the matrix
	do i = 1, n
	diag(i) = a(i, i)
	end do
c Test nonnegativity of diagonal
	do i = 1, n
	if(diag(i).le.0) then
	ier = 1
	return
	endif
	end do

	a(1,1) = dsqrt(a(1, 1))
	if(n.eq.1)return

	do 1 i = 2, n
	do 2 j = 1 ,i-1
	z = 0.d0
	if(j.gt.1)then
	do 3 k=1,j-1
3	z = z + a(i,k) * a(j,k)
	endif
	a(i,j) = (a(i,j) - z)/a(j,j)
2	continue
	z = 0.d0
	do 4 j=1,i-1
4	z = z + a(i,j)**2
	temp = a(i, i) - z

c   Test positive definiteness
	if(temp.le.0.d0) then
	ier = i
c   Restore the diagonal
	do ii = 1, n
	a(ii, ii) = diag(ii)
	end do
c   Restore lower triangular part
	do ii = 2, n
	do j = 1, ii-1
	a(ii, j) = a(j, ii)
	end do
	end do
	return
	endif

	a(i,i) = dsqrt(temp)
1	continue
	return
	end


	subroutine sicho (n, a, x, b, aux, nlin)
	implicit double precision (a-h,o-z)
	dimension a(nlin, n),x(n),b(n),aux(n)
	aux(1) = b(1)/a(1,1)
	if(n.gt.1)then
	do 1 i=2,n
	z = 0.d0
	do 2 j=1,i-1
2	z = z + a(i,j)*aux(j)
	aux(i) = (b(i) - z) / a(i,i)
1	continue
	endif
	x(n) = aux(n)/a(n,n)
	if(n.eq.1)return
	do 3 i=n-1,1,-1
	z = 0.d0
	do 4 j=i+1,n
4	z = z + a(j,i)*x(j)
	x(i) = (aux(i) - z)/a(i,i)
3	continue
	return
	end

