c Subroutine printm
c
c	Imprime matrizes de n x m

	subroutine printm(a,n,m,p)

      include 'sizes.i'

      integer n, i, j
	double precision a(maxd,maxd)
      character*1 p

      if(p.eq.'s') then
        jl = 0
        nseq = (m - mod(m,6)) / 6
        do ii = 1, nseq 
          jf = jl + 1  
          jl = jf + 5
          write(*,2) (j, j = jf, jl)
          do i = 1, n
            write(*,3) i, (a(i,j), j = jf, jl)
          end do
        end do
        jf = jl + 1
        jl = jf + mod(m,6) - 1
        write(*,2) (j, j = jf, jl)
        do i = 1, n
          write(*,3) i, (a(i,j), j = jf, jl)
        end do
      end if

      if(p.eq.'d') then
        jl = 0
        nseq = (m - mod(m,4)) / 4
        do ii = 1, nseq 
          jf = jl + 1  
          jl = jf + 3
          write(*,4) (j, j = jf, jl)
          do i = 1, n
            write(*,5) i, (a(i,j), j = jf, jl)
          end do
        end do
        jf = jl + 1
        jl = jf + mod(m,4) - 1
        write(*,4) (j, j = jf, jl)
        do i = 1, n
          write(*,5) i, (a(i,j), j = jf, jl)
        end do
      end if

  2   format(tr6,100(tr2,i10))
  3	format(i5,tr1,100(tr2,f10.4))
  4   format(tr6,100(tr1,i17))
  5   format(i5,tr1,100(tr1,D17.10))	

      write(*,*)
	return
	end
 

c Subroutine guess: creates an initial point
c
c
      subroutine guess(ndim,k,c,iguess)
      
      include 'sizes.i'
      include 'files.i'
      
      integer i, j, k, ndim, iguess
      double precision c(maxd,maxd)
      integer hms(3)

c zero density matrix

      if(iguess.eq.0) then
        do i = 1, k
          do j = 1, k
            c(i,j) = 0.d0
          end do 
        end do
      end if

c Guessing a feasible initial point 

      if(iguess.eq.2) then
        call itime(hms)
        isem = hms(1)*hms(2)*hms(3)
        write(*,*) ' Random seed: ', isem
        do i = 1, k
          do j = i, k
            c(i,j) = rnd(isem)
          end do
        end do
        call gs(c,k,k,maxd)
      end if

c Getting guess from file

      if(iguess.eq.3) then
        open(10, file = densfile, status = 'old')
        do i = 1, k
          read(10,40) (c(i,j), j = 1, ndim)
        end do
        close(10)
      end if
   40 format(300(tr2,d17.10))

      return
      end 


c Subroutine diag: chama a subrotina jacobi para diagonalizar as matrizes
c de n x n deste problema.
c	
c     input:
c	d matriz que vai ser diagonalizada
c     k dimensao da matriz (kdim x kdim)
c
c	return:
c	d matriz (diagonal) de autovalores de d(input)
c	v matriz de autovetores (ortonormais) de d(input)
c

	subroutine diag(d,v,kdim)	

	implicit double precision (a-h,o-z)

      include 'sizes.i'
	dimension d(maxd,maxd), v(maxd,maxd)
      dimension dl(maxd,maxd), vl(maxd,maxd) 

	call jacobi(d,v,kdim)

      do i = 1, kdim
        do j = 1, kdim
          if(i.ne.j) d(i,j) = 0.d0
        end do
      end do

	return
	end

c   Subroutine pmatrix: makes matrix products a . b = c
c 
c	input:
c	n number of lines of matrix a 
c	m number of columns of matrix a (equal to the number of lines of matrix b)  
c	l number of columns of matrix b
c	a matrix a(n,m)
c	b matrix b(m,l)
c
c	return:
c	c matrix c(n,l) (product of a and b)
c

      subroutine pmatrix(n,m,l,a,b,c)

	implicit double precision (a-h,o-z)
      include 'sizes.i'

	dimension a(maxd,maxd), b(maxd,maxd), c(maxd,maxd)
      dimension al(maxd,maxd), bl(maxd,maxd), cl(maxd,maxd)

      do j = 1, m
        do i = 1, n
          al(i,j) = a(i,j)
        end do
        do i = 1, l
          bl(j,i) = b(j,i)
        end do
      end do 	
     
	do ila = 1, n
        do icb = 1, l
          cl(ila,icb) = 0.d0
          do i = 1, m
            cl(ila, icb) = cl(ila,icb) + a(ila,i) * b(i,icb)
          end do
        end do
	end do

      do i = 1, n
        do j = 1, l
          c(i,j) = cl(i,j)
        end do
      end do


      return
	end
 
c
c   Subroutine chari
c   Given an integer i less than 9999, returns the corresponding number as characters.

        subroutine chari(i, char)

        character*4 char
        character digit(0:9)
        integer i, d3, d2, d1, d0
        data digit /'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/

        char = '    '

        d3 = aint(dfloat(i/1000))
        d2 = aint(dfloat(i-d3*1000)/100.d0)
        d1 = aint(dfloat(i-d3*1000-d2*100)/10.d0)
        d0 = aint(dfloat(i-d3*1000-d2*100-d1*10))

        if(i.ge.1000)
     *         char = digit(d3)//digit(d2)//digit(d1)//digit(d0)
        if(i.lt.1000.and.i.ge.100)
     *         char = digit(d2)//digit(d1)//digit(d0)
        if(i.lt.100.and.i.ge.10)
     *         char = digit(d1)//digit(d0)
        if(i.lt.10)
     *         char = digit(d0)

        return
        end
c
c  Subroutine Graham-Schmidt: Performs vector ortonormalization
c     using the naive Graham-Schmidt algorithm.
c 
c  Input:    a    square matrix (maybe overdimensioned) containing the 
c                 vectors as columns 1 to n
c            m    number of elements of each vector
c            n    number of vectors
c            maxd dimension of the matrix a
c
c  Return:   a    matrix containing the ortonormal vector in columns
c                 1 to n (with m elements each vector).
c
      subroutine gs(a,m,n,maxd)

      implicit double precision (a-h,o-z)

      dimension a(maxd,maxd)

      do j = 1, n
        anorm = 0.d0
        do i = 1, m
          anorm = anorm + a(i,j)**2
        end do
        anorm = dsqrt(anorm)
        do i = 1, m
          a(i,j) = a(i,j) / anorm
        end do
      end do

      do j1 = 1, n - 1
        do j2 = j1 + 1, n
          p12 = 0.d0
          do i = 1, m
            p12 = p12 + a(i,j1)*a(i,j2)
          end do
          anorm = 0.d0
          do i = 1, m          
            a(i,j2) = a(i,j2) - p12 * a(i,j1)
            anorm = anorm + a(i,j2)**2
          end do
          anorm = dsqrt(anorm)
          do i = 1, m
            a(i,j2) = a(i,j2) / anorm
          end do
        end do
      end do

      return
      end

c
c  Subroutine check: checks if a point is feasible
c

      subroutine check(x,s,kdim,ndm)

      implicit double precision (a-h,o-z)
      include 'sizes.i'

      dimension x(maxd,maxd), s(maxd,maxd), 
     *          xt(maxd,maxd), sx(maxd,maxd)
      logical feasible

      feasible = .true.

      tol = 1.d-8

      call pmatrix(kdim,kdim,ndm,s,x,sx)

      do i = 1, kdim
        do j = 1, ndm
          xt(j,i) = x(i,j)
        end do
      end do

      call pmatrix(ndm,kdim,ndm,xt,sx,xt)

      do i = 1, ndm
        do j = 1, ndm
          if(i.eq.j) then
            test = abs(xt(i,j) - 1.d0)
          else
            test = abs(xt(i,j))
          end if  
          if(test.gt.tol) then
            feasible = .false.
          end if           
        end do
      end do

      if(.not.feasible) then
        write(*,10)
10      format(/,' POINT NOT FEASIBLE: ',/)
        call printm(xt,ndm,ndm,'d')
      else
        write(*,30)
30      format(/,'  FEASIBILITY: OK ',/)
      end if

      return
      end
 
c
c The following subroutines were not writen by me.
c

       
      SUBROUTINE JACOBI(A,V,N) 
C*********************************************************************** 
C 
C     DIAGONALISATION OF REAL SYMMETRIC MATICES BY JACOBI METHOD 
C 
C*********************************************************************** 
      IMPLICIT double precision(A-H,O-Z) 
      include 'sizes.i'       
      DIMENSION A(maxd,maxd),V(maxd,maxd)
      RHO=1.0E-10 
      TES=0.0 
      SCL=0.0 
      DO 10 I=1,N 
   10 SCL=SCL+A(I,I)**2 
      SCL=DSQRT(SCL)/DFLOAT(N) 
      DO 20 I=1,N 
      DO 20 J=1,N 
   20 A(I,J)=A(I,J)/SCL 
      DO 30 I=1,N 
      DO 30 J=1,N 
      V(I,J)=0.0 
      IF(I.EQ.J)V(I,J)=1.0 
   30 CONTINUE 
      DO 100 I=2,N 
      DO 100 J=1,I-1 
  100 TES=TES+2.0*A(I,J)*A(I,J) 
      TES=DSQRT(TES) 
      M=0 
  105 TES=TES/DFLOAT(N) 
      IF(TES.LT.RHO)TES=RHO 
  110 DO 165 I=2,N 
      DO 165 J=1,I-1 
      IF(DABS(A(I,J))-TES)165,115,115 
  115 M=1 
      V1=A(J,J) 
      V2=A(I,J) 
      V3=A(I,I) 
      U=0.5*(V1-V3) 
      IF(DABS(U)-RHO)120,125,125 
  120 OMG=-1.0 
      GO TO 130 
  125 OMG=-V2/DSQRT(V2*V2+U*U) 
      IF(U.LT.0.0)OMG=-OMG 
  130 S=OMG/DSQRT(2.0*(1.0+DSQRT(1.0-OMG*OMG))) 
      C=DSQRT(1.0-S*S) 
      DO 160 K=1,N 
      IF(K-I)140,135,135 
  135 TEM=A(K,J)*C-A(K,I)*S 
      A(K,I)=A(K,J)*S+A(K,I)*C 
      A(K,J)=TEM 
      GO TO 155 
  140 IF(K-J)145,150,150 
  145 TEM=A(J,K)*C-A(I,K)*S 
      A(I,K)=A(J,K)*S+A(I,K)*C 
      A(J,K)=TEM 
      GO TO 155 
  150 TEM=A(K,J)*C-A(I,K)*S 
      A(I,K)=A(K,J)*S+A(I,K)*C 
      A(K,J)=TEM 
  155 TEM=V(K,J)*C-V(K,I)*S 
      V(K,I)=V(K,J)*S+V(K,I)*C 
      V(K,J)=TEM 
  160 CONTINUE 
      A(J,J)=V1*C*C+V3*S*S-2.0*V2*S*C 
      A(I,I)=V1*S*S+V3*C*C+2.0*V2*S*C 
      A(I,J)=(V1-V3)*S*C+V2*(C*C-S*S) 
  165 CONTINUE 
      IF(M-1)175,170,170 
  170 M=0 
      GO TO 110 
  175 IF(TES-RHO)180,180,105 
  180 DO 190 I=1,N 
      DO 190 J=1,N 
  190 A(I,J)=SCL*A(I,J) 
      RETURN 
      END 

C---------------------------------------------------------------------- 
C 
	real*8 function rnd(sem) 
C 
	integer sem, mult 
C 
	sem = mod( mult( sem, 3141581) + 1, 100000000) 
	rnd = sem/100000000.0d0 
	return 
	end 
C 
C---------------------------------------------------------------------- 
 

c 
C---------------------------------------------------------------------- 
C 
      	integer function mult( p, q) 
C 
	Integer p, q, p0, p1, q0, q1 
C 
	p1 = p/10000 
	p0 = mod(p,10000) 
	q1 = q/10000 
	q0 = mod(q,10000) 
	mult = mod( mod( p0*q1+p1*q0,10000)*10000+p0*q0,100000000) 
	return 
	end 
C 
C---------------------------------------------------------------------- 
C 

