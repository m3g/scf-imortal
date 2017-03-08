
      SUBROUTINE integral
C*********************************************************
C Compute all required molecular integrals for Cartesian               
C Gaussian Type Functions.
c  
c This subroutine was originaly writen by Luciano Nassif Vidal
c and merely adapted to the actual package by Leandro Martínez.
c
C*********************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

      include 'sizes.i'
      include 'arrays.i'
      include 'integrals.i'
      include 'pi.i'

      double precision an(maxatom,maxl,maxl,maxl,maxbt)
      common / norm / an

c      real*8 ABCDtp(maxd,maxd,maxd,maxd), Htp(maxd,maxd), Stp(maxd,maxd)
      real*8 Htp(maxd,maxd), Stp(maxd,maxd)
    
      LOGICAL WTODO
      logical uncont, efield

c ?   COMMON /HF/ P(maxbase,maxbase),RNORM(maxbase),EHF,NSCF

      efield = .false.

c Abrindo o arquivo no qual as integrais sao escritas:

      open(20,file="intvalues.dat",access="direct",
     *        form="unformatted",recl=17) 


c      x(natoms+1) = 0.d0
c      y(natoms+1) = 0.d0
c      z(natoms+1) = 0.d0


C
C Norm of CGTF's
C
      DO I = 1,natoms
         DO J = 1,NANG(I)
            DO K = 1,NCGTF(I,J)

               DO LZ = 0,J-1
               DO LY = 0,J-1
               DO LX = 0,J-1
                  LMN = LX + LY + LZ
                  IF ( LMN .EQ. J-1 ) THEN

      e1   = DFLOAT(2*LMN) + 1.5D0
      e2   = DFLOAT(  LMN) + 1.5D0
      FACX = FACFAC(2*LX-1)
      FACY = FACFAC(2*LY-1)
      FACZ = FACFAC(2*LZ-1)

      AN(I,LX,LY,LZ,K) = ( 2.D0**e1 ) * ( ALPHA(I,J,K)**e2 )
      AN(I,LX,LY,LZ,K) = AN(I,LX,LY,LZ,K) /
     !                   ( FACX * FACY * FACZ * PI**1.5D0 )
      AN(I,LX,LY,LZ,K) = DSQRT( AN(I,LX,LY,LZ,K) )

                  END IF
               END DO
               END DO
               END DO
            END DO
         END DO
      END DO
C
C
C One Electron Integrals
C
      NONEINT = 0
      NTWOINT = 0
      I = 0
      J = 0
      K = 0 
      L = 0
      DO I1 = 1,natoms
         DO J1 = 1,NANG(I1)
            DO K1 = 1,NCGTF(I1,J1)

      DO LZ1 = 0,J1-1
      DO LY1 = 0,J1-1
      DO LX1 = 0,J1-1
         LMN1 = LX1 + LY1 + LZ1
         IF ( LMN1 .EQ. J1-1 ) THEN
              I = I + 1

      DO I2 = 1,natoms
         DO J2 = 1,NANG(I2)
            DO K2 = 1,NCGTF(I2,J2)

      DO LZ2 = 0,J2-1
      DO LY2 = 0,J2-1
      DO LX2 = 0,J2-1
         LMN2 = LX2 + LY2 + LZ2
         IF ( LMN2 .EQ. J2-1 ) THEN
              J  = J + 1

      A1 = ALPHA(I1,J1,K1)
      A2 = ALPHA(I2,J2,K2)

      IF ( J.GE.I ) THEN
C
C Overlap Integrals  
C
      S(I,J) = SAB(A1,A2,I1,J1,K1,I2,J2,K2,LX1,LY1,LZ1,LX2,LY2,LZ2,'T')
      NONEINT = NONEINT + 1

C
C Kinetic Energy Integrals
C
      Hcore(I,J) = T(A1,A2,I1,J1,K1,I2,J2,K2,LX1,LY1,LZ1,LX2,LY2,LZ2)
      NONEINT = NONEINT + 1

C
C Finite Electric Field Integrals
C
c      IF ( EFIELD ) THEN 
c
c           Hcore(I,J) = Hcore(I,J) + 
c     !DFIELD(1) * 
c     !DIPINT(LX1,LY1,LZ1,I1,J1,K1,A1,LX2,LY2,LZ2,I2,J2,K2,A2,1,0,0) +
c     !DFIELD(2) *
c     !DIPINT(LX1,LY1,LZ1,I1,J1,K1,A1,LX2,LY2,LZ2,I2,J2,K2,A2,0,1,0) +
c     !DFIELD(3) *
c     !DIPINT(LX1,LY1,LZ1,I1,J1,K1,A1,LX2,LY2,LZ2,I2,J2,K2,A2,0,0,1)
c
c      END IF
c
C
C Potential Energy Integrals
C
      DO I3 = 1,natoms
         Hcore(I,J) = Hcore(I,J) + V(A1,A2,I1,J1,K1,I2,J2,K2,
     !                               LX1,LY1,LZ1,LX2,LY2,LZ2,I3)
         NONEINT = NONEINT + 1
      END DO

      ELSE

      S(I,J) = S(J,I)
      Hcore(I,J) = Hcore(J,I)

      END IF

      Stp(I,J) = S(I,J)
      Htp(I,J) = Hcore(I,J)

C
C Two Electrons Integrals
C
      DO I3 = 1,natoms
         DO J3 = 1,NANG(I3)
            DO K3 = 1,NCGTF(I3,J3)

      DO LZ3 = 0,J3-1
      DO LY3 = 0,J3-1
      DO LX3 = 0,J3-1
         LMN3 = LX3 + LY3 + LZ3
         IF ( LMN3 .EQ. J3-1 ) THEN
              K = K + 1

      DO I4 = 1,natoms
         DO J4 = 1,NANG(I4)
            DO K4 = 1,NCGTF(I4,J4)

      DO LZ4 = 0,J4-1
      DO LY4 = 0,J4-1
      DO LX4 = 0,J4-1
         LMN4 = LX4 + LY4 + LZ4
         IF ( LMN4 .EQ. J4-1 ) THEN
              L  = L + 1

      CALL JUMP(I,J,K,L,NI,NJ,NK,NL,WTODO)
      IF ( WTODO ) THEN

      A3 = ALPHA(I3,J3,K3)
      A4 = ALPHA(I4,J4,K4)

c Adicionado para escrever as integrais em arquivo:


      ABCDtemp = Two(A1,A2,A3,A4,I1,J1,K1,I2,J2,K2,I3,J3,K3,
     !I4,J4,K4,LX1,LY1,LZ1,LX2,LY2,LZ2,LX3,LY3,LZ3,LX4,LY4,LZ4)

      index = maxd*maxd*maxd*(l-1) + maxd*maxd*(k-1) + maxd*(j-1) + i
      write(20,rec=index) ABCDtemp


c      ABCD(I,J,K,L) = Two(A1,A2,A3,A4,I1,J1,K1,I2,J2,K2,I3,J3,K3,
c     !I4,J4,K4,LX1,LY1,LZ1,LX2,LY2,LZ2,LX3,LY3,LZ3,LX4,LY4,LZ4)
c      ABCDtp(I,J,K,L) = ABCD(I,J,K,L) voltar
c      NTWOINT = NTWOINT + 1

c Fim da mudanca, 17/11/03

      END IF

         END IF
       END DO
       END DO
       END DO
             END DO
          END DO
       END DO

         END IF
       L = 0
       END DO
       END DO
       END DO
             END DO
          END DO
       END DO

         END IF
       K = 0
       END DO
       END DO
       END DO
             END DO
          END DO
       END DO

         END IF
       J = 0
       END DO
       END DO
       END DO
             END DO
          END DO
       END DO

      close(20)
      return

C
C Transformation of One and Two Electrons Integrals from
C uncontrated basis set to contrated basis set.
C
c      IF ( UNCONT ) RETURN
c      DO I = 1,NDIM
c         DO J = 1,NDIM
c            Stp(I,J) = 0.D0
c            Htp(I,J) = 0.D0
c            DO K = 1,NDIM
c               DO L = 1,NDIM
c                  ABCDtp(I,J,K,L) = 0.D0
c               END DO
c            END DO
c         END DO
c      END DO
c
c      I  = 0 
c      J  = 0 
c      K  = 0 
c      L  = 0
c      IC = 0 
c      JC = 0 
c      KC = 0 
c      LC = 0
c      IREF  = 1 
c      JREF  = 1 
c      KREF  = 1 
c      LREF  = 1
c      ICREF = 1 
c      JCREF = 1 
c      KCREF = 1 
c      LCREF = 1
c      DO I1 = 1,natoms
c         DO J1 = 1,NANG(I1)
c            DO L1 = 1,NCONT(I1,J1)
c
c               IF (L1.EQ.1) THEN
c                   IREF = I
c               ELSE
c                   I = IREF
c               END IF
c
c            DO K1 = 1,NCGTF(I1,J1)
c 
c               IF (K1.EQ.1) THEN
c                   ICREF = IC
c               ELSE
c                   IC = ICREF
c               END IF
c
c      DO LZ1 = 0,J1-1
c      DO LY1 = 0,J1-1
c      DO LX1 = 0,J1-1
c         LMN1 = LX1 + LY1 + LZ1
c         IF ( LMN1 .EQ. J1-1 ) THEN
c              I  = I  + 1 
c              IC = IC + 1 
c               
c      DO I2 = 1,natoms
c         DO J2 = 1,NANG(I2)
c            DO L2 = 1,NCONT(I2,J2)
c
c               IF (L2.EQ.1) THEN
c                   JREF = J
c               ELSE
c                   J = JREF
c               END IF
c
c            DO K2 = 1,NCGTF(I2,J2)
c 
c               IF (K2.EQ.1) THEN
c                   JCREF = JC
c               ELSE
c                   JC = JCREF
c               END IF
c
c      DO LZ2 = 0,J2-1
c      DO LY2 = 0,J2-1
c      DO LX2 = 0,J2-1
c         LMN2 = LX2 + LY2 + LZ2
c         IF ( LMN2 .EQ. J2-1 ) THEN
c              J   = J + 1
c              JC  = JC + 1
c
c      IF ( JC.GE.IC ) THEN
c
c      Stp(IC,JC) = Stp(IC,JC) + 
c     !             S(I,J)     * CONT(I1,J1,K1,L1) * CONT(I2,J2,K2,L2)
c      Htp(IC,JC) = Htp(IC,JC) + 
c     !             Hcore(I,J) * CONT(I1,J1,K1,L1) * CONT(I2,J2,K2,L2)
c
c      ELSE
c
c         Stp(IC,JC) = Stp(JC,IC)
c         Htp(IC,JC) = Htp(JC,IC)
c
c      END IF
c
c      DO I3 = 1,natoms
c         DO J3 = 1,NANG(I3)
c            DO L3 = 1,NCONT(I3,J3)
c
c               IF (L3.EQ.1) THEN
c                   KREF = K
c               ELSE
c                   K = KREF
c               END IF
c
c            DO K3 = 1,NCGTF(I3,J3)
c         
c               IF (K3.EQ.1) THEN
c                   KCREF = KC
c               ELSE
c                   KC = KCREF
c               END IF
c
c      DO LZ3 = 0,J3-1
c      DO LY3 = 0,J3-1
c      DO LX3 = 0,J3-1
c         LMN3 = LX3 + LY3 + LZ3
c         IF ( LMN3 .EQ. J3-1 ) THEN
c              K = K + 1
c              KC = KC + 1
c
c      DO I4 = 1,natoms
c         DO J4 = 1,NANG(I4)
c            DO L4 = 1,NCONT(I4,J4)
c
c               IF (L4.EQ.1) THEN
c                   LREF = L
c               ELSE
c                   L = LREF
c               END IF
c
c            DO K4 = 1,NCGTF(I4,J4)
c
c               IF (K4.EQ.1) THEN
c                   LCREF = LC
c               ELSE
c                   LC = LCREF
c               END IF
c
c      DO LZ4 = 0,J4-1
c      DO LY4 = 0,J4-1
c      DO LX4 = 0,J4-1
c         LMN4 = LX4 + LY4 + LZ4
c         IF ( LMN4 .EQ. J4-1 ) THEN
c              L  = L + 1
c              LC  = LC + 1
c
c      CALL JUMP(I,J,K,L,NI,NJ,NK,NL,WTODO)
c      CALL JUMP(IC,JC,KC,LC,NIC,NJC,NKC,NLC,WTODO)
c
c      IF (WTODO) THEN
c
c      ABCDtp(NIC,NJC,NKC,NLC) = ABCDtp(NIC,NJC,NKC,NLC) + 
c     !                          ABCD(NI,NJ,NK,NL) 
c     !                  * CONT(I1,J1,K1,L1) * CONT(I2,J2,K2,L2)
c     !                  * CONT(I3,J3,K3,L3) * CONT(I4,J4,K4,L4)
c
c      END IF
c
c         END IF
c       END DO
c       END DO
c       END DO
c             END DO  
c             END DO  
c          END DO
c       END DO
c
c         END IF
c       L = 0
c       LC = 0
c       END DO
c       END DO
c       END DO
c             END DO
c             END DO  
c          END DO
c       END DO
c
c         END IF
c       K = 0
c       KC = 0
c       END DO
c       END DO
c       END DO
c             END DO  
c             END DO  
c          END DO
c       END DO
c
c         END IF
c       J = 0
c       JC = 0
c       END DO
c       END DO
c       END DO
c             END DO
c             END DO  
c          END DO
c       END DO
c
cC
cC Renormalization of Integrals.
cC
c      DO I = 1,NDIM
c         RNORM(I) = DSQRT(DABS(Stp(I,I)))
c      END DO
c      DO I = 1,NDIM
c         DO J = 1,NDIM
c            Stp(I,J) = Stp(I,J) / RNORM(I) / RNORM(J)
c            Htp(I,J) = Htp(I,J) / RNORM(I) / RNORM(J)
c            DO K = 1,NDIM
c               DO L = 1,NDIM
c                  CALL JUMP(I,J,K,L,NI,NJ,NK,NL,WTODO)
c                  IF (WTODO) THEN
c                  ABCDtp(I,J,K,L) = ABCDtp(I,J,K,L) / RNORM(I) / 
c     !                              RNORM(J) / RNORM(K) / RNORM(L)
c                  END IF
c               END DO
c            END DO
c         END DO
c      END DO
c
c      RETURN
      END  

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Two Electrons Integrals
C
      DOUBLE PRECISION FUNCTION Two(A1,A2,A3,A4,I1,J1,K1,I2,J2,K2,
     !I3,J3,K3,I4,J4,K4,LX1,LY1,LZ1,LX2,LY2,LZ2,LX3,LY3,LZ3,LX4,LY4,LZ4)
      IMPLICIT REAL*8 (A-H,O-Z)  

      include 'sizes.i'
      include 'pi.i'
      include 'integrals.i'

      double precision an(maxatom,maxl,maxl,maxl,maxbt)
      common / norm / an

      Two = 0.D0

      P1 = AN(I1,LX1,LY1,LZ1,K1) * AN(I2,LX2,LY2,LZ2,K2) * 
     !     AN(I3,LX3,LY3,LZ3,K3) * AN(I4,LX4,LY4,LZ4,K4)

      P2 = 2.D0 * PI**2.5D0
      P2 = P2 / ( (A1+A2) * (A3+A4) * DSQRT(A1+A2+A3+A4) )

      GAMMA1 = (A1*A2) * AB2(I1,I2) / (A1+A2)
      GAMMA2 = (A3*A4) * AB2(I3,I4) / (A3+A4)
      P3 = DEXP(-(GAMMA1 + GAMMA2))

      PAx = PAi(I1,J1,K1,I2,J2,K2,I1,'X')
      PBx = PAi(I1,J1,K1,I2,J2,K2,I2,'X')
      PAy = PAi(I1,J1,K1,I2,J2,K2,I1,'Y')
      PBy = PAi(I1,J1,K1,I2,J2,K2,I2,'Y')
      PAz = PAi(I1,J1,K1,I2,J2,K2,I1,'Z')
      PBz = PAi(I1,J1,K1,I2,J2,K2,I2,'Z')

      QCx = PAi(I3,J3,K3,I4,J4,K4,I3,'X')
      QDx = PAi(I3,J3,K3,I4,J4,K4,I4,'X')
      QCy = PAi(I3,J3,K3,I4,J4,K4,I3,'Y')
      QDy = PAi(I3,J3,K3,I4,J4,K4,I4,'Y')
      QCz = PAi(I3,J3,K3,I4,J4,K4,I3,'Z')
      QDz = PAi(I3,J3,K3,I4,J4,K4,I4,'Z') 

      P4 = SIGMAB(A1,I1,J1,K1,LX1,LY1,LZ1,PAx,PAy,PAz,
     !            A2,I2,J2,K2,LX2,LY2,LZ2,PBx,PBy,PBz,
     !            A3,I3,J3,K3,LX3,LY3,LZ3,QCx,QCy,QCz,
     !            A4,I4,J4,K4,LX4,LY4,LZ4,QDx,QDy,QDz)

      Two = P1 * P2 * P3 * P4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Compute Function SIGMA B
C
      DOUBLE PRECISION FUNCTION SIGMAB (
     !                         A1,I1,J1,K1,LX1,LY1,LZ1,PAx,PAy,PAz,
     !                         A2,I2,J2,K2,LX2,LY2,LZ2,PBx,PBy,PBz,
     !                         A3,I3,J3,K3,LX3,LY3,LZ3,QCx,QCy,QCz,
     !                         A4,I4,J4,K4,LX4,LY4,LZ4,QDx,QDy,QDz
     !                                 )
      IMPLICIT REAL*8 (A-H,O-Z)

      include 'sizes.i'
      include 'arrays.i'

      DIMENSION PQ(3)

      SIGMAB = 0.D0

      PQ2 = 0.D0
      DO I = 1,3
      if(I.eq.1) then
         C1 = x(I1)
         C2 = x(I2)
         C3 = x(I3)
         C4 = x(I4)
      else if(I.eq.2) then
         C1 = y(I1)
         C2 = y(I2)
         C3 = y(I3)
         C4 = y(I4)
      else if(I.eq.3) then
         C1 = z(I1)
         C2 = z(I2)
         C3 = z(I3)
         C4 = z(I4)
      end if
      PQ(I) = ( -(C1*A1 + C2*A2)*(A3 + A4)+(C3*A3 + C4*A4)*(A1 + A2) ) / 
     !        ( (A1 + A2) * (A3 + A4) )
         PQ2 = PQ2 + PQ(I)**2
       END DO

      tau = PQ2 * (A1 + A2) * (A3 + A4) / (A1 + A2 + A3 + A4) 

      DO IAX = 0,LX1+LX2
         DO IBX = 0,LX3+LX4
            DO ICX = 0,INT(DFLOAT(IAX)/2.D0)
               DO IDX = 0,INT(DFLOAT(IBX)/2.D0)
                  DO IEX = 0,INT(DFLOAT(IAX + IBX -2*(ICX + IDX) )/2.D0)

      DO IAY = 0,LY1+LY2
         DO IBY = 0,LY3+LY4
            DO ICY = 0,INT(DFLOAT(IAY)/2.D0)
               DO IDY = 0,INT(DFLOAT(IBY)/2.D0)
                  DO IEY = 0,INT(DFLOAT(IAY + IBY -2*(ICY + IDY) )/2.D0)
                 
      DO IAZ = 0,LZ1+LZ2
         DO IBZ = 0,LZ3+LZ4
            DO ICZ = 0,INT(DFLOAT(IAZ)/2.D0)
               DO IDZ = 0,INT(DFLOAT(IBZ)/2.D0)
                  DO IEZ = 0,INT(DFLOAT(IAZ + IBZ -2*(ICZ + IDZ) )/2.D0)

      nu  = IAX +  IBX + IAY + IBY + IAZ + IBZ 
     !      -2  * (ICX + IDX + ICY + IDY + ICZ + IDZ) 
     !           -(IEX + IEY + IEZ)

      SIGMAB = SIGMAB + Fn(nu,tau) *
     !         Bx(IAX,IBX,ICX,IDX,IEX,
     !            LX1,LX2,LX3,LX4,PAx,PBx,QCx,QDx,A1,A2,A3,A4,PQ(1)) *
     !         Bx(IAY,IBY,ICY,IDY,IEY,
     !            LY1,LY2,LY3,LY4,PAy,PBy,QCy,QDy,A1,A2,A3,A4,PQ(2)) *
     !         Bx(IAZ,IBZ,ICZ,IDZ,IEZ,
     !            LZ1,LZ2,LZ3,LZ4,PAz,PBz,QCz,QDz,A1,A2,A3,A4,PQ(3))

                  END DO
               END DO
            END DO
         END DO
      END DO


                  END DO
               END DO
            END DO
         END DO
      END DO


                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Compute Function Bx
C
      DOUBLE PRECISION FUNCTION 
     !              Bx(IA,IB,IC,ID,IE,
     !                 L1,L2,L3,L4,PA,PB,QC,QD,A1,A2,A3,A4,PQ) 
      IMPLICIT REAL*8 (A-H,O-Z)

      Bx = 0.D0

      del = .25D0 * ( A1 + A2 + A3 + A4 ) / ( (A1 + A2) * (A3 + A4) )

      P1 = DFLOAT((-1)**(IB+IE)) * fj(IA,L1,L2,PA,PB)
     !                           * fj(IB,L3,L4,QC,QD)

      P2 = FACT(IA) * FACT(IB) / ( ( 4.D0 * (A1 + A2) )**IA ) 
     !                         / ( ( 4.D0 * (A3 + A4) )**IB )
     !   / ( del**(IA + IB -2*(IC + ID) -IE ) )

      P3 = ( ( 4.D0 * (A1 + A2) )**IC ) *
     !     ( ( 4.D0 * (A3 + A4) )**ID )
     !     / FACT(IC) / FACT(ID) / FACT(IA-2*IC) 
     !     / FACT(IB-2*ID) 

      P4 = FACT(IA + IB -2*(IC + ID))

      P5 = ( PQ**(IA+IB-2*(IC+ID+IE)) )
     !     / FACT(IE) / FACT(IA+IB-2*(IC+ID+IE))

      Bx = P1 * P2 * P3 * P4 * P5

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Potential Energy Integrals
C
      DOUBLE PRECISION FUNCTION V(A1,A2,I1,J1,K1,I2,J2,K2,
     !                            LX1,LY1,LZ1,LX2,LY2,LZ2,I3)
      IMPLICIT REAL*8 (A-H,O-Z)
 
      include 'sizes.i'
      include 'arrays.i'
      include 'integrals.i'
      include 'pi.i'

      double precision an(maxatom,maxl,maxl,maxl,maxbt)
      common / norm / an

      V = 0.D0

      PAx = PAi(I1,J1,K1,I2,J2,K2,I1,'X')
      PBx = PAi(I1,J1,K1,I2,J2,K2,I2,'X')
      PAy = PAi(I1,J1,K1,I2,J2,K2,I1,'Y')
      PBy = PAi(I1,J1,K1,I2,J2,K2,I2,'Y')
      PAz = PAi(I1,J1,K1,I2,J2,K2,I1,'Z')
      PBz = PAi(I1,J1,K1,I2,J2,K2,I2,'Z')
      
      P0 = -2.D0 * PI * ZATOM(I3) / (A1 + A2)

      P1 = AN(I1,LX1,LY1,LZ1,K1) * AN(I2,LX2,LY2,LZ2,K2)
         
      P2 = DEXP( -( A1 * A2 ) * AB2(I1,I2) / ( A1 + A2 ) )
      
      P3 = SIGMAA(LX1,LX2,LY1,LY2,LZ1,LZ2,
     !PAx,PBx,PAy,PBy,PAz,PBz,A1,A2,I1,J1,K1,I2,J2,K2,I3)

      V = P0 * P1 * P2 * P3

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Compute Function SIGMA A
C
      DOUBLE PRECISION FUNCTION SIGMAA(LX1,LX2,LY1,LY2,LZ1,LZ2,
     !PAx,PBx,PAy,PBy,PAz,PBz,A1,A2,IA,JA,KA,IB,JB,KB,IC)
      IMPLICIT REAL*8 (A-H,O-Z)

      SIGMAA = 0.D0
      
      eps = 0.25D0 / ( A1 + A2 )

      CPx = -PAi(IA,JA,KA,IB,JB,KB,IC,'X')
      CPy = -PAi(IA,JA,KA,IB,JB,KB,IC,'Y')
      CPz = -PAi(IA,JA,KA,IB,JB,KB,IC,'Z')
      CP2 = (CPx**2) + (CPy**2) + (CPz**2)

      tau = ( A1 + A2 ) * CP2

      NIX = LX1 + LX2
      NIY = LY1 + LY2
      NIZ = LZ1 + LZ2

      DO IX = 0,NIX
         NJX = INT(DFLOAT(IX)/2.D0)
         DO JX = 0,NJX
            NKX = INT( DFLOAT(IX-2*JX) / 2.D0 )
            DO KX = 0,NKX

      DO IY = 0,NIY
         NJY = INT(DFLOAT(IY)/2.D0)
         DO JY = 0,NJY
            NKY = INT( DFLOAT(IY-2*JY) / 2.D0 )
            DO KY = 0,NKY

      DO IZ = 0,NIZ
         NJZ = INT(DFLOAT(IZ)/2.D0)
         DO JZ = 0,NJZ
            NKZ = INT( DFLOAT(IZ-2*JZ) / 2.D0 )
            DO KZ = 0,NKZ
     
      nu = IX + IY + IZ -2*(JX + JY + JZ) -(KX + KY + KZ)

      SIGMAA = SIGMAA + Fn(nu,tau) *
     !         Ax(IX,JX,KX,LX1,LX2,PAx,PBx,CPx,eps) *
     !         Ax(IY,JY,KY,LY1,LY2,PAy,PBy,CPy,eps) *
     !         Ax(IZ,JZ,KZ,LZ1,LZ2,PAz,PBz,CPz,eps)

            END DO
         END DO
      END DO


            END DO
         END DO
      END DO


            END DO
         END DO
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Compute Function Ax
C
      DOUBLE PRECISION FUNCTION Ax(IX,JX,KX,LX1,LX2,PAx,PBx,CPx,eps)
      IMPLICIT REAL*8 (A-H,O-Z)

      Ax1 = DFLOAT( (-1)**KX ) * fj(IX,LX1,LX2,PAx,PBx)
      Ax2 = FACT(IX) * ( CPx**( IX -2*(JX + KX) ) ) * 
     !      ( eps**(JX + KX) )
      Ax3 = FACT(JX) * FACT(KX) * FACT( IX -2*JX -2*KX )
      Ax  = Ax1 * Ax2 / Ax3

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C F_nu(tau) Function
C
       DOUBLE PRECISION FUNCTION Fn(nu,tau)
       IMPLICIT REAL*8 (A-H,O-Z)

      include 'pi.i'

       Fn = 0.D0

       IF (nu.EQ.0) THEN
           Fn = 1.D0
           IF (tau .LT. 1.D-10) RETURN
           Fn = .5D0 * DSQRT(PI/tau) * DERF(DSQRT(tau))
           RETURN
       ELSE IF (tau .LT. 1.D-10) THEN
           Fn = (DFLOAT( 2*nu+1 ))**(-1)
           RETURN
       ELSE IF (tau .GT. 20.D0) THEN
           P1 = FACFAC( 2*nu-1 )
           P2 = DSQRT(PI / tau**(2*nu+1) )
           P3 = 2D0**(nu+1)
           Fn = FACFAC( 2*nu-1 ) * DSQRT(PI / (tau**(2*nu+1)) ) /
     !          ( 2.D0**(nu + 1) )
           RETURN
       ELSE IF (tau .LT. .1D0) THEN
           I = 0
           FnBack = 1.
           DO WHILE (DABS(Fn-FnBack) .GT. 1.d-15)
              FnBack = Fn
              Fn = Fn + ((-tau)**I) / (FACT(I) * DFLOAT(2*nu+2*I+1))
              I = I + 1
           END DO
           RETURN
       ELSE 
           Fn = .5D0 * DSQRT(PI/tau) * DERF(DSQRT(tau))
           DO i = 1,nu
              Fn = ( Fn * DFLOAT(2*i-1) - DEXP(-tau) ) / 2.D0 / tau
           END DO
           RETURN
       END IF

       I = 0
 01    CONTINUE
       FnBack = Fn
       Fn = Fn + DEXP(-tau) * FACFAC(2*nu-1) * ((2.D0*tau)**I) / 
     !           FACFAC(2*nu+2*I+1)
       I = I + 1
       IF (DABS(Fn-FnBack) .GT. 1.d-15) GO TO 01

       RETURN
       END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Kinetic Energy Integrals
C
      DOUBLE PRECISION FUNCTION T(A1,A2,I1,J1,K1,I2,J2,K2,
     !                            LX1,LY1,LZ1,LX2,LY2,LZ2)
      IMPLICIT REAL*8 (A-H,O-Z)
 
      include 'sizes.i'
      include 'integrals.i'

      double precision an(maxatom,maxl,maxl,maxl,maxbt)
      common / norm / an

      T = 0.D0

      LMN2 = LX2 + LY2 + LZ2
      P1 = A2 * DFLOAT( 2 * LMN2 + 3 ) * 
     !     SAB(A1,A2,I1,J1,K1,I2,J2,K2,LX1,LY1,LZ1,LX2,LY2,LZ2,'F')
   
      S1 = SAB(A1,A2,I1,J1,K1,I2,J2,K2,LX1,LY1,LZ1,LX2+2,LY2,LZ2,'F')
      S2 = SAB(A1,A2,I1,J1,K1,I2,J2,K2,LX1,LY1,LZ1,LX2,LY2+2,LZ2,'F')
      S3 = SAB(A1,A2,I1,J1,K1,I2,J2,K2,LX1,LY1,LZ1,LX2,LY2,LZ2+2,'F')
      P2 = -2.D0 * (A2**2) * ( S1 + S2 + S3 )

      S1     = 0.D0
      S2     = 0.D0
      S3     = 0.D0
      P3     = 0.D0

      IF ( LX2-2 .GE. 0 ) THEN
           S1 = -.5D0 * DFLOAT(LX2) * (DFLOAT(LX2-1)) * 
     ! SAB(A1,A2,I1,J1,K1,I2,J2,K2,LX1,LY1,LZ1,LX2-2,LY2,LZ2,'F')
      END IF

      IF ( LY2-2 .GE. 0 ) THEN
           S2 = -.5D0 * DFLOAT(LY2) * (DFLOAT(LY2-1)) * 
     ! SAB(A1,A2,I1,J1,K1,I2,J2,K2,LX1,LY1,LZ1,LX2,LY2-2,LZ2,'F')
      END IF

      IF ( LZ2-2 .GE. 0 ) THEN
           S3 = -.5D0 * DFLOAT(LZ2) * (DFLOAT(LZ2-1)) * 
     ! SAB(A1,A2,I1,J1,K1,I2,J2,K2,LX1,LY1,LZ1,LX2,LY2,LZ2-2,'F')
      END IF

      P3 = S1 + S2 + S3

      P4 = AN(I1,LX1,LY1,LZ1,K1) * AN(I2,LX2,LY2,LZ2,K2)

      T = ( P1 + P2 + P3 ) * P4

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Overlap Integrals
C
      DOUBLE PRECISION FUNCTION SAB(A1,A2,I1,J1,K1,I2,J2,K2,
     !                              LX1,LY1,LZ1,LX2,LY2,LZ2,NORM)
      IMPLICIT REAL*8 (A-H,O-Z)
      
      include 'sizes.i'
      include 'pi.i'
      include 'integrals.i'

      double precision an(maxatom,maxl,maxl,maxl,maxbt)
      common / norm / an

      CHARACTER*1 NORM

      SAB = 0.D0

      PAx = PAi(I1,J1,K1,I2,J2,K2,I1,'X')
      PBx = PAi(I1,J1,K1,I2,J2,K2,I2,'X')
      PAy = PAi(I1,J1,K1,I2,J2,K2,I1,'Y')
      PBy = PAi(I1,J1,K1,I2,J2,K2,I2,'Y')
      PAz = PAi(I1,J1,K1,I2,J2,K2,I1,'Z')
      PBz = PAi(I1,J1,K1,I2,J2,K2,I2,'Z')

      P0 = ( PI/(A1+A2) )**1.5D0

      IF (NORM .EQ. 'F') THEN
          P1 = 1.D0
      ELSE
          P1 = AN(I1,LX1,LY1,LZ1,K1) * AN(I2,LX2,LY2,LZ2,K2)
      END IF

      P2 = DEXP( -( A1 * A2 ) * AB2(I1,I2) / ( A1 + A2 ) )

      P3 = Sx(A1,LX1,PAx,A2,LX2,PBx) * Sx(A1,LY1,PAy,A2,LY2,PBy) *
     !     Sx(A1,LZ1,PAz,A2,LZ2,PBz)

      SAB = P0 * P1 * P2 * P3

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Compute Function Ix (Overlap Integrals)
C
      DOUBLE PRECISION FUNCTION Sx(A1,LX1,PAx,A2,LX2,PBx)
      IMPLICIT REAL*8 (A-H,O-Z)

      include 'pi.i'

      Sx = 0.D0

      DO IX = 0,INT(DFLOAT(LX1+LX2)/2.D0)

      Sx = Sx + fj(2*IX,LX1,LX2,PAx,PBx) * FACFAC(2*IX-1)
     !        / ((2.D0 * (A1+A2) )**IX)

      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Compute Function C(n,r)
C
      DOUBLE PRECISION FUNCTION Cnr(M,N)
      IMPLICIT REAL*8 (A-H,O-Z)

      Cnr = 1.D0
      IF (N.GT.M) STOP ' Error in function Cnr(M,N) '
      Cnr = FACT(M) / ( FACT(M-N) * FACT(N) )
      
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Compute Factorial: n!
C
      DOUBLE PRECISION FUNCTION FACT(M)
      IMPLICIT REAL*8 (A-H,O-Z)

      FACT = 1.D0
      IF (M.LE.-1) STOP ' Error in function FACT(M) '
      IF (M.EQ.1 .OR. M.EQ.0) RETURN

      DO I = 1,M
         FACT = FACT * DFLOAT(I)
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Compute Function ABi = (r_A - r_B) dot vec(e_i)
C
      DOUBLE PRECISION FUNCTION PAi(IA,JA,KA,IB,JB,KB,IC,DIREC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*1 DIREC

      include 'sizes.i'
      include 'arrays.i'
      include 'integrals.i'

      double precision an(maxatom,maxl,maxl,maxl,maxbt)
      common / norm / an

      PAi = 0.D0
      IF (IA.EQ.IB .AND. IB.EQ.IC) RETURN

      if(direc.eq.'X') then
        C1 = x(IA)
        C2 = x(IB)
        C3 = x(IC)
      else if(direc.eq.'Y') then
        C1 = y(IA)
        C2 = y(IB)
        C3 = y(IC)
      else if(direc.eq.'Z') then
        C1 = z(IA)
        C2 = z(IB)
        C3 = z(IC)
      end if

      A1 = ALPHA(IA,JA,KA)
      A2 = ALPHA(IB,JB,KB)

      PAi = ( A1*(C1-C3) + A2*(C2-C3) ) / (A1+A2)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Compute Function (2*k-1)!!
C
C L = 2*k-1
C
      DOUBLE PRECISION FUNCTION FACFAC(L)
      IMPLICIT REAL*8 (A-H,O-Z)
   
      FACFAC = 1.D0
      IF (L .EQ. 0 .OR. L .EQ. -1) RETURN
      IF (L .LT. -1) STOP 'Error in FACFAC: L .LT. -1' 

      IODD = 1
      DO WHILE (IODD.LE.L)
         FACFAC = FACFAC * DFLOAT(IODD)
         IODD = IODD + 2
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Compute Function AB2 = |r_A - r_B|^2
C
      DOUBLE PRECISION FUNCTION AB2(IA,IB)
      IMPLICIT REAL*8 (A-H,O-Z)

      include 'sizes.i'
      include 'arrays.i'

      AB2 = 0.D0
      IF (IA.EQ.IB) RETURN

      AB2 = ( x(IA) - x(IB) )**2 + 
     *      ( y(IA) - y(IB) )**2 +
     *      ( z(IA) - z(IB) )**2
      
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Compute Function fj(l,m,a,b)
C
      DOUBLE PRECISION FUNCTION fj(J,L,M,A,B)
      IMPLICIT REAL*8 (A-H,O-Z)
      
      fj = 0.D0

      Ist = MAX0(0,J-M)
      Ied = MIN0(J,L)
 
      DO I = Ist,Ied

      fj = fj + Cnr(L,I) * Cnr(M,J-I) *
     !        ( A**(L-I) ) * ( B**(M+I-J) )

      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Decide if ERI must be evaluated
C
      SUBROUTINE JUMP(I,J,K,L,NI,NJ,NK,NL,WTODO)
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL WTODO

      WTODO = .FALSE.
      NI = I 
      NJ = J 
      NK = K 
      NL = L 
      
      IF ( I.GE.J .AND. K.GE.L .AND. IJUMP(I,J).GE.IJUMP(K,L) ) THEN
           WTODO = .TRUE.
           RETURN
      END IF

      IF ( NI.LT.NJ ) THEN
           NI = J 
           NJ = I
      END IF
      IF ( NK.LT.NL ) THEN
           NK = L 
           NL = K
      END IF

      IF ( IJUMP(NI,NJ).LT.IJUMP(NK,NL) ) THEN

           NITMP = NI 
           NJTMP = NJ 
           NKTMP = NK 
           NLTMP = NL
           NI = NKTMP 
           NJ = NLTMP 
           NK = NITMP 
           NL = NJTMP

      END IF

      RETURN 
      END

      INTEGER FUNCTION IJUMP(I,J)

      IJUMP = ( ( I * (I-1) )/2 ) + J

      RETURN
      END

