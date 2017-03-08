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
c   Subroutine fmm: solves the SCF problem using the
c       algorithm of Francisco, Martinez and Martinez
c   
c

      subroutine fmm

      implicit double precision (a-h,o-z)      

      include 'sizes.i'
      include 'arrays.i'
      include 'options.i'
      include 'files.i'

      integer nc(maxd)
      double precision ss(maxd,maxd)
      double precision u(maxd,maxd), ut(maxd,maxd)
      double precision temp(maxd,maxd), ftemp(maxd,maxd),
     *                 fl(maxd,maxd), temp2(maxd,maxd)
      double precision xx(maxd,maxd), xlast(maxd,maxd), xy(maxd,maxd),
     *                 fx(maxd,maxd), fy(maxd,maxd),
     *                 gx(maxd,maxd), gxlast(maxd,maxd),
     *                 ssrp(maxd,maxd), ssrm(maxd,maxd),
     *                 zmat(maxd,maxd), vvec(maxd,maxd), uvec(maxd,maxd)
      character*5 source
      character*80 current

c Checking if the total number of electron is even

      if(mod(nel,2).ne.0) then
        write(*,*) ' ERROR: Odd number of electrons '
        stop
      else
        ndm = nel / 2
      end if

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
          temp2(j,i) = temp(i,j)
        end do
      end do
      
      call pmatrix(kdim, kdim, kdim, ssrp, temp2, ssrp)
      call pmatrix(kdim, kdim, kdim, temp, ssrp, ssrp)

      call pmatrix(kdim, kdim, kdim, ssrm, temp2, ssrm)
      call pmatrix(kdim, kdim, kdim, temp, ssrm, ssrm)

c Obtaining U and Ut such that Ut S U = I
   
      do i = 1, kdim
        ss(i,i) = 1.d0 / dsqrt(ss(i,i))
      end do
 
      call pmatrix(kdim, kdim, kdim, temp, ss, u)

      do i = 1, kdim
        do j = 1, kdim
          ut(i,j) = u(j,i)
        end do
      end do

c Getting feasible initial point (iguess = 2 or 3)

999   continue
      call guess(ndm,kdim,xx,iguess)
      if(iguess.eq.2) call pmatrix(kdim,kdim,kdim,ssrm,xx,xx)

c Testing if the point is feasible

      call pmatrix(kdim, kdim, ndm, s, xx, temp)
      do i = 1, kdim    
        do j = 1, ndm   
          temp2(j,i) = xx(i,j)
        end do
      end do
      call pmatrix(ndm, kdim, ndm, temp2, temp, temp2)
      
      do i = 1, ndm
        do j = 1, ndm
          if(i.eq.j) then
            test = abs(temp2(i,j) - 1.d0)
          else
            test = abs(temp2(i,j))
          end if  
          if(test.gt.1.d-8) then
            write(*,*) ' FATAL ERROR: Initial Point is not feasible '
            write(*,*) ' I will try again... '
            goto 999
            stop  
          end if
        end do
      end do

c Naming current point file

      call setlength(xyzfile,length)
      i = length
      do while(xyzfile(i:i).ne.'/')
        i = i - 1
      end do
      current = 'CURRENT_'//xyzfile(i+1:length)
   40 format(200(tr2,d17.10))

c Calling fixed point algorithm

      do i = 1, kdim
        do j = 1, ndm
          temp(i,j) = xx(i,j)
        end do
      end do
      call fixpoint(temp)

c
c Starting fmm algorithm
c

      nit = 1
      dehf = 1.d0

      sigmak = 0.5
      alpha = 1.d-4

c Computing the first fock matrix and corresponding energy

      call computefe(xx,fx,ex)

	write(*,10) ex
   10 format(/,72('-'),/,t5,' HF Energy at initial point = ', e17.10,/)  

c Starting iterative procedure

      t = 1.d0
      do while (dehf.gt.ecrit.and.nit.le.maxitfix)

c
c Performing a fixed-point iteration
c

c Computing F' = Ut F U

        k = kdim
        call pmatrix(k, k, k, fx, u, temp)
        call pmatrix(k, k, k, ut, temp, fl)         

c Diagonalize F'

        call diag(fl,xy,kdim)

c Ordering eigenvectors and eigenvalues from lower to higher

        do i = 1, kdim
          nc(i) = 1
          do j = 1, kdim
            if(fl(i,i).gt.fl(j,j)) nc(i) = nc(i) + 1
            if(fl(i,i).eq.fl(j,j).and.i.gt.j) nc(i) = nc(i) + 1
          end do
        end do

        do i = 1, kdim
          ftemp(i,i) = fl(i,i)
          do j = 1, kdim
            temp(i,j) = xy(i,j)
          end do
        end do

        do i = 1, kdim
          fl(nc(i),nc(i)) = ftemp(i,i)
          do j = 1, kdim
            xy(j,nc(i)) = temp(j,i)
          end do
        end do

c Computing the matrix Y from the eigenvectors of F'

        k = kdim
        call pmatrix(k, k, k, u, xy, temp) 

        do i = 1, kdim
          do j = 1, kdim
            xy(i,j) = temp(i,j)
          end do
        end do 

c Computing the new fock and density matrix
      
        call computefe(xy,fy,ey)

c Checking if convergence was achieved

        ared = ex - ey
        dehf = abs(ared)
 
        write(*,20) nit, ex, ey, dehf
20      format(/,
     *         i4,' HF Energy at current point = ', e17.10,/ 
     *         t5,' HF Energy at trial point   = ', e17.10,/ 
     *        t22,' Variation = ', f17.10,/ )  
 
c If not converged, test current point

        if(dehf.gt.ecrit) then

c Computing gx(i) = 4F(X)X(i)

          call pmatrix(kdim,kdim,ndm,fx,xx,gx)
          do i = 1, kdim
            do j = 1, ndm
              gx(i,j) = 4.d0 * gx(i,j) 
            end do
          end do

c Computing the new sigma_k

          if(nit.gt.0) then
            anumer = 0.d0
            do i = 1, kdim
              do j = 1, ndm
                temp(i,j) = xx(i,j) - xlast(i,j) 
                anumer = anumer + temp(i,j) * (gx(i,j) - gxlast(i,j))
              end do
            end do

            call pmatrix(kdim,kdim,ndm,s,temp,temp2)

            adenon = 0.d0
            do i = 1, kdim
              do k = 1, ndm
                adenon = adenon + temp(i,j) * temp2(i,j)
              end do
            end do          

            acheck = anumer / adenon          
        
            sigmak = dmax1(1.d-2,dmin1(100.d0,acheck)) 
          end if 

c Testing if the trial point c is good enough (Step 2)

          neasy = 0
          if(t.gt.4.d0) t = t / 4.d0

          do while(ared.lt.(alpha*pred).or.neasy.eq.0)
    
            ared = ex - ey

            do i = 1, kdim
              do j = 1, ndm
                temp(i,j) = xy(i,j) - xx(i,j)
              end do 
            end do

            pred = 0.d0
            if(neasy.eq.0) then
              call pmatrix(kdim,kdim,ndm,fx,temp,temp2)
              do i = 1, kdim
                do j = 1, ndm
                  pred = pred - gx(i,j) * temp(i,j)
     *                 -0.5 * temp(i,j) * 4.d0 * temp2(i,j)
                end do
              end do
            else
              call pmatrix(kdim,kdim,ndm,s,temp,temp2)
              do i = 1, kdim
                do j = 1, ndm
                  pred = pred - gx(i,j) * temp(i,j)
     *                  -0.5 * temp(i,j) * sigmak * temp2(i,j)
                end do
              end do
            end if

c Printing iteration's data

            if(neasy.eq.0) then
              source = 'FIXED'
            else
              if(neasy.eq.1) write(*,*) 
     *                       '    Entering EASY solver with t = ', t 
              source = 'EASY '
              t = 2.d0 * t
              if(t.gt.1.d14) write(*,*)
     *                     ' WARNING: t > 1.d14: ',
     *                     ' Probably working under the',
     *                     ' floating point accuracy. '
            end if
           
            write(*,21) neasy, source, ey, pred, ared 
21          format(i4,tr1,a5,
     *             ' E = ',e17.10,
     *             ' PRed. = ', e17.10, 
     *             ' ARed. =  ', e17.10)

c If pred < 0.d0, stop
 
            if(pred.lt.0.d0) then

c Testing the feasibility of the current point
              call pmatrix(kdim, kdim, ndm, s, xx, temp)
              do i = 1, kdim    
                do j = 1, ndm   
                  temp2(j,i) = xx(i,j)
                end do
              end do
              call pmatrix(ndm, kdim, ndm, temp2, temp, temp2)
              do i = 1, ndm
                do j = 1, ndm
                  if(i.eq.j) then
                    test = abs(temp2(i,j) - 1.d0)
                  else
                    test = abs(temp2(i,j))
                  end if  
                  if(test.gt.1.d-4) then
                    write(*,*) ' FATAL ERROR: ',
     *                         'Current Point is not feasible '
                    current = 'CURRENTNOTFEASIBLE_'//current(9:80)
                    open(10,file=current)
                    do ii = 1, kdim
                      write(10,40) (xx(ii,jj), jj = 1, ndm)
                    end do
                    close(10)
                    stop  
                  end if
                end do
              end do

c Testing the feasibilty of the trial point
              call pmatrix(kdim, kdim, ndm, s, xy, temp)
              do i = 1, kdim    
                do j = 1, ndm   
                  temp2(j,i) = xy(i,j)
                end do
              end do
              call pmatrix(ndm, kdim, ndm, temp2, temp, temp2)
              do i = 1, ndm
                do j = 1, ndm
                  if(i.eq.j) then
                    test = abs(temp2(i,j) - 1.d0)
                  else
                    test = abs(temp2(i,j))
                  end if  
                  if(test.gt.1.d-4) then
                    write(*,*) ' FATAL ERROR: ',
     *                         'Trial Point is not feasible '
                    current = 'TRIALNOTFEASIBLE_'//current(9:80)
                    open(10,file=current)
                    do ii = 1, kdim
                      write(10,40) (xy(ii,jj), jj = 1, ndm)
                    end do
                    close(10)
                    stop  
                  end if
                end do
              end do

              if(pred.gt.1.d-8) then
                write(*,*) '  FATAL ERROR: ',
     *                     'Current and Trial points are feasible, ',
     *                     'but PRed is too large: something is wrong.'
              else
                write(*,*) '  BOTH the Current and the Trial points',
     *                     ' area feasible and Pred < 0: '
                write(*,*) '  Probably the current point is a solution.'
              end if
              stop

            end if 
 
c If the trial point is not good enough, solve the Easy problem (Step 3)

            neasy = neasy + 1
            if(ared.lt.alpha*pred) then

              call pmatrix(kdim,kdim,ndm,ssrm,gx,temp2)
              call pmatrix(kdim,kdim,ndm,ssrp,xx,temp)

              sigt = 1.d0 / (sigmak * t)
              do i = 1, kdim
                do j = 1, ndm
                  zmat(i,j) = temp(i,j) - sigt * temp2(i,j)
                end do 
              end do

              do i = 1, kdim
                do j = 1, ndm
                  temp(j,i) = zmat(i,j)
                end do
              end do

              call pmatrix(ndm,kdim,ndm,temp,zmat,temp2)
              call diag(temp2,vvec,ndm)
              call pmatrix(kdim,ndm,ndm,zmat,vvec,uvec)

              do j = 1, ndm
                vecnorm = 0.d0
                do i = 1, kdim
                  vecnorm = vecnorm + uvec(i,j)**2
                end do
                vecnorm = dsqrt(vecnorm)
                do i = 1, kdim
                  uvec(i,j) = uvec(i,j) / vecnorm
                end do
              end do

              do i = 1, ndm
                do j = 1, ndm
                  temp(j,i) = vvec(i,j)
                end do
              end do
           
              call pmatrix(kdim,ndm,ndm,uvec,temp,temp2)
              call pmatrix(kdim,kdim,ndm,ssrm,temp2,xy)
              call computefe(xy,fy,ey)

            end if

          end do

        end if


c Saving gxlast <- gx and accepting point

        dehf = abs(ex - ey)
        ex = ey

        do i = 1, kdim
          do j = 1, ndm
            gxlast(i,j) = gx(i,j)
            xlast(i,j) = xx(i,j)
            xx(i,j) = xy(i,j)
          end do
        end do

        do i = 1, kdim
          do j = 1, kdim
            fx(i,j) = fy(i,j)
          end do
        end do

        nit = nit + 1
  
C Saving current point

        open(10,file=current)
        do i = 1, kdim
          write(10,40) (xx(i,j), j = 1, ndm)
        end do
        close(10)
 
	end do  

      if(neasy.eq.1) then
        write(*,50) ex
50      format(/,72('-'),/,/,
     *         ' Convergence achieved in fixed point iteration',/,
     *         ' Hartree-Fock Energy at the solution: ', e17.10,/,/,
     *         ' The solution satisfies de AUFBAU principle! ',
     *         /,/,72('-'),/)
      else
        write(*,51) ex
51      format(/,72('-'),/,/,
     *         ' Convergence achieved ',/,
     *         ' Hartree-Fock Energy at the solution: ', e17.10,/,/,
     *         ' The solution does NOT satisfy de AUFBAU principle.',
     *         /,/,72('-'),/)
      end if      

      return
      end



