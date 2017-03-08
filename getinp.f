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
c   Subroutine getinp: reads the input file and sets the 
c   values of each keyword found.
c
      subroutine getinp

      include 'kill.i'
      include 'sizes.i'
      include 'files.i'
      include 'options.i'
      include 'arrays.i'

      logical empty
      integer i, j, k
      integer nline, ilast
      character*80 record, keyword, keyvalue
      character*80 inpline(maxline)

c Reading input file, line by line

      nline = 0
      do while (.true.)
        read(5,15,err=10,end=10) record

   15 format(a80)

c Checking if line is empty

        i = 0
        empty = .true.
        do while (empty.and.i.lt.80)
          i = i + 1
          if(record(i:i).gt.' ') empty = .false.
        end do
        
c If line is not empty nor is a comentarie, save it in array inpline

        if(.not.empty) then
          if(record(1:1).ne.'#') then
            nline = nline + 1
            inpline(nline) = record(1:80)
          end if
        end if

      end do
   10 continue
        
c Searching for keywords in the inpline array

      do i = 1, nline

c Reading first word in input line (must be the keyword)

        j = 0
        record = inpline(i)
        do k = 1, 80
          if(record(k:k).gt.' ') then
            j = j + 1
            keyword(j:j) = record(k:k)
          end if
          if(record(k:k).le.' '.and.j.ne.0) goto 30
        end do
   30   continue
        ilast = j

c Reading second word in input line (must be keyvalue)

        j = 0
        do k = ilast + 1, 80
          if(record(k:k).gt.' ') then
            j = j + 1
            keyvalue(j:j) = record(k:k)
          end if
          if(record(k:k).le.' '.and.j.ne.0) goto 40
        end do
   40   continue

c Checking keywords assigning key values

        if(keyword(1:ilast).eq.'xyzfile') then
          if(j.eq.0) goto 60
          xyzfile = keyvalue(1:j)

        else if(keyword(1:ilast).eq.'charge') then
          if(j.eq.0) goto 60
          read(keyvalue(1:j),*) icharge
          nel = - icharge
      
        else if(keyword(1:ilast).eq.'basefile') then
          if(j.eq.0) goto 60
          basefile = keyvalue(1:j)

        else if(keyword(1:ilast).eq.'contract') then
          if(j.eq.0) goto 60
          if(keyvalue(1:j).eq.'yes') then 
            contract = .true.
          else if(keyvalue(1:j).ne.'no') then
            write(*,42)
            kill = .true.
   42       format(' ERROR: Invalid contract value ')
          end if
          
        else if(keyword(1:ilast).eq.'guess') then
          if(j.eq.0) goto 60
          if(keyvalue(1:j).eq.'zero') then
            iguess = 0
          else if(keyvalue(1:j).eq.'huckel') then
            iguess = 1
          else if(keyvalue(1:j).eq.'feasible') then
            iguess = 2
          else
            iguess = 3
            densfile = keyvalue(1:j)
          end if

        else if(keyword(1:ilast).eq.'maxitfix') then
          if(j.eq.0) goto 60
          read(keyvalue(1:j),*,err=60) maxitfix

        else if(keyword(1:ilast).eq.'ecrit') then
          if(j.eq.0) goto 60
          read(keyvalue(1:j),*,err=60) ecrit  

        else if(keyword(1:ilast).eq.'printlevel') then
          if(j.eq.0) goto 60
          read(keyvalue(1:1),*,err=44) printlevel
   44     if(keyvalue(1:1).ne.'0'.and.
     &       keyvalue(1:1).ne.'1'.and.
     &       keyvalue(1:1).ne.'2'.and.
     &       keyvalue(1:1).ne.'3') then
               write(*,45)
               kill = .true.
   45          format(' ERROR: Invalid printlevel value ')
          end if

c Stop if errors are encountered in the input file

        else
          write(*,50) keyword(1:ilast)
   50     format(' ERROR: unrecongnized keyword found: ', a15)
          kill = .true.
        end if

   60   if(j.eq.0) then 
          write(*,70) keyword(1:ilast)
   70     format(' ERROR: keyword without value or invalid value: ',a15)
          kill = .true.
        end if

      end do

c Checking if there are missing necessary parameters

      if(xyzfile(1:1).eq.'#') then
        write(*,80)
   80   format(' ERROR: Cartesian coordinate file not specified')
        killnow = .true.
      end if

      if(basefile(1:1).eq.'#') then
        write(*,90)
   90   format(' ERROR: Base file not specified')
        killnow = .true.
      end if 

      return
      end

c
c  Subroutine writedata: Writes the specificities of the problem
c                        being run
c

      subroutine writedata
 
      include 'kill.i'
      include 'sizes.i'
      include 'files.i'
      include 'options.i'
      include 'arrays.i'

      write(*,*)       

      call setlength(xyzfile,length)
      write(*,*) ' Cartesian coordiante file: ', xyzfile(1:length)

      call setlength(basefile,length)
      write(*,*) ' Basis set file: ', basefile(1:length)
      
      write(*,*) ' Charge = ', icharge

      write(*,*) ' Energy convergence criterium = ', ecrit

      write(*,*) ' Maximum number of iterations = ', maxitfix

      write(*,*)       

      return
      end

      subroutine setlength(file,length)

      implicit double precision (a-h,o-z)
      character*80 file

      length = 1
      do while(file(length:length).gt.' ')
        length = length + 1
      end do
      length = length - 1

      return
      end
        


 
      













