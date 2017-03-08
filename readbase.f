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
c   Subroutine readbase: reads the basis data file
c   
c

      subroutine readbase

      include 'kill.i'
      include 'sizes.i'         
      include 'files.i'
      include 'options.i'
      include 'arrays.i'

      integer i

c Setting the total number of functions to 0

      ntotf = 0

c Reading basis for each element in the coordinate cartesian file:
c Uncontracted basis sets:

      if(.not.contract) then
        do i = 1, natoms
          if(kill) goto 100
          if(ele(i).eq.'H') then
            zatom(i) = 1.d0
            nel = nel + int(zatom(i))
            call asubase('HYDROGEN',i)
          else if(ele(i).eq.'He') then
            zatom(i) = 2.d0
            nel = nel + int(zatom(i))
            call asubase('HELIUM',i)
          else if(ele(i).eq.'Li') then
            zatom(i) = 3
            nel = nel + int(zatom(i))
            call asubase('LITHIUM',i)
          else if(ele(i).eq.'N') then
            zatom(i) = 7.d0
            nel = nel + int(zatom(i))
            call asubase('NITROGEN',i)
          else if(ele(i).eq.'C') then
            zatom(i) = 6.d0
            nel = nel + int(zatom(i))
            call asubase('CARBON',i)
          else if(ele(i).eq.'O') then
            zatom(i) = 8.d0
            nel = nel + int(zatom(i))
            call asubase('OXYGEN',i)
          else if(ele(i).eq.'S') then
            zatom(i) = 16.d0
            nel = nel + int(zatom(i))
            call asubase('SULFUR',i)
          else if(ele(i).eq.'Cr') then
            zatom(i) = 24.d0
            nel = nel + int(zatom(i))
            call asubase('CHROMIUM',i)
          end if
        end do
      end if

c Contracted basis sets: 

      if(contract) then
         kcont = 0
         do i = 1, natoms
          if(kill) goto 100
          if(ele(i).eq.'H') then
            zatom(i) = 1.d0
            nel = nel + int(zatom(i))
            call ascbase('HYDROGEN',i,k)
          else if(ele(i).eq.'He') then
            zatom(i) = 2.d0
            nel = nel + int(zatom(i))
            call ascbase('HELIUM',i,k)
          else if(ele(i).eq.'Li') then
            zatom(i) = 3
            nel = nel + int(zatom(i))
            call ascbase('LITHIUM',i,k)
          else if(ele(i).eq.'N') then
            zatom(i) = 7.d0
            nel = nel + int(zatom(i))
            call ascbase('NITROGEN',i,k)
          else if(ele(i).eq.'C') then
            zatom(i) = 6.d0
            nel = nel + int(zatom(i))
            call ascbase('CARBON',i,k)
          else if(ele(i).eq.'O') then
            zatom(i) = 8.d0
            nel = nel + int(zatom(i))
            call ascbase('OXYGEN',i,k)
          else if(ele(i).eq.'S') then
            zatom(i) = 16.d0
            nel = nel + int(zatom(i))
            call ascbase('SULFUR',i,k)
          else if(ele(i).eq.'Cr') then
            zatom(i) = 24.d0
            nel = nel + int(zatom(i))
            call ascbase('CHROMIUM',i,k)
          end if
          kcont = kcont + k
        end do
      end if

  100 return
      end

c
c ####################################################
c # Subroutine asubase: Author: Leandro Martínez     #
c ####################################################
c
c  Subroutine asubase: Assign the exponents and functions types 
c  for uncontracted basis sets
c

      subroutine asubase(element,iatom)

      include 'kill.i'
      include 'files.i'
      include 'sizes.i'
      include 'arrays.i'
      include 'integrals.i'

      integer i
      integer iatom
      integer nchar, nfunc, nftype
      character*(*) element
      character*80 record, string
      character*1 ftype
        
c Defining the length of the string associated with each element.

      if(element.eq.'HYDROGEN') nchar = 8
      if(element.eq.'HELIUM') nchar = 6
      if(element.eq.'LITHIUM') nchar = 7
      if(element.eq.'NITROGEN') nchar = 8
      if(element.eq.'CARBON') nchar = 6
      if(element.eq.'OXYGEN') nchar = 6
      if(element.eq.'SULFUR') nchar = 6
      if(element.eq.'CHROMIUM') nchar = 8

c Opening the base file

      open(10, file = basefile, err=70 , status='old')
      goto 90

   70 write(*,80)
   80 format(' ERROR: Could not find base file')
      kill = .true.
      goto 200 

c Searching for the correct element base:
      
   90 string = '#'
      do while (string.ne.element)
        read(10,100,err=140,end=140) record
        read(record,*,err=120,end=120) string
  100   format(a80)
  120 end do


c Reading the number of functions of all types in this base

      read(record,*) string, nfunc

c Reading exponents (zeta) for function types S, P, D, F

      is = 0
      ip = 0
      id = 0
      if = 0
      ifunc = 0
      do while(ifunc.lt.nfunc)
        read(10,*,err=130) ftype, nftype

          i = 0
          if(ftype.eq.'S') then
            do while(i.lt.nftype)
              ifunc = ifunc + 1 
              is = is + 1
              call checkbase(is)
              read(10,*,err=130) i, alpha(iatom,1,is)
            end do

          i = 0
          else if(ftype.eq.'P') then
            do while(i.lt.nftype)
              ifunc = ifunc + 1 
              ip = ip + 1
              call checkbase(ip)
              read(10,*,err=130) i, alpha(iatom,1,ip) 
            end do
      
          i = 0
          else if(ftype.eq.'L') then
            do while(i.lt.nftype)
              ifunc = ifunc + 2 
              is = is + 1
              ip = ip + 1
              call checkbase(is)
              call checkbase(ip)
              read(10,*,err=130) i, alpha(iatom,1,is)
              alpha(iatom,2,ip) = alpha(iatom,1,is)
            end do

          i = 0
          else if(ftype.eq.'D') then
            do while(i.lt.nftype)
              ifunc = ifunc + 1 
              id = id + 1
              call checkbase(id)
              read(10,*,err=130) i, alpha(iatom,3,id)
            end do

          i = 0
          else if(ftype.eq.'F') then
            do while(i.lt.nftype)
              ifunc = ifunc + 1 
              if = if + 1
              call checkbase(if)
              read(10,*,err=130) i, alpha(iatom,4,if)
            end do
          end if
         
      end do

c Assigning the number of different ang. mom. functions of each atom

      if(is.ne.0) nang(iatom) = 1
      if(ip.ne.0) nang(iatom) = 2
      if(id.ne.0) nang(iatom) = 3
      if(if.ne.0) nang(iatom) = 4

c Assigning the number of functions of each type (array ncgtf)

      ncgtf(iatom,1) = is
      ncgtf(iatom,2) = ip
      ncgtf(iatom,3) = id
      ncgtf(iatom,4) = if

c Adding nfunc to ntotf 

      ntotf = ntotf + nfunc

c Defining the dimension of the problem:

      kdim = kdim + is + 3*ip + 6*id + 10*if

      close(10)
      return 

c Stopping if some necessary element was not found in base:

  130 write(*,135) element(1:nchar)
  135 format(' ERROR: Format error in base data for element ', a15)
      kill = .true.
      goto 200

  140 write(*,150) element(1:nchar)
  150 format(' ERROR: Could not find base for ', a15)
      close(10)
      kill = .true.

  200 return
      end

c
c ####################################################
c # Subroutine ascbase: Author: Leandro Martínez     #
c ####################################################
c
c  Subroutine ascbase: Assign the exponents, functions types 
c  and contractions coefficients for contracted basis sets
c

      subroutine ascbase(element,iatom,k)

      include 'kill.i'
      include 'files.i'
      include 'sizes.i'
      include 'arrays.i'
      include 'integrals.i'

      integer i
      integer iatom
      integer nchar, nfunc, nftype
      character*(*) element
      character*80 record, string
      character*1 ftype
        
c Defining the length of the string associated with each element.

      if(element.eq.'HYDROGEN') nchar = 8
      if(element.eq.'HELIUM') nchar = 6
      if(element.eq.'LITHIUM') nchar = 7
      if(element.eq.'NITROGEN') nchar = 8
      if(element.eq.'CARBON') nchar = 6
      if(element.eq.'OXYGEN') nchar = 6
      if(element.eq.'SULFUR') nchar = 6
      if(element.eq.'CHROMIUM') nchar = 8

c Opening the base file

      open(10, file = basefile, err=70 , status='old')
      goto 90

   70 write(*,80)
   80 format(' ERROR: Could not find base file')
      kill = .true.
      goto 200 

c Searching for the correct element base:
      
   90 string = '#'
      do while (string.ne.element)
        read(10,100,err=140,end=140) record
        read(record,*,err=120,end=120) string
  100   format(a80)
  120 end do


c Reading the number of functions of all types in this base

      read(record,*) string, nfunc

c Reading exponents (zeta) and contraction coefficients
c         for function types S, P, D, F

      is = 0
      ip = 0
      id = 0
      if = 0
      ifunc = 0
      k = 0
      ngcf(iatom,1) = 0
      do while(ifunc.lt.nfunc)
        read(10,*,err=130) ftype, nftype
          k = k + 1
      
          i = 0
          if(ftype.eq.'S') then
            ig = ngcf(iatom,1) + 1
            ngcf(iatom,1) = ig
            ncont(iatom,1,ig) = nftype
            do while(i.lt.nftype)
              ifunc = ifunc + 1 
              is = is + 1
              call checkbase(is)
              read(10,*,err=130) i, alpha(iatom,1,is), 
     *                              cont(iatom,1,ig,i)
           end do

          i = 0
          else if(ftype.eq.'P') then
            do while(i.lt.nftype)
              ifunc = ifunc + 1 
              ip = ip + 1
              call checkbase(ip)
              read(10,*,err=130) i, alpha(iatom,2,ip), 
     *                              cont(iatom,2,ip,i)
              ncont(iatom,2,ip) = i
            end do
      
          i = 0
          else if(ftype.eq.'L') then
            do while(i.lt.nftype)
              ifunc = ifunc + 2 
              is = is + 1
              ip = ip + 1
              call checkbase(is)
              call checkbase(ip)
              read(10,*,err=130) i, alpha(iatom,1,is),
     *                              cont(iatom,1,is,i),
     *                              cont(iatom,2,ip,i)
              alpha(iatom,2,ip) = alpha(iatom,1,is)
              ncont(iatom,1,is) = i
              ncont(iatom,2,ip) = i
            end do

          i = 0
          else if(ftype.eq.'D') then
            do while(i.lt.nftype)
              ifunc = ifunc + 1 
              id = id + 1
              call checkbase(id)
              read(10,*,err=130) i, alpha(iatom,3,id),
     *                              cont(iatom,3,id,i)
              ncont(iatom,3,id) = i
            end do

          i = 0
          else if(ftype.eq.'F') then
            do while(i.lt.nftype)
              ifunc = ifunc + 1 
              if = if + 1
              call checkbase(if)
              read(10,*,err=130) i, alpha(iatom,4,if),
     *                              cont(iatom,4,if,i)
              ncont(iatom,4,if) = i
            end do
          end if
         
      end do

c Assigning the number of different ang. mom. functions of each atom

      if(is.ne.0) nang(iatom) = 1
      if(ip.ne.0) nang(iatom) = 2
      if(id.ne.0) nang(iatom) = 3
      if(if.ne.0) nang(iatom) = 4

c Assigning the number of functions of each type (array ncgtf)

      ncgtf(iatom,1) = is
      ncgtf(iatom,2) = ip
      ncgtf(iatom,3) = id
      ncgtf(iatom,4) = if

c Adding nfunc to ntotf 

      ntotf = ntotf + nfunc

c Defining the dimension of the problem:

      kdim = kdim + is + 3*ip + 6*id + 10*if

      close(10)
      return 

c Stopping if some necessary element was not found in base:

  130 write(*,135) element(1:nchar)
  135 format(' ERROR: Format error in base data for element ', a15)
      kill = .true.
      goto 200

  140 write(*,150) element(1:nchar)
  150 format(' ERROR: Could not find base for ', a15)
      close(10)
      kill = .true.

  200 return
      end
 
c
c Subroutine checkbase: checks if the number of basis has 
c exceeded maxbase (see the definition of maxbase in sizes.i)
c

      subroutine checkbase(i)

      include 'sizes.i'
      integer i

      if(i.gt.maxbase) then
        write(*,10)
   10   format(' ERROR: Number of base functions exceeded maxbase.')
        stop
      end if

      return
      end






