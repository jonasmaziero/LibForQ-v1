!###################################################################################################################################
!                                       Random number generators (uniform in [0,1], or integers)
!###################################################################################################################################
!                                      Merseene-Twister RNG (the module is in a separated file)
subroutine sub_rand_seed(seed) ! Uses the computer time to generate a random seed for the initialization of the random number generators
implicit none
character(8) :: date
character(10) :: time
character(5) :: zone
integer, dimension(8) :: values
integer(4) :: seed

call date_and_time(date,time,zone,values) ;   seed = 1.d3*values(8)
  
end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine init_mt() ! Initializes the Merseene Twister pseudo-random number generator
use mt95
implicit none
integer(4) :: seed
real(8) :: rn

call sub_rand_seed(seed) ;   call genrand_init(seed) ;   call genrand_real1(rn)
 
end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rng_mt(d, rn) ! Calls the Mersenne Twister pseudo-random number generator
use mt95
implicit none
integer :: d  ! Dimension of the vector of random numbers
real(8) :: rn(1:d) ! Vector whose components are random numbers uniformly distributed in [0,1]

call genrand_real1(rn)

end
!------------------------------------------------------------------------------------------------------------------------------------
!                                                    Gnu's RNG
subroutine init_gnu()  ! Initialization for the GNU RNG RANDOM_NUMBER()
IMPLICIT NONE
INTEGER, ALLOCATABLE :: seed(:)
INTEGER :: i, n, un, istat, dt(8), pid, t(2), s
INTEGER(8) :: count, tms
  
CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))
! First try if the OS provides a random number generator
OPEN(newunit=un, file="/dev/urandom", access="stream", form="unformatted", action="read", status="old", iostat=istat)
IF (istat == 0) THEN
  read(un) seed
  close(un)
ELSE
  ! Fallback to XOR:ing the current time and pid. The PID is useful in case one launches multiple instances of the same program in parallel.
  CALL SYSTEM_CLOCK(count)
  IF (count /= 0) THEN
    t = TRANSFER(count, t)
  ELSE
    CALL DATE_AND_TIME(values=dt)
    tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 + dt(2) * 31_8 * 24 * 60 * 60 * 1000 + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 + dt(6) * 60 * 1000 + dt(7) * 1000 + dt(8)
    t = TRANSFER(tms, t)
  ENDIF
  s = ieor(t(1), t(2))
  pid = getpid() + 1099279 ! Add a prime
  s = ieor(s, pid)
  IF (n >= 3) THEN
     seed(1) = t(1) + 36269
     seed(2) = t(2) + 72551
     seed(3) = pid
     IF (n > 3) THEN
        seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
     ENDIF
  ELSE
     seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
  ENDIF
ENDIF
CALL RANDOM_SEED(put=seed)
          
end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rng_gnu(d, rn) ! Calls the Gnu "standard" random number generator. See https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fNUMBER.html 
! The runtime-library implements George Marsaglia's KISS (Keep It Simple Stupid) random number generator.
implicit none
integer :: d  ! Dimension of the vector of random numbers
real(8) :: rn(1:d) ! Vector whose components are random numbers uniformly distributed in [0,1)

call random_number(rn)

end
!------------------------------------------------------------------------------------------------------------------------------------
!                                                    Netlib's RNG
subroutine init_netlib() ! Initializes Netlib's pseudo-random number generator zufall
implicit none
integer :: seed

call sub_rand_seed(seed) ;   call zufalli(seed)
 
end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rng_netlib(d, rn) ! Calls Netlib's pseudo-random number generator zufall
implicit none
integer :: d  ! Dimension of the vector of random numbers
real(8) :: rn(1:d) ! Vector whose components are random numbers uniformly distributed in [0,1]

call zufall(d,rn)

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine zufall(n,a)  ! Netlib's RNG (this sequence of 4 subroutines is needed)
! README for zufall random number package
! ------ --- ------ ------ ------ -------
! This package contains a portable random number generator set
! for: uniform (u in [0,1)), normal (<g> = 0, <g^2> = 1), and
! Poisson distributions. The basic module, the uniform generator,
! uses a lagged Fibonacci series generator:
! 
!               t    = u(n-273) + u(n-607)
!               u(n) = t - float(int(t))
! 
! where each number generated, u(k), is floating point. Since
! the numbers are floating point, the left end boundary of the
! range contains zero. This package is nearly portable except
! for the following. (1) It is written in lower case, (2) the
! test package contains a timer (second) which is not portable,
! and (3) there are cycle times (in seconds) in data statements 
! for NEC SX-3, Fujitsu VP2200, and Cray Y-MP. Select your 
! favorite and comment out the others. Replacement functions 
! for 'second' are included - comment out the others. Otherwise 
! the package is portable and returns the same set of floating 
! point numbers up to word precision on any machine. There are 
! compiler directives ($cdir for Cray, *vdir for SX-3, and VOCL 
! for Fujitsu VP2200) which should be otherwise ignored.
! 
! To compile this beast, note that all floating point numbers
! are declared 'double precision'. On Cray X-MP, Y-MP, and C-90
! machines, use the cft77 (cf77) option -dp to run this in 64
! bit mode (not 128 bit double).
! 
! External documentation, "Lagged Fibonacci Random Number Generators
! for the NEC SX-3," is to be published in the International
! Journal of High Speed Computing (1994). Otherwise, ask the
! author: 
! 
!          W. P. Petersen 
!          IPS, RZ F-5
!          ETHZ
!          CH 8092, Zurich
!          Switzerland
! 
! e-mail:  wpp@ips.ethz.ch.
! 
! The package contains the following routines:
! 
! ------------------------------------------------------
! UNIFORM generator routines:
! 
!       subroutine zufalli(seed)
!       integer seed
! c initializes common block containing seeds. if seed=0,
! c the default value is 1802.
! 
!       subroutine zufall(n,u)
!       integer n
!       double precision u(n)
! c returns set of n uniforms u(1), ..., u(n).
! 
!       subroutine zufallsv(zusave)
!       double precision zusave(608)
! c saves buffer and pointer in zusave, for later restarts
! 
!       subroutine zufallrs(zusave)
!       double precision zusave(608)
! c restores seed buffer and pointer from zusave
! ------------------------------------------------------
      implicit none
!
! portable lagged Fibonacci series uniform random number
! generator with "lags" -273 und -607:
!
!       t    = u(i-273)+buff(i-607)  (floating pt.)
!       u(i) = t - float(int(t))
!
! W.P. Petersen, IPS, ETH Zurich, 19 Mar. 92
!
      double precision a(*)
      double precision buff(607)
      double precision t
      integer i,k,ptr,VL,k273,k607
      integer buffsz,nn,n,left,q,qq
      integer aptr,aptr0,bptr

      common /klotz0/buff,ptr
      data buffsz/607/

      aptr = 0
      nn   = n

1     continue

      if(nn .le. 0) return

! factor nn = q*607 + r

      q    = (nn-1)/607
      left = buffsz - ptr

      if(q .le. 1) then

! only one or fewer full segments

         if(nn .lt. left) then
            do 2 i=1,nn
               a(i+aptr) = buff(ptr+i)
2           continue
            ptr  = ptr + nn
            return
         else
            do 3 i=1,left
               a(i+aptr) = buff(ptr+i)
3           continue
            ptr  = 0
            aptr = aptr + left
            nn   = nn - left
!  buff -> buff case
            VL   = 273
            k273 = 334
            k607 = 0
            do 4 k=1,3
! cdir$ ivdep
! vdir nodep
! VOCL LOOP, TEMP(t), NOVREC(buff)
               do 5 i=1,VL
                  t            = buff(k273+i) + buff(k607+i)
                  buff(k607+i) = t - float(int(t))
5              continue
               k607 = k607 + VL
               k273 = k273 + VL
               VL   = 167
               if(k.eq.1) k273 = 0
4           continue

            goto 1
         endif
      else

! more than 1 full segment
 
          do 6 i=1,left
             a(i+aptr) = buff(ptr+i)
6         continue
          nn   = nn - left
          ptr  = 0
          aptr = aptr+left
 
! buff -> a(aptr0)
 
          VL   = 273
          k273 = 334
          k607 = 0
          do 7 k=1,3
             if(k.eq.1)then
! VOCL LOOP, TEMP(t)
                do 8 i=1,VL
                   t         = buff(k273+i) + buff(k607+i)
                   a(aptr+i) = t - float(int(t))
8               continue
                k273 = aptr
                k607 = k607 + VL
                aptr = aptr + VL
                VL   = 167
             else
!cdir$ ivdep
!vdir nodep
!VOCL LOOP, TEMP(t)
                do 9 i=1,VL
                   t         = a(k273+i) + buff(k607+i)
                   a(aptr+i) = t - float(int(t))
9               continue
                k607 = k607 + VL
                k273 = k273 + VL
                aptr = aptr + VL
             endif
7         continue
          nn = nn - 607

! a(aptr-607) -> a(aptr) for last of the q-1 segments

          aptr0 = aptr - 607
          VL    = 607

!vdir novector
          do 10 qq=1,q-2
             k273 = 334 + aptr0
! cdir$ ivdep
!vdir nodep
!VOCL LOOP, TEMP(t), NOVREC(a)
             do 11 i=1,VL
                t         = a(k273+i) + a(aptr0+i)
                a(aptr+i) = t - float(int(t))
11           continue
             nn    = nn - 607
             aptr  = aptr + VL
             aptr0 = aptr0 + VL
10        continue

! a(aptr0) -> buff, last segment before residual

          VL   = 273
          k273 = 334 + aptr0
          k607 = aptr0
          bptr = 0
          do 12 k=1,3
             if(k.eq.1) then
! VOCL LOOP, TEMP(t)
                do 13 i=1,VL
                   t            = a(k273+i) + a(k607+i)
                   buff(bptr+i) = t - float(int(t))
13              continue
                k273 = 0
                k607 = k607 + VL
                bptr = bptr + VL
                VL   = 167
             else
!cdir$ ivdep
! vdir nodep
! VOCL LOOP, TEMP(t), NOVREC(buff)
                do 14 i=1,VL
                   t            = buff(k273+i) + a(k607+i)
                   buff(bptr+i) = t - float(int(t))
14              continue
                k607 = k607 + VL
                k273 = k273 + VL
                bptr = bptr + VL
             endif
12        continue
          goto 1
      endif
      
end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine zufalli(seed)
      implicit none

!  generates initial seed buffer by linear congruential
!  method. Taken from Marsaglia, FSU report FSU-SCRI-87-50
!  variable seed should be 0 < seed <31328

      integer seed
      integer ptr
      double precision s,t
      double precision buff(607)
      integer ij,kl,i,ii,j,jj,k,l,m
      common /klotz0/buff,ptr
      data ij/1802/,kl/9373/

      if(seed.ne.0) ij = seed

      i = mod(ij/177,177) + 2
      j = mod(ij,177) + 2
      k = mod(kl/169,178) + 1
      l = mod(kl,169)
      do ii=1,607
         s = 0.0
         t = 0.5
         do jj=1,24
            m = mod(mod(i*j,179)*k,179)
            i = j
            j = k
            k = m
            l = mod(53*l+1,169)
            if(mod(l*m,64).ge.32) s = s+t
            t = .5*t
         enddo
         buff(ii) = s
      enddo
      
end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine zufallsv(svblk)
      implicit none
!  saves common blocks klotz0, containing seeds and 
!  pointer to position in seed block. IMPORTANT: svblk must be
!  dimensioned at least 608 in driver. The entire contents
!  of klotz0 (pointer in buff, and buff) must be saved.
      double precision buff(607)
      integer ptr,i
      double precision svblk(*)
      common /klotz0/buff,ptr
      
      svblk(1) = ptr
      do i=1,607
         svblk(i+1) = buff(i)
      enddo

      
end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine zufallrs(svblk)
      implicit none
!  restores common block klotz0, containing seeds and pointer
!  to position in seed block. IMPORTANT: svblk must be
!  dimensioned at least 608 in driver. The entire contents
!  of klotz0 must be restored.
!
      double precision buff(607)
      integer i,ptr
      double precision svblk(*)
      common /klotz0/buff,ptr

      ptr = svblk(1)
      do i=1,607
         buff(i) = svblk(i+1)
      enddo
 
      
end
!###################################################################################################################################
!                                           RN with different DISTRIBUTIONS and/or DOMAINS)
!###################################################################################################################################
subroutine rng_gauss(d, grn) ! Returns a vector of random (real) numbers with Gaussian (standard normal) probability distributions
! Ref: Katzgraber, H. G. (2010). Random numbers in scientific computing: An introduction, arXiv:1005.4117.
implicit none
integer :: d  ! Dimension of the vectors of random numbers
real(8), allocatable :: rn(:)  ! Vector of uniformly distributed random numbers
real(8) :: grn(1:d)  ! Vector of Gaussianily distributed random numbers
real(8), parameter :: pi = 4.d0*atan(1.d0)
real(8) :: logterm, angle  ! Auxiliary variables
integer :: j  ! Auxiliary variable for counters

allocate( rn(1:2*d) ) ;   call rng(2*d, rn)  ! Allocate memory for and get the uniformly distributed random numbers

do j = 1, 2*d, 2   ! Gets the Gaussianily distributed random numbers
  if ( rn(j+1) < 1.d-15 ) rn(j+1) = 1.d-15 ;   logterm = sqrt(-2.d0*log(rn(j+1))) ;   angle = 2.d0*pi*rn(j)
  grn((j+1)/2) = logterm*cos(angle)
enddo

deallocate( rn )

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rng_exp(d, ern) ! Returns a vector of random (real) numbers with exponential probability distributions
! Ref: Katzgraber, H. G. (2010). Random numbers in scientific computing: An introduction, arXiv:1005.4117.
implicit none
integer :: d  ! Dimension of the vectors of random numbers
real(8), allocatable :: rn(:)  ! Vector of uniformly distributed random numbers
real(8) :: ern(1:d)  ! Vector of exponentianilly distributed random numbers
integer :: j  ! Auxiliary variable for counters

allocate( rn(1:d) ) ;   call rng(d, rn)  ! Allocate memory for and get the uniformly distributed random numbers
 
! Gets the exponentianily distributed random numbers 
do j = 1, d ;   if ( rn(j) < 1.d-15 ) rn(j+1) = 1.d-15 ;   ern(j) = -log(rn(j)) ;   enddo

deallocate( rn )

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rng_unif(d, a, b, rn) ! Returns a vector of random (real) numbers with each component uniformly distributed in the interval [a,b]
! Ref: Katzgraber, H. G. (2010). Random numbers in scientific computing: An introduction, arXiv:1005.4117.
implicit none
integer :: d  ! Dimension of the vectors of random numbers
real(8) :: rn(1:d)  ! Vector of random numbers
real(8) :: a, b  ! New lower and upper limits to the domain of the random numbers

call rng(d, rn) ;   rn = a + (b-a)*rn  ! Gets rn in [0,1] and sets the new limits [a,b]

end
!###################################################################################################################################
!                                                  Calling subroutines - RNG
!###################################################################################################################################
subroutine rng_init() ! Calls the initializer subroutine for the choosed random number generator
use meths
implicit none
integer :: seed

     if ( opt_rng == "mt" ) then ;   call init_mt() ! Initializes the Mersenne-Twister rng
else if ( opt_rng == "gnu" ) then ;   call init_gnu() ! Initializes the Gnu's rng
else if ( opt_rng == "netlib" ) then ;   call init_netlib() ! Initializes the Netlib's rng
     end if

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rng(d, rn) ! Calls the choosed random number generator
use meths
implicit none
integer :: d  ! Dimension of the vector of random numbers
real(8) :: rn(1:d)  ! Vector of random numbers

     if (opt_rng == "mt") then ;   call rng_mt(d, rn)  ! Calls the Mersenne-Twister rng
else if (opt_rng == "gnu") then ;   call rng_gnu(d, rn)  ! ! Calls the Gnu's rng
else if (opt_rng == "netlib") then ;   call rng_netlib(d, rn)  ! Calls the Netlib's rng
     endif
  
end
!###################################################################################################################################