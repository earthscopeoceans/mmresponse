program getdeltat

! Compile
! on my Mac: g7 getdeltat getpz2subs sac bin
! elsewhere:
! gfortran -o ~/bin/getdeltat getdeltat.f90 getpzsubs.f sacio.a

! Run this program on a long test signal to estimate the true sampling
! interval on the CAN before running getpz2, so that there is no
! misalignment at later time because of time drift

! Use editcsv.f90 to create de-noised input file from the Excel files
! provided by OSEAN

! Note: we assume the CSV file has the Tektronics sampling interval
! correct, while the CAN may be sampled with a (slightly) different dt.
! Routines pz2y and y2pz differ from their namesakes in getpz2.

! Input from screen: file name and estimate of time to half relaxation


! fin and fobs are read from the same can file and
! must start at the same time. Therefore do not cut data at the start, 
! since a slight difference in dt may then cause start times to differ.

! The structure of the data file is:
! line0: comment header
! line1: n (number of data input)
! line2: dttek (sampling of Tektronics, in seconds)
! lines3: n data (free format, units Pascal)

implicit none

! we allow for up to NPAR (complex) poles and zeroes
integer, parameter :: ND=32768, M2=15, NPAR=40         ! ND=2**M2
integer, parameter :: NY=2*NPAR

character*80 :: fname,fnout,pzfile,header,line80,sacfile
character*1 :: h
character*8 :: pz,kstnm
character :: date*8, time*10
complex(kind=4) :: fcmplx(ND)   ! FT of the input in Pa, scale to energy 1
complex(kind=4) :: ft(ND)               ! temporary fourier trf
complex(kind=4) :: pdflt(5)             ! default set of poles for startup
complex(kind=4) :: pole(NPAR)           ! response poles
complex(kind=4) :: r                    ! response
complex(kind=4) :: zdflt(4)             ! default set of zeroes for startup`
complex(kind=4) :: zero(NPAR)           ! response zeroes
integer :: i,ios,i1,i2,j,k,m,n,nerr
integer :: ifl1,ifl2            ! indices of flat spectrum limits
integer :: iter                 ! Number of iterations needed by Powel
integer :: itry                 ! counts trials with increasing nr of poles
integer :: jbug                 ! flags if bug output is wanted
integer :: jday                 ! Julian day (for sac files)
integer :: kontinue             ! 1 if continuing with more poles/zeroes
integer :: kplot                ! flags if plot files are made
integer :: kpow                 ! power of 2 for all FT's
integer :: n2                   ! length of fft (n2=2**kpow)
integer :: ndim                 ! nr of parameters druing trial itry
integer :: nfile                ! number of files with response data
integer :: np,nq,nz             ! length of pole(),zero(),q()
integer :: npts                 ! length of each fin,fobs
integer :: nz0                  ! no of zero() being kept 0.
integer :: t0(6)                ! origin time (for SAC files)
real*4 :: alfa                  ! inverse of relaxation time (exp -alfa*t)
real*4 :: ampDC                 ! small amplitude considered close to DC
real*4 :: ampmax                ! maximum fobs()
real*4 :: ampmin                ! minimum fobs()
real*4 :: ampobs                ! scaling for fobs()
real*4 :: ampred                ! scaling for fout()
real*4 :: df                    ! FFT bin in Herz
real*4 :: dtcan                 ! true sampling interval of the CAN
real*4 :: dttek                 ! sampling interval of the Tektronics
real*4 :: dum,x,yy,z
real*4 :: dw                    ! step in circle frequency w
real*4 :: f2                    ! Low pass frequency during optimization
real*4 :: fin(ND)               ! input in Pa, after scaling/tapering
real*4 :: finit                 ! initial misfit (value of func)
real*4 :: fl1,fl2               ! limits of the ideally flat response
real*4 :: fnyquist              ! Nyquist frequency in Hertz
real*4 :: fobs(ND)              ! output (scaled, tapered
real*4 :: fret                  ! value of func returned by Powel
real*4 :: ftol                  ! tolerance for convergence of Powel
real*4 :: func                  ! function func computes the misfit
real*4 :: pi=3.14159265, twopi=6.2831853
real*4 :: q(NY)                 ! real mapping of complex poles, zero
real*4 :: t                     ! time
real*4 :: thalf                 ! time between 0 and relaxation to max/2
real*4 :: w                     ! circle frequency
real*4 :: xi(NY,NY)             ! search direction in Powel's algorithm
real*4 :: Xre,Xim               ! real, imaginary component of pole, zero
real*4 :: y(NY)                 ! mapping of poles & zeroes to parameters y

data pdflt/(-0.525,0.),(-0.86,1.20),(-0.86,-1.20),(-0.71E-01,-0.35E-04), &
  (-0.71E-01,0.35E-04)/
data zdflt/(-1.68,0.165E-02),(-1.68,-0.165E-02),(-0.423E-01,0.),  &
  (-0.397E-01,0.)/

common ampDC,ampred,dttek,dw,fcmplx,fl1,fl2,fin,fobs,i2,ifl1,ifl2,jbug,&
  kplot,kpow,n2,ndim,nfile,np,npts,nq,nz,nz0,pole,q,zero,dtcan

! initialize
fl1=4.0                         ! flat spectrum start
fl2=8.0                         ! flat spectrum test limit
ftol=0.0001                     ! tolerance for convergence
open(3,file='out.getdeltat',action='write')
open(4,file='diagnostic.getdeltat',action='write')

ampobs=0.
ampmax=-1.0e30
ampmin=+1.0e30
fobs=0.
write(3,'(a)') '  n    in   obs      dt     ampmax     ampmin     ampobs'
nfile=1
print *,'Input signal file name (e.g. test1.can):'
read(5,*,iostat=ios) fname
print *,'Opening ',trim(fname)

open(1,file=fname,iostat=ios,action='read')
if(ios.ne.0) stop 'Cannot open input file'
read(1,'(a)') header
read(1,*) npts
if(2*n>ND) stop 'Signal too long'
read(1,*) dttek
write(13,*) 'dttek=',dttek
! initialize dtcan to the Tektronics sampling interval
dtcan=dttek
! read pressure (fin) and response (fobs)
fin=0.
fobs=0.
do i=1,npts
  read(1,*,iostat=ios) fin(i),fobs(i)
  if(ios.ne.0) stop 'Error reading input file'
enddo  
close(1)
write(13,*) 'fobs read as:'
do i=1,20
  write(13,*) i,fobs(i)
enddo  

! compute energy of response
do i=1,npts
  ampobs=ampobs+fobs(i)**2
enddo
write(13,*) 'ampobs=',ampobs
write(3,'(i3,i6,f8.3,2i11,e11.2,1x,a)') nfile, &
  npts,dttek,nint(maxval(fobs(1:npts))), &
  nint(minval(fobs(1:npts))),sqrt(ampobs),trim(fname)

! scale total fobs energy to 1. This scaling is needed to have some
! reference level for the misfit level
ampobs=sqrt(ampobs)
write(3,'(a,e12.3)') 'Observed response fobs is scaled by ampobs= ',ampobs
fobs=fobs/ampobs
ampmax=maxval(fobs)
ampmin=minval(fobs)
ampDC=0.1*max(abs(ampmin),ampmax)
write(3,'(a,2g13.5)') 'After scaling ampmin,max are:',ampmin,ampmax

! find nearest power of 2 
kpow=M2
do while(2**kpow .ge. npts)
  kpow=kpow-1
enddo
kpow=kpow+1
n2=2**kpow
df=1.0/(n2*dttek)            ! frequency interval in Hz
dw=twopi*df

fnyquist=0.5/dttek           ! Nyquist frequency
write(3,'(a,i10,a)') 'All signals are zero-padded to ',n2,' samples'
write(3,'(a,f10.6,a)') 'FFT frequency step is ',df,' Hz'
write(3,'(a,f10.3)') 'Nyquist frequency (Hz): ',fnyquist
flush(3)

! write SAC files of unfiltered signals (fin scaled to physical units)
call date_and_time(date,time)
read(date,'(i4,3i2)') (t0(i),i=1,3)
read(time,'(3i2)') (t0(i),i=4,6)
call newhdr
call setnhv('npts',n2,nerr)
call setfhv('b',0.,nerr)
call setfhv('e',(n2-1)*dttek,nerr)
call setfhv('delta',dttek,nerr)
call setnhv('nzyear',t0(1),nerr)
call julian(t0(1),t0(2),t0(3),jday)
call setnhv('nzjday',jday,nerr)
call setnhv('nzhour',t0(4),nerr)
call setnhv('nzmin',t0(5),nerr)
call setnhv('nzsec',t0(6),nerr)
call setnhv('nzmsec',0,nerr)
call setfhv('o',0.,nerr)
call setihv('iztype','IB',nerr)
! Plot files for fin and fobs (input and response after scaling)
! also output of input and observed signal to GMT plot format
write(kstnm,'(a5,i2.2)') 'test_',n
call setkhv('kstnm',kstnm,nerr)
kstnm(1:5)='input'
call setkhv('kstnm',kstnm,nerr)
call wsac0(trim(kstnm)//'.sac',dum,fin,nerr)
open(2,file=trim(kstnm)//'.xy')
do i=1,npts
  write(2,'(2f10.3)') (i-1)*dttek,fin(i)
enddo
close(2)
kstnm(1:5)='outpt'
call setkhv('kstnm',kstnm,nerr)
call wsac0(trim(kstnm)//'.sac',dum,fobs,nerr)
open(2,file=trim(kstnm)//'.xy')
do i=1,npts
  write(2,'(2f10.3)') (i-1)*dtcan,fobs(i)
enddo
close(2)

! we always assume the spectrum is flat at high frequency (but see the
! commented section below this segment)
! set limits of deviation from flat spectrum penalty
fl1=4.0                         ! lowest frequency of flat part
fl2=8.0                         ! and highest frequency
if(fl2>0.) then
  ifl1=fl1/df+1
  ifl2=fl2/df+1
else
  ifl1=n2
  ifl2=0
endif  
 
write(3,'(a,2f8.1)') 'We want a flat amplitude between (Hz):',fl1,fl2

fl1=twopi*fl1
fl2=twopi*fl2

print *,'Give time to half relaxation:'
read(5,*) thalf

! we work with 3 poles, 2 nonzero zeros and one zero kept (0,0).
alfa=0.69/thalf
nz0=1
nz=2
np=3
pole(1)=cmplx(-alfa,0.)
pole(2)=cmplx(0.05,0.50)
pole(3)=cmplx(0.05,-0.05)
zero(1)=cmplx(0.05,0.05)
zero(2)=cmplx(0.05,-0.05)
write(3,'(a,f8.2)') 'Initial set derived from thalf=',thalf
write(4,'(a,f10.4)') 'alfa=',alfa
    
write(3,'(/,a,i2,a)') 'Starting with ',np,' poles:'
write(3,'(i3,2f10.4)') (i,pole(i),i=1,np)
write(3,'(a,i2,a)') 'and ',nz,' zeros:'
write(3,'(i3,2f10.4)') (i,zero(i),i=1,nz)
write(3,'(i3,a)') nz0,' are being kept fixed to (0,0)'
if(np.ne.nz+nz0) then
  print *,'WARNING! np is not equal to nz+nz0, the response at high'
  print *,'frequency is therefore not flat'
  write(3,'(a)') 'WARNING! np is not equal to nz+nz0, the response at high'
  write(3,'(a)') 'frequency is therefore not flat'
endif  
flush(3)

! FFT using dttek  
do n=1,nfile
  fcmplx(1:n2)=fin(1:n2)
  call clogp(kpow,fcmplx,-1.0,dttek)      ! compute FT of fin (pressure)
enddo

! from here on fin and fobs contain the scaled (filtered) in- and output
! and fcmplx contains the FFT of fin

! Now start nonlinear search for poles and zeroes and optimal dtcan
itry=1
kplot=0                 ! causes func to write fout0.sac for initial fit
kstnm(1:5)='dtcan'
call setkhv('kstnm',kstnm,nerr)
call pz2y(pole,np,zero,nz,y,ndim,q,nq,dtcan)
if(ndim.ne.6) stop 'BUG ndim != 6'
finit=func(y)
kstnm(1:5)='init_'
call setkhv('kstnm',kstnm,nerr)
kplot=-1                ! suppress plots while optimizing
write(6,'(a,g12.3)') 'Initial misfit:',finit
write(6,'(a,f6.3)') 'Relative misfit: ',1.0

open(7,file='inpz')     ! open file for resulting poles/zeroes

! set up starting parameter vector 
call pz2y(pole,np,zero,nz,y,ndim,q,nq,dtcan)

! at this point y contains poles, zeroes (except zero(0)=0)

! set up the directions for routine powell which will search
! iteratively in these directions for a better y vector
! (starting with latest additions as first directions)
xi=0.
do i=1,ndim-1
  xi(i,i)=0.02
enddo
xi(ndim,ndim)=0.000002        ! sampling interval

! optimize poles and zeroes  (in y) using powell
print *,'Now optimizing...have patience'
print *,'Calling powell with ndim=',ndim
write(4,'(6x,2a)') 'pflat       ysom       asom      worst  n        ', &
  'func  y'
call powell(y,xi,ndim,NY,ftol,iter,fret)
flush(4)

kplot=min(9,itry)   ! SAC file for this iteration
t=func(y)           ! compute penalty and plot results
write(6,'(a,2g10.3)') 'Final misfit:',fret
write(6,'(a,f6.3)') 'Relative misfit: ',fret/finit
write(3,'(a,i4,e11.3,f9.3,e12.3)') 'Results from iter: ',itry,fret,  &
  fret/finit,ampred

! map back from y to poles/zeroes
call y2pz(pole,np,zero,nz,y,ndim,q,nq,dtcan)
write(6,'(a,f10.6,a,f10.6,a)') 'dtcan=',dtcan,' (cf dttek=',dttek,')'

! write inpz file with the new solution
write(7,'(/,a,i3)') '------------- iteration: ',itry
write(7,'(a,i5)') 'POLES ',np
do i=1,np
  write(7,'(2g14.5)') pole(i)
enddo  
write(7,'(a,i5)') 'ZEROS ',nz+nz0
do i=1,nz
  write(7,'(2g14.5)') zero(i)
enddo  
do i=1,nz0
  write(7,'(a)') '   0.   0.'
enddo  

write(3,'(/,a,f10.6)') 'Optimal sampling interval (s): ',dtcan

! The unscaled fobs is fobs(scaled)*ampobs, the unscaled fout is 
! fout(scaled)/ampred. The amplificationis thus fobs/fout=ampobs*ampred
write(7,'(a,g14.5)') 'CONSTANT',ampobs*ampred

! tell me what you have
write(3,'(/,a,i3,a)') 'Results of iteration ',itry,':'
write(3,'(a,e11.3,a,f9.3)') 'Absolute misfit:',fret,', Relative:',  &
  fret/finit
write(3,'(a,g14.5,a,e11.3)') 'A0=',ampobs*ampred,', scaling for fout:', &
  ampred
write(3,*) 'Poles:'
do i=1,np
  write(3,'(i3,2g14.5)') i,real(pole(i)),aimag(pole(i))
enddo  
write(3,*) 'Zeroes:'
if(nz>0) then
  do i=1,nz
    write(3,'(i3,2g14.5)') i,zero(i)
  enddo  
endif
do i=1,nz0
  write(3,'(i3,2g14.5)') i+nz,0.,0.
enddo

end

subroutine resp(omega,ip,p,iz,z,nz0,r)

!	gives response at complex frequency s=i*omega=i*2*pi*f
!       but ignoring amplification factor
!	ip poles stored in p, iz zeroes in z
!	output r is response
!       for complex conjugate pairs, only one is given in p,z
!       nz0 extra zeroes z=0 are added 

implicit none
integer, parameter :: NPAR=40
complex(kind=4), intent(in) :: p(NPAR),z(NPAR)
real*4, intent(in) :: omega
integer, intent(in) :: ip,iz,nz0
complex(kind=4), intent(out) :: r
real*4 :: eps=1.0e-20
complex(kind=4) :: s
integer :: j

s=cmplx(0.,omega)               ! in SEED, s = +i*omega

r = cmplx(1.,0.)                ! startup of response r

! multiply the factors in the nominator (zeros)
do j = 1,iz  
  r = r*(s-z(j))
enddo

! multiply by fixed 0's
do j=1,nz0
  r=r*s
enddo

! do the denominator (poles)
do j = 1,ip 
  r = r/(s-p(j))
enddo

return
end subroutine resp

! Subroutines pz2y and y2pz are needed because complex poles and
! zeroes exist in conjugate pairs, and only one of each pair needs
! to be optimized, while for real values only the real part of
! the (complex) variable needs to be optimized. The optimization
! therefore deals not with pole() and zero(), but with y().

subroutine pz2y(p,np,z,nz,y,ndim,q,nq,dtcan)

! maps magnification poles and zeroes into y via intermediate
! q (q remembers 0 imaginary parts of real poles of zeroes)
! q needs to be saved for the inverse mapping with y2pz

implicit none
integer, parameter :: NPAR=40
complex(kind=4), intent(in) :: p(NPAR), z(NPAR)
real*4, intent(in) :: dtcan
integer, intent(in) :: np,nz
real*4, intent(out) :: y(2*NPAR),q(2*NPAR)
integer, intent(out) :: ndim,nq
integer :: i,j


! first map p and z into q, ignoring conjugates
j=0
i=0
do while(i<np)       ! loop over poles
  i=i+1
  j=j+2
  q(j-1)=real(p(i))
  q(j)=aimag(p(i))
  if(abs(q(j))>0.) i=i+1        ! conjugate not needed in q
enddo  

i=0
do while(i<nz)       ! loop over zeroes
  i=i+1
  j=j+2
  q(j-1)=real(z(i))
  q(j)=aimag(z(i))
  if(abs(q(j))>0.) i=i+1        ! conjugate not needed in q
enddo  
nq=j
if(2*(nq/2).ne.nq) stop 'nq is odd'     ! DEBUG

! now map q into y, ignoring zero imaginary parts
j=0
do i=1,nq
  if(abs(q(i))>0.) then
    j=j+1
    y(j)=q(i)
  endif
enddo
y(j+1)=dtcan
ndim=j+1

return 
end subroutine pz2y


subroutine y2pz(p,np,z,nz,y,ndim,q,nq,dtcan)

! maps y back into q, poles and zeroes
! existing real/imaginary nature of pole/zero is used to determine whether
! they have a complex conjugate, so arrays pole/zero are both in- and output
! q is needed on input and may be changed on output.

implicit none
integer, parameter :: NPAR=40
real*4, intent(inout) :: q(2*NPAR),dtcan
complex(kind=4), intent(out) :: p(NPAR), z(NPAR)
integer, intent(in) :: np,nz
real*4, intent(in) :: y(2*NPAR)
integer, intent(in) :: ndim
integer, intent(out) :: nq
integer :: i,j,k

! map back from y to poles/zeroes
j=0
do i=1,nq,2
  j=j+1
  q(i)=y(j)
  if(abs(q(i+1))>0.) then       ! change nonzero q as well (=imaginary)
    j=j+1
    q(i+1)=y(j)
  endif  
enddo
dtcan=y(j+1)

j=0
k=0
do i=1,nq,2
  if(j<np) then
    j=j+1
    p(j)=cmplx(q(i),q(i+1))
    if(abs(q(i+1))>0.) then
      j=j+1
      p(j)=cmplx(q(i),-q(i+1))
    endif
  else  
    k=k+1
    z(k)=cmplx(q(i),q(i+1))
    if(abs(q(i+1))>0.) then
      k=k+1
      z(k)=cmplx(q(i),-q(i+1))
    endif
  endif
enddo

return
end subroutine y2pz

real function func(y)

! given parameter vector y with np poles and nz zeroes, compute the
! instrument response and set func equal to quadratic misfit with data d()

! data etc are passed on via common

! Note: this func differs from the one in getpz2.f90

implicit none

integer, parameter :: ND=32768, M2=15, NPAR=40    ! ND=2**M2
real*4, intent(in) :: y(2*NPAR)
real*4 :: dttek,dw,twopi=6.283185,w,ampred,t,som,q(2*NPAR),ysom,asom,aver
real*8 :: somx,somxy
complex(kind=4) :: pole(NPAR),zero(NPAR),fcmplx(ND),s(ND)
complex(kind=4) :: cc,r
real*4 :: dtcan,dum,eps,fobs(ND),fin(ND),fout(ND),fl1,fl2,pflat
integer :: n,np,nz,nz0,i,i1,i2,j,kpow,n2,iter,itry,ndim,nq,nerr,kplot
integer ::  jbug,ifl1,ifl2,npts,nfile,nworst
real*4 :: ampDC,worst

common ampDC,ampred,dttek,dw,fcmplx,fl1,fl2,fin,fobs,i2,ifl1,ifl2,jbug,&
  kplot,kpow,n2,ndim,nfile,np,npts,nq,nz,nz0,pole,q,zero,dtcan

eps=0.005              ! damping term for poles and zeroes (in y)
cc=cmplx(0.,1.)

! map back from y to poles/zeros since resp works with them
call y2pz(pole,np,zero,nz,y,ndim,q,nq,dtcan)
write(13,*) 'in func with ndim, dtcan =',ndim,dtcan
write(13,'(a,5f7.4,f10.6)') 'y=',(y(i),i=1,6)
write(13,'(a,6f7.4)') 'z=',(zero(i),i=1,3)
write(13,'(a,6f7.4)') 'p=',(pole(i),i=1,3)

! compute the penalty due to rms difference from flat response
pflat=0.
do i=ifl1,ifl2
  w=(i-1)*dw
  call resp(w,np,pole,nz,zero,nz0,r)
  pflat=pflat+(abs(r)-1.0)**2
enddo
write(13,*) 'last resp=',r
pflat=0.1*sqrt(pflat/(ifl2-ifl1+1))

! create spectrum of the response to each fin, store in s() and compute the
! penalty to deviations from each fobs
! note that fcmplx has the spectrum of fin
somx=0.d0
somxy=0.d0

s(1:n2)=fcmplx(1:n2)
w=0.
do i=1,n2/2-1
  call resp(w,np,pole,nz,zero,nz0,r)
  s(i)=r*s(i)
  w=w+dw
enddo
s(n2/2)=0.
! use dtcan as sampling interval for inverse FFT, so that the theoretical
!response is synchronized with the observed fobs, which has dtcan
call ftinv(s,kpow,dtcan,fout)

! find optimal amplification fobs/fout (y/x) with least squares:
do i=1,npts
  if(abs(fobs(i))<ampDC) cycle      ! do not include DC offsets
  somx=somx+fout(i)**2
  somxy=somxy+fobs(i)*fout(i)
enddo

write(13,*) 'somx,somxy=',somx,somxy
do i=1,30
  write(13,*) i,fout(i),fobs(i)
enddo
flush(13)

! we minimize |fobs - ampred*fout|^2 before computing waveform misfit
ampred=somxy/somx
fout=fout*ampred       ! scales energy of fout to 1
if(kplot.eq.0) write(3,'(2a,e12.3,a)') 'Initial response fout is', &
  ' multiplied by ',ampred, ' to minimize misfit'

! compute misfit between observed and predicted signals
som=0.                  ! misfit for individual signals
asom=0.                 ! total for al n signals
worst=0.
do n=1,nfile
  if(kplot.ge.0) call plts(n,fout,kplot,npts,dtcan)
  som=0.
  do i=1,npts
    som=som+(fobs(i)-fout(i))**2
  enddo
  if(som>worst) then
    worst=som
    nworst=n
  endif  
  asom=asom+som
enddo

! damping penalty
ysom=0.
do i=1,ndim-1
  ysom=ysom+y(i)**2
enddo 
ysom=eps*ysom/(ndim-1)

func=pflat+asom+ysom+abs(dttek-dtcan)
if(kplot>0) write(4,'(a,i2)') 'func called only for plot=',kplot
write(4,'(4e11.3,i3,f12.6,5f7.3,f10.6)') pflat,ysom,asom,worst,nworst,  &
  func,(y(i),i=1,ndim)

if(isnan(func)) stop 'ERROR: func returns NaN for misfit'

return
end function func

subroutine plts(n,fout,kplot,npts,dt)

! writes pole-zero response (fout) in both SAC and GMT formats

implicit none
integer, intent(in) :: n,kplot,npts
real*4, intent(in) :: fout(npts),dt
real*4 :: dum
integer :: i,nerr
character*40 :: fname
character*8 :: kstnm

write(fname,'(a4,i1,a1,i2.2,a4)') 'fout',kplot,'_',n,'.sac'
write(kstnm,'(a4,i1,a1,i2.2)') 'iter',kplot,'_',n
call setnhv('npts',npts,nerr)
call setfhv('b',0.,nerr)
call setfhv('e',(npts-1)*dt,nerr)
call setfhv('delta',dt,nerr)
call setfhv('o',0.,nerr)
call setihv('iztype','IB',nerr)
call setkhv('kstnm',kstnm,nerr)
call wsac0(fname,dum,fout,nerr)

! also output of predicted signal to GMT
fname=fname(1:9)//'xy'
open(2,file=fname)
do i=1,npts
  write(2,'(2f10.3)') (i-1)*dt,fout(i)
enddo
close(2)

return
end

subroutine julian(year,month,day,jday)

! returns Julian day if called with year, month, day

implicit none
integer :: year,month,day,jday

integer :: mday(12),kday(12),i,k,m,n,leap
data mday/31,28,31,30,31,30,31,31,30,31,30,31/

leap=0
mday(2)=28
if(mod(year,4).eq.0.and.(mod(year,100).ne.0.or.mod(year,400).eq.0)) leap=1
if(leap.eq.1) mday(2)=29
kday(1)=0
do i=2,12
  kday(i)=kday(i-1)+mday(i-1)
enddo
jday=kday(month)+day

end

