program getpz2

! Compile
! on my Mac: g7 getpz2 getpz2subs sac bin
! elsewhere:
! gfortran -o ~/bin/getpz2 getpz2.f90 getpzsubs.f sacio.a

! Given one of more input signals fin (in Pascal) and observed responses
! fobs (in Counts), find the Mermaid response in terms of poles and zeroes

! Use editcsv.f90 to create de-noised input file(s) from the Excel files
! provided by OSEAN

! Input:
! file in.getpz2 with input file names like test1.can 
! from screen: low pass frequency and commands to continue or stop

! format of input files is as written by programs editcsh and readcsh
! array fin has the input (in Pascal), fobs the output (in counts)
! of the tests.

! The program can also be run without using the sacio.a library
! if you comment all calls to routine wsac. You can plot the
! *.xy files with GMT or any other plotting program.

! From a digitized input signal fin(t), and a digitized output
! signal fobs(t), find optimal poles, zeroes and the amplification
! of an instrument. The only other input needed is an estimate of
! the high-pass frequency, which is used to start the iterations,
! a specification of the low-pass that defines the frequency band
! over which we fit in- and output, and a maximum allowable misfit 
! (in percent) to stop iterations.

! For a block-type signal, you can create a proper input file
! using makeio.f90

! fin and fobs must start at the same time and have equal sampling
! interval dt. To alleviate edge effects it is good to start and
! end fin and fobs with zero segments (no input, response quieted down).

! We assume that the response has at least one zero equal to 0 (i.e.
! there is zero amplification of a DC signal).

! The structure of the data files is:
! line0: comment header
! line1: n (number of data input)
! line2: dt (in seconds)
! lines3: n data (free format, units Pascal)

! Output:
! out.getpz2 with the result
! sacpz with a file to be used in sac (trans from polezero s sacpz to ...)
! inpz intermediate pz results for SAC
! diagnostics.getpz allows you to follow the optimization
! inputnn.sac SAC file of input nn (Pa, tapered but not scaled)
! outptnn.sac SAC file of the observed response to nn (energy scaled to 1)
! fout0.sac response predicted with the starting set of poles and zeroes
! foutk_nn.sac response to nn predicted after (outer) iteration rkn

implicit none

! we allow for up to NPAR (complex) poles and zeroes
integer, parameter :: ND=32768, M2=15, NPAR=40, NREC=100         ! ND=2**M2
integer, parameter :: NY=2*NPAR

character*80 :: fname(NREC),fnout(NREC),pzfile,header,line80,sacfile
character*1 :: h
character*8 :: pz,kstnm
character :: date*8, time*10
complex(kind=4) :: fcmplx(ND,NREC) ! FT of the input in Pa, with energy 1
complex(kind=4) :: ft(ND)          ! temporary fourier trf
complex(kind=4) :: pdflt(5)        ! default set of poles for startup
complex(kind=4) :: pole(NPAR)      ! response poles
complex(kind=4) :: r               ! response
complex(kind=4) :: zdflt(4)        ! default set of zeroes for startup`
complex(kind=4) :: zero(NPAR)      ! response zeroes
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
integer :: nmax                 ! length of largest fobs()
integer :: np,nq,nz             ! length of pole(),zero(),q()
integer :: npts(NREC)           ! length of each fin,fobs
integer :: nz0                  ! no of zero() being kept 0.
integer :: t0(6)                ! origin time (for SAC files)
real*4 :: alfa                  ! inverse of relaxation time (exp -alfa*t)
real*4 :: ampDC                 ! small amplitude considered close to DC
real*4 :: ampmax                ! maximum fobs()
real*4 :: ampmin                ! minimum fobs()
real*4 :: ampobs                ! scaling for fobs()
real*4 :: ampred                ! scaling for fout()
real*4 :: df                    ! FFT bin in Herz
real*4 :: dt                    ! sampling interval in seconds (Tektronics)
real*4 :: dtcan                 ! actual interval for the CAN (see below)
real*4 :: dum,x,yy,z
real*4 :: dw                    ! step in circle frequency w
real*4 :: f2                    ! Low pass frequency during optimization
real*4 :: fin(ND,NREC)          ! input in Pa, after scaling/tapering
real*4 :: finit                 ! initial misfit (value of func)
real*4 :: fl1,fl2               ! limits of the ideally flat response
real*4 :: fnyquist              ! Nyquist frequency in Hertz
real*4 :: fobs(ND,NREC)         ! output (scaled, tapered
real*4 :: fret                  ! value of func returned by Powel
real*4 :: ftol                  ! tolerance for convergence of Powel
real*4 :: func                  ! function func computes the misfit
real*4 :: pi=3.14159265, twopi=6.2831853
real*4 :: q(NY)                 ! real mapping of complex poles, zero
real*4 :: t                     ! time
real*4 :: thalf                 ! time between 0 and relaxation to max/2
real*4 :: tshift                ! time until first jump
real*4 :: w                     ! circle frequency
real*4 :: xi(NY,NY)             ! search direction in Powel's algorithm
real*4 :: Xre,Xim               ! real, imaginary component of pole, zero
real*4 :: y(NY)                 ! mapping of poles & zeroes to parameters y

data pdflt/(0.49343E-01,-0.53623E-02), (0.49343E-01, 0.53623E-02),  &
  (-0.72407, 0.0000), (-0.58644E-01,-0.15366E-03),  &
  (-0.58644E-01, 0.15366E-03)/
data zdflt/(0.55110E-01, 0.45550E-01), (0.55110E-01,-0.45550E-01), &
  (-0.23670E-01, 0.41947E-01), (-0.23670E-01,-0.41947E-01)/

common ampDC,ampred,dt,dw,fcmplx,fl1,fl2,fin,fobs,i2,ifl1,ifl2,jbug,&
  kplot,kpow,n2,ndim,nfile,np,npts,nq,nz,nz0,pole,q,zero,dtcan

! initialize
fl1=4.0                         ! flat spectrum start
fl2=8.0                         ! flat spectrum test limit
ftol=0.0001                     ! tolerance for convergence
open(3,file='out.getpz2',action='write')
open(4,file='diagnostic.getpz2',action='write')

! We assume the two columns in the CSV files are perfectly synchronic.
! However, if doubts exist, run getdeltat.f90 on one of the files to
! find the sampling interval of the CAN (dtcan) & uncomment these lines:
!print *,'We assume you have found the correct sampling interval for'
!print *,'the Mermaid card (CAN) using program getdeltat. If you'
!print *,'answer 0, it will assume CAN and Tektronics have the same dt.'
!print *,'Give the true sampling interval for the CAN (in s):'
!read(5,*) dtcan

! Read input signals
open(15,file='in.getpz2',action='read')
nfile=0
nmax=0
ampobs=0.
ampmax=-1.0e30
ampmin=+1.0e30
fobs=0.                 ! ensures zero padding whatever true length of fobs
write(3,'(a)') ' n    obs      dt tshift     ampmax     ampmin     ampobs'
do
  nfile=nfile+1
  if(nfile>NREC) stop 'Too many files, increase array sizes'
  read(15,*,iostat=ios) fname(nfile)
  if(is_iostat_end(ios) .or. fname(nfile)(1:4).eq.'stop') then
    nfile=nfile-1
    exit
  endif  
  print *,'Opening ',trim(fname(nfile))

  open(1,file=fname(nfile),iostat=ios,action='read')
  if(ios.ne.0) stop 'Cannot open input file'
  read(1,'(a)') header
  read(1,*) npts(nfile)
  if(2*n>ND) stop 'Signal too long'
  read(1,*) x
  if(nfile.eq.1) then
    dt=x
  else
    if(abs(dt-x)>0.0001) stop 'Sampling intervals not equal'
  endif
  if(abs(dtcan).le.0.) dtcan=dt
  do i=1,npts(nfile)
    read(1,*,iostat=ios) fin(i,nfile),fobs(i,nfile)
    if(ios.ne.0) stop 'Error reading input file'
  enddo  
  close(1)

  ! If necessary, we resample fobs do match the dt of the Tektronics 
  if(abs(dt-dtcan)>0.) then
    kpow=M2
    do while(2**kpow .ge. npts(nfile))
      kpow=kpow-1
    enddo
    kpow=kpow+1
    n2=2**kpow
    call resample(fobs,n2,kpow,dtcan,dt)
  endif  

  ! remove start of signal until first pulse
  m=npts(nfile)
  call shift0(fin(1,nfile),fobs(1,nfile),npts(nfile))
  tshift=(npts(nfile)-m)*dt
  m=npts(nfile)
  do i=1,m
    ampobs=ampobs+fobs(i,nfile)**2
  enddo
  write(3,'(i3,i6,2f8.3,2i11,e11.2,1x,a)') nfile, &
    m,dt,tshift,nint(maxval(fobs(1:m,nfile))), &
    nint(minval(fobs(1:m,nfile))),sqrt(ampobs),trim(fname(nfile))
  if(npts(nfile)<400) then
    print *,'Warning: file ',trim(fname(nfile)),' is too short and ignored'
    nfile=nfile-1
    cycle
  endif  
  nmax=max(nmax,npts(nfile))          ! find largest file size
enddo  
close(15)

! scale total fobs energy to 1. This scaling is needed to have some
! reference for mthe misfit level
ampobs=sqrt(ampobs)
write(3,'(a,e12.3)') 'Observed response fobs is scaled by ampobs= ',ampobs
fobs=fobs/ampobs
ampmax=maxval(fobs)
ampmin=minval(fobs)
ampDC=0.1*max(abs(ampmin),ampmax)
write(3,'(a,2g13.5)') 'After scaling ampmin,max are:',ampmin,ampmax

! find nearest power of 2 
kpow=M2
do while(2**kpow .ge. nmax)
  kpow=kpow-1
enddo
kpow=kpow+1
n2=2**kpow
df=1.0/(n2*dt)            ! frequency interval in Hz
dw=twopi*df

fnyquist=0.5/dt           ! Nyquist frequency
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
call setfhv('e',(n2-1)*dt,nerr)
call setfhv('delta',dt,nerr)
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
do n=1,nfile
  write(kstnm,'(a5,i2.2)') 'test_',n
  call setkhv('kstnm',kstnm,nerr)
  kstnm(1:5)='input'
  call setkhv('kstnm',kstnm,nerr)
  call wsac0(trim(kstnm)//'.sac',dum,fin(1,n),nerr)
  open(2,file=trim(kstnm)//'.xy')
  do i=1,npts(n)
    write(2,*) (i-1)*dt,fin(i,n)
  enddo
  close(2)
  kstnm(1:5)='outpt'
  call setkhv('kstnm',kstnm,nerr)
  call wsac0(trim(kstnm)//'.sac',dum,fobs(1,n),nerr)
  open(2,file=trim(kstnm)//'.xy')
  do i=1,npts(n)
    write(2,*) (i-1)*dt,fobs(i,n)
  enddo
  close(2)
enddo  

! we always assume the spectrum is flat at high frequency (but see the
! commented section below this segment)
! set limits of deviation from flat spectrum penalty
fl1=4.0                         ! lowest frequency of flat part
fl2=8.0                         ! and highest frequency
! Uncomment this if you want more flexibility on fl1,fl2:
! print *,'You may impose flat amplitude between fl1 and fl2 Hz'
! print *,'Give fl1 and fl2 (or 0,0):'
! read(5,*) fl1, fl2
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

print *,'There are three options for the starting set of poles and'
print *,'zeroes: when nothing else is known, give the time (in s) at which'
print *,'the response has relaxed to half its maximum (measured from t=0)'
print *,'If this is a Mermaid from a known series, you can answer'
print *,'<default> (without the <>) and it will start from the standard'
print *,'response. If you have a result from an earlier iteration saved'
print *,'in a file, give that file name. Getpz2 saves the latest result'
print *,'in file <inpz>.'

print *,'Give time to half relaxation, or file name, or <default>:'
read(5,'(a)') pzfile
read(pzfile,*,iostat=ios) thalf

if(ios.eq.0) then
  alfa=0.69/thalf
  nz0=1
  nz=2
  np=3
  ! start with one set of complex conjugates.
  pole(1)=cmplx(-alfa,0.)
  pole(2)=cmplx(0.05,0.05)
  pole(3)=cmplx(0.05,-0.05)
  zero(1)=cmplx(0.05,0.05)
  zero(2)=cmplx(0.05,-0.05)
  write(3,'(a,f8.2)') 'Initial set derived from thalf=',thalf
  write(4,'(a,f10.4)') 'alfa=',alfa
else if(pzfile.eq.'default') then
  nz0=1
  nz=4
  np=5
  pole(1:5)=pdflt
  zero(1:4)=zdflt
  write(3,'(a)') 'Initial set is the default set.'
else
  write(3,*) 'Reading starting values from pz file: ',trim(pzfile)
  open(2,file=pzfile,action='read')
  ! format should as in SAC: first POLES, then ZEROS, then CONSTANT
  ! set nz equal to *all* zeros and include the ones that are fixed (0,0).
  ! (all zeroes that have starting value (0,0) will be kept fixed).
  ! avoid any comments in the pzfile
  read(2,*,iostat=ios) pz,np
  if(ios.ne.0) stop 'Error in first line of pole-zero input'
  do i=1,np
    read(2,*,iostat=ios) Xre,Xim
    pole(i)=cmplx(Xre,Xim)
    if(ios.ne.0) stop 'Error in reading poles'
  enddo
  read(2,*,iostat=ios) pz,nz
  if(ios.ne.0) stop 'Error in reading ZEROS line'
  do i=1,nz
    read(2,*,iostat=ios) Xre,Xim
    zero(i)=cmplx(Xre,Xim)
    if(ios.ne.0) stop 'Error in reading zeroes'
  enddo
  ! remove any zeroes equal to (0,0) from the optimizable set
  k=nz
  nz0=0
  do i=1,k
    if(abs(zero(i)).le.0.) then         ! remove and move up the rest
      do j=i+1,nz
        zero(j-1)=zero(j)
      enddo  
      nz=nz-1
      nz0=nz0+1
    endif
  enddo  
  
  if(nz.ne.k) print *,nz0,' starting zeros are (0,0) & will be kept (0,0)'
  if(nz0.eq.0) print *,'WARNING! DC response is not guaranteed zero'
  close(2)              ! we ignore the CONSTANT
endif  
    
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

! Optional low pass filtering to exclude spurious noise at high frequency
print *,'Give low pass frequency, or 0 for no filter (e.g. 4 Hz):'
read(5,*) f2            ! f2=lowpass
i1=1
i2=n2
if(f2>0.) then
  if(f2>fnyquist) stop 'lowpass frequency exceeds Nyquist frequency'
  write(3,'(a,2f10.3)') 'Low pass frequency: ',f2
  ! taper 10% to f2
  i1=0.9*nint(f2/df)-1      
  if(i1<1) stop 'lowpass frequency too low'
  i2=nint(f2/df)
  yy=i2-i1+1
else
  write(3,*) 'No low pass applied'
endif  
flush(3)
  
! lowpass observed and input signals
if(f2>0.) then  
  do n=1,nfile
    ! lowpass fobs, using taper, use ft for fobs spectrum
    ft(1:n2)=fobs(1:n2,n)
    call clogp(kpow,ft,-1.0,dt)      ! compute FT
    do i=i1,i2
      x=float(i-i1)/yy
      z=0.5*(1.0+cos(pi*x))
      ft(i)=z*ft(i)
    enddo  
    ft(i2:n2)=0.            ! negative spectrum is added in ftinv
    call ftinv(ft,kpow,dt,fobs(1,n))      ! and back to time domain
    ! DEBUG
    write(sacfile,'(a,i1,a)') 'debugfobs',n,'.sac'
    call wsac0(sacfile,dum,fobs(1,n),nerr)
  
    ! do the same with fin, and from now on preserve fcmplx for each fin
    fcmplx(1:n2,n)=fin(1:n2,n)
    call clogp(kpow,fcmplx(1,n),-1.0,dt)      ! compute FT
    do i=i1,i2
      x=float(i-i1)/yy
      z=0.5*(1.0+cos(pi*x))
      fcmplx(i,n)=z*fcmplx(i,n)
    enddo  
    fcmplx(i2:n2,n)=0.            ! negative spectrum is added in ftinv
  enddo
  
else            ! if no filtering get FT of fin
  do n=1,nfile
    fcmplx(1:n2,n)=fin(1:n2,n)
    call clogp(kpow,fcmplx(1,n),-1.0,dt)      ! compute FT of fin
  enddo
endif  

! from here on fin and fobs contain the scaled (filtered) in- and output
! and fcmplx contains the FFT of fin

! Now start nonlinear search for poles and zeroes, starting with
! np poles and nz+nz0 zeroes. Note that the spectrum can only be flat
! at high frequency if np=nz+nz0
itry=0
kplot=0                 ! causes func to write fout0.sac for initial fit
kstnm(1:5)='iter0'
call setkhv('kstnm',kstnm,nerr)
call pz2y(pole,np,zero,nz,y,ndim,q,nq)
finit=func(y)
kstnm(1:5)='init_'
call setkhv('kstnm',kstnm,nerr)
kplot=-1                ! suppress plots while optimizing
write(6,'(a,g12.3)') 'Initial misfit:',finit
write(6,'(a,f6.3)') 'Relative misfit: ',1.0

open(7,file='inpz')     ! open file for resulting poles/zeroes

do

  itry=itry+1
  
  ! set up starting parameter vector (note: real input stays real)
  call pz2y(pole,np,zero,nz,y,ndim,q,nq)

  ! at this point y contains poles, zeroes (except zero(0)=0)
  
  ! set up the directions for routine powell which will search
  ! iteratively in these directions for a better y vector
  ! (starting with latest additions as first directions)
  xi=0.
  do i=1,ndim
    xi(i,i)=0.02
  enddo
  
  ! optimize poles and zeroes  (in y) using powell
  print *,'Now optimizing...have patience'
  print *,'Calling powell with ndim=',ndim
  write(4,'(6x,2a)') 'pflat       ysom       asom      worst  n        ', &
    'func  y'
  kplot=-1
  call powell(y,xi,ndim,NY,ftol,iter,fret)
  flush(4)
  
  kplot=min(9,itry)   ! SAC file for this iteration
  t=func(y)           ! compute penalty and plot results
  write(6,'(a,2g10.3)') 'Final misfit:',fret
  write(6,'(a,f6.3)') 'Relative misfit: ',fret/finit
  write(3,'(a,i4,e11.3,f9.3,e12.3)') 'Results from iter: ',itry,fret,  &
    fret/finit,ampred
  
  ! map back from y to poles/zeroes
  call y2pz(pole,np,zero,nz,y,ndim,q,nq)
  
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

  if(np+nz<NPAR-3) then
    print *,'Add more poles/zeroes (1) or stop (0)?'
    read(5,*)  kontinue
    if(kontinue.le.0) exit
  else
    exit
  endif  

!  np=np+1
!  nz=nz+1
!  pole(np)=cmplx(0.05,0.05)
!  zero(np)=cmplx(0.05,0.05)
!  np=np+1
!  nz=nz+1
!  pole(np)=cmplx(0.05,-0.05)
!  zero(nz)=cmplx(0.05,-0.05)

  ! add 2 poles and 2 zeros at the start of the parameter list
  do i=np,1,-1
    pole(i+2)=pole(i)
  enddo
  do i=nz,1,-1
    zero(i+2)=zero(i)
  enddo  
  pole(1)=cmplx(0.05,0.05)
  zero(1)=cmplx(0.05,0.05)
  pole(2)=cmplx(0.05,-0.05)
  zero(2)=cmplx(0.05,-0.05)
  np=np+2
  nz=nz+2

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

subroutine pz2y(p,np,z,nz,y,ndim,q,nq)

! maps magnification poles and zeroes into y via intermediate
! q (q remembers 0 imaginary parts of real poles of zeroes)
! q needs to be saved for the inverse mapping with y2pz

implicit none
integer, parameter :: NPAR=40
complex(kind=4), intent(in) :: p(NPAR), z(NPAR)
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
ndim=j

return 
end subroutine pz2y


subroutine y2pz(p,np,z,nz,y,ndim,q,nq)

! maps y back into q, poles and zeroes
! existing real/imaginary nature of pole/zero is used to determine whether
! they have a complex conjugate, so arrays pole/zero are both in- and output
! q is needed on input and may be changed on output.

implicit none
integer, parameter :: NPAR=40
real*4, intent(inout) :: q(2*NPAR)
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

! a full history of iterations is written to powell.xy

implicit none

integer, parameter :: ND=32768, M2=15, NPAR=40, NREC=100   ! ND=2**M2
real*4, intent(in) :: y(2*NPAR)
real*4 :: dt,dtcan,dw,twopi=6.283185,w,ampred,t,q(2*NPAR)
real*4 :: som,ysom,asom,aver
real*8 :: somx,somxy
complex(kind=4) :: pole(NPAR),zero(NPAR),fcmplx(ND,NREC),s(ND)
complex(kind=4) :: cc,r
real*4 :: dum,eps,fobs(ND,NREC),fin(ND,NREC),fout(ND,NREC),fl1,fl2,pflat
integer :: n,np,nz,nz0,i,i1,i2,j,kpow,n2,iter,itry,ndim,nq,nerr,kplot
integer ::  jbug,ifl1,ifl2,npts(NREC),nfile,nworst
real*4 :: ampDC,worst

common ampDC,ampred,dt,dw,fcmplx,fl1,fl2,fin,fobs,i2,ifl1,ifl2,jbug,&
  kplot,kpow,n2,ndim,nfile,np,npts,nq,nz,nz0,pole,q,zero,dtcan

eps=0.005              ! damping term for poles and zeroes (in y)
cc=cmplx(0.,1.)

! map back from y to poles/zeros since resp works with them
call y2pz(pole,np,zero,nz,y,ndim,q,nq)

! compute the penalty due to rms difference from flat response
pflat=0.
do i=ifl1,ifl2
  w=(i-1)*dw
  call resp(w,np,pole,nz,zero,nz0,r)
  pflat=pflat+(abs(r)-1.0)**2
enddo
pflat=0.1*sqrt(pflat/(ifl2-ifl1+1))

! create spectrum of the response to each fin, store in s() and compute the
! penalty to deviations from each fobs
! note that fcmplx has the spectrum of fin
somx=0.d0
somxy=0.d0
do n=1,nfile
  s(1:n2)=fcmplx(1:n2,n)
  w=0.
  do i=1,n2/2-1
    call resp(w,np,pole,nz,zero,nz0,r)
    s(i)=r*s(i)
    w=w+dw
  enddo
  s(n2/2)=0.
  call ftinv(s,kpow,dt,fout(1,n))

  ! find optimal amplification fobs/fout (y/x) with least squares:
  do i=1,npts(n)
    if(abs(fobs(i,n))<ampDC) cycle      ! do not include DC offsets
    somx=somx+fout(i,n)**2
    somxy=somxy+fobs(i,n)*fout(i,n)
  enddo
enddo  

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
  if(kplot.ge.0) call plts(n,fout(1,n),kplot,npts(n),dt)
  som=0.
  do i=1,npts(n)
    som=som+(fobs(i,n)-fout(i,n))**2
  enddo
  if(som>worst) then
    worst=som
    nworst=n
  endif  
  asom=asom+som
enddo

! damping penalty
ysom=0.
do i=1,ndim
  ysom=ysom+y(i)**2
enddo 
ysom=eps*ysom/ndim 

func=pflat+asom+ysom          ! add all penalty terms
if(kplot>0) write(4,'(a,i2)') 'func called only for plot=',kplot
write(4,'(4e11.3,i3,e12.4,40f7.3)') pflat,ysom,asom,worst,nworst,  &
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
  write(2,*) (i-1)*dt,fout(i)
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

subroutine shift0(fin,fobs,npts)

! shift zero time to onset of first pulse (first fin not 0)
integer, parameter :: ND=32768
real*4 :: fin(ND),fobs(ND)
integer :: i,j,npts,n

n=1
do while (abs(fin(n)).le.0.)
  n=n+1
  if(n>npts) stop 'Error in shift0'
enddo
if(n.le.1) return
fin(1)=fin(n)
fobs(1)=0.
n=n-1
do i=2,npts-n
  fin(i)=fin(i+n)
  fobs(i)=fobs(i+n)
enddo
npts=npts-n
return
end

subroutine resample(y,n,npow,dtold,dtnew)

real*4 :: y(n)
complex(kind=4) :: ft(n)

ft=y
call clogp(npow,ft,-1.0,dtold)
call ftinv(ft,npow,dtnew,y)

return
end
