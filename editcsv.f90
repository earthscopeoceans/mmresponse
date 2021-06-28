program editcsv

! Compile: g7 editcsv bin
! or gfortran -o editcsv editcsv.f90

! This program is used to edit the data obtained with instrument repsonse
! tests in the basin of OSEAN - in which the hydrophone is lowered by
! a few cm (variable: depth)
! For input, transform Excel files into CSV format with four columns: 
! time (s), Potentiometer (V), hydrophone (V) and CAN (counts)
! Its skips the first (header) linea in the csv file.

! It writes input for program getpz2. 

! Usage: editcsv filename depth

! Note: the CSV files received from OSEAN have the fourth column (CAN)
! ahead by one sample with respect to the second column (potentiometer).
! This is corrected by this program immediately after reading the CSV file.

implicit none

integer, parameter :: ND=8192           ! max signal length 2^13 samples

real*4 :: average
real*4, allocatable, dimension(:,:) :: a        ! matrix for LSQR
real*4 :: can(ND)                       ! CAN reading in counts
real*4 :: can0(ND)                      ! CAN reading before editing
real*4 :: density=997.77                ! water density at room temperature
real*4 :: depth                         ! depth range of hydrophone moves
real*4 :: dt,dt1,dt2,halfdt             ! sampling interval, dt/2
real*4 :: dum
real*4 :: dw                            ! frequency step
real*4 :: dy(ND),dy2(ND)                ! 1st and 2nd derivates of can
real*4 :: fmax=2.0                      ! Hz, max frequency used in fitting
real*4 :: g=9.82                        ! local acceleration of gravity
real*4 :: hyd(ND)                       ! hydrophone voltage
integer :: itmax=100                    ! max nr of iterations of LSQR
integer :: histogram(200)               ! potentiometer voltage histogram
integer :: jmp                          ! numer of sudden signal changes
integer :: kbreak(10)                   ! points where to split the signal
integer :: kjmp(200),kjmp1(200)         ! jump indices
integer :: nf                           ! nr of frequencies used in fitting
integer :: nsmp                         ! length of total signal
integer :: m                            ! number of coefficients
real*4 :: pot(ND)                       ! potentiometer V, converted to Pa
real*4 :: pDC                           ! assumed rest pressure
real*4 :: pmid,pmin,pmax                ! mid,min,max pressure or Voltage
real*4 :: PperV                         ! Pa per Volt conversion factor
real*4 :: r                             ! residual of fit by LSQR
integer :: relax=160                    ! nr of samples for relaxation 
real*4 :: som,som1,som2                 ! used for summing/averaging
real*4 :: time(ND)                      ! time axis for can,pot,hyd
real*4 :: tjmp(200)                     ! time of jumps (start of moves)
real*4 :: t                             ! time variable
real*4 :: tlen                          ! length of segment to edit in sec
real*4 :: twopi=6.283185
real*4 :: ty(ND)                        ! time axis for y()
real*4 :: w                             ! circle frequency
real*4 :: x(257)                        ! fitting coefficients
real*4 :: xlim                          ! noise level estimate
real*4 :: ylow,yhigh,diff               ! limits of can in segments
real*4 :: y(ND)                         ! edited signal segments
character*1 :: plot(100)                ! for editable print plot
character*1 :: flag,q                   ! for editing
character*5 :: ctime
character*80 :: fname
character*124 :: line
integer :: i,i1,i2,ib,ios,iout,j,k,k0,n,nmin,nmax,ny,npow

open(4,file='out.editcsv',action='write')
open(10,file='diagnostics.editcsv')
n=command_argument_count()
if(n<2) then
  print *, "Usage: editcsv file.csv depth (in cm)"
  stop
endif  
call get_command_argument(1,fname)
open(1,file=fname,action='read')
read(1,*)                       ! skip header
write(4,'(2a)') 'CSV file: ',trim(fname)
call get_command_argument(2,line)
read(line,*) depth
write(4,'(a,f8.2,a)') 'Depth range: ',depth,' cm'
depth=depth*0.01                ! convert to m

j=0
k=0
pmax=-1e20
pmin=1e20
histogram=0
i=0
do
  i=i+1
  read(1,*,iostat=ios) time(i),pot(i),hyd(i),can(i)
  if(is_iostat_end(ios)) exit
  if(i.eq.2) dt1=time(2)-time(1)
  if(i.eq.3) dt2=time(3)-time(2)
  if(i.eq.3) write(10,'(a,2f10.4)') 'dt1, dt2=',dt1,dt2
  if(i.eq.4) then               ! get first estimate of dt
    if(abs(dt1-dt2)>0.002) stop 'Check on dt in first 3 samples'
    dt=0.5*(dt1+dt2)
    dt1=0.
  endif  
  ! repair missing data
  if(i>4) dt1=time(i)-time(i-1)-dt
  if (i>3.and.dt1>0.002) then
    j=nint(time(i)-time(i-1))/dt
    if(j>1) stop 'Too many consecutive samples are missing'
    time(i+1)=time(i)
    pot(i+1)=pot(i)
    hyd(i+1)=hyd(i)
    can(i+1)=can(i)
    time(i)=0.5*(time(i-1)+time(i+1))
    pot(i)=0.5*(pot(i-1)+pot(i+1))
    hyd(i)=0.5*(hyd(i-1)+hyd(i+1))
    can(i)=0.5*(can(i-1)+can(i+1))
    k=k+1
    i=i+1
    dt1=0.
    if(k.eq.1) write(10,'(a,/,a)') 'Interpolated missing data',  &
      '    i         t        can'
    write(10,'(i5,f10.3,i11)') i,time(i),nint(can(i))
  endif  
  pmax=max(pot(i),pmax)
  pmin=min(pot(i),pmin)
  j=nint((pot(i)+4.0)*25.0)
  histogram(j)=histogram(j)+1
enddo
close(1)
nsmp=i-1

! the CAN runs one sample ahead of the pot, we throw out that first sample
do i=1,nsmp-1
  can(i)=can(i+1)
enddo
nsmp=nsmp-1  

can0=can                ! save unedited version for plots

! find upper- and lower constant levels for potentiometer
pmid=0.5*(pmax+pmin)
write(10,*) 'min, max pot=',pmin,pmax
nmin=1
nmax=1
do i=1,nsmp
  j=nint((pot(i)+4.0)*25.0)
  if(pot(i)<pmid) then
    if(histogram(j)>histogram(nmin)) nmin=j
  else
    if(histogram(j)>histogram(nmax)) nmax=j
  endif
enddo  
pmax=nmax/25.-4.0
pmin=nmin/25.-4.0
write(10,'(a,2f8.2)') 'Histogram extrema: ',pmin,pmax

! find average on both upper and lower plateau
pmid=0.
nmax=0
nmin=0
som1=0.
som2=0.
do i=1,nsmp
  if(abs(pot(i)-pmax)<0.1) then
    nmax=nmax+1
    som2=som2+pot(i)
  else if(abs(pot(i)-pmin)<0.1) then
    nmin=nmin+1
    som1=som1+pot(i)
  endif
enddo
pmax=som2/nmax          ! average high constant level
pmin=som1/nmin          ! average low constant level
write(10,'(a,2f8.2)') 'Average constant voltage at ',pmin,pmax

! remove noise
! we try to detect spikes of amplitude up to 0.5V
do i=2,nsmp-1
  i1=max(1,i-20)
  i2=min(nsmp,i+20)
  ! do nothing if we are near jump
  if(abs(pot(i1)-pot(i2))>1.5) cycle
  if(abs(pot(i-1)-pmax)<0.5 .and. abs(pot(i)-pmax)<0.5 .and.  &
     abs(pot(i+1)-pmax)<0.5) pot(i)=pmax
  if(abs(pot(i-1)-pmin)<0.5 .and. abs(pot(i)-pmin)<0.5 .and.  &
     abs(pot(i+1)-pmin)<0.5) pot(i)=pmin
enddo
pDC=pmin
if(abs(pot(1)-pmax)<abs(pot(1)-pmin)) pDC=pmax
pot(1)=pDC
pot(nsmp)=pDC

PperV=depth*g*density/(pmax-pmin)
write(10,'(a,f12.2)') 'Pascal per Volt =',PperV
do i=1,nsmp
  dum=pot(i)
  pot(i)=(pot(i)-pDC)*PperV            ! Convert V to Pa
enddo  
dt=(time(nsmp)-time(1))/(nsmp-1.0)
halfdt=0.5*dt

write(4,'(a,i7,a)') 'The CVS file has ',nsmp,' samples'
write(4,'(a,f8.4,a)') 'Sampling interval:',dt,' sec'
write(4,'(a,f8.1,a,f8.1)') 'Times range from ',time(1),' to ',time(nsmp)
if(k>0) write(4,'(a,i5,a)') 'WARNING: the file missed ',k,' samples'
if(k>0) write(6,'(a,i5,a,/,a)') 'WARNING: the file missed ',k,  &
  ' samples','Missing data have been linearly interpolated, see diagnostics'

! get first and second derivatives of can
dy(1)=can(2)-can(i)
do i=2,nsmp-1
  dy(i)=can(i+1)-can(i-1)
enddo
dy(nsmp)=can(nsmp)-can(nsmp-1)

dy2(1)=dy(2)-dy(1)
do i=2,nsmp-1
  dy2(i)=dy(i+1)-dy(i-1)
enddo
dy2(nsmp)=dy(nsmp)-dy(nsmp-1)

! write gmt files
open(2,file='dy.xy',action='write')
open(3,file='dy2.xy',action='write')
open(7,file='can0.xy',action='write')
open(8,file='Pa.xy',action='write')
open(9,file='hyd.xy',action='write')
do i=1,nsmp
  write(2,'(f10.3,e12.3)') time(i),dy(i)
  write(3,'(f10.3,e12.3)') time(i),dy2(i)
  write(7,'(f10.3,i11)') time(i),nint(can0(i))
  write(8,'(f10.3,f10.1)') time(i),pot(i)
  write(9,'(f10.3,f10.3)') time(i),hyd(i)
enddo  
close(2)
close(3)
close(7)
close(8)
close(9)

! find jumps, assume approximate relaxation in relax(=160) samples to
! define the minimum length of segment to edit
! jumps may represent hydrophone up/down moves, but also noise bursts
jmp=0
i=2
xlim=0.                 ! noise level estimate from first few samples
write(10,'(a,/,a)') 'Jumps','jmp  kjmp      tjmp        can        dy2'
do while(i<nsmp-relax)  
  i=i+1
  if(jmp.eq.0) xlim=xlim+abs(can(i))
  if(abs(can(i))>1.0e7) then            ! if amplitude indicated a jump
    jmp=jmp+1
    if(jmp>200) stop 'jmp>200:increase max nr of jumps'
    j=max(1,i-4)
    kjmp(jmp)=j
    if(jmp.eq.1) xlim=xlim/(i-2)
    tjmp(jmp)=time(j)
    write(10,'(i3,i6,f10.3,2i11)') jmp,i,time(i),nint(can(i)),nint(dy2(i))
    i=min(nsmp,i+relax)
    do while (abs(can(i))>xlim.or.abs(dy2(i))>1.0e7) 
      i=i+1
      if(i.ge.nsmp) exit
    enddo
    kjmp1(jmp)=i
  endif  
enddo  

write(4,'(/,i3,a,e12.3)') jmp,' jumps exceeding ',xlim
write(4,'(a)') ' i       t    k1    k2'
do i=1,jmp
  write(4,'(i2,f8.3,2i6)') i,tjmp(i),kjmp(i),kjmp1(i)
enddo  

! Now write the segments that may need editing to file out.editcsv
write(4,'(/,a)') 'Edit segments below where needed. # identifies outliers'
write(4,'(a)') 'or clips. A dot signals a probably OK value (scaled by 1e7)'
write(4,'(a)') 'Replace # or . by R to require rejection & interpolation'
write(4,'(a)') 'Replace # or . by S to split the data at this point so that'
write(4,'(a,/)') 'a segment can be rejected. Otherwise leave # or . intact'

do j=1,jmp
  ylow=1.0e20
  yhigh=-1.0e20
  do i=kjmp(j),kjmp1(j)
    ylow=min(ylow,can(i))
    yhigh=max(yhigh,can(i))
  enddo  
  write(4,'(a,i3,a,2e12.3,a)') 'jump ',j,'  ylow,yhigh=',ylow*1.0e-7,  &
    yhigh*1.0e-7,'(x 1.0e7)'
  diff=99./(yhigh-ylow)
  k0=(0-ylow)*diff+1
  do i=kjmp(j),kjmp1(j)
    plot=' '
    if(mod(i,10).eq.0) plot='-'
    k=(can(i)-ylow)*diff+1
    plot(k)='*'
    if(k0.ne.k.and.k0>1.and.k0<100) plot(k0)=':'
    flag='.'
    if(abs(dy2(i))>1.0e7) flag='#'
    if(abs(can(i))>1.34e8) flag='#'     ! clipping
    q='|'
    if(mod(i,5).eq.0) q='+'
    ctime=''
    if(abs(time(i)-nint(10*time(i))*0.1)<halfdt) write(ctime,'(f5.1)') &
      time(i)
    write(4,'(i5,f8.3,1x,a5,2x,a1,a2,101a1)') i,can(i)*1.0e-7,ctime, &
      flag,q,(plot(k),k=1,100),q
  enddo
enddo  

! close out.editcsv and call vi to edit it
close(4)

call system('vi out.editcsv')

! reopen the edited file and read it
open(4,file='out.editcsv',action='read')

! skip the headers
line=''
do while(line(1:4).ne.'jump')
  read(4,'(a124)',iostat=ios) line
enddo  

! now loop over segments
ib=1
kbreak(1)=1
do
  if(is_iostat_end(ios).or.len_trim(line).eq.0) exit
  if(line(1:5).ne.'break') read(line,*,iostat=ios) ctime,jmp
  if(ios.ne.0) then
    write(6,*) 'ios error  in line',ios
    write(6,*) trim(line)
    stop
  endif  

  j=0           ! counts CAN readings in this segment
  do 
    read(4,'(a124)',iostat=ios) line
    if(line(1:4).eq.'jump' .or. is_iostat_end(ios)) exit
    if(len_trim(line).eq.0) exit
    read(line,'(i5,16x,a1)',iostat=ios) i,flag
    if(i.le.0.or.ios.ne.0) then
      write(6,*) 'i or ios error ',i,ios
      write(6,*) trim(line)
      stop
    endif  
    if(j.eq.0) i1=i             ! start of segment to interpolate
    if(flag.eq.'.' .or. flag.eq.'#') then       ! accept
      j=j+1
      y(j)=can(i)
      ty(j)=time(i)
    else if(flag.eq.'R') then   ! mark datum for editing
      can(i)=1.1e30
    else if(flag.eq.'S') then
      line(1:5)='break'
      exit
    endif
  enddo
  if(line(1:5).eq.'break') then
    ib=ib+1
    if(ib>9) stop 'Too many splits. Increase dimension of kbreak'
    kbreak(ib)=i
    write(10,'(a,i7)') 'Break at sample: ',i
  endif
  n=j           ! nr of samples of segment
  i2=i          ! end of segment to interpolate
  
  ! interpolate y by fitting harmonics < fmax Hz
  tlen=(i2-i1)*dt
  if(tlen.le.0.) then
    print *,'Failure in harmonic interpolation between i1,i2=',i1,i2
    print *,'j,y(j),ty(j)=',j,y(j),ty(j)
    print *,'Last line read: ',trim(line)
  endif
  dw=twopi/tlen
  nf=fmax*tlen+1
  m=2*nf+1              ! number of coefficients
  write(10,'(a,f8.0,2i8)') 'tlen, nf, m=',tlen,nf,m
  if(m>257) stop 'Segment too long. Split or increase dimension of array x'

  ! construct matrix
  allocate(a(n,m))
  do i=1,n
    t=ty(i)
    w=0.
    j=1
    a(i,1)=1.0          ! zero frequency
    do j=2,m-1,2
      w=w+dw
      a(i,j)=cos(w*t)
      a(i,j+1)=sin(w*t)
    enddo  
  enddo  
  ! solve A*x=y for coefficients x
  ! if no pick-up, this is an orthogonal systemand and converges in 1 
  itmax=100               
  iout=10               
  write(10,*) 'LSQR:'
  call flsqr(a,n,n,m,x,y,itmax,iout,r)    
  deallocate(a)

  ! identify outliers that need editing
  write(10,'(a,i5,a,i5)') 'Editing from',i1,' to ',i2
  write(10,'(a)') '     i         t         can0          can'
  do i=i1,i2
    if(can(i)<1.0e30) cycle
    som=x(1)
    j=1
    w=0.
    t=time(i)
    do k=2,m,2
      w=w+dw
      som=som+x(k)*cos(w*t)+x(k+1)*sin(w*t)
    enddo
    can(i)=som          ! replace can(i) by sum of cos,sin
    write(10,'(i6,f10.3,2f13.0)') i,t,can0(i),can(i)
  enddo
enddo

! write edited file to can.xy
open(7,file='can.xy',action='write')
do i=1,nsmp
  write(7,'(f10.3,i11)') time(i),nint(can(i))
enddo
close(7)

! Now write the (edited) can(s) in format suitable for getpz2
ib=ib+1
kbreak(ib)=nsmp
k=index(fname,'.csv')
if(k>0) then
  fname(k:k+3)='.can'
endif
k=len_trim(fname)
if(ib.eq.2) then
  write(10,*) 'No breaks detected'
else    
  write(10,'(a,/,a)') 'Split signal','ib   j-1     j        n'
endif

! Write input for getpz2. Split if there are breaks (ib>2).
do j=2,ib
  n=kbreak(j)-kbreak(j-1)+1
  if(ib.eq.2) then
    open(2,file=trim(fname),action='write')
    write(2,'(2a)') 'Edited hydrophone data from ',trim(fname)
  else
    write(fname(k+1:k+1),'(i1)') j-1
    open(2,file=trim(fname),action='write')
    write(2,'(3a,i2)') 'Edited hydrophone data from ',trim(fname),  &
      ', split',j-1
    write(10,'(i2,3i6)') j-1,kbreak(j-1),kbreak(j),n
  endif
  write(2,*) n
  write(2,*) dt
  do i=kbreak(j-1),kbreak(j)
    write(2,'(f10.1,i11)') pot(i),nint(can(i))
  enddo  
  close(2)
enddo  

write(6,'(a)') 'Succesful end of editcsv.'
write(6,'(a)') 'Plot interpolations with gmtcan, check results'
write(6,'(a)') '  in file diagnostics.editcsh.'
if(ib.eq.2) then
  write(6,'(2a)') 'input getpz2: ',trim(fname)
else
  write(6,'(3a)') 'input getpz2: ',trim(fname),' etc.'
endif  


end

! SUBROUTINES

subroutine flsqr(a,nm,n,m,x,rhs,itmax,iout,r)
 
!   subroutine to solve the linear tomographic problem Ax=u using the
!   lsqr algorithm. As lsqr but for full matrix A, and rhs is preserved.
 
!   reference: C.C.Paige and M.A.Saunders, ACM Trans.Math.Softw. 8, 
!       43-71, 1982
!   and ACM Trans.Math.Softw. 8, 195-209, 1982. 
!   See also A.v.d.Sluis and H.v.d.Vorst in: G. Nolet (ed.), Seismic 
!       Tomography, Reidel, 1987.
 
!   Input: 
!   a(nm,m) is the matrix
!   n is the number of data (rows), 
!   m the number of unknowns (colunms),
!   rhs contains the data, 
!   itmax is the number of iterations allowed
!   iout is the output channel (set iout=6 for screen output, 0 for none)

!   Output: x is the solution

!   If you set itmax very large, the program will quit as soon as the
!   misfit improves less than 1.0e-8 (relative) between two iterations.

!   Scratch (allocatable arrays): arrays v(m) and w(m)

!   Subroutines used (and included with the source):
!   avpu computes u=u+A*v for given input u,v (overwrites u)
!   atupv computes v=v+A(transpose)*u for given u,v )(overwrites v)
!   normlz normalizes the vector u or v
 
implicit none

integer, intent(in) :: nm,n,m           ! matrix size, rows, columns
integer, intent(in) :: iout             ! output unit (or 0 for none)
integer, intent(in) :: itmax            ! max nr of iterations
real*4, intent(in) :: a(nm,m)           ! matrix
real*4, intent(in) :: rhs(n)            ! right hand side (data)
real*4, intent(out) :: x(m)             ! solution
real*4, intent(out) :: r                ! relative misfit

real*4, allocatable, dimension(:) :: u,v,w
real*4 :: alfa,beta,b1,c,rho,rhobar,phi,phibar,rold,s,t1,t2,teta
integer :: iter

allocate(u(n))
allocate(v(m))
allocate(w(m))
x=0.
v=0.
u=rhs

! initialize
call normlz(n,u,beta)
b1=beta
call atupv(a,nm,n,m,u,v)
call normlz(m,v,alfa)
rhobar=alfa
phibar=beta
w=v
if(iout>0) write(iout,fmt='(a,/,i5,2g12.4,f10.6)')   &
  ' iter        x(1)        beta         r',0,x(1),beta,1.0
rold=1.0e30
do iter =1,itmax  
! repeat for fixed nr of iterations
  u=-alfa*u
! bidiagonalization
  call avpu(a,nm,n,m,u,v)
  call normlz(n,u,beta)
  v=-beta*v
  call atupv(a,nm,n,m,u,v)
  call normlz(m,v,alfa)
  rho=sqrt(rhobar*rhobar+beta*beta)	! modified QR factorization
  c=rhobar/rho
  s=beta/rho
  teta=s*alfa
  rhobar=-c*alfa
  phi=c*phibar
  phibar=s*phibar
  t1=phi/rho
  t2=-teta/rho
  x=t1*w+x
  w=t2*w+v
  r=phibar/b1
  if(iout>0) write(iout,fmt='(i5,2g12.4,f10.6)') iter,x(1),phibar,r
  if(abs(r-rold)<1.0e-8) exit
  rold=r
end do 
deallocate(u)
deallocate(v)
deallocate(w)
return
end

subroutine normlz(n,x,s) 
dimension x(n)
s=0.
do i=1,n  
  s=s+x(i)**2
end do 
s=sqrt(s)
ss=1./s
do i=1,n  
  x(i)=x(i)*ss
end do 
return
end  


subroutine atupv(a,mn,m,n,u,v)

real*4 :: a(mn,n),u(m),v(n)

do j=1,n
  do i=1,m
    v(j)=v(j)+a(i,j)*u(i)
  enddo
enddo
return
end

subroutine avpu(a,mn,m,n,u,v)

real*4 :: a(mn,n),u(m),v(n)

do j=1,n
  do i=1,m
    u(i)=u(i)+a(i,j)*v(j)
  enddo
enddo

return
end
