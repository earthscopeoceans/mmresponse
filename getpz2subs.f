c     FFT routines for real-valued time signals. Not the most efficient,
c     but these have the correct scaling to mimic the Fourier *integral*
c     using the following conventions:

c     H(f) = \int h(t) exp(-2*pi*i*f*t) dt
c     h(t) = \int H(f) exp(+2*pi*i*f*t) df

c     You still need a complex time signal (x) and call clogp with
c     zzsign=-1 to go from time to frequency
c     But if you have a positive spectrum in the first half of (c0mplex)
c     array s, you can use ftinv to go back to a real-valued signal (r)
c     in the time domain.

      subroutine ftinv(s,npow,dt,r)
 
      parameter (ND=32768)
      dimension r(ND)
      complex s(ND)
 
c     Combines construction of negative spectrum, inverse fft, and
c     storage in a real array. It is possible to equivalence s and
c     r in the main program.
c     s(1)...s(nhalf) contain positive complex spectrum of real signal
c     nsmp=2*nhalf=2**npow, dt sampling interval, r destination array
c     WARNING: uses the sign convention of clogp not clogc

      nsmp=2**npow
      nhalf=nsmp/2
      call rspec(s,nhalf)
      call clogp(npow,s,+1.0,dt)
      do 10 i=1,nsmp
10    r(i)=real(s(i))
      return
      end

      subroutine rspec(s,np2)
 
c     Constructs negative spectrum
c     input: positive spectrum in s(1)...s(np2).
c     output: complete spectrum in s(1)...s(2*np2)
 
      parameter (ND=32768)
      complex s(ND)
      n=2*np2
      n1=np2+1
c     s(n1)=0.                  ! comment to keep Nyquist frequency
c     s(1)=0.                   ! comment to keep DC offset
      do 20 i=1,np2
20    s(np2+i)=conjg(s(np2+2-i))
      return
      end

      subroutine clogp(n,x,zzign,dt)
c--- performs fft on signals with length 2**n and sampling interval
c--- of dt seconds (if in the time domain; notice that dt*df=1/2**n).
c--- the signal is stored in x. it may be complex.
c--- the spectrum is returned in x. it is almost always complex.
c--- a time-to-frequency transform is done with zign=-1. (conform
c--- the convention adopted in SEED - the alternative
c--- convention may be obtained by taking complex conjugates after
c--- the call to clogp).
c--- the normalization factor 1./twopi occurs in the frequency-to
c--- time transform.
c--- normalization is such that physical dimensions are respected.
c--- thus, if the time signal is dimensioned in meters, the
c--- resulting spectral density in x is in meters/hz. for example,
c--- if the time signal is the unit sinc function of width dt, centered
c--- at t=0, the spectral density is dt for all values of the frequency.
c
c--- array locations: if x contains the spectrum, it has the spectrum
c--- for positive frequencies in the first 2**n/2+1 elements, such that
c--- x(1) is at 0 hz, x(2) at df hertz, and x(2**n/2+1) at the nyquist,
c--- where df=1./(2**n*dt) and the nyquist is 1./(2*dt) hz.
c--- the second half of x contains the spectrum for negative frequencies
c--- such that x(2**n) is at -df, x(2**n-1) at -2*df hz etcetera.
c--- if x contains the time signal, x(1) is at time 0, x(2)
c--- at time dt etc.
c
      parameter (ND=32768)
      dimension x(ND),m(21)
      complex x,wk,hold,q
      zign=zzign
      if(zign.ge.0.) then
        zign=1.
      else
        zign=-1.
      endif
      lx=2**n
      do 1 i=1,n
    1 m(i)=2**(n-i)
      do 4 l=1,n
      nblock=2**(l-1)
      lblock=lx/nblock
      lbhalf=lblock/2
      k=0
      do 4 iblock=1,nblock
      fk=k
      flx=lx
      v=zign*6.283185308*fk/flx
      wk=cmplx(cos(v),sin(v))
      istart=lblock*(iblock-1)
      do 2 i=1,lbhalf
      j=istart+i
      jh=j+lbhalf
      q=x(jh)*wk
      x(jh)=x(j)-q
      x(j)=x(j)+q
    2 continue
      do 3 i=2,n
      ii=i
      if(k.lt.m(i)) go to 4
    3 k=k-m(i)
    4 k=k+m(ii)
      k=0
      do 7 j=1,lx
      if(k.lt.j) go to 5
      hold=x(j)
      x(j)=x(k+1)
      x(k+1)=hold
    5 do 6 i=1,n
      ii=i
      if(k.lt.m(i)) go to 7
    6 k=k-m(i)
    7 k=k+m(ii)
      if(zign.lt.0.) go to 9
      flx=flx*dt
      do 8 i=1,lx
    8 x(i)=x(i)/flx
      return
    9 do 10 i=1,lx
   10 x(i)=x(i)*dt
      return
      end


c  optimization routines from Numerical Recipes, adjusted where needed

      FUNCTION BRENT(AX,BX,CX,F,TOL,XMIN)
      PARAMETER (ITMAX=200,CGOLD=.3819660,ZEPS=1.0E-10)
      A=MIN(AX,CX)
      B=MAX(AX,CX)
      V=BX
      W=V
      X=V
      E=0.
      FX=F(X)
      FV=FX
      FW=FX
      DO 11 ITER=1,ITMAX
        XM=0.5*(A+B)
        TOL1=TOL*ABS(X)+ZEPS
        TOL2=2.*TOL1
        IF(ABS(X-XM).LE.(TOL2-.5*(B-A))) GOTO 3
        IF(ABS(E).GT.TOL1) THEN
          R=(X-W)*(FX-FV)
          Q=(X-V)*(FX-FW)
          P=(X-V)*Q-(X-W)*R
          Q=2.*(Q-R)
          IF(Q.GT.0.) P=-P
          Q=ABS(Q)
          ETEMP=E
          E=D
          IF(ABS(P).GE.ABS(.5*Q*ETEMP).OR.P.LE.Q*(A-X).OR. 
     *        P.GE.Q*(B-X)) GOTO 1
          D=P/Q
          U=X+D
          IF(U-A.LT.TOL2 .OR. B-U.LT.TOL2) D=SIGN(TOL1,XM-X)
          GOTO 2
        ENDIF
1       IF(X.GE.XM) THEN
          E=A-X
        ELSE
          E=B-X
        ENDIF
        D=CGOLD*E
2       IF(ABS(D).GE.TOL1) THEN
          U=X+D
        ELSE
          U=X+SIGN(TOL1,D)
        ENDIF
        FU=F(U)
        IF(FU.LE.FX) THEN
          IF(U.GE.X) THEN
            A=X
          ELSE
            B=X
          ENDIF
          V=W
          FV=FW
          W=X
          FW=FX
          X=U
          FX=FU
        ELSE
          IF(U.LT.X) THEN
            A=U
          ELSE
            B=U
          ENDIF
          IF(FU.LE.FW .OR. abs(W-X)<zeps) THEN
            V=W
            FV=FW
            W=U
            FW=FU
          ELSE IF(FU.LE.FV .OR. abs(V-X)<zeps .OR. abs(V-W)<zeps) THEN
            V=U
            FV=FU
          ENDIF
        ENDIF
11    CONTINUE
      print *, 'Brent exceeds maximum iterations.'
      stop
3     XMIN=X
      BRENT=FX
      RETURN
      END
      FUNCTION F1DIM(X)
      PARAMETER (NMAX=50)
      COMMON /F1COM/ NCOM,PCOM(NMAX),XICOM(NMAX)
      DIMENSION XT(NMAX)
      DO 11 J=1,NCOM
        XT(J)=PCOM(J)+X*XICOM(J)
11    CONTINUE
      F1DIM=FUNC(XT)
      RETURN
      END
      SUBROUTINE LINMIN(P,XI,N,FRET)
      PARAMETER (NMAX=50,TOL=1.E-4)
      EXTERNAL F1DIM
      DIMENSION P(N),XI(N)
      COMMON /F1COM/ NCOM,PCOM(NMAX),XICOM(NMAX)
      !write(13,*) 'calling linmin, n,p=',n,(p(i),i=1,N)
      !write(13,*) 'xi=',(xi(i),i=1,n)
      NCOM=N
      DO 11 J=1,N
        PCOM(J)=P(J)
        XICOM(J)=XI(J)
11    CONTINUE
      AX=0.
      XX=1.
      BX=2.
      CALL MNBRAK(AX,XX,BX,FA,FX,FB,F1DIM)
      FRET=BRENT(AX,XX,BX,F1DIM,TOL,XMIN)
      !write(13,*) 'brent returns fret,xmin=',fret,xmin
      !flush(13)
      DO 12 J=1,N
        XI(J)=XMIN*XI(J)
        P(J)=P(J)+XI(J)
12    CONTINUE
      RETURN
      END
      SUBROUTINE MNBRAK(AX,BX,CX,FA,FB,FC,FUNC)
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.E-20)
      FA=FUNC(AX)
      FB=FUNC(BX)
      IF(FB.GT.FA)THEN
        DUM=AX
        AX=BX
        BX=DUM
        DUM=FB
        FB=FA
        FA=DUM
      ENDIF
      CX=BX+GOLD*(BX-AX)
      FC=FUNC(CX)
1     IF(FB.GE.FC)THEN
        R=(BX-AX)*(FB-FC)
        Q=(BX-CX)*(FB-FA)
        U=BX-((BX-CX)*Q-(BX-AX)*R)/(2.*SIGN(MAX(ABS(Q-R),TINY),Q-R))
        ULIM=BX+GLIMIT*(CX-BX)
        IF((BX-U)*(U-CX).GT.0.)THEN
          FU=FUNC(U)
          IF(FU.LT.FC)THEN
            AX=BX
            FA=FB
            BX=U
            FB=FU
            GO TO 1
          ELSE IF(FU.GT.FB)THEN
            CX=U
            FC=FU
            GO TO 1
          ENDIF
          U=CX+GOLD*(CX-BX)
          FU=FUNC(U)
        ELSE IF((CX-U)*(U-ULIM).GT.0.)THEN
          FU=FUNC(U)
          IF(FU.LT.FC)THEN
            BX=CX
            CX=U
            U=CX+GOLD*(CX-BX)
            FB=FC
            FC=FU
            FU=FUNC(U)
          ENDIF
        ELSE IF((U-ULIM)*(ULIM-CX).GE.0.)THEN
          U=ULIM
          FU=FUNC(U)
        ELSE
          U=CX+GOLD*(CX-BX)
          FU=FUNC(U)
        ENDIF
        AX=BX
        BX=CX
        CX=U
        FA=FB
        FB=FC
        FC=FU
        GO TO 1
      ENDIF
      RETURN
      END
      SUBROUTINE POWELL(P,XI,N,NP,FTOL,ITER,FRET)
      PARAMETER (NMAX=40,ITMAX=1500)
      DIMENSION P(NP),XI(NP,NP),PT(NMAX),PTT(NMAX),XIT(NMAX)
      FRET=FUNC(P)
      DO 11 J=1,N
        PT(J)=P(J)
11    CONTINUE
      ITER=0
1     ITER=ITER+1
      FP=FRET
      IBIG=0
      DEL=0.
      DO 13 I=1,N
        DO 12 J=1,N
          XIT(J)=XI(J,I)
12      CONTINUE
        CALL LINMIN(P,XIT,N,FRET)
        IF(ABS(FP-FRET).GT.DEL)THEN
          DEL=ABS(FP-FRET)
          IBIG=I
        ENDIF
13    CONTINUE
      IF(2.*ABS(FP-FRET).LE.FTOL*(ABS(FP)+ABS(FRET))) then
        RETURN
      endif
      IF(ITER.EQ.ITMAX) print *, 'Powell exceeding maximum iterations.'
      IF(ITER.EQ.ITMAX) stop
      DO 14 J=1,N
        PTT(J)=2.*P(J)-PT(J)
        XIT(J)=P(J)-PT(J)
        PT(J)=P(J)
14    CONTINUE
      FPTT=FUNC(PTT)
      IF(FPTT.GE.FP)GO TO 1
      T=2.*(FP-2.*FRET+FPTT)*(FP-FRET-DEL)**2-DEL*(FP-FPTT)**2
      IF(T.GE.0.)GO TO 1
      CALL LINMIN(P,XIT,N,FRET)
      DO 15 J=1,N
        XI(J,IBIG)=XIT(J)
15    CONTINUE
      GO TO 1
      END

