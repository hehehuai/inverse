c
c          Program ----- teletmg1.f -------
c
c          by Dapeng Zhao     January 2, 1993
c          updated by weiwei  
c          for adaptive tomography
c
      COMMON/CONTRL/NSTS,NEQS,NITLOC,RMSCUT,DVMAX,VDAMP,
     &              NHITCT,NITMAX,RMSTOP,STEPL,TLIM,NITPB
      COMMON/ITER/NIT
      COMMON/ITER1/NITRLQ
      OPEN(1,FILE="datacp")
      OPEN(2,FILE="station.txt")
      OPEN(3,FILE="base_grid.txt")
      OPEN(4,FILE="aftloc")
      OPEN(8,FILE="resjp1")
      OPEN(9,FILE="resjp2")
      open(10,file="out")
	open(11,file="des")
	open(12,file="hit")
      open(20,file="raypath")
      open(30,file="relocated")
	open(41,file="latmin")
	open(42,file="latmax")
	open(43,file="lonmin")
	open(44,file="lonmax")
	open(45,file="depmin")
	open(46,file="depmax")
c	open(50,file="wrong")

      PRINT*," Determine Tomography:"
      PRINT*," P- (1) or S-wave (2) or not (0)....."
c      READ(*,*) MPS
      MPS=1
      PRINT*," Relocate hypocenter (1) or not(0)..."
C      READ(*,*) IHYPO
      IHYPO=0
      PRINT*," See every residual in the scren:"
      PRINT*," Yes (1) or not(0) ....."
C      READ(*,*) ISEE
      ISEE=0
      PRINT*," Use local & regional events only (1)..6."
      PRINT*," Use teleseismic events only (2)........"
      PRINT*," Use all events (3)....................."
C      READ(*,*) IEVENT
      IEVENT=1
      CALL INPUT1
      CALL INPUT2
      CALL INPUT3(MPS)
      CALL INPUT4
      CALL INPUT5
	CALL INPUT6
     
      NIT     = 0
1     CONTINUE
      NIT     = NIT+1
      NITRLQ  = NIT
      IF(NIT.gt.1) Vdamp=100.0
      CALL STRT(NIT)
c   use local and regional events 
      IF(IEVENT.NE.2)    THEN
c     PRINT*,"INPUT NEQ1 & NEQ2 ......"
c     READ(*,*) NEQ1,NEQ2
c     DO 2 IE = NEQ1,NEQ2
      DO 2 IE = 1,2347

      PRINT*," Neq. = ",IE
      NE      = IE
      IF(IHYPO.EQ.1) CALL LOCEQK(NE,ISEE)
      IF(MPS.NE.0)   CALL FORWRD(NE,MPS,ISEE)
2     CONTINUE
      END IF
c   use teleseismic events
      IF(IEVENT.NE.1)    THEN
c     PRINT*,"INPUT NEQ1 & NEQ2 ......"
c     READ(*,*) NEQ1,NEQ2
c     DO 3 IE = NEQ1,NEQ2
      DO 3 IE = 1001,NEQS
      PRINT*," Neq. = ",IE
      NE      = IE
c      CALL FORTELE(NE,MPS,ISEE)
3     CONTINUE
      END IF
c      CALL SMOOTH
      CALL VELADJ
      CALL OUTADJ
      IF(NIT.GE.NITMAX)  GO TO 9
      GO TO 1
9     STOP
      END

      SUBROUTINE STRT(NIT)
      INCLUDE "PARAM"
      COMMON/CONTRL/NSTS,NEQS,NITLOC,RMSCUT,DVMAX,VDAMP,
     &              NHITCT,NITMAX,RMSTOP,STEPL,TLIM,NITPB
      COMMON/MODINV/A(MG),JA(MG),NA(MD),RDS(MD),KHIT(MU),NEQ4,
     &              JNDEX(MU),NOD,NOU,NOUEQ,ITOT,RNORM,XNORM
	COMMON/NEPICENT/N1,N2,N3,N4,N5,N6,N7
c   Option to take source parameters as unknowns or not.
c   If NEQ4 is zero, NO;   If NEQ4 is NEQS*4, YES.
c     NEQ4    = NEQS*4
      N1 = 0 
	N2 = 0
	N3 = 0
	N4 = 0
	N5 = 0
	N6 = 0
	N7 = 0
	
	NEQ4    = 0
      NOD     = 0
      ITOT    = 0
	
      DO 1 I  = 1,MU
1     KHIT(I) = 0
      WRITE(8,100) NIT
100   FORMAT(/" ITERATION STEP:",I3/)
      RETURN
      END

      SUBROUTINE STRTRAY
      INCLUDE "PARAM"
      COMMON/STATIN/PHIS(MST),RAMS(MST),HIGS(MST),STC(3,MST)
      COMMON/HLASTA/HLSS(MST,2)
      COMMON/RAYPATH/IWV(8,5),R0,PIDEG,EPS
      DIMENSION JWV(8,5)
      DATA JWV/1,2,2,6*1,2,2,6*1,2,2,1,2,4*1,
     &         2,1,2,2,4*1,2,1,2,2,3*1/
      R0       = 6371.0
      PIDEG    = 0.017453
      EPS      = 1.0E-6
      DO  1  I = 1,5
      DO  1  J = 1,8
1     IWV(J,I) = JWV(J,I)
      DO  2  I = 1,MST
      PS       = STC(1,I)
      RS       = STC(2,I)
      DO  2  J = 1,2
      K        = J
      CALL HLAY(PS,RS,HK,K)
      HLSS(I,K)= HK
2     CONTINUE
      RETURN
      END

      SUBROUTINE INPUT1
      COMMON/CONTRL/NSTS,NEQS,NITLOC,RMSCUT,DVMAX,VDAMP,
     &              NHITCT,NITMAX,RMSTOP,STEPL,TLIM,NITPB
      READ(1,100)   NSTS,NEQS,NITLOC,RMSCUT,DVMAX,VDAMP
      READ(1,200)   NHITCT,NITMAX,RMSTOP,STEPL,TLIM,NITPB
      WRITE(8,300)
      WRITE(8,100)  NSTS,NEQS,NITLOC,RMSCUT,DVMAX,VDAMP
      WRITE(8,400)
      WRITE(8,200)  NHITCT,NITMAX,RMSTOP,STEPL,TLIM,NITPB
100   FORMAT(i4,1x,I5,i4,3F7.3)
200   FORMAT(2I4,3F7.3,I3)
300   FORMAT(" CONTRL: NSTS,NEQS,NITLOC,RMSCUT,DVMAX,VDAMP")
400   FORMAT(" CONTRL: NHITCT,NITMAX,RMSTOP,STEPL,TLIM,NITPB")
      RETURN
      END

      SUBROUTINE INPUT2
      INCLUDE "PARAM" 
      COMMON/CONTRL/NSTS,NEQS,NITLOC,RMSCUT,DVMAX,VDAMP,
     &              NHITCT,NITMAX,RMSTOP,STEPL,TLIM,NITPB
      COMMON/STATIN/PHIS(MST),RAMS(MST),HIGS(MST),STC(3,MST)
      PIDEG    = 0.017453
      DO  1  N = 1,NSTS
      READ(2,100) PS,RS,HS
      STC(1,N) = (90.0-PS)*PIDEG
      STC(2,N) = RS*PIDEG
      STC(3,N) =-HS
      PHIS(N)  = PS
      RAMS(N)  = RS
      HIGS(N)  =-HS
1     CONTINUE
100   FORMAT(3F10.5)
      RETURN
      END

      SUBROUTINE INPUT3(MPS)
      INCLUDE "PARAM" 
      COMMON/VMOD3D/NPA,NRA,NHA,PNA(MPA),RNA(MRA),HNA(MHA),
     &              DVAP(MPA,MRA,MHA),VELAP(MPA,MRA,MHA) 
      COMMON/GRINET/NPA2,NPRA2,NODESA2,NODETOT
      CALL INPUTJB
100   FORMAT(3I3)
      READ(3,100)  NPA,NRA,NHA
      WRITE(8,100) NPA,NRA,NHA
      CALL PUT3(NPA,NRA,NHA,PNA,RNA,HNA,DVAP, 
     &          VELAP,NPA2,NPRA2,NODESA2,MPS)
      NODETOT= NODESA2
      CALL BLDMAP
      RETURN
      END

      SUBROUTINE PUT3(NPX,NRX,NHX,PNX,RNX,HNX,DVXP,
     &                VELXP,NPX2,NPRX2,NODESX2,MPS)
      DIMENSION  VELXP(NPX,NRX,NHX),DVXP(NPX,NRX,NHX),
     &           PNX(NPX),RNX(NRX),HNX(NHX),hlae(2)
      READ(3,110) (PNX(I),I=1,NPX)
      READ(3,110) (RNX(I),I=1,NRX)
      READ(3,110) (HNX(I),I=1,NHX)
      DO 1 K = 1,NHX
      DO 1 I = 1,NPX
      READ(3,110)  (DVXP(I,J,K),J=1,NRX)
1     CONTINUE
      DO 2 K = 1,NHX
      DO 2 J = 1,NRX
      DO 2 I = 1,NPX
      PP     = PNX(I)
      RR     = RNX(J)
      HH     = HNX(K)
      CALL LAYHE(PP,RR,HH,HLAE,LAY)
      CALL VEL1D(HH,V1,LAY,MPS)
      WRITE(6,200) HH,V1
      DV     = 0.01*DVXP(I,J,K)
      VELXP(I,J,K) = V1*(1.0+DV)
2     CONTINUE
      NPX2   = NPX-2
      NPRX2  = NPX2*(NRX-2)
      NODESX2= (NHX-2)*NPRX2
110   FORMAT(10(9F7.2/))
200   FORMAT(" HH = ",F6.1," V1 = ",F5.2)
      WRITE(8,110) (PNX(I),I=1,NPX)
      WRITE(8,110) (RNX(I),I=1,NRX)
      WRITE(8,110) (HNX(I),I=1,NHX)
      RETURN
      END



      SUBROUTINE INPUT4
      INCLUDE "PARAM" 
      COMMON/CONTRL/NSTS,NEQS,NITLOC,RMSCUT,DVMAX,VDAMP,
     &              NHITCT,NITMAX,RMSTOP,STEPL,TLIM,NITPB
      COMMON/EVENTS/MINO(MEQ),SECO(MEQ),PHIE(MEQ),RAME(MEQ),
     &              DEPE(MEQ),EVC(3,MEQ),KOBS(MEQ)
      COMMON/STODAT/IYMSTO(MEQ),IDSTO(MEQ),IHRSTO(MEQ)
      COMMON/STATIN/PHIS(MST),RAMS(MST),HIGS(MST),STC(3,MST)
      COMMON/OBSERV/ISTO(MAX,MEQ),SECT(MAX,MEQ),KWV(MAX,MEQ)
       DIMENSION     NSTN(4),TT(4),IKPS(4),IW(4),IUD(4)
      PIDEG     = 0.017453
      DO 10  N  = 1,NEQS
      READ(4,100) IY,IM,ID,IH,IMIN,SECO(N),DT,PHIE(N),DPHI,
     &            RAME(N),DRAM,DEPE(N),DDEP,NSTM,FMG,RMS
      IYMSTO(N) = IY*100+IM
      IDSTO(N)  = ID
      IHRSTO(N) = IH
      MINO(N)   = IMIN
      EVC(1,N)  = (90.0-PHIE(N))*PIDEG
      EVC(2,N)  = RAME(N)*PIDEG
      EVC(3,N)  = DEPE(N)
      L         = 0
      XNUM      = FLOAT(NSTM)/4.0
      NUMB      = IFIX(XNUM)
      IF(XNUM.GT.FLOAT(NUMB))  NUMB = NUMB+1
      DO 20  K  = 1,NUMB
      READ(4,200)  (NSTN(J),TT(J),IW(J),IKPS(J),IUD(J),J=1,4)
      DO 30 J   = 1,4
      L         = L+1
      ISTO(L,N) = NSTN(J)
      SECT(L,N) = TT(J)
      KWV(L,N)  = IKPS(J)
30    CONTINUE
20    CONTINUE
      KOBS(N)   = NSTM
10    CONTINUE
100   FORMAT(5I2,F6.2,F5.2,F7.3,F6.2,F9.3,
     &       2F6.2,F5.2,I3,F4.1,F5.2)
200   FORMAT(4(I4,1x,F6.2,3I2))
      RETURN
      END

      SUBROUTINE INPUT5
c      CALL READUBPP
      CALL ELLIPSE
      CALL STRTRAY
      RETURN
      END

      SUBROUTINE INPUT6
	INCLUDE "PARAM" 
	COMMON/NCELL/PEMIN(NKA),PEMAX(NKA),REMIN(NKA),REMAX(NKA),
     &             HEMIN(NKA),HEMAX(NKA)
	M=NC
	DO 1 I = 1,NC
	READ(43,*) RLONMIN
	REMIN(I)=RLONMIN
    1 CONTINUE

      DO 2 I = 1,NC
	READ(44,*) RLONMAX
	REMAX(I)=RLONMAX
    2 CONTINUE

      DO 3 I = 1,NC
	READ(41,*) RLATMIN
	PEMIN(I)=RLATMIN
   3  CONTINUE
  
	DO 4 I = 1,NC
	READ(42,*) RLATMAX
	PEMAX(I)=RLATMAX
   4  CONTINUE

      DO 5 I = 1,NC
	READ(45,*) DEPMIN
	HEMIN(I)=DEPMIN
   5  CONTINUE

      DO 6 I = 1,NC
	READ(46,*) DEPMAX
	HEMAX(I)=DEPMAX
   6  CONTINUE

	END
      SUBROUTINE LOCEQK(NE,ISEE)
      INCLUDE "PARAM" 
      COMMON/CONTRL/NSTS,NEQS,NITLOC,RMSCUT,DVMAX,VDAMP,
     &              NHITCT,NITMAX,RMSTOP,STEPL,TLIM,NITPB
      COMMON/STODAT/IYMSTO(MEQ),IDSTO(MEQ),IHRSTO(MEQ)
      COMMON/EVENTS/MINO(MEQ),SECO(MEQ),PHIE(MEQ),RAME(MEQ),
     &              DEPE(MEQ),EVC(3,MEQ),KOBS(MEQ)
      COMMON/STATIN/PHIS(MST),RAMS(MST),HIGS(MST),STC(3,MST)
      COMMON/OBSERV/ISTO(MAX,MEQ),SECT(MAX,MEQ),KWV(MAX,MEQ)
      COMMON/HYPINV/RES(MAX),DTH(MAX,4),W(MAX)
      COMMON/ITER1/NITRLQ
      DIMENSION     ADJ(4),ERR(4),GF(5,11),ANS(5,6),HLAE(2)
      DATA PID,C0,HCOR,HMAX/0.017453,111.19,200.0,100.0/
      IUN       = 8
      IF(ISEE.EQ.1)  IUN = 6
      NOBS      = KOBS(NE)
      WRITE(IUN,100) NE,SECO(NE),PHIE(NE),RAME(NE),DEPE(NE),NOBS
      OSTRT     = SECO(NE)
      PSTRT     = EVC(1,NE)
      RSTRT     = EVC(2,NE)
      HSTRT     = EVC(3,NE)
      NIT       = 0
10    NIT       = NIT+1
      PE        = EVC(1,NE)
      RE        = EVC(2,NE)
      HE        = EVC(3,NE)
      CALL LAYHE(PE,RE,HE,HLAE,LAY)
      DO 1  MO  = 1,NOBS
      NO        = MO
      KPS       = KWV(NO,NE)
      IF(KPS.LT.1.OR.KPS.GT.8)   GO TO 1
c	IF(KPS.EQ.2)              GO TO 1
      NS        = ISTO(NO,NE)
      IF(NS.LT.1.OR.NS.GT.NSTS)  GO TO 1
      PS        = STC(1,NS)
      RS        = STC(2,NS)
      HS        = STC(3,NS)
      PO        = SECT(NO,NE)-SECO(NE)
	Num        = 0       !  0  location
	                    !  1  forwrd
      CALL TRAVT(PE,RE,HE,PS,RS,HS,KPS,NS,HLAE,LAY,
     &           TLIM,NITPB,PO,DEL,AZ,VE,FSTIME,Num)
      DT        = 0.0
      DELC      = DEL/C0
      IF(DEL.GT.HCOR)  CALL CORREC(HE,DELC,KPS,PE,AZ,DT)
      FSTIME    = FSTIME+DT
      IF(ISEE.NE.0)              THEN
      RESS      = PO-FSTIME
      WRITE(IUN,150)  NE,NOBS,NO,KPS,DEL,VE,RESS
      END IF
      CALL TTMDER(NE,NO,KPS,VE,FSTIME,0,RRES,DEL)
1     CONTINUE
      CALL WTHYP(NOBS,GF,RMSWT,NWR)
      IF(NWR.LE.4)               GO TO 6
      CALL LINER(4,9,GF,ANS)
      DO  2  I  = 1,4
      ADJ(I)    = ANS(I,1)
      IF(I.EQ.2.OR.I.EQ.3) ADJ(I) = ADJ(I)*0.01
2     CONTINUE
      ANS(5,1)  =-1.0
      SG        = 0.0
      DO 3  K   = 1,5
      DO 3  J   = 1,5
      SG        = SG+GF(K,J)*ANS(K,1)*ANS(J,1)
3     CONTINUE
      DT        = SQRT(ABS(SG)/(NWR-4))
      DO 4  K   = 1,4
      J         = K+1
      ERR(K)    = SQRT(ABS(ANS(K,J)))*DT
      IF(K.EQ.2.OR.K.EQ.3) ERR(K) = ERR(K)*0.01/PID
4     CONTINUE
      HADJ      = ADJ(4)+EVC(3,NE)
      IF(HADJ.LT.0.0)    ADJ(4) = 0.0
      SECO(NE)  = SECO(NE)+ADJ(1)
      DO 5  I   = 1,3
      J         = I+1
      EVC(I,NE) = EVC(I,NE)+ADJ(J)
5     CONTINUE
      IF(EVC(3,NE).GT.HMAX)  EVC(3,NE) = HMAX
6     SNEW      = SECO(NE)
      PNEW      = 90.0-EVC(1,NE)/PID
      RNEW      = EVC(2,NE)/PID
      HNEW      = EVC(3,NE)
      WRITE(IUN,200) NIT,SNEW,ERR(1),PNEW,ERR(2),RNEW,  
     &               ERR(3),HNEW,ERR(4),NWR,RMSWT
      IF(NIT.LT.NITLOC.AND.RMSWT.GT.RMSCUT)  GO TO 10
      DPE       = (EVC(1,NE)-PSTRT)/PID
      DRE       = (EVC(2,NE)-RSTRT)/PID
      DHE       =  EVC(3,NE)-HSTRT
      DOT       =  SECO(NE) -OSTRT
      WRITE(IUN,300) DOT,DPE,DRE,DHE
      WRITE(9,400)   NE,SNEW,ERR(1),PNEW,ERR(2),RNEW, 
     &               ERR(3),HNEW,ERR(4),NWR,RMSWT
      write(30,110) PNEW,RNEW,HNEW
110   format(3f9.4)
100   FORMAT(/" IEQ:",I4," initial value----- ", 
     &         F7.2,F8.3,F9.3,F8.2,I3)
150   FORMAT(" NE,NOBS,NO,KPS,DEL,VE,RES = ",4I5,F8.1,F6.1,F7.2)
200   FORMAT(I2,F6.2,"(",F5.2,")",F7.3,"(",F5.3,")",F8.3,
     &       "(",F5.3,")",F7.2,"(",F5.2,")",I4,F6.2)
300   FORMAT(F6.2,2F6.3,F8.2)
400   FORMAT(I7,F8.2,F7.2,F9.3,F8.3,F10.3,F8.3,F9.2,F7.2,I6,F7.2)
      RETURN
      END

      SUBROUTINE LAYHE(PE,RE,HE,HLAE,LAY)
      DIMENSION  HLAE(2)
      LAY     = 1
      DO  1 K = 1,2
      I       = K
      CALL HLAY(PE,RE,HI,I)
      HLAE(I) = HI
      IF(HE.GT.HI) LAY = LAY+1
1     CONTINUE
      RETURN
      END

      SUBROUTINE LINER(M,N,GF,ANS)
      DIMENSION  GF(5,11),C(5,11),ANS(5,6)
      DO 1  K  = 1,5
      DO 1  L  = 1,11
1     C(K,L)   = GF(K,L)
      DO 2  K  = 1,M
      W        = C(K,K)
      IF(ABS(W).LT.1.0E-5)   W = 1.0
      KK       = K+1
      DO 3  J  = KK,N
3     C(K,J)   = C(K,J)/W
      DO 4  I  = 1,M
      IF(I.EQ.K)             GO TO 4
      W        = C(I,K)
      DO 5  J  = KK,N
5     C(I,J)   = C(I,J)-W*C(K,J)
4     CONTINUE
2     CONTINUE
      J        = N-M
      DO 6  I  = 1,J
      JJ       = M+I
      DO 6  K  = 1,M
6     ANS(K,I) = C(K,JJ)
      RETURN
      END

      SUBROUTINE WTHYP(NOBS,GF,RMSWT,NWR)
      INCLUDE "PARAM" 
      COMMON/HYPINV/RES(MAX),DTH(MAX,4),W(MAX)
      DIMENSION     AQ(5),GF(5,11)
      EPC     = 1.0E-5
      DO  1 K = 1,5
      DO  2 J = 1,9
2     GF(K,J) = 0.0
      M       = K+5
1     GF(K,M) = 1.0
      NWR     = 0
      WNORM   = 0.0
      RMSWT   = 0.0
      DO  4 I = 1,NOBS
      WI      = W(I)
      WNORM   = WNORM+WI*WI
      RMSW    = RES(I)*WI
      RMSWT   = RMSWT+RMSW*RMSW
      IF(WI.GT.EPC) NWR = NWR+1
      DO  5 J = 1,4
5     AQ(J)   = DTH(I,J)*W(I)
      AQ(5)   = RES(I)*W(I)
      DO  6 K = 1,5
      DO  6 J = 1,5
6     GF(K,J) = GF(K,J)+AQ(K)*AQ(J)
4     CONTINUE
      IF(WNORM.GT.EPC) RMSWT = SQRT(RMSWT/WNORM)
      RETURN
      END

      SUBROUTINE FORWRD(NE,MPS,ISEE)
      INCLUDE "PARAM" 
      COMMON/CONTRL/NSTS,NEQS,NITLOC,RMSCUT,DVMAX,VDAMP,
     &              NHITCT,NITMAX,RMSTOP,STEPL,TLIM,NITPB
      COMMON/EVENTS/MINO(MEQ),SECO(MEQ),PHIE(MEQ),RAME(MEQ),
     &              DEPE(MEQ),EVC(3,MEQ),KOBS(MEQ)
      COMMON/STATIN/PHIS(MST),RAMS(MST),HIGS(MST),STC(3,MST)
      COMMON/OBSERV/ISTO(MAX,MEQ),SECT(MAX,MEQ),KWV(MAX,MEQ)
      COMMON/RAYLOC/NRP,RP(3,MDT),IVK(MDT),IWK(MDT)
      COMMON/ITER/NIT


      DIMENSION HLAE(2)
      C0      = 111.19
	PIDEG   = 0.017453
      NOBS    = KOBS(NE)
      PE      = EVC(1,NE)
      RE      = EVC(2,NE)
      HE      = EVC(3,NE)
      CALL LAYHE(PE,RE,HE,HLAE,LAY)
      IF(ISEE.EQ.1)
     &   WRITE(6,100) NE,PHIE(NE),RAME(NE),DEPE(NE),NOBS
      DO 1 MO = 1,NOBS
      NO      = MO
      KPS     = KWV(NO,NE)
      IF(KPS.LT.1.OR.KPS.GT.8)  GO TO 1
      IF(KPS.EQ.2)              GO TO 1
      NS      = ISTO(NO,NE)
      IF(NS.LT.1.OR.NS.GT.NSTS) GO TO 1
      PS      = STC(1,NS)
      RS      = STC(2,NS)
      HS      = STC(3,NS)
      PO      = SECT(NO,NE)-SECO(NE)
	Num      = 1 !  (0 IS FOR LOCATION, 1 IS FOR FORWRD)
      CALL TRAVT(PE,RE,HE,PS,RS,HS,KPS,NS,HLAE,LAY,
     &           TLIM,NITPB,PO,DEL,AZ,VE,FSTIME,Num)
        
      IF(NIT.EQ.NITMAX) THEN
      DO J=NRP,1,-1
	WRITE(20,30) 90-RP(1,J)/PIDEG,RP(2,J)/PIDEG,RP(3,J),DEL
      enddo
      write(20,'(">")')
      endif
30    format(4f10.3)
      DT      = 0.0
      DELC    = DEL/C0
      IF(DEL.GT.250.0) CALL CORREC(HE,DELC,KPS,PE,AZ,DT)
      FSTIME  = FSTIME+DT
      IF(ISEE.NE.1)             GO TO 2
      RESS    = PO-FSTIME
      WRITE(6,200) NE,NOBS,NO,KPS,DEL,VE,RESS
2     CALL TTMDER(NE,NO,KPS,VE,FSTIME,1,RRES,DEL)
c      IF(NIT.EQ.NITMAX.AND.DEL.GT.500) THEN
        IF(NIT.EQ.NITMAX)  THEN
	WRITE(11,*)  DEL,RRES
	ELSE
	ENDIF
1     CONTINUE
100   FORMAT(" NEQ,PHI,RAM,DEP,NOBS = ",I5,F6.2,F7.2,F6.1,I4)
200   FORMAT(" NE,NOBS,NO,KPS,DEL,VE,RES = ",4I5,F8.1,F6.1,F7.2)
      RETURN
      END

      
      SUBROUTINE ELLIPSE
      DIMENSION  QD0T0(5),QD0T1(5),QD0T2(5),QD3T0(5),QD3T1(5), 
     &           QD3T2(5),QD6T0(5),QD6T1(5),QD6T2(5),TD0T0(5),  
     &           TD0T1(5),TD0T2(5),TD3T0(5),TD3T1(5),TD3T2(5),  
     &           TD6T0(5),TD6T1(5),TD6T2(5),DELDA(5)
      COMMON/TOYUAN/PD0T0(5),PD0T1(5),PD0T2(5),PD3T0(5),PD3T1(5), 
     &              PD3T2(5),PD6T0(5),PD6T1(5),PD6T2(5),SD0T0(5),  
     &              SD0T1(5),SD0T2(5),SD3T0(5),SD3T1(5),SD3T2(5),  
     &              SD6T0(5),SD6T1(5),SD6T2(5),DELRA(5)
      DATA QD0T0/ 0.00,-0.18,-0.34,-0.48,-0.57/
      DATA QD0T1/ 0.00,-0.02,-0.05,-0.10,-0.16/
      DATA QD0T2/ 0.00, 0.00,-0.01,-0.03,-0.06/
      DATA QD3T0/-0.12,-0.17,-0.29,-0.38,-0.46/
      DATA QD3T1/ 0.00,-0.11,-0.14,-0.19,-0.24/
      DATA QD3T2/ 0.00, 0.00,-0.01,-0.03,-0.07/
      DATA QD6T0/-0.26,-0.27,-0.31,-0.37,-0.41/
      DATA QD6T1/ 0.00,-0.20,-0.28,-0.34,-0.37/
      DATA QD6T2/ 0.00, 0.00,-0.02,-0.04,-0.08/
      DATA TD0T0/ 0.00,-0.33,-0.57,-0.82,-1.00/
      DATA TD0T1/ 0.00,-0.02,-0.10,-0.19,-0.30/
      DATA TD0T2/ 0.00, 0.00,-0.02,-0.05,-0.12/
      DATA TD3T0/-0.22,-0.31,-0.52,-0.68,-0.84/
      DATA TD3T1/ 0.00,-0.20,-0.25,-0.34,-0.43/
      DATA TD3T2/ 0.00, 0.00,-0.02,-0.05,-0.13/
      DATA TD6T0/-0.47,-0.49,-0.56,-0.66,-0.72/
      DATA TD6T1/ 0.00,-0.36,-0.50,-0.61,-0.67/
      DATA TD6T2/ 0.00, 0.00,-0.03,-0.07,-0.14/
      DATA DELDA/ 0.00, 5.00,10.00,15.00,20.00/
      DO 1  L  = 1,5
      PD0T0(L) = QD0T0(L)
      PD0T1(L) = QD0T1(L)
      PD0T2(L) = QD0T2(L)
      PD3T0(L) = QD3T0(L)
      PD3T1(L) = QD3T1(L)
      PD3T2(L) = QD3T2(L)
      PD6T0(L) = QD6T0(L)
      PD6T1(L) = QD6T1(L)
      PD6T2(L) = QD6T2(L)
      SD0T0(L) = TD0T0(L)
      SD0T1(L) = TD0T1(L)
      SD0T2(L) = TD0T2(L)
      SD3T0(L) = TD3T0(L)
      SD3T1(L) = TD3T1(L)
      SD3T2(L) = TD3T2(L)
      SD6T0(L) = TD6T0(L)
      SD6T1(L) = TD6T1(L)
      SD6T2(L) = TD6T2(L)
      DELRA(L) = DELDA(L)
1     CONTINUE
      RETURN
      END

      SUBROUTINE CORREC(DEP,DEL,KPS,PHI,AZM,DT)
      COMMON/TOYUAN/PD0T0(5),PD0T1(5),PD0T2(5),PD3T0(5),PD3T1(5), 
     &              PD3T2(5),PD6T0(5),PD6T1(5),PD6T2(5),SD0T0(5),  
     &              SD0T1(5),SD0T2(5),SD3T0(5),SD3T1(5),SD3T2(5),  
     &              SD6T0(5),SD6T1(5),SD6T2(5),DELRA(5)
      COMMON/ELLIPS/DELTB(5),IDEL,IDEL1,T0,T1,T2
      DO 10 L = 1,5
      DELTB(L) = DELRA(L)
10    CONTINUE
      IPS    = 2
      IF(KPS.NE.2) IPS = 1
      IDEL   = 0
      DO 1 L = 1,5
      IF(DEL.GE.DELTB(L)) IDEL = IDEL+1
1     CONTINUE
      IDEL1  = IDEL+1
      IF(IDEL1.GE.5) IDEL1 = 5
      IF(IPS.EQ.1)                            THEN
       IF(DEP.GE.0.0.AND.DEP.LT.300.0)        THEN
        CALL T012(PD0T0,PD3T0,PD0T1,PD3T1,PD0T2,PD3T2,
     &            0.0,300.0,DEP,DEL)
       ELSE IF(DEP.GE.300.0.AND.DEP.LT.650.0) THEN
        CALL T012(PD3T0,PD6T0,PD3T1,PD6T1,PD3T2,PD6T2,
     &            300.0,650.0,DEP,DEL)
       ELSE
        CALL T012(PD6T0,PD6T0,PD6T1,PD6T1,PD6T2,PD6T2,
     &            650.0,651.0,DEP,DEL)
       END IF
      ELSE IF(IPS.EQ.2)                       THEN
       IF(DEP.GE.0.0.AND.DEP.LT.300.0)        THEN
        CALL T012(SD0T0,SD3T0,SD0T1,SD3T1,SD0T2,SD3T2,
     &            0.0,300.0,DEP,DEL)
       ELSE IF(DEP.GE.300.0.AND.DEP.LT.650.0) THEN
        CALL T012(SD3T0,SD6T0,SD3T1,SD6T1,SD3T2,SD6T2,
     &            300.0,650.0,DEP,DEL)
       ELSE
        CALL T012(SD6T0,SD6T0,SD6T1,SD6T1,SD6T2,SD6T2,
     &            650.0,651.0,DEP,DEL)
       END IF
      ELSE
       T0    = 0.0
       T1    = 0.0
       T2    = 0.0
      END IF
      PHI2   = 2.0*PHI
      SIN2P  = SIN(PHI2)
      COS2P  = COS(PHI2)
      COSAZ  = COS(AZM)
      COS2AZ = 2.0*COSAZ*COSAZ-1.0
      DT     = (0.25+0.75*COS2P)*T0+0.866*SIN2P*COSAZ*T1
     &         +0.433*(1.0-COS2P)*COS2AZ*T2
      RETURN
      END

      SUBROUTINE T012(PD0T0,PD3T0,PD0T1,PD3T1,PD0T2,PD3T2,
     &                DEP1,DEP2,DEP,DEL)
      DIMENSION PD0T0(5),PD3T0(5),PD0T1(5),PD3T1(5),PD0T2(5),PD3T2(5)
      COMMON/ELLIPS/DELTB(5),IDEL,IDEL1,T0,T1,T2
      DEPP   = (DEP-DEP1)/(DEP2-DEP1)
      DELT   = (DEL-DELTB(IDEL))/5.0
      T0A    = (PD0T0(IDEL1)-PD0T0(IDEL))*DELT+PD0T0(IDEL)
      T0B    = (PD3T0(IDEL1)-PD3T0(IDEL))*DELT+PD3T0(IDEL)
      T1A    = (PD0T1(IDEL1)-PD0T1(IDEL))*DELT+PD0T1(IDEL)
      T1B    = (PD3T1(IDEL1)-PD3T1(IDEL))*DELT+PD3T1(IDEL)
      T2A    = (PD0T2(IDEL1)-PD0T2(IDEL))*DELT+PD0T2(IDEL)
      T2B    = (PD3T2(IDEL1)-PD3T2(IDEL))*DELT+PD3T2(IDEL)
      T0     = T0A+(T0B-T0A)*DEPP
      T1     = T1A+(T1B-T1A)*DEPP
      T2     = T2A+(T2B-T2A)*DEPP
      RETURN
      END

      SUBROUTINE TTMDER(NE,NO,KPS,VE,TTIME,KOPT,RRES,DEL)
      INCLUDE "PARAM" 
      COMMON/CONTRL/NSTS,NEQS,NITLOC,RMSCUT,DVMAX,VDAMP,
     &              NHITCT,NITMAX,RMSTOP,STEPL,TLIM,NITPB
      COMMON/VMOD3D/NPA,NRA,NHA,PNA(MPA),RNA(MRA),HNA(MHA),
     &              DVAP(MPA,MRA,MHA),VELAP(MPA,MRA,MHA)
      COMMON/MODINV/A(MG),JA(MG),NA(MD),RDS(MD),KHIT(MU),NEQ4,
     &              JNDEX(MU),NOD,NOU,NOUEQ,ITOT,RNORM,XNORM
      COMMON/GRINET/NPA2,NPRA2,NODESA2,NODETOT
      COMMON/EVENTS/MINO(MEQ),SECO(MEQ),PHIE(MEQ),RAME(MEQ),
     &              DEPE(MEQ),EVC(3,MEQ),KOBS(MEQ)
      COMMON/OBSERV/ISTO(MAX,MEQ),SECT(MAX,MEQ),KWV(MAX,MEQ)
      COMMON/HYPINV/RES(MAX),DTH(MAX,4),W(MAX)
      COMMON/WEIGHT/WV(8),IP,JP,KP,IP1,JP1,KP1
      COMMON/NPRH12/NP1,NR1,NH1,NP2,NPR2
      COMMON/RAYPATH/IWV(8,5),R0,PIDEG,EPS
      COMMON/RAYLOC/NRP,RP(3,MDT),IVK(MDT),IWK(MDT)
      COMMON/ONED/DTM(300),NDM(300),NUMB
	COMMON/ITER/NIT
	COMMON/HITS/MHIT(MNUM),NCELL
      DIMENSION     IN(8),LX(24)
      R0        = 6371.0
	PIDEG     = 0.017453
      EPC       = 1.0E-5
      PE        = RP(1,NRP)
      RE        = RP(2,NRP)
      HE        = RP(3,NRP)
      NRP1      = NRP-1
      PE1       = RP(1,NRP1)
      RE1       = RP(2,NRP1)
      HE1       = RP(3,NRP1)
      CALL CONAZ(RE,PE,RE1,PE1,DES,AZ)
      DES       = DES*(1.0-HE/R0)
      DHE       = HE1-HE
      ADH       = ABS(DHE)
      IF(DES.LT.EPC.AND.ADH.LT.EPC) DHE = 1.0
      TH        = ATAN2(DES,DHE)
      DTDDEL    = SIN(TH)*(R0-HE)/VE
      DTDR      = COS(TH)/VE
      DTH(NO,1) = 1.0
      DTH(NO,2) = DTDDEL*COS(AZ)*0.01
      DTH(NO,3) =-DTDDEL*SIN(AZ)*SIN(PE)*0.01
      DTH(NO,4) =-DTDR
      RES(NO)   = SECT(NO,NE)-SECO(NE)-TTIME
      RRES      = RES(NO)
      AR        = ABS(RES(NO))
      W(NO)     = 0.0
      RESM      = 2.5
	IF(DEL.GE.500.AND.DEL.LT.600)  RESM=3.0
	IF(DEL.GE.600.AND.DEL.LT.700)  RESM=3.5
	IF(DEL.GE.700.AND.DEL.LT.800)  RESM=4.0
	IF(DEL.GE.800.AND.DEL.LT.900)  RESM=4.5
      IF(DEL.GE.900.AND.DEL.LT.1000)  RESM=5.0
	
      IF(AR.LE.RESM)    W(NO) = 1.0
      IF(KPS.EQ.2)      W(NO) = W(NO)*0.5
      IF(KOPT.EQ.0)     RETURN
      IF(W(NO).LT.EPC)  RETURN
	
c  travel time derivatives with respect to velocity parameters
c   loop over all segments comprising the ray path
      HALF    = 0.5
      NOD     = NOD+1
      NUMB    = 0
      DO 10 I = 1,300
      DTM(I)  = 0.0
      NDM(I)  = 0
10    CONTINUE
      DO  1 I = 1,NRP1
      LAY     = IVK(I)
      IPS     = IWV(KPS,IWK(I))
      CALL IJMNX(LAY,IJMN)
      P       = RP(1,I)
      R       = RP(2,I)
      H       = RP(3,I)
      I1      = I+1
      P1      = RP(1,I1)
      R1      = RP(2,I1)
      H1      = RP(3,I1)
      DP      = P1-P
      DR      = R1-R
      DH      = H1-H
c  compute segment length
      CALL LENGTH(P,R,H,P1,R1,H1,SL)
c  decide the number of subsegments and compute length
      NSEG    = IFIX(SL/STEPL)+1
      FNSEGI  = 1.0/FLOAT(NSEG)
      SSL     = SL*FNSEGI
      DPS     = DP*FNSEGI
      DRS     = DR*FNSEGI
      DHS     = DH*FNSEGI
      PX      = P-HALF*DPS
      RX      = R-HALF*DRS
      HX      = H-HALF*DHS
c  loop over all subsegments
      DO 2 IS = 1,NSEG
      PX      = PX+DPS
      RX      = RX+DRS
      HX      = HX+DHS
c----------------------------------------------------------
c      calculate  the  hit number in cell
      RLAT = 90-PX/PIDEG
	RLON = RX/PIDEG
	RH   = HX
C	CALL CELL(RLAT,RLON,RH)
	CALL HIT(PX,RX,HX)

      CALL VEL33(PX,RX,HX,VX,LAY,IPS)
      DT      = SSL/VX
      CALL NODO(NP1,NR1,NH1,IP,JP,KP,IP1,JP1,KP1,LX,IND)
c  nodes with non-zero weight
C      IN(1)   = IP-1+NP2*(JP-2)+NPR2*(KP-2)
C      IN(2)   = IN(1)+1
C      IN(3)   = IN(1)+NP2
C      IN(4)   = IN(3)+1
C      IN(5)   = IN(1)+NPR2
C      IN(6)   = IN(5)+1
C      IN(7)   = IN(5)+NP2
C      IN(8)   = IN(7)+1

	IN(1)   =  IP-1+NP2*(JP-2)+NPR2*(KP-2)
      IN(2)   = IP1-1+NP2*(JP-2)+NPR2*(KP-2)
      IN(3)   =  IP-1+NP2*(JP1-2)+NPR2*(KP-2)
      IN(4)   = IP1-1+NP2*(JP1-2)+NPR2*(KP-2)
      IN(5)   =  IP-1+NP2*(JP-2)+NPR2*(KP1-2)
      IN(6)   = IP1-1+NP2*(JP-2)+NPR2*(KP1-2)
      IN(7)   =  IP-1+NP2*(JP1-2)+NPR2*(KP1-2)
      IN(8)   = IP1-1+NP2*(JP1-2)+NPR2*(KP1-2)
c  accumulate model partial derivatives
      DO 3 KK = 1,2
      KK1     = KK-1
      DO 3 JJ = 1,2
      JJ1     = JJ-1
      DO 3 II = 1,2
      IJK     = II+2*JJ1+4*KK1
      IF(IND.EQ.0)        GO TO 5
      DO 4 KZ = 1,24
      IF(IJK.EQ.LX(KZ))   GO TO 3
4     CONTINUE
5     II1 = II-1
      INI = IJMN+IN(IJK)
	IF(II1.EQ.0) THEN
       L  = IP
	ELSE
	 L  = IP1
	ENDIF
      
	IF(JJ1.EQ.0) THEN
       M  = JP
	ELSE
	 M  = JP1
	ENDIF

	IF(KK1.EQ.0) THEN
       N  = KP
	ELSE
	 N  = KP1
	ENDIF
     
      VIJK    = VELAP(L,M,N)
      DTV     = DT*WV(IJK)*VIJK/VX
      CALL ONERAY(INI,DTV)
3     CONTINUE
2     CONTINUE
1     CONTINUE
      WN      = W(NO)
      NUMC    = 0
      DO 6  I = 1,NUMB
      DTMI    = DTM(I)
      IF(DTMI.LT.EPC)     GO TO 6
      NUMC    = NUMC+1
      K       = NDM(I)+NEQ4
      KHIT(K) = KHIT(K)+1
      ITOT    = ITOT+1
      A(ITOT) = DTMI*WN
      JA(ITOT)= K
6     CONTINUE
      IF(NEQ4.EQ.0)       GO TO 8
      DO 7 I  = 1,4
      K       = (NE-1)*4+I
      KHIT(K) = KHIT(K)+1
      ITOT    = ITOT+1
      A(ITOT) = DTH(NO,I)*WN
      JA(ITOT)= K
7     CONTINUE
8     IHYPO   = NEQ4/NEQS
      NA(NOD) = NUMC+IHYPO
      RDS(NOD)= RES(NO)*WN
      RETURN
      END

      

      SUBROUTINE ONERAY(INI,DTV)
      COMMON/ONED/DTM(300),NDM(300),NUMB
      NB       = MAX0(1,NUMB)
      DO 1  I  = 1,NB
      IF(INI.EQ.NDM(I))  GO TO 2
1     CONTINUE
      NUMB     = NUMB+1
      NDM(NUMB)= INI
      DTM(NUMB)= DTV
      RETURN
2     DTM(I)   = DTM(I)+DTV
      RETURN
      END

      SUBROUTINE IJMNX(LAY,IJMN)
      INCLUDE "PARAM" 
      COMMON/VMOD3D/NPA,NRA,NHA,PNA(MPA),RNA(MRA),HNA(MHA),
     &              DVAP(MPA,MRA,MHA),VELAP(MPA,MRA,MHA)
      COMMON/GRINET/NPA2,NPRA2,NODESA2,NODETOT
      IJMN  = 0
      CALL NW(NPA,NRA,NHA,NPA2,NPRA2)
      RETURN
      END

      SUBROUTINE NW(NPA,NRA,NHA,NPA2,NPRA2)
      COMMON/NPRH12/NP1,NR1,NH1,NP2,NPR2
      NP1      = NPA-1
      NR1      = NRA-1
      NH1      = NHA-1
      NP2      = NPA2
      NPR2     = NPRA2
      RETURN
      END

      SUBROUTINE NODO(NP1,NR1,NH1,IP,JP,KP,IP1,JP1,KP1,LX,IND)
      DIMENSION  LXD(24),LX(24)
      DATA LXD/1,3,5,7,2,4,6,8,1,2,5,6, 
     &         3,4,7,8,1,2,3,4,5,6,7,8/
      DO 1 I = 1,24
1     LX(I)  = 0
      IF(IP.EQ.1.OR.IP1.EQ.NP1+1)  GO TO 10
      IF(JP.EQ.1.OR.JP1.EQ.NR1+1)  GO TO 30
      IF(KP.EQ.1.OR.KP1.EQ.NH1+1)  GO TO 50
      IND    = 0
      RETURN
10    IF(IP.NE.1)               GO TO 20
      DO 2 I = 1,4
2     LX(I)  = LXD(I)
20    IF(IP1.NE.NP1+1)             GO TO 30
      DO 3 I = 5,8
3     LX(I)  = LXD(I)
30    IF(JP.NE.1)               GO TO 40
      DO 4 I = 9,12
4     LX(I)  = LXD(I)
40    IF(JP1.NE.NR1+1)             GO TO 50
      DO 5 I = 13,16
5     LX(I)  = LXD(I)
50    IF(KP.NE.1)               GO TO 60
      DO 6 I = 17,20
6     LX(I)  = LXD(I)
60    IF(KP1.NE.NH1+1)             GO TO 70
      DO 7 I = 21,24
7     LX(I)  = LXD(I)
70    IND    = 1
      RETURN
      END

      SUBROUTINE BLDMAP
      INCLUDE "PARAM" 
      COMMON/VMOD3D/NPA,NRA,NHA,PNA(MPA),RNA(MRA),HNA(MHA),
     &              DVAP(MPA,MRA,MHA),VELAP(MPA,MRA,MHA)
      COMMON/LOCATE/PLA,RLA,HLA,IPLOCA(MKA),IRLOCA(MKA),IHLOCA(MKA)
      CALL LOCX(PNA,RNA,HNA,NPA,NRA,NHA,MKA,
     &     PLA,RLA,HLA,IPLOCA,IRLOCA,IHLOCA)
      RETURN
      END

      SUBROUTINE LOCX(PNX,RNX,HNX,NPX,NRX,NHX,MKX,
     &           PLX,RLX,HLX,IPLOCX,IRLOCX,IHLOCX)
      DIMENSION  IPLOCX(MKX),IRLOCX(MKX),IHLOCX(MKX),
     &           PNX(NPX),RNX(NRX),HNX(NHX)
      PLX      = 1.0-PNX(1)*100.0
      IPMAX    = IFIX(PNX(NPX)*100.0+PLX)
      IP       = 1
      DO 10 I  = 1,IPMAX
      IP1      = IP+1
      PNOW     = (FLOAT(I)-PLX)/100.0
      IF(PNOW.GE.PNX(IP1))   IP = IP1
      IPLOCX(I)= IP
10    CONTINUE
      RLX      = 1.0-RNX(1)*100.0
      IRMAX    = IFIX(RNX(NRX)*100.0+RLX)
      IR       = 1
      DO 20 I  = 1,IRMAX
      IR1      = IR+1
      RNOW     = (FLOAT(I)-RLX)/100.0
      IF(RNOW.GE.RNX(IR1))   IR = IR1
      IRLOCX(I)= IR
20    CONTINUE
      HLX      = 1.0-HNX(1)
      IHMAX    = IFIX(HNX(NHX)+HLX)
      IH       = 1
      DO 30 I  = 1,IHMAX
      IH1      = IH+1
      HNOW     = FLOAT(I)-HLX
      IF(HNOW.GE.HNX(IH1))   IH = IH1
      IHLOCX(I)= IH
30    CONTINUE
      RETURN
      END

      SUBROUTINE READUBPP
      COMMON/DISCON/PN(51),RN(63),DEPA(51,63), 
     &              DEPB(51,63),DEPC(51,63)
      READ(3,100)  NP,NR
      READ(3,110) (PN(I),I=1,NP)
      READ(3,120) (RN(I),I=1,NR)
      DO 1  I = NP,1,-1
      READ(3,130) (DEPA(I,J),J=1,NR)
1     CONTINUE
      DO 2  I = NP,1,-1
      READ(3,130) (DEPB(I,J),J=1,NR)
2     CONTINUE
      DO 3  I = NP,1,-1
      READ(3,130) (DEPC(I,J),J=1,NR)
3     CONTINUE
100   FORMAT(2I6)
110   FORMAT(5(10F7.2/),F7.2)
120   FORMAT(6(10F7.2/),3F7.2)
130   FORMAT(6(10F7.1/),3F7.1)
      RETURN
      END



      SUBROUTINE HLAY(PE,RE,HE,IJK)
      IF(IJK.EQ.1)      THEN
         HE  = 20.0
      ELSE 
         HE  = 50.0
      END IF
      RETURN
      END
 
      SUBROUTINE VEL3(PE,RE,HE,V,LAY,IPS)
      CALL VEL1D(HE,V,LAY,IPS)
c     CALL VEL33(PE,RE,HE,V,LAY,IPS)
      RETURN
      END

      SUBROUTINE VEL33(PE,RE,HE,V,LAY,IPS)
      INCLUDE "PARAM" 
      COMMON/VMOD3D/NPA,NRA,NHA,PNA(MPA),RNA(MRA),HNA(MHA),
     &              DVAP(MPA,MRA,MHA),VELAP(MPA,MRA,MHA)
      COMMON/LOCATE/PLA,RLA,HLA,IPLOCA(MKA),IRLOCA(MKA),IHLOCA(MKA)
      COMMON/WEIGHT/WV(8),IP,JP,KP,IP1,JP1,KP1
      COMMON/PRHFD/P,R,H,PF,RF,HF,PF1,RF1,HF1,PD,RD,HD
      PIDEG = 0.017453
      P     = 90.0-PE/PIDEG
      R     = RE/PIDEG
      H     = HE
      CALL PRHF(IPLOCA,IRLOCA,IHLOCA,PLA,RLA,HLA,
     &          PNA,RNA,HNA,MPA,MRA,MHA,MKA)
      WV(1) = PF1*RF1*HF1
      WV(2) = PF*RF1*HF1
      WV(3) = PF1*RF*HF1
      WV(4) = PF*RF*HF1
      WV(5) = PF1*RF1*HF
      WV(6) = PF*RF1*HF
      WV(7) = PF1*RF*HF
      WV(8) = PF*RF*HF
      CALL VABPS(MPA,MRA,MHA,DVAP,DV)
      CALL VEL1D(HE,V1,LAY,IPS)
      V     = V1*(1.0+DV*0.01)
c      IF(IPS.EQ.2) V = V1
      RETURN
      END

      SUBROUTINE VABPS(MP,MR,MH,V,VEL)
      COMMON/WEIGHT/WV(8),IP,JP,KP,IP1,JP1,KP1
      DIMENSION  V(MP,MR,MH)
      VEL = WV(1)*V(IP,JP,KP)  + WV(2)*V(IP1,JP,KP)
     &    + WV(3)*V(IP,JP1,KP) + WV(4)*V(IP1,JP1,KP)
     &    + WV(5)*V(IP,JP,KP1) + WV(6)*V(IP1,JP,KP1)
     &    + WV(7)*V(IP,JP1,KP1)+ WV(8)*V(IP1,JP1,KP1)
      RETURN
      END
      
	
	SUBROUTINE VEL1D(HE,V,LAY,IPS)
      IF(LAY.EQ.1)      THEN
	  IF(HE.LE.0) V = 4.60
	  IF((HE.GT.0).AND.(HE.LE.3))  V = 4.60+(HE-0)*(5.75-4.60)/3
        IF((HE.GT.3).AND.(HE.LE.20)) V = 5.75+(HE-3)*(6.20-5.75)/17
        IF(IPS.EQ.2) V = V/1.73
     
     
      ELSE IF(LAY.EQ.2) THEN
        V=6.60+(HE-20)*(7.0-6.60)/30
	  IF(IPS.EQ.2) V=  V/1.73
      ELSE
       
        IF(HE.LE.85) V=7.8+(HE-50)*(7.94-7.80)/35
        IF((HE.GT.85).AND.(HE.LE.120))  V=7.94+(HE-85)*(8.08-7.94)/35
	  IF((HE.GT.120).AND.(HE.LE.200)) V=8.08+(HE-120)*(8.39-8.08)/80
        IF((HE.GT.200).AND.(HE.LE.300)) V=8.39+(HE-200)*(8.79-8.39)/100
	  IF((HE.GT.300).AND.(HE.LE.450)) V=8.79+(HE-300)*(9.38-8.79)/150
        IF(IPS.EQ.2) V = V/1.73
      END IF
      RETURN
      END
     



      SUBROUTINE INPUTJB
      COMMON/JBMODV/VP(29),VS(29),RA(29),DEPJ(29)
      DIMENSION VP1(29),VS1(29),RA1(29)
      DATA VP1/7.75, 7.94, 8.13, 8.33, 8.54, 8.75, 8.97,
     &         9.50, 9.91,10.26,10.55,10.99,11.29,11.50, 
     &        11.67,11.85,12.03,12.20,12.37,12.54,12.71,  
     &        12.87,13.02,13.16,13.32,13.46,13.60,13.64,13.64/
      DATA VS1/4.353,4.444,4.539,4.638,4.741,4.850,4.962, 
     &         5.227,5.463,5.670,5.850,6.125,6.295,6.395,  
     &         6.483,6.564,6.637,6.706,6.770,6.833,6.893,  
     &         6.953,7.012,7.074,7.137,7.199,7.258,7.314,7.304/
      DATA RA1/1.00,0.99,0.98,0.97,0.96,0.95,0.94,0.93, 
     &         0.92,0.91,0.90,0.88,0.86,0.84,0.82,0.80,  
     &         0.78,0.76,0.74,0.72,0.70,0.68,0.66,0.64,  
     &         0.62,0.60,0.58,0.56,0.55/
      DO 1 L  = 1,29
      VP(L)   = VP1(L)
      VS(L)   = VS1(L)
      RA(L)   = RA1(L)
      DEPJ(L) = 40.0+6325.59*(1.0-RA1(L))
1     CONTINUE
      RETURN
      END

      SUBROUTINE JBMODEL(IPS,H,V)
      COMMON/JBMODV/VP(29),VS(29),RA(29),DEPJ(29)
      DO 2 K = 1,28
      K1     = K+1
      H1     = DEPJ(K)
      H2     = DEPJ(K1)
      IF(H.GE.H1.AND.H.LT.H2) GO TO 3
2     CONTINUE
3     CONTINUE
      H12    = (H-H1)/(H2-H1)
      IF(IPS.EQ.1)  THEN
         V   = (VP(K1)-VP(K))*H12+VP(K)
      ELSE
         V   = (VS(K1)-VS(K))*H12+VS(K)
      END IF
      RETURN
      END

      SUBROUTINE INTMAP(R,IRLOC,NR,RL,IR)
      DIMENSION IRLOC(NR)
      IS      = IFIX(R+RL)
      IR      = IRLOC(IS)
      RETURN
      END

      SUBROUTINE PRHF(IPLOCX,IRLOCX,IHLOCX,PLX,RLX,HLX,
     &                PNX,RNX,HNX,MPX,MRX,MHX,MKX)
	INCLUDE "PARAM" 
      COMMON/WEIGHT/WV(8),IP,JP,KP,IP1,JP1,KP1
      COMMON/PRHFD/P,R,H,PF,RF,HF,PF1,RF1,HF1,PD,RD,HD
	COMMON/NCELL/PEMIN(NKA),PEMAX(NKA),REMIN(NKA),REMAX(NKA),
     &             HEMIN(NKA),HEMAX(NKA)
	DIMENSION  IPLOCX(MKX),IRLOCX(MKX),IHLOCX(MKX),
     &           PNX(MPX),RNX(MRX),HNX(MHX)
c	ntest=0
	DO 10 I = 1,NC
	NCC=NC
	IF((P.GE.PEMIN(I)).AND.(P.LT.PEMAX(I))) THEN
	  IF((R.GE.REMIN(I)).AND.(R.LT.REMAX(I))) THEN
	     IF((H.GE.HEMIN(I)).AND.(H.LT.HEMAX(I))) THEN
	     M=I
c      ntest=ntest+1
c	if(ntest.ge.2) then
c	 write(*,*) 'erro'
c	WRITE(50,*) P,R,H
c	endif
	     ELSE
	     ENDIF
	  ELSE
	  ENDIF
	ELSE
	ENDIF	
   10 CONTINUE
      RLAT1=PEMIN(M)
	RLAT2=PEMAX(M)
	RLON1=REMIN(M)
	RLON2=REMAX(M)
	RDEP1=HEMIN(M)
	RDEP2=HEMAX(M)

      CALL LIMIT(PNX(1),PNX(MPX),RLAT1)
      CALL LIMIT(RNX(1),RNX(MRX),RLON1)
      CALL LIMIT(HNX(1),HNX(MHX),RDEP1)
	CALL LIMIT(PNX(1),PNX(MPX),RLAT2)
      CALL LIMIT(RNX(1),RNX(MRX),RLON2)
      CALL LIMIT(HNX(1),HNX(MHX),RDEP2)
      CALL INTMAP(RLAT1*100.0,IPLOCX,MKX,PLX,IP)
      CALL INTMAP(RLON1*100.0,IRLOCX,MKX,RLX,JP)
      CALL INTMAP(RDEP1,IHLOCX,MKX,HLX,KP)
	CALL INTMAP(RLAT2*100.0,IPLOCX,MKX,PLX,IP1)
      CALL INTMAP(RLON2*100.0,IRLOCX,MKX,RLX,JP1)
      CALL INTMAP(RDEP2,IHLOCX,MKX,HLX,KP1)
      PD    = PNX(IP1)-PNX(IP)
      RD    = RNX(JP1)-RNX(JP)
      HD    = HNX(KP1)-HNX(KP)
      PF    = (P-PNX(IP))/PD
      RF    = (R-RNX(JP))/RD
      HF    = (H-HNX(KP))/HD
      PF1   = 1.0-PF
      RF1   = 1.0-RF
      HF1   = 1.0-HF
      RETURN
      END

      SUBROUTINE AZIN(P,R,DTH,DAZ,FAZ,FIN)
      EPS  =  1.0E-6
      AZM  = -1.570796-DAZ
      SDTH =  SIN(DTH)
      X    =  SDTH*COS(AZM)
      Y    =  SDTH*SIN(AZM)
      Z    =  COS(DTH)
      SR   =  SIN(R)
      CR   =  COS(R)
      SP   =  SIN(P)
      CP   =  COS(P)
      XG   =  SR*X+CR*CP*Y+CR*SP*Z
      YG   = -CR*X+SR*CP*Y+SR*SP*Z
      ZG   = -SP*Y+CP*Z
      QG   =  SQRT(XG*XG+YG*YG)
      FIN  =  ATAN2(QG,ZG)
      FAZ  =  0.0
      IF(QG.GT.EPS)   FAZ = ATAN2(YG,XG)
      RETURN
      END

      SUBROUTINE CONAZ(RE,PE,RS,PS,DEL,AZ)
      DATA R0,EPS,PI15/6371.05,1.0E-7,4.712389/
      SPS  = SIN(PS)
      CPS  = COS(PS)
      SPE  = SIN(PE)
      CPE  = COS(PE)
      SES  = SIN(RE-RS)
      CES  = COS(RE-RS)
      X    = SPS*SES
      Y    = CPE*SPS*CES - SPE*CPS
      S    = SQRT(X*X + Y*Y)
      C    = SPE*SPS*CES + CPE*CPS
      DEL  = ATAN2(S,C)*R0
      AZ   = 0.0
      AX   = ABS(X)
      AY   = ABS(Y)
      IF(AX.LT.EPS.AND.AY.LT.EPS) RETURN
      AZ   = PI15-ATAN2(Y,X)
      RETURN
      END

      SUBROUTINE VELADJ
      INCLUDE "PARAM" 
      COMMON/CONTRL/NSTS,NEQS,NITLOC,RMSCUT,DVMAX,VDAMP,
     &              NHITCT,NITMAX,RMSTOP,STEPL,TLIM,NITPB
      COMMON/MODINV/A(MG),JA(MG),NA(MD),RDS(MD),KHIT(MU),NEQ4,
     &              JNDEX(MU),NOD,NOU,NOUEQ,ITOT,RNORM,XNORM
      COMMON/GRINET/NPA2,NPRA2,NODESA2,NODETOT
      COMMON/VMOD3D/NPA,NRA,NHA,PNA(MPA),RNA(MRA),HNA(MHA),
     &              DVAP(MPA,MRA,MHA),VELAP(MPA,MRA,MHA)
      COMMON/STDERR/SEVAP(MPA,MRA,MHA) 
      DIMENSION     X(MU),SE(MU),V(MU),W(MU)
c   reorder the elements in A,JA,NA,RDS
      CALL REORDER
c   using CG-type algorithm LSQR to carry out the inversion
      ATOL    = 0.0
      BTOL    = 0.0
      CONLIM  = 0.0
      DAMP    = VDAMP
      ITMAX   = 250
      NOUT    = 8
      CALL LSQR(NOD,NOU,DAMP,RDS,V,W,X,SE,ATOL,BTOL,CONLIM,ITMAX,
     &          NOUT,ISTOP,ANORM,ACOND,RNORM,ARNORM,XNORM)
c   inversion complete
c   apply adjustments to velocity model
      RESSQ   = 0.0
      N1      = NEQ4+1
      N2      = NEQ4+NODETOT
      DO 10 N = N1,N2
      KTN     = KHIT(N)
      IF(KTN.LT.NHITCT)     GO TO 10
      JNN     = JNDEX(N)
c   X is slowness perturbation vector, RH=ds/s=-dv/v
      RH      = X(JNN)
c   standard error of velocity perturbation in %
      SV      = SE(JNN)*100.0
c   calculate P,R and H indeces of velocity grid
      CALL NPNR(NPA2,NPRA2,NEQ4,NP2,NPR2,NVEL)
      NQ      = N-NVEL
      K       = (NQ-1)/NPR2+2
      J       = 2+(NQ-1+(2-K)*NPR2)/NP2
      I       = 1+NQ+NP2*(2-J)+NPR2*(2-K)
      VIJK    = VELAP(I,J,K)
      DVN     = DVAP(I,J,K)
      RH1     = RH/(1.0+RH)
      RESSQ   = RESSQ+RH1*RH1
      DELM    =-VIJK*RH1
      DELMA   = ABS(DELM)
c   place upper bound on velocity adjustment
      IF(DELMA.GT.DVMAX)  DELM = DVMAX*DELM/DELMA
      VIJKD   = VIJK+DELM
      DVND    = DVN+100.0*DELM/VIJKD
c   apply adjustment to model
      DVAP(I,J,K)  = DVND
      VELAP(I,J,K) = VIJKD
      SEVAP(I,J,K) = SV
10    CONTINUE
      RESRMS = XNORM/SQRT(FLOAT(NOU))
      WRITE(8,100)  XNORM,NOU,RESRMS
100   FORMAT(" Norm of solution(vel model) (%)= ",F15.4/,
     &       " Total number of model parameters = ",I5/,
     &       " Total model rms(%) =",F10.2)
      RETURN
      END

      SUBROUTINE NPNR(NPX2,NPRX2,NVX,NP2,NPR2,NVEL)
      NP2    = NPX2
      NPR2   = NPRX2
      NVEL   = NVX
      RETURN
      END

      SUBROUTINE OUTADJ
      INCLUDE "PARAM" 
      COMMON/CONTRL/NSTS,NEQS,NITLOC,RMSCUT,DVMAX,VDAMP,
     &              NHITCT,NITMAX,RMSTOP,STEPL,TLIM,NITPB
      COMMON/VMOD3D/NPA,NRA,NHA,PNA(MPA),RNA(MRA),HNA(MHA),
     &              DVAP(MPA,MRA,MHA),VELAP(MPA,MRA,MHA)
      COMMON/MODINV/A(MG),JA(MG),NA(MD),RDS(MD),KHIT(MU),NEQ4,
     &              JNDEX(MU),NOD,NOU,NOUEQ,ITOT,RNORM,XNORM
      COMMON/STDERR/SEVAP(MPA,MRA,MHA) 
      COMMON/GRINET/NPA2,NPRA2,NODESA2,NODETOT
      COMMON/EVENTS/MINO(MEQ),SECO(MEQ),PHIE(MEQ),RAME(MEQ),
     &              DEPE(MEQ),EVC(3,MEQ),KOBS(MEQ)
      COMMON/STODAT/IYMSTO(MEQ),IDSTO(MEQ),IHRSTO(MEQ)
	COMMON/HITS/MHIT(MNUM),NCELL
      PID    = 0.017453
      SSQR   = RNORM*RNORM
      RMSRES = SQRT(SSQR/FLOAT(NOD))
      WRITE(8,110) SSQR,NOD,RMSRES
c   output the inversion result
      WRITE(8,120)
      WRITE(8,130)
      WRITE(8,140)
      WRITE(8,150)
      NEQ41  = NEQ4-1
      JM     = NEQ41
      CALL OUTA(NPA,NRA,NHA,PNA,RNA,HNA,DVAP,VELAP,SEVAP,JM)
C------------------------------------------------------------
       DO 10  I = 1,NCELL

       WRITE(12,*) MHIT(I)
   10  CONTINUE
C
110   FORMAT(/" Sum of Squared Residuals (sec*sec)=",F12.3/,
     &        " Total number of observations =",I6/,
     &        " Total rms travel time residual (sec)=",F8.5)
120   FORMAT(/" Updated velocity model")
130   FORMAT(" Velocity model difference (%)"/)
140   FORMAT(" Standard error of velocity (%)"/)
150   FORMAT(" Hit conuter "/)
      RETURN
      END

      SUBROUTINE OUTA(NPX,NRX,NHX,PNX,RNX,HNX,DVP,VELP,SEV,JM)
      INCLUDE "PARAM" 
      COMMON/CONTRL/NSTS,NEQS,NITLOC,RMSCUT,DVMAX,VDAMP,
     &              NHITCT,NITMAX,RMSTOP,STEPL,TLIM,NITPB
      COMMON/MODINV/A(MG),JA(MG),NA(MD),RDS(MD),KHIT(MU),NEQ4,
     &              JNDEX(MU),NOD,NOU,NOUEQ,ITOT,RNORM,XNORM
      COMMON/ITER/NIT
	COMMON/NEPICENT/N1,N2,N3,N4,N5,N6,N7
      DIMENSION PNX(NPX),RNX(NRX),HNX(NHX),VELP(NPX,NRX,NHX),
     &          SEV(NPX,NRX,NHX),KITP(MPA,MRA,MHA),DVP(MPA,MRA,MHA)
      NPX1    = NPX-1
      NRX1    = NRX-1
      NHX1    = NHX-1
      NPX2    = NPX-2
      NPRX2   = NPX2*(NRX-2)
      NODESX2 = NPRX2*(NHX-2)
      DO 1  K = 1,NHX
      DO 1  J = 1,NRX
      DO 1  I = 1,NPX
      KITP(I,J,K) = 0
1     CONTINUE
      DO 2  K = 2,NHX1
      DO 2  J = 2,NRX1
      JJ      = JM+(K-2)*NPRX2+(J-2)*NPX2
      DO 2  I = 2,NPX1
      II      = JJ+I
      KITP(I,J,K) =  KHIT(II)
2     CONTINUE
      DO 3  K = 2,NHX1
      WRITE(9,100)  K,HNX(K),(RNX(J),J=1,NRX)
      DO 3  I = NPX,1,-1
      WRITE(9,200)   PNX(I),(VELP(I,J,K),J=1,NRX)
      WRITE(9,300)  (DVP(I,J,K), J=1,NRX)
      WRITE(9,300)  (SEV(I,J,K), J=1,NRX)
      WRITE(9,400)  (KITP(I,J,K),J=1,NRX)
3     CONTINUE
      if(nit.eq.nitmax) then
        write(10,500)  npx,nrx,nhx
        write(10,600) (pnx(i),i=1,npx)
        write(10,600) (rnx(i),i=1,nrx)
        write(10,600) (hnx(i),i=1,nhx)
        do  k=1,nhx
      	do  i=1,npx
      	write(10,300) (dvp(i,j,k),j=1,nrx)
      	enddo
        enddo
      endif

	write(*,*)  n1,n2,n3,n4,n5,n6,n7
 
100   FORMAT("Layer",I2," Depth =",F7.2,"KM"/6(10F7.2/))
200   FORMAT(F5.2/6(10F7.2/))
300   FORMAT(6(10F7.2/))
400   FORMAT(6(10I7/))
500   format(3i3)
600   format(10(10f7.2/))
      RETURN
      END

      SUBROUTINE REORDER
      INCLUDE "PARAM" 
      COMMON/CONTRL/NSTS,NEQS,NITLOC,RMSCUT,DVMAX,VDAMP,
     &              NHITCT,NITMAX,RMSTOP,STEPL,TLIM,NITPB
      COMMON/MODINV/A(MG),JA(MG),NA(MD),RDS(MD),KHIT(MU),NEQ4,
     &              JNDEX(MU),NOD,NOU,NOUEQ,ITOT,RNORM,XNORM
      COMMON/GRINET/NPA2,NPRA2,NODESA2,NODETOT
      DIMENSION     B(MG),JB(MG),NB(MD),RBS(MD)
c     NOD, true number of data;
c     NOU, true number of unknowns;
c     ITOT, true number of elements (all non-zero) in A;
c     JNDEX(i), i-th grid point correspond to
c               which number of the unknowns
      ADIT      = FLOAT(ITOT)/FLOAT(NOD)
      WRITE(8,100) NOD,ITOT,ADIT
100   FORMAT(" Before REORDER: "/,
     &       "   Number of data : ",I8/,
     &       "   Number of non-zero elements in A : ",I15/,
     &       "   Number of non-zero elements/per ray : ",F7.2)
c     To find events with data less than 5, and
c     grids with hitcount smaller NHITCT, these events
c     and grids are not combined in inversion
      NOUEQ     = 0
      IF(NEQ4.EQ.0)            GO TO 11
      DO 10  J  = 1,NEQ4
      KTJ       = KHIT(J)
      IF(KTJ.GE.5)        THEN
        NOUEQ   = NOUEQ+1
        JNDEX(J)= NOUEQ
      ELSE
        JNDEX(J)= 0
      END IF
10    CONTINUE
11    J1        = NEQ4+1
      J2        = NEQ4+NODETOT
      NOU       = NOUEQ
      DO 20  J  = J1,J2
      KTJ       = KHIT(J)
      IF(KTJ.GE.NHITCT)  THEN
        NOU     = NOU+1
        JNDEX(J)= NOU
      ELSE
        JNDEX(J)= 0
      END IF
20    CONTINUE
      NOUGD     = NOU-NOUEQ
c     To winnow those elements in A that correspond to
c     the grids with hitcount smaller than NHITCT
      ID        = 0
      JU        = 0
      L2        = 0
      DO  2  I  = 1,NOD
      L1        = L2+1
      L2        = L2+NA(I)
      JJ        = 0
      DO  3  L  = L1,L2
      JL        = JA(L)
      JNK       = JNDEX(JL)
      IF(JNK.EQ.0)            GO TO 3
      JJ        = JJ+1
      JU        = JU+1
      JB(JU)    = JNK
      B(JU)     = A(L)
3     CONTINUE
      IF(JJ.EQ.0)             GO TO 2
      ID        = ID+1
      NB(ID)    = JJ
      RBS(ID)   = RDS(I)
2     CONTINUE
      NOD       = ID
      ITOT      = JU
      DO  4  I  = 1,NOD
      NA(I)     = NB(I)
      RDS(I)    = RBS(I)
4     CONTINUE
      DO  5  I  = 1,ITOT
      JA(I)     = JB(I)
      A(I)      = B(I)
5     CONTINUE
      NEQ       = NOUEQ/4
      WRITE(8,200) NOD,NOUEQ,NEQ,NOUGD,ITOT
200   FORMAT(" After REORDER :"/,
     &       "   Number of data :",I8/,
     &       "   Number of Eq.  unknowns : ",I8/,
     &       "   Number of Eq.  in A     : ",I5/,
     &       "   Number of grid unknowns : ",I8/,
     &       "   Number of non-zero elements in A : ",I15)
      RETURN
      END

      SUBROUTINE LSQR(M,N,DAMP,U,V,W,X,SE,ATOL,BTOL,CONLIM,ITMAX,
     &                NOUT,ISTOP,ANORM,ACOND,RNORM,ARNORM,XNORM)
      DIMENSION  U(M),V(N),W(N),X(N),SE(N)
      IF(NOUT.GT.0) WRITE(NOUT,1000) M,N,DAMP,ATOL,CONLIM,BTOL,ITMAX
c     initialization
      ZERO    = 0.0
      ONE     = 1.0
      CTOL    = ZERO
      IF(CONLIM.GT.ZERO) CTOL = ONE/CONLIM
      DAMPSQ  = DAMP*DAMP
      ANORM   = ZERO
      ACOND   = ZERO
      BBNORM  = ZERO
      DDNORM  = ZERO
      RES2    = ZERO
      XNORM   = ZERO
      XXNORM  = ZERO
      CS2     =-ONE
      SN2     = ZERO
      Z       = ZERO
      ITN     = 0
      ISTOP   = 0
      NSTOP   = 0
      DO 10 I = 1, N
      V(I)    = ZERO
      X(I)    = ZERO
      SE(I)   = ZERO
10    CONTINUE
      CALL NORMLZ(M,U,BETA)
      CALL APROD(2,M,N,V,U)
      CALL NORMLZ(N,V,ALFA)
      DO 20 I = 1,N
20    W(I)    = V(I)
      RHOBAR  = ALFA
      PHIBAR  = BETA
      BNORM   = BETA
      RNORM   = BETA
      ARNORM  = ALFA*BETA
      IF(ARNORM.LE.ZERO)       GO TO 800
      IF(NOUT.LE.0)            GO TO 100
      IF(DAMPSQ.LE.ZERO)       WRITE(NOUT, 1200)
      IF(DAMPSQ.GT.ZERO)       WRITE(NOUT, 1300)
      TEST1   = ONE
      TEST2   = ALFA/BETA
      WRITE(NOUT,1500) ITN,X(1),RNORM,TEST1,TEST2
      WRITE(NOUT,1600)
c     main iteration loop.
100   ITN     = ITN + 1
c     bidiagonalization
      AF      =-ALFA
      DO 30 I = 1,M
30    U(I)    = AF*U(I)
      CALL APROD(1,M,N,V,U)
      CALL NORMLZ(M,U,BETA)
      BBNORM  = BBNORM+ALFA*ALFA+BETA*BETA+DAMPSQ
      BT      =-BETA
      DO 40 I = 1,N
40    V(I)    = BT*V(I)
      CALL APROD(2,M,N,V,U)
      CALL NORMLZ(N,V,ALFA)
c     modified QR
      RHBAR2  = RHOBAR*RHOBAR+DAMPSQ
      RHBAR1  = SQRT(RHBAR2)
      CS1     = RHOBAR/RHBAR1
      SN1     = DAMP/RHBAR1
      PSI     = SN1*PHIBAR
      PHIBAR  = CS1*PHIBAR
      RHO     = SQRT(RHBAR2+BETA*BETA)
      CS      = RHBAR1/RHO
      SN      = BETA/RHO
      THETA   = SN*ALFA
      RHOBAR  =-CS*ALFA
      PHI     = CS*PHIBAR
      PHIBAR  = SN*PHIBAR
      TAU     = SN*PHI
c     update X, W and the standard error estimates.
      T1      = PHI/RHO
      T2      =-THETA/RHO
      T3      = ONE/RHO
      DO 50 I = 1, N
      T       = W(I)
      X(I)    = T1*T+X(I)
      W(I)    = T2*T+V(I)
      T3T     = T3*T
      T       = T3T*T3T
      SE(I)   = T+SE(I)
      DDNORM  = T+DDNORM
50    CONTINUE
      DELTA   = SN2*RHO
      GAMBAR  =-CS2*RHO
      RHS     = PHI-DELTA*Z
      ZBAR    = RHS/GAMBAR
      XNORM   = SQRT(XXNORM+ZBAR*ZBAR)
      GAMMA   = SQRT(GAMBAR*GAMBAR+THETA*THETA)
      CS2     = GAMBAR/GAMMA
      SN2     = THETA/GAMMA
      Z       = RHS/GAMMA
      XXNORM  = XXNORM + Z*Z
      ANORM   = SQRT(BBNORM)
      ACOND   = ANORM*SQRT(DDNORM)
      RES1    = PHIBAR*PHIBAR
      RES2    = RES2+PSI*PSI
      RNORM   = SQRT(RES1+RES2)
      ARNORM  = ALFA*ABS(TAU)
      TEST1   = RNORM/BNORM
      TEST2   = ARNORM/(ANORM*RNORM)
      TEST3   = ONE/ACOND
      T1      = TEST1/(ONE+ANORM*XNORM/BNORM)
      RTOL    = BTOL+ATOL*ANORM*XNORM/BNORM
      T3      = ONE+TEST3
      T2      = ONE+TEST2
      T1      = ONE+T1
      WRITE(8,1250) ITN,(X(J),J=1,4),RNORM
      write(6,1250) ITN,(X(J),J=1,4),RNORM
1250  FORMAT(I4,4F12.6,F10.3)
      IF(ITN.LT.60)             GO TO 100
      IF(ITN.GE.ITMAX)          ISTOP = 7
      IF(T3.LE.ONE)             ISTOP = 6
      IF(T2.LE.ONE)             ISTOP = 5
      IF(T1.LE.ONE)             ISTOP = 4
      IF(TEST3.LE.CTOL)         ISTOP = 3
      IF(TEST2.LE.ATOL)         ISTOP = 2
      IF(TEST1.LE.RTOL)         ISTOP = 1
      IF(NOUT.LE.0)             GO TO 600
      IF(M.LE.40.OR.N.LE.40)    GO TO 400
      IF(ITN.LE.10)             GO TO 400
      IF(ITN.GE.ITMAX-10)       GO TO 400
      IF(MOD(ITN,10).EQ.0)      GO TO 400
      IF(TEST3.LE.2.0*CTOL)     GO TO 400
      IF(TEST2.LE.10.0*ATOL)    GO TO 400
      IF(TEST1.LE.10.0*RTOL)    GO TO 400
      GO TO 600
400   CONTINUE
c     WRITE(NOUT,1500) ITN,X(1),RNORM,TEST1,TEST2,ANORM,ACOND
      IF(MOD(ITN,10).EQ.0)      WRITE(NOUT, 1600)
600   IF(ISTOP.EQ.0)            NSTOP = 0
      IF(ISTOP.EQ.0)            GO TO 100
      NCONV   = 1
      NSTOP   = NSTOP+1
      IF(NSTOP.LT.NCONV.AND.ITN.LT.ITMAX)  ISTOP = 0
      IF(ISTOP.EQ.0)            GO TO 100
c     end of iteration loop
      T       = ONE
      IF(M.GT.N)                T = M - N
      IF(DAMPSQ.GT.ZERO)        T = M
      T       = RNORM/SQRT(T)
      DO 70 I = 1, N
      SE(I)   = T*SQRT(SE(I))
70    CONTINUE
800   IF(NOUT.LE.0)             GO TO 900
      WRITE(NOUT,1900)          ITN,ISTOP
      IF(ISTOP.EQ.0)            WRITE(NOUT,2000)
      IF(ISTOP.EQ.1)            WRITE(NOUT,2100)
      IF(ISTOP.EQ.2)            WRITE(NOUT,2200)
      IF(ISTOP.EQ.3)            WRITE(NOUT,2300)
      IF(ISTOP.EQ.4)            WRITE(NOUT,2400)
      IF(ISTOP.EQ.5)            WRITE(NOUT,2500)
      IF(ISTOP.EQ.6)            WRITE(NOUT,2600)
      IF(ISTOP.EQ.7)            WRITE(NOUT,2700)
 1000 FORMAT(// 25X,"LSQR   --   Least-squares solution of  A*X = B"
     &   // 25X,"The matrix  A  has ", I6,"rows and ", I6," cols"
     &    / 25X,"The damping parameter is DAMP = ", 1PE10.2
     &   // 25X,"ATOL  = ", 1PE10.2, 10X, " CONLIM = ", 1PE10.2
     &    / 25X,"BTOL  = ", 1PE10.2, 10X, " ITMAX  = ", I10)
 1200 FORMAT(// 3X, 3HITN, 9X, 4HX(1), 14X, 8HFUNCTION, 7X,
     &   "COMPATIBLE INCOMPATIBLE    NORM(A)  COND(A)"/)
 1300 FORMAT(// 3X, 3HITN, 9X, 4HX(1), 14X, 8HFUNCTION, 7X,
     &   "COMPATIBLE INCOMPATIBLE NORM(ABAR) COND(ABAR)"/)
 1500 FORMAT(I6, 1PE20.10, 1PE19.10, 1P2E13.3, 1P2E11.2)
 1600 FORMAT(1X)
 1900 FORMAT(/" No. of iterations =",I6,8X," Stopping condition =",I3)
 2000 FORMAT(/" The exact solution is  X = 0.")
 2100 FORMAT(/" A*X - B  is small enough, given  ATOL, BTOL")
 2200 FORMAT(/" The least-sqrs solu is good enough, given  ATOL")
 2300 FORMAT(/" The estimate of  COND(ABAR)  has exceeded  CONLIM")
 2400 FORMAT(/" A*X - B  is small enough for this machine")
 2500 FORMAT(/" The least-sqrs solu is good enough for this machine")
 2600 FORMAT(/" COND(ABAR)  seems to be too large for this machine")
 2700 FORMAT(/" The iteration limit has been reached")
 900  RETURN
      END

      SUBROUTINE NORMLZ(N,X,S)
      DIMENSION  X(N)
c     normalizes vector X
      EPS    = 1.0E-10
      S      = 0.0
      DO 1 I = 1,N
      XI     = X(I)
1     S      = S+XI*XI
      S      = SQRT(S)
      SS     = 0.0
      IF(S.GT.EPS) SS = 1.0/S
      DO 2 I = 1,N
2     X(I)   = X(I)*SS
      RETURN
      END

      SUBROUTINE APROD(MODE,M,N,X,Y)
      INCLUDE "PARAM" 
      COMMON/MODINV/A(MG),JA(MG),NA(MD),RDS(MD),KHIT(MU),NEQ4,
     &              JNDEX(MU),NOD,NOU,NOUEQ,ITOT,RNORM,XNORM
      DIMENSION  X(N),Y(M)
      ZERO   = 0.0
      L2     = 0
      IF(MODE.NE.1)  GO TO 4
c     MODE = 1 -- SET  Y = Y+A*X.
      DO 1 I = 1,M
      SUM    = ZERO
      L1     = L2+1
      L2     = L2+NA(I)
      DO 2 L = L1,L2
      J      = JA(L)
      SUM    = SUM+A(L)*X(J)
2     CONTINUE
      Y(I)   = Y(I)+SUM
1     CONTINUE
      RETURN
4     CONTINUE
c     MODE = 2 -- SET  X = X+AT*Y
      DO 5 I = 1,M
      YI     = Y(I)
      L1     = L2+1
      L2     = L2+NA(I)
      DO 6 L = L1,L2
      J      = JA(L)
      X(J)   = X(J)+A(L)*YI
6     CONTINUE
5     CONTINUE
      RETURN
      END
      SUBROUTINE TRAVT(PE,RE,DEP,PS,RS,HS,KPS,NSTN,HLAE, 
     &                 LAY,TACC,NITPB,PO,DEL,AZ,VE,FST,Num)
      DIMENSION P(5,200),R(5,200),H(5,200),V(5,200),NS(5),IH(5),IV(5), 
     &   PA(5,200),RA(5,200),HA(5,200),VA(5,200),NSA(5),IHA(5),IVA(5),
     &   PB(5,200),RB(5,200),HB(5,200),VB(5,200),NSB(5),IHB(5),IVB(5),
     &   HLAE(2),HLAS(2)
      INCLUDE "PARAM"
      COMMON/HLASTA/HLSS(MST,2)
      COMMON/RAYPATH/IWV(8,5),R0,PIDEG,EPS
      COMMON/RAYLOC/NRP,RP(3,MDT),IVK(MDT),IWK(MDT)
	COMMON/NEPICENT/N1,N2,N3,N4,N5,N6,N7
      DO  1 I = 1,2
      HLAS(I) = HLSS(NSTN,I)
1     CONTINUE
      CALL CONAZ(RE,PE,RS,PS,DEL,AZ)

	IF(Num.EQ.0)  THEN

	 IF(DEL.GT. 0 .AND.DEL.LE.100)  N1=N1+1
	 IF(DEL.GT.100.AND.DEL.LE.200)  N2=N2+1
	 IF(DEL.GT.200.AND.DEL.LE.300)  N3=N3+1
	 IF(DEL.GT.300.AND.DEL.LE.400)  N4=N4+1
	 IF(DEL.GT.400.AND.DEL.LE.500)  N5=N5+1
	 IF(DEL.GT.500.AND.DEL.LE.600)  N6=N6+1
	 IF(DEL.GT.600)                 N7=N7+1
	 

	ELSE

	ENDIF
c      IF(DEL.GT.800.0) RETURN
      IF(LAY.EQ.1)                    THEN
        CALL BEELINE(PE,RE,DEP,PS,RS,HS,KPS,VE,FST0)
        CALL FIND(PE,RE,DEP,PS,RS,HS,KPS,LAY,1,HLAE,HLAS,
     &            NUMSG,NS,P,R,H,V,IH,IV,TACC,NITPB,FST)
        CALL FIND(PE,RE,DEP,PS,RS,HS,KPS,LAY,2,HLAE,HLAS,NUMSGA, 
     &            NSA,PA,RA,HA,VA,IHA,IVA,TACC,NITPB,FST2)
        
	FSTM  = AMIN1(FST,FST0,FST2)
	IF(FST0.LE.FSTM) THEN
	  FST = FST0
	  RETURN
        END IF
        IF(FST2.LE.FSTM) CALL EQ3(FST2,NUMSGA,NSA,PA,RA,HA,VA,IHA,
     &                            IVA,FST,NUMSG,NS,P,R,H,V,IH,IV)
       
      ELSE IF(LAY.EQ.2)               THEN
        CALL FIND(PE,RE,DEP,PS,RS,HS,KPS,LAY,0,HLAE,HLAS,
     &            NUMSG,NS,P,R,H,V,IH,IV,TACC,NITPB,FST)
        CALL FIND(PE,RE,DEP,PS,RS,HS,KPS,LAY,2,HLAE,HLAS,NUMSGA, 
     &            NSA,PA,RA,HA,VA,IHA,IVA,TACC,NITPB,FST2)
        
	FSTM  = AMIN1(FST,FST2)
        IF(FST2.LE.FSTM) CALL EQ3(FST2,NUMSGA,NSA,PA,RA,HA,VA,IHA,
     &                            IVA,FST,NUMSG,NS,P,R,H,V,IH,IV)
       

      ELSE IF(LAY.EQ.3)               THEN
          CALL FIND(PE,RE,DEP,PS,RS,HS,KPS,LAY,0,HLAE,HLAS,
     &              NUMSG,NS,P,R,H,V,IH,IV,TACC,NITPB,FST)
      ELSE
      END IF
      NA        = NS(NUMSG)
      VE        = V(NUMSG,NA)
      NRP       = 0
      DO  3  I  = 1,NUMSG
      II        = NS(I)-1
      IF(I.EQ.NUMSG)      II = NS(I)
      DO  3  J  = 1,II
      NRP       = NRP+1
      RP(1,NRP) = P(I,J)
      RP(2,NRP) = R(I,J)
      RP(3,NRP) = H(I,J)
      IVK(NRP)  = IV(I)
      IWK(NRP)  = I
3     CONTINUE
      RETURN
      END
      

      

      SUBROUTINE FIND(PE,RE,DEP,PS,RS,HS,KPS,LAY,IHEAD,HLAE,HLAS,
     &                NUMSG,NS,P,R,H,V,IH,IV,TACC,NITPB,FST)
      DIMENSION PG(5,200),RG(5,200),HG(5,200),VG(5,200),
     &          P(5,200),R(5,200),H(5,200),V(5,200),NS(5),
     &          IH(5),IV(5),HLAE(2),HLAS(2)
      CALL INIPATH(PE,RE,DEP,PS,RS,HS,KPS,LAY,IHEAD,
     &             HLAE,HLAS,NUMSG,PG,RG,HG,VG,IH,IV)
      CALL MINIMA(NUMSG,PG,RG,HG,VG,NS,P,R,H,V,IH,IV,
     &             KPS,TACC,NITPB,FST)
      RETURN
      END

      SUBROUTINE INIPATH(PE,RE,DEP,PS,RS,HS,KPS,LAY,IHEAD,
     &                   HLAE,HLAS,NUMSG,P,R,H,V,IH,IV)
      DIMENSION P(5,200),R(5,200),H(5,200),V(5,200),
     &          IH(5),IV(5),HLAE(2),HLAS(2)
      COMMON/RAYPATH/IWV(8,5),R0,PIDEG,EPS
      IF(IHEAD.EQ.0)        THEN
         NUMSG  = LAY
         DO 1 J = 1,NUMSG
         I      = J-1
         IF(J.EQ.1)     H1 = HS
         IF(J.GT.1)     H1 = HLAS(I)
         IF(J.LT.NUMSG) H2 = HLAS(J)
         IF(J.EQ.NUMSG) H2 = DEP
         IPS    = IWV(KPS,J)
         JJ     = J
         CALL VEL3(PS,RS,H1,V1,JJ,IPS)
         CALL VEL3(PS,RS,H2,V2,JJ,IPS)
         CALL EQ1(PS,RS,H1,V1,P(J,1),R(J,1),H(J,1),V(J,1))
         CALL EQ1(PS,RS,H2,V2,P(J,2),R(J,2),H(J,2),V(J,2))
         IF(J.EQ.NUMSG)  THEN
         CALL VEL3(PE,RE,DEP,V2,JJ,IPS)
         CALL EQ1(PE,RE,DEP,V2,P(J,2),R(J,2),H(J,2),V(J,2))
         END IF
         IH(J)  = J
         IV(J)  = J
1        CONTINUE
      ELSE
         NUMSG  = 2*IHEAD-LAY+2
         DO 2 J = 1,IHEAD
         I      = J-1
         IF(J.EQ.1) H1 = HS
         IF(J.GT.1) H1 = HLAS(I)
         H2     = HLAS(J)
         IPS    = IWV(KPS,J)
         JJ     = J
         CALL VEL3(PS,RS,H1,V1,JJ,IPS)
         CALL VEL3(PS,RS,H2,V2,JJ,IPS)
         CALL EQ1(PS,RS,H1,V1,P(J,1),R(J,1),H(J,1),V(J,1))
         CALL EQ1(PS,RS,H2,V2,P(J,2),R(J,2),H(J,2),V(J,2))
         IH(J)  = J
         IV(J)  = J
2        CONTINUE
         IHE    = IHEAD+1
         H1     = HLAS(IHEAD)
         H2     = HLAE(IHEAD)
         IPS    = IWV(KPS,IHE)
         CALL VEL3(PS,RS,H1,V1,IHE,IPS)
         CALL VEL3(PS,RS,H2,V2,IHE,IPS)
         CALL EQ1(PS,RS,H1,V1,P(IHE,1),R(IHE,1),H(IHE,1),V(IHE,1))
         CALL EQ1(PS,RS,H2,V2,P(IHE,2),R(IHE,2),H(IHE,2),V(IHE,2))
         IV(IHE)= IHE
         IH(IHE)= IHEAD
         IHEA   = NUMSG-IHE
         DO 3 J = 1,IHEA
         J1     = IHEAD-J+1
         J2     = J1-1
         J3     = IHE+J
         H1     = HLAE(J1)
         IF(J2.EQ.0)   H2 = DEP
         IF(J2.GT.0)   H2 = HLAE(J2)
         IF(J.EQ.IHEA) H2 = DEP
         IPS    = IWV(KPS,J3)
         CALL VEL3(PE,RE,H1,V1,J1,IPS)
         CALL VEL3(PE,RE,H2,V2,J1,IPS)
         CALL EQ1(PE,RE,H1,V1,P(J3,1),R(J3,1),H(J3,1),V(J3,1))
         CALL EQ1(PE,RE,H2,V2,P(J3,2),R(J3,2),H(J3,2),V(J3,2))
         IF(J2.EQ.0) J2 = 1
         IV(J3) = J1
         IH(J3) = J2
3        CONTINUE
      END IF
      RETURN
      END

      SUBROUTINE MINIMA(NUMSG,PG,RG,HG,VG,NS,P,R,H,V,IH,IV,
     &                  KPS,TACC,NITPB,FSTIME)
      DIMENSION P(5,200),R(5,200),H(5,200),V(5,200),PG(5,200),RG(5,200),
     &          HG(5,200),VG(5,200),NSG(5),NS(5),IH(5),IV(5),SLENG(5)
      COMMON/RAYPATH/IWV(8,5),R0,PIDEG,EPS
      DO 10 I = 1,NUMSG
      NS(I)   = 2
      NSG(I)  = 2
10    CONTINUE
      CALL EQ2(NUMSG,NS,PG,RG,HG,VG,P,R,H,V)
      CALL TRAVEL(NUMSG,NS,P,R,H,V,SLENG,FSTIME)
      LMN     = NUMSG*2
      KMN     = 2
      ITERA   = 0
      GO TO 2
c     -------- double the points in the raypath ------
1     CONTINUE
      ITERA   = ITERA+1
      IF(ITERA.GE.NITPB)                GO TO 9
      FST     = FSTIME
      DO 15 I = 1,2
      DO 15 J = 1,2
      CALL EQ1(PG(I,J),RG(I,J),HG(I,J),VG(I,J),
     &          P(I,J),R(I,J),H(I,J),V(I,J))
15    CONTINUE
      DO 20 J = 1,NUMSG
      NSG(J)  = NS(J)
      IF(IV(J).LE.1)                    GO TO 20
      IF(SLENG(J).LE.5.0)               GO TO 20
      NS(J)   = 2*NSG(J)-1
      LP      = NSG(J)
      LQ      = NS(J)
      CALL EQ1(PG(J,LP),RG(J,LP),HG(J,LP),VG(J,LP),
     &          P(J,LQ), R(J,LQ), H(J,LQ), V(J,LQ))
      IVJ     = IV(J)
      IWJ     = IWV(KPS,J)
      DO 30 I = 1,LP-1
      II      = I+1
      I2      = 2*I
      I3      = 2*I-1
      CALL EQ1(PG(J,I),RG(J,I),HG(J,I),VG(J,I),
     &         P(J,I3),R(J,I3),H(J,I3),V(J,I3))
      CALL MIDPOT(PG(J,I),RG(J,I),HG(J,I),PG(J,II),RG(J,II),
     &            HG(J,II),P2,R2,H2)
      CALL VEL3(P2,R2,H2,V2,IVJ,IWJ)
      CALL EQ1(P2,R2,H2,V2,P(J,I2),R(J,I2),H(J,I2),V(J,I2))
30    CONTINUE
20    CONTINUE
      CALL TRAVEL(NUMSG,NS,P,R,H,V,SLENG,FSTIME)
      DT      = ABS(FST-FSTIME)
      LMN     = 0
      DO 35 I = 1,NUMSG
35    LMN     = LMN+NS(I)
      IF(LMN.GE.200.OR.LMN.GT.50.AND.DT.LE.TACC) GO TO 9
      CALL EQ2(NUMSG,NS,P,R,H,V,PG,RG,HG,VG)
c     --------- perturbation ---------
2     CONTINUE
      KMN     = LMN/NUMSG
      DO 40 J = 1,6
      FST     = FSTIME
      DO 50 I = 1,NUMSG
      IA      = I
      M       = I+1
      II      = NS(I)
      DO 60 K = 2,II
      IF(I.EQ.NUMSG.AND.K.EQ.II)   GO TO 60
      KK      = K-1
      KJ      = K+1
      IF(K.LT.II)           THEN
c     for "continuous points"
      CALL BEND(P(I,KK),R(I,KK),H(I,KK),V(I,KK),P(I,KJ),R(I,KJ),H(I,KJ),
     &          V(I,KJ),P(I,K),R(I,K),H(I,K),V(I,K),IV(I),IA,KPS)
      ELSE
c     for "discontinuous points"
      IVI     = IV(I)
      IVM     = IV(M)
      IF(IVI.GT.IVM)      THEN
         CALL SNELL(P(I,KK),R(I,KK),H(I,KK),V(I,KK),IV(I),IA,
     &              P(M,2),R(M,2),H(M,2),V(M,2),IV(M),M,IH(I),KPS,
     &              P(I,K),R(I,K),H(I,K),V(I,K),VM)
      ELSE IF(IVI.LT.IVM) THEN
         CALL SNELL(P(M,2),R(M,2),H(M,2),V(M,2),IV(M),M,
     &              P(I,KK),R(I,KK),H(I,KK),V(I,KK),IV(I),IA,IH(I),KPS,
     &              P(I,K),R(I,K),H(I,K),VM,V(I,K))
      ELSE
      END IF
      CALL EQ1(P(I,K),R(I,K),H(I,K),VM,P(M,1),R(M,1),H(M,1),V(M,1))
      END IF
60    CONTINUE
50    CONTINUE
      CALL TRAVEL(NUMSG,NS,P,R,H,V,SLENG,FSTIME)
      DT      = FST-FSTIME
      IF(KMN.LE.2.AND.J.LE.5.OR.DT.GT.0.0)  GO TO 3
      FSTIME  = FST
      GO TO 1
3     CALL EQ2(NUMSG,NS,P,R,H,V,PG,RG,HG,VG)
      IF(KMN.LE.2.AND.J.LE.5)    GO TO 40
      IF(DT.LT.TACC)             GO TO 1
40    CONTINUE
      GO TO 1
9     CONTINUE
      RETURN
      END

      SUBROUTINE BEND(A1,B1,H1,V1,A3,B3,H3,V3,A2,B2,H2,V2,IV,IW,KPS)
c   Pseudo-bending for a spherical coordinates.  A1,A2,A3: colatitude
c   in radian; B1,B2,B3: longitude in radian;  H1,H2,H3: depth in km.
      COMMON/RAYPATH/IWV(8,5),R0,PIDEG,EPS
      XFAC  = 1.15
      HMAXM = 1000.0
      IWW  = IWV(KPS,IW)
      R1   = R0-H1
      R3   = R0-H3
      SIA1 = SIN(A1)
      SIA3 = SIN(A3)
      X1   = R1*SIA1*COS(B1)                           
      Y1   = R1*SIA1*SIN(B1)                           
      Z1   = R1*COS(A1)                                   
      X3   = R3*SIA3*COS(B3)                           
      Y3   = R3*SIA3*SIN(B3)                           
      Z3   = R3*COS(A3)                                   
      DX   = X3 - X1                                      
      DY   = Y3 - Y1                                      
      DZ   = Z3 - Z1                                      
      X2   = X1 + DX/2.0
      Y2   = Y1 + DY/2.0
      Z2   = Z1 + DZ/2.0
      R2   = SQRT(X2*X2 + Y2*Y2 + Z2*Z2)
      A2   = ACOS(Z2/R2)                                  
      SINA = SIN(A2)                                    
      B2   = ATAN2(Y2,X2)
      H2   = R0-R2
      CALL VEL3(A2,B2,H2,V2,IV,IWW)
      CALL VELD(A2,B2,H2,VA,VB,VH,IV,IWW)
      VR   = -VH
      DN   = DX*DX + DY*DY + DZ*DZ                    
      DDN  = SQRT(DN)                                    
      DR   = (R3-R1) / DDN                                
      DA   = (A3-A1) / DDN                                
      DB   = (B3-B1) / DDN
      pr   = dr
      pa   = r2 * da
      pb   = r2 * sina * db
      vrd  = pr*vr + pa*va + pb*vb
      rvr  = vr - vrd*pr
      rva  = va - vrd*pa
      rvb  = vb - vrd*pb
      RVS  = SQRT(rvr*rvr + rva*rva + rvb*rvb)
      IF(RVS .EQ. 0.)            GO TO 9
      RVr  = RVr / RVS                                 
      RVa  = RVa / RVS                                 
      RVb  = RVb / RVS                                 
      CC   = (1/V1+1/V3)/2.0
      RCUR = Vr*RVr + Va*RVa + Vb*RVb               
      RCUR = (CC*V2+1) / (4*CC*RCUR)                
      RCUR = -RCUR + SQRT(RCUR*RCUR+DN/(8*CC*V2))     
      rdr  = rvr * rcur
      rda  = rva * rcur
      rdb  = rvb * rcur
      rp   = r2 + rdr
      ap   = a2 + rda/r2
      bp   = b2 + rdb/(r2*sina)
      r2   = (rp-r2)*xfac + r2
      a2   = (ap-a2)*xfac + a2
      b2   = (bp-b2)*xfac + b2
      H2   = R0-R2
      CALL LIMIT(0.0,HMAXM,H2)
      CALL VEL3(A2,B2,H2,V2,IV,IWW)
9     RETURN
      END

      subroutine LENGTH(pe,re,he,ps,rs,hs,ds)
c   spatial distance between two points (pe,re,he) and (ps,rs,hs)
      r0     = 6371.0
      r0e    = r0-he
      r0s    = r0-hs
      sinpe  = sin(pe)
      sinps  = sin(ps)
      xe     = r0e*sinpe*cos(re)
      ye     = r0e*sinpe*sin(re)
      ze     = r0e*cos(pe)
      xs     = r0s*sinps*cos(rs)
      ys     = r0s*sinps*sin(rs)
      zs     = r0s*cos(ps)
      dx     = xe-xs
      dy     = ye-ys
      dz     = ze-zs
      ds     = dx*dx + dy*dy + dz*dz
      ds     = sqrt(ds)
      return
      end

      SUBROUTINE VELD(PE,RE,HE,VP,VR,VH,LAY,IPS)
c     DATA R0,PID,SL1D/6371.0,0.017453,111.19/
      VP    = 0.0
      VR    = 0.0
      VH    = 0.0
      IF(LAY.LE.2)         RETURN
      DH    = 5.0
      H1    = HE-DH
      H2    = HE+DH
      CALL VEL3(PE,RE,H1,V1,LAY,IPS)
      CALL VEL3(PE,RE,H2,V2,LAY,IPS)
      VH    = (V2-V1)/(H2-H1)
c     RHE   = R0-HE
c     DP    = PID*DH/SL1D
c     P1    = PE-DP
c     P2    = PE+DP
c     CALL VEL3(P1,RE,HE,V1,LAY,IPS)
c     CALL VEL3(P2,RE,HE,V2,LAY,IPS)
c     VP    = (V2-V1)/(P2-P1)
c     VP    = VP/RHE
c     DR    = DP
c     R1    = RE-DR
c     R2    = RE+DR
c     CALL VEL3(PE,R1,HE,V1,LAY,IPS)
c     CALL VEL3(PE,R2,HE,V2,LAY,IPS)
c     VR    = (V2-V1)/(R2-R1)
c     VR    = VR/RHE/SIN(PE)
      RETURN
      END

      SUBROUTINE SNELL(PE,RE,HE,VE,IVE,IWE,PS,RS,HS,VS,IVS,IWS,
     &                 IH,KPS,PL,RL,HL,VEL,VSL)
      COMMON/RAYPATH/IWV(8,5),R0,PIDEG,EPS
      PI     = 3.14159265
      CALL MIDPOT(PE,RE,HE,PS,RS,HS,PL,RL,H0)
      CALL HLAY(PL,RL,HL,IH)
      IWVE   = IWV(KPS,IWE)
      IWVS   = IWV(KPS,IWS)
      CALL VEL3(PL,RL,HL,VEL,IVE,IWVE)
      CALL VEL3(PL,RL,HL,VSL,IVS,IWVS)
      CALL CONAZ(RE,PE,RS,PS,DEL,AZ)
      IF(DEL.LT.EPS)       RETURN
      CALL HLAY(PE,RE,HEL,IH)
      CALL HLAY(PS,RS,HSL,IH)
      DH     = HEL-HSL
      BZ     = AZ+PI
      AZS    = AZ
      AZE    = BZ
      IF(DH.GT.0.0)         GO TO 1
      AZS    = BZ
      AZE    = AZ
1     TANP   = ABS(DH)*(1.0+0.5*HL/R0)/DEL
c     dip angle of the interface
      THP    = ATAN(TANP)
      DE     = ABS(HE-HEL)*TANP/(R0-HEL)
      DS     = ABS(HSL-HS)*TANP/(R0-HSL)
c     pedals from two ends to the interface
      CALL AZIN(PE,RE,DE,AZE,REE,PEE)
      CALL AZIN(PS,RS,DS,AZS,RSS,PSS)
      CALL HLAY(PEE,REE,HEE,IH)
      CALL HLAY(PSS,RSS,HSS,IH)
      CALL LENGTH(PEE,REE,HEE,PSS,RSS,HSS,D)
c     lengths of pedal lines from two ends to the interface
      AZES   = AZ
      CTHP   = COS(THP)
      S1     = ABS(HE-HEL)*CTHP
      S2     = ABS(HSL-HS)*CTHP
      SX     = 0.001
      DO 2 K = 1,3
      V1     = (VE+VEL)*0.5
      V2     = (VS+VSL)*0.5
      CALL FASPO(S1,S2,D,V1,V2,X)
      DX     = 0.2*D
      AX     = AMIN1(SX,DX)
      IF(X.LE.AX)  X = AX
      XC     = X*CTHP/(R0-HL)
      CALL AZIN(PEE,REE,XC,AZES,RL,PL)
      CALL LIMIT(PEE,PSS,PL)
      CALL LIMIT(REE,RSS,RL)
      CALL HLAY(PL,RL,HL,IH)
      CALL VEL3(PL,RL,HL,VSL,IVS,IWVS)
      CALL VEL3(PL,RL,HL,VEL,IVE,IWVE)
2     CONTINUE
      RETURN
      END

      subroutine MIDPOT(p1,r1,h1,p3,r3,h3,p2,r2,h2)
      r0    = 6371.0
      call CONAZ(r1,p1,r3,p3,del,az)
      del2  = 0.5*del/r0
      call AZIN(p1,r1,del2,az,r2,p2)
      h2    = (h1+h3)*0.5
      dh    = (r0-h2)*(1.0-cos(del2))
      h2    = h2+dh
      call LIMIT(0.0,2888.0,h2)
      return
      end

      SUBROUTINE LIMIT(C1,C2,C)
      A1   = AMIN1(C1,C2)
      A2   = AMAX1(C1,C2)
      IF(C.LT.A1)  C = A1
      IF(C.GT.A2)  C = A2
      RETURN
      END

      SUBROUTINE TRAVEL(NUMSG,NS,P,R,H,V,SLENG,FST)
      DIMENSION  NS(5),P(5,200),R(5,200),H(5,200),V(5,200),SLENG(5)
      FST      = 0.0
      DO 1  I  = 1,NUMSG
      SLENG(I) = 0.0
      II       = NS(I)
      DO 2  J  = 2,II
      K        = J-1
      CALL LENGTH(P(I,J),R(I,J),H(I,J),P(I,K),R(I,K),H(I,K),DS)
      VIJ      = V(I,J)
      VIK      = V(I,K)
      FST      = FST+DS/VIJ+DS/VIK
      SLENG(I) = SLENG(I)+DS
2     CONTINUE
      SLENG(I) = SLENG(I)/(FLOAT(II)-1.0)
1     CONTINUE
      FST      = FST*0.5
      IF(FST.LT.0.0)  FST = 999.9
      RETURN
      END

      SUBROUTINE FASPO(S1,S2,D,V1,V2,X)
      EPS  = 1.0E-7
      IF(ABS(S1*S2).LT.EPS)  GO TO 5
c     using dichomatical method to solve the equation "fun"
      X1   = 0.0
      X2   = D
      CALL FUN(S1,S2,D,V1,V2,X1,Y)
      NK   = 0
1     NK   = NK+1
      XM   = (X1+X2)*0.5
      CALL FUN(S1,S2,D,V1,V2,XM,F)
      IF(Y*F.LT.0.0)         GO TO 2
      X1   = XM
      GO TO 3
2     X2   = XM
3     DX   = ABS(X2-X1)
      IF(NK.GT.30)           GO TO 4
      IF(DX.GT.0.05)         GO TO 1
4     X    = XM
      RETURN
5     V12  = V1/V2
      IF(V12.GT.1.0)         V12 = 1.0/V12
      DV   = V12/SQRT(1.0-V12*V12)
      IF(ABS(S1).LT.EPS)     X   = D-S2*DV
      IF(ABS(S2).LT.EPS)     X   = S1*DV
      RETURN
      END

      SUBROUTINE FUN(S1,S2,D,V1,V2,X,Y)
c     according to Snell's law:
c     x/sqrt(s1*s1+x*x)/v1 = dx/sqrt(s2*s2+dx*dx)/v2
      DX   = D-X
      D2   = DX*DX
      X2   = X*X
      Y1   = V1*V1*D2*(X2+S1*S1)
      Y2   = V2*V2*X2*(D2+S2*S2)
      Y    = Y1-Y2
      RETURN
      END

      SUBROUTINE BEELINE(PE,RE,HE,PS,RS,HS,KPS,VE,FST)
      INCLUDE "PARAM" 
      COMMON/RAYLOC/NRP,RP(3,MDT),IVK(MDT),IWK(MDT)
      CALL LENGTH(PE,RE,HE,PS,RS,HS,SL)
      STEPL   = 5.0
      NRP     = IFIX(SL/STEPL)+2
      NRP1    = NRP-1
      AB      = 1.0/FLOAT(NRP1)
      DP      = (PE-PS)*AB
      DR      = (RE-RS)*AB
      DH      = (HE-HS)*AB
      FST     = 0.0
      DO 1 K1 = 0,NRP1
      K       = K1+1
      C       = FLOAT(K1)
      RP(1,K) = PS+DP*C
      RP(2,K) = RS+DR*C
      RP(3,K) = HS+DH*C
      IVK(K)  = 1
      IWK(K)  = 1
      CALL VEL3(RP(1,K),RP(2,K),RP(3,K),V,1,KPS)
      IF(K.EQ.1.OR.K.EQ.NRP) V = V*2.0
      FST     = FST+1.0/V
1     CONTINUE
      FST     = FST*SL*AB
      IF(FST.LT.0.0)  FST = 999.9
      VE      = V*0.5
      RETURN
      END

      SUBROUTINE EQ1(PG,RG,HG,VG,P,R,H,V)
      P       = PG
      R       = RG
      H       = HG
      V       = VG
      RETURN
      END

      SUBROUTINE EQ2(NUMSG,NS,PG,RG,HG,VG,P,R,H,V)
      DIMENSION NS(5),P(5,200),R(5,200),H(5,200),V(5,200),
     &          PG(5,200),RG(5,200),HG(5,200),VG(5,200)
      DO 1 I  = 1,NUMSG
      DO 1 J  = 1,NS(I)
      P(I,J)  = PG(I,J)
      R(I,J)  = RG(I,J)
      H(I,J)  = HG(I,J)
      V(I,J)  = VG(I,J)
1     CONTINUE
      RETURN
      END

      SUBROUTINE EQ3(FSTA,NUMSGA,NSA,PA,RA,HA,VA,IHA,IVA,
     &               FST,NUMSG,NS,P,R,H,V,IH,IV)
      DIMENSION PA(5,200),RA(5,200),HA(5,200),VA(5,200),
     &          P(5,200),R(5,200),H(5,200),V(5,200),
     &          NS(5),IH(5),IV(5),NSA(5),IHA(5),IVA(5)
      FST     = FSTA
      NUMSG   = NUMSGA
      DO 1 K  = 1,NUMSG
      NS(K)   = NSA(K)
      IV(K)   = IVA(K)
      IH(K)   = IHA(K)
1     CONTINUE
      CALL EQ2(NUMSG,NS,PA,RA,HA,VA,P,R,H,V)
      RETURN
      END

      SUBROUTINE SMOOTH
c   set a priori parameters to the inversion matricies for models
c   (1st order smoothness)
      INCLUDE "PARAM"
      COMMON/CONTRL/NSTS,NEQS,NITLOC,RMSCUT,DVMAX,VDAMP,
     &              NHITCT,NITMAX,RMSTOP,STEPL,TLIM,NITPB
      COMMON/VMOD3D/NPA,NRA,NHA,PNA(MPA),RNA(MRA),HNA(MHA),
     &             DVAP(MPA,MRA,MHA),VELAP(MPA,MRA,MHA)
    
      COMMON/MODINV/A(MG),JA(MG),NA(MD),RDS(MD),KHIT(MU),NEQ4,
     &              JNDEX(MU),NOD,NOU,NOUEQ,ITOT,RNORM,XNORM
      COMMON/GRINET/NPA2,NPRA2,NODESA2,NODETOT
      IFSMO   =  1
      SIGMP0  =  1
      XSMOHP  =  0.02
      XSMOVP  =  0.02
      PIDEG = 0.01745329
c   c0 = 111.1949;  length of 1 degree in km.
      C0      = 111.1949
c      write(9,100) NOD
      write(6,100) NOD
100   format("Before SMOOTH:  Real number of data used = ",I9/)
      do 10 K = 2,NHA-1
            do 10 I = 2,NPA-1
                   PPN     = (90.0-PNA(I))*PIDEG
                   SWT     = 1.0/sqrt(sin(PPN))
                      do 10 j = 2,nra-1
                          ali1    = 1.0/abs(pna(i)-pna(i-1))/c0
                          ali2    = 1.0/abs(pna(i+1)-pna(i))/c0
                          alj1    = 1.0/abs(rna(j)-rna(j-1))/c0
                          alj2    = 1.0/abs(rna(j+1)-rna(j))/c0
                          alk1    = 1.0/abs(hna(k)-hna(k-1))
                          alk2    = 1.0/abs(hna(k+1)-hna(k))
                          alp     = ali1+ali2
                          alr     = alj1+alj2
                          alh     = alk1+alk2
                          npara    = 0
                          do 20 ii = i-1,i+1
                          do 20 jj = j-1,j+1
                          do 20 kk = k-1,k+1
                             if(ii.eq.1.or.ii.eq.npa)  go to 20
                             if(jj.eq.1.or.jj.eq.nra)  go to 20
                             if(kk.eq.1.or.kk.eq.nha)  go to 20
                             ijk      = (ii-1)+(jj-2)*npa2+(kk-2)*npra2
                             n1       = neq4+ijk
                             if(ii.eq.i.and.jj.eq.j)      then
                                if(kk.eq.k)           then
                                   itot    = itot+1
                            a(itot) =-((alp+alr*swt)/xsmohp+alh/xsmovp)
                                   npara   = npara+1
                                  else if(kk.eq.(k-1))  then
                                  itot    = itot+1
                                  a(itot) = alk1/xsmovp
                                  npara   = npara+1
                                  else if(kk.eq.(k+1))  then
                                  itot    = itot+1
                                   a(itot) = alk2/xsmovp
                                   npara   = npara+1
                                    else
                                      go to 20
                                   end if
                                    else if(jj.eq.j.and.kk.eq.k) then
                                   if(ii.eq.(i-1))       then
                                   itot    = itot+1
                                   a(itot) = ali1/xsmohp
                                   npara   = npara+1
                                   else if(ii.eq.(i+1))  then
                                   itot    = itot+1
                                   a(itot) = ali2/xsmohp
                                   npara   = npara+1
                                    else
                                     go to 20
                                    end if
                                    else if(ii.eq.i.and.kk.eq.k) then
                                    if(jj.eq.(j-1))       then
                                     itot    = itot+1
                                      a(itot) = alj1*swt/xsmohp
                                      npara   = npara+1
                                      else if(jj.eq.(j+1))  then
                                            itot    = itot+1
                                            a(itot) = alj2*swt/xsmohp
                                         npara   = npara+1
                                      else
                                      go to 20
                                      end if
                                     else
                                      go to 20
                                       end if
                                     ja(itot)    = n1
20              continue
                nod         = nod+1
                na(nod)     = npara
               rds(nod)    = 0.0
10    continue
      write(6,200) nod
200   format("After smooth1: Real number of data used = ",I9/)
      return
      end

	subroutine iaspi (He,VEL,IPS)
	
	logical table
      character*1 answer,name*50
      common /tab/ table
      common /file/ luout
      data luout, iopt /6,0/
      data table /.false./
      lout=6

  
c      write
c     >  (6,'(/," Enter depth [km] (negative depth stops program): ",$)')
c         read(5,*) z
c         if(z.lt.0) stop "program finished"
         r=6371.-he
c         if(z.gt.6371.) then
c            write(6,'(" Unrealistic value for depth: ",
c     >                "give new value: ",$)')
c            read(5,*) z
c            if(z.lt.0) stop "program finished"
cc            r=6371.-z
c         endif
         
c         write(luout,'(" radius [km]  depth [km]  v-p [km/sec]",
c     >          "  v-s [km/sec]     rho")')
      xr=r/6371.                  
     
      call mod2(r,xr,a,b,p)
c      call print(r,a,b,p)
c     call print(r,a,b)
c      go to 1000
      IF(IPS.EQ.1) THEN
      vel=a
      ELSE
	VEL=b
	endif
      end
     
     
      subroutine mod2(r,xr,a,b,p)

c iasp91 model: 
c B.L.N. Kennett and E. R. Engdahl, Geophys. J. In., 105, 429-466, 1991
      
      dimension rad(0:11)
      common /file/ luout
      data rad /0.0,1217.1,3482.0,3631.0,5611.0,5711.0,5961.0,
     >          6161.0,6251.0,6336.0,6351.0,6371.0/


      if(r.eq.0.0) then
          call value(1,p,a,b,xr,r,rad(1))
          return
          endif
      do 1 ii=0,11
1     if(r.gt.rad(ii).and.r.le.rad(ii+1)) 
     >    call value(ii+1,p,a,b,xr,r,rad(ii+1))
      return
      end
      
       subroutine value(i,p,a,b,xr,r,r1)

      common /file/ luout
      common /tab/ table
 
      logical table,log(11),first
      data first /.true./

      if(first) then
        do 100 j=1,11
100     log(j)=.true.
        first=.false.
      endif
	

c      if(table.and.log(i)) call discon(i,r,r1)
      log(i)=.false.
      
         p=rho2(i,xr)
         a=pvel2(i,xr)
         b=svel2(i,xr)
      return
      end
     
     
     
      function pvel2(i,xr)

c computes P-velocities for IASP90 model

      goto (1,2,3,4,5,6,7,8,9,10,11) i
1     pvel2=11.24094-(4.09689*xr**2)
      return
2     pvel2=10.03904+(3.75665*xr)-(13.67046*xr**2)
      return
3     pvel2=14.49470-(1.47089*xr)
      return
4     pvel2=25.1486 -(41.1538*xr)+(51.9932*xr**2)-(26.6083*xr**3)
      return
5     pvel2=25.96984-(16.93412*xr)
      return
6     pvel2=29.38896-(21.40656*xr)
      return
7     pvel2=30.78765-(23.25415*xr)
      return
8     pvel2=25.41389-(17.69722*xr)
      return
9     pvel2= 8.78541-(0.74953*xr)
      return
10    pvel2= 6.50
      return
11    pvel2= 5.80
      return
      end
c
      function svel2(i,xr)

c computes S-velocities for iasp91 model:

      goto (1,2,3,4,5,6,7,8,9,10,11) i
1     svel2= 3.56454-(3.45241*xr**2)
      return
2     svel2= 0.
      return
3     svel2= 8.16616-(1.58206*xr)
      return
4     svel2=12.9303 -(21.2590*xr)+(27.8988*xr**2)-(14.1080*xr**3)
      return
5     svel2=20.76890-(16.53147*xr)
      return
6     svel2=17.70732-(13.50652*xr)
      return
7     svel2=15.24213-(11.08552*xr)
      return
8     svel2= 5.75020-(1.27420*xr)
      return
9     svel2= 6.706231-(2.248585*xr)
      return
10    svel2= 3.75
      return
11    svel2= 3.36
      return
      end
      
      function rho2(i,xr)

c computation of iasP91 densities not supported in present version

      rho2=0
      return

      goto (1,2,3,4,5,6,7,8,9,10,11) i
1     rho2=13.01219-(8.45292*xr**2)
      return
2     rho2=12.58416-(1.69929*xr)-(1.94128*xr**2)-(7.11215*xr**3)
      return
3     rho2= 6.8143 -(1.66273*xr)-(1.18531*xr**2)
      return
4     rho2= 6.8143 -(1.66273*xr)-(1.18531*xr**2)
      return
5     rho2= 6.8143 -(1.66273*xr)-(1.18531*xr**2)
      return
6     rho2=11.11978-( 7.87054*xr)
      return
7     rho2= 7.15855-( 3.85999*xr)
      return
8     rho2= 7.15855-( 3.85999*xr)
      return
9     rho2= 7.15855-(3.85999*xr)
      return
10    rho2= 2.92
      return
11    rho2= 2.72
      return
      end
      
c      subroutine print(r,a,b,p)
	subroutine print(r,a,b)

      common /file/ luout
      write(luout,10) r,6371.-r,a,b
      return

10    format(2(3x,f8.1),2(5x,f9.4)) 
   
      end
c
c
      SUBROUTINE HIT(A,B,C)
	INCLUDE "PARAM"
      COMMON/VMOD3D/NPA,NRA,NHA,PNA(MPA),RNA(MRA),HNA(MHA),
     &              DVAP(MPA,MRA,MHA),VELAP(MPA,MRA,MHA)
	COMMON/HITS/MHIT(MNUM),NCELL
	ILAT = NPA-3
	ILON = NRA-3
	IDEP = NHA-3
	NCELL = ILAT*ILON*IDEP
      PIDEG = 0.017453
      RLAT = 90-A/PIDEG
	RLON = B/PIDEG

	DO 1 I = 2,NPA-2
      IF((RLAT.GE.PNA(I)).AND.(RLAT.LT.PNA(I+1))) THEN
	M = I-1
	ENDIF
   1  CONTINUE

      DO 2 J = 2,NRA-2
      IF((RLON.GE.RNA(J)).AND.(RLON.LT.RNA(J+1))) THEN
	N = J-1
	ENDIF
   2  CONTINUE

      DO 3 K = 2,NHA-2
      IF((C.GE.HNA(K)).AND.(C.LT.HNA(K+1))) THEN
	L = K-1
	ENDIF
   3  CONTINUE
     

	JJ=M+(N-1)*ILAT+(L-1)*ILAT*ILON
	MHIT(JJ)=MHIT(JJ)+1

	END

	SUBROUTINE CELL(RLAT,RLON,RH)
	INCLUDE "PARAM"
	COMMON/NCELL/PEMIN(NKA),PEMAX(NKA),REMIN(NKA),REMAX(NKA),
     &             HEMIN(NKA),HEMAX(NKA)
	COMMON/VMOD3D/NPA,NRA,NHA,PNA(MPA),RNA(MRA),HNA(MHA),
     &              DVAP(MPA,MRA,MHA),VELAP(MPA,MRA,MHA)
      COMMON/LOCATE/PLA,RLA,HLA,IPLOCA(MKA),IRLOCA(MKA),IHLOCA(MKA)
	DO 10 I = 1,NC
	IF((RLAT.GE.PEMIN(I)).AND.(RLAT.LT.PEMAX(I))) THEN
	  IF((RLON.GE.REMIN(I)).AND.(RLON.LT.REMAX(I))) THEN
	     IF((RH.GE.HEMIN(I)).AND.(RH.LT.HEMAX(I))) THEN
	     M=I
	     ELSE
	     ENDIF
	  ELSE
	  ENDIF
	ELSE
	ENDIF	
   10 CONTINUE
      RLAT1=PEMIN(M)
	RLAT2=PEMAX(M)
	RLON1=REMIN(M)
	RLON2=REMAX(M)
	RDEP1=HEMIN(M)
	RDEP2=HEMAX(M)
	END