C***************************************************************
C 20-12-01
C D=4 SU(NCOL) code with  Q,  J=0,2   P,C=+,+  glueballs 
C N-ality = 1,2,3 flux loops
C
C*****
C
C blocking uses products of 2 smeared links :
C itype different parameter sets
C 
C WARNING: POLY VACS ARE REAL -- NO GOOD FOR FINITE T!
C
C***************************************************************
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER(LX1=64,LX2=24,LX3=24,LX4=24)
        PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
        PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
        PARAMETER(ICALLG=1,NITER=200,NMEASL=NITER/ICALLG)
        PARAMETER(ICMIN=1,ICMAX=NITER+ICMIN-1)
c        PARAMETER(ICMIN=1,ICMAX=1)
        PARAMETER(IBING=8,NUMBIN=25)
        PARAMETER(IBLOK=5,ITYPE=1)
        PARAMETER(IDIAG=1)
        PARAMETER(MAXDTLS=12,MAXDTG=12,MAXDTLL=2)
        PARAMETER(NF=2)
        PARAMETER(IENDIAN=1,ISWAPTX=1)
        PARAMETER(IIMAX=1)
C
C        PARAMETER(ICALLQ=50,NCOOL=20)
C        PARAMETER(QCUT=0.02,RRCUT=4.0,MAX=999,MAXP=99)
C
        COMMON/ARRAYS/U11(NCOL2*LSIZE*4)
        COMMON/PARAM/PARBS(ITYPE),PARBDS(ITYPE)
        REAL*8 rndnum
        COMPLEX U11
C
        CALL STARTF
        ISEED=19752
c        REWIND(40)
c        READ(40) U11,ISEED
        CALL RINI(ISEED)
C
        open(file='kappa.in',status='old',unit=11)
        READ (11,*) dkappa
        close(11)

        CALL SETUP
        CALL MOMENTM
C
        PARBS(1)=0.40
        PARBDS(1)=0.16 
c        PARBS(1)=0.30
c        PARBDS(1)=0.12
c        PARBS(2)=0.375
c        PARBDS(2)=0.15
c        PARBS(3)=0.45
c        PARBDS(3)=0.18
C
        BETAG= 2.25d0
C
        RATIO=1.0
c        RATIO=0.25
        USTADPOLE=1.0d0
        UTTADPOLE=1.0d0
        UST2=USTADPOLE*USTADPOLE 
        UTT2=UTTADPOLE*UTTADPOLE
        SRAT=RATIO/(UST2*UST2)
        TRAT=1./(RATIO*UST2*UTT2)
C
c        write (*,*) SRAT, TRAT
c        stop 111
C
        IHEAT=0
c        IHEAT=10000
        ITOT=NITER
C
        WRITE(6,90)
90      FORMAT(' ****************************************************')
        WRITE(6,91)
91      FORMAT(' *')
        WRITE(6,92) LX1,LX2,LX3,LX4,BETAG,NCOL
92      FORMAT('  LX = ',4I6,'    BETA =',F8.4,'   NCOLOR = ',I4)
        WRITE(6,91)
        WRITE(6,96) RATIO
96      FORMAT('  at/as  tree level  =',F8.4)
        WRITE(6,91)
        WRITE(6,93) IHEAT,ITOT
93      FORMAT('  NUM HEATS =',I6,'   NUM ITER =',I6)
        WRITE(6,91)
        WRITE(6,94) ICALLG,IBING
94      FORMAT('  SWEEPS PER MEAS =',I4,'   MEAS PER BIN =',I4)
        WRITE(6,91)
        WRITE(6,951) PARBS
951     FORMAT('  PAR BLSMEAR:  AXIS  =',8F8.4)
        WRITE(6,952) PARBDS
952     FORMAT('  PAR BLSMEAR:  DIAG  =',8F8.4)
        WRITE(6,91)
        IF(IDIAG.NE.1) WRITE(6,953)
953     FORMAT('  DIAG SMEARING NOT USED ')
        WRITE(6,91)
        WRITE(6,90)
C
c        ACTLL=0.0d0
c        IBIN=IHEAT/10
c        DO 10 ITER=1,IHEAT
c        IBIAS=2
c        IF(ITER.LE.200) IBIAS=1
c        IF(ITER.EQ.(ITER/5 )*5 ) IBIAS=1
c        CALL UPDATG(IBIAS,ACTL,BETAG,SRAT,TRAT)
c        CALL RENORM
c        ACTLL=ACTLL+ACTL
c        IF(ITER.EQ.(ITER/IBIN)*IBIN)THEN
c        NBIN=ITER/IBIN
c        ACTLL=ACTLL/IBIN
c        WRITE(6,75) NBIN,ACTLL
c75      FORMAT('  NBIN =',I4,' ACTLL HEATS =',F12.6)
c        ACTLL=0.0d0
c        ENDIF
c10      CONTINUE
C        goto999
C
 211    format('[FM][0]Check plaq = ', f8.6)
        IFILE = ICMIN - 1
        DO 20 ITER=1,ITOT
           IFILE = IFILE + 1
c        IBIAS=2
c        IF(ITER.EQ.(ITER/5 )*5 ) IBIAS=1
c        CALL UPDATG(IBIAS,ACTL,BETAG,SRAT,TRAT)
c        IF(ITER.EQ.5*(ITER/5))CALL RENORM
c           CALL fortran_read_old_gauge(u11,LX4,LX1,LX2,LX3,NCOL,NF,BETAG
c     &          ,DKAPPA,IFILE,0)
           CALL fortran_read_new_gauge(u11,LX4,LX1,LX2,LX3,NCOL,NF,BETAG
     &          ,DKAPPA,IFILE,IENDIAN,ISWAPTX)
           CALL ACTION(0,ITER,TOTACT)
           write (*,211) TOTACT 
           CALL POLYLOOP()
c           if (iter.eq.(icmin+30)) goto 998
           write (*,*) ' ' 
           CALL MEASURE(ITER,ITOT)
C        IF(ITER.EQ.(ITER/ICALLQ)*ICALLQ) CALL QCOOL
20      CONTINUE
c
           CALL ACTION(1,0,TOTACT)
           CALL TODISK
c999     ISEED = INT(rndnum()*float(259200))
c        REWIND(40)
c        WRITE(40) U11,ISEED
 998    continue
C
        STOP
        END
C*********************************************************************
C  actn(nb,ip) : average of trace(plaq)/2 for meas in bin 'nb'
C                for plaqs in plane  ip=6-mu-nu
C
C acorl(nb, mom, dist, bl , bl):  corrln  therm lines
C avacf(nb, bl ) :                averages therm lines
C*********************************************************************
        SUBROUTINE TODISK
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER(IBLOK=5,ITYPE=1)
        PARAMETER(NUMBIN=25)
        PARAMETER(LX1=64,LX2=24,LX3=24,LX4=24)
        PARAMETER(LMAX=LX4/2+1)
C
        COMMON/PBLOCK/ABCORL1(NUMBIN,7,LMAX,IBLOK*IBLOK,ITYPE)
     &  ,ABVACL1(NUMBIN,IBLOK,ITYPE)
     &  ,ABCORL2(NUMBIN,7,2*2,LMAX,IBLOK*IBLOK,ITYPE)
     &  ,ABVACL2(NUMBIN,2,IBLOK,ITYPE)
     &  ,ABCORL3(NUMBIN,7,3*3,LMAX,IBLOK*IBLOK,ITYPE)
     &  ,ABVACL3(NUMBIN,3,IBLOK,ITYPE)
C
        COMMON/XSTORE/AVAC0(NUMBIN,2,IBLOK,ITYPE)
     &  ,ACOR0(NUMBIN,2*2,LMAX,IBLOK*IBLOK,ITYPE)
     &  ,ACORE(NUMBIN,6*6,LMAX,IBLOK*IBLOK,ITYPE)
     &  ,ACORT(NUMBIN,3*3,LMAX,IBLOK*IBLOK,ITYPE)
C
        COMMON/ASTORE/ACTN(NUMBIN,6)
C
        COMPLEX*16 ABCORL1,ABCORL2,ABCORL3
C
        REWIND(50)
        WRITE(50) ACTN
        WRITE(50) ABCORL1,ABVACL1
        WRITE(50) ABCORL2,ABVACL2
C        WRITE(50) ABCORL3,ABVACL3
        WRITE(50) ACOR0,AVAC0,ACORE,ACORT
C
        RETURN
        END
C***************************************************************
C***************************************************************
        SUBROUTINE QCOOL
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER(LX1=64,LX2=24,LX3=24,LX4=24)
        PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
        PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
        PARAMETER(ICALLQ=50,NCOOL=20)
        PARAMETER(QCUT=0.02,RRCUT=4.0,MAX=999,MAXP=99)
C
c        COMMON/ARRAYS/U11(NCOL2*LSIZE*4)
c        DIMENSION V11(NCOL2*LSIZE*4)
c        COMPLEX U11,V11
C
c        DO 10 NS=1,NCOL2*LSIZE*4
c        V11(NS)=U11(NS)
c10      CONTINUE
C
c        IBIAS=3
c        DO 20 IC=1,NCOOL
c        CALL UPDATG(IBIAS,ACTL,BETAG,SRAT,TRAT)
c        CALL FFDACT
c        IF(IC.EQ.NCOOL) CALL SEARCH
c20      CONTINUE
C
c        DO 11 NS=1,NCOL2*LSIZE*4
c        U11(NS)=V11(NS)
c11      CONTINUE
C
        RETURN
        END
C******************************************************************
C
C*********************************************************************
        SUBROUTINE MEASURE(ITER,NTOT)
        IMPLICIT REAL*8 (A-H,O-Z)
C
        PARAMETER(ICALLG=1,IBING=8,NUMBIN=25)
C
        COMMON/ITEM/ITERG,JBING,NTOTG
C
c        write (*,*) 'in measure'
        ITERG=ITER/ICALLG
        JBING=IBING
        NTOTG=IBING*NTOT/(ICALLG*IBING)
C
        IF(ITER.EQ.(ITERG*ICALLG).AND.ITERG.LE.NTOTG) THEN
c           write (*,*) 'measure:1'
        CALL SETUPB
        CALL BLOCK
        ENDIF
C
        IF(ITER.EQ.NTOTG*ICALLG) THEN
        CALL MBPLAN
        CALL MBLINE
        ENDIF
C
        RETURN
        END
C*********************************************************************
C Print correlations and energies from THERML
C*********************************************************************
        SUBROUTINE MBLINE
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER(LX1=64,LX2=24,LX3=24,LX4=24,IBLOK=5,ITYPE=1)
        PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
        PARAMETER(LMAX=LX4/2+1,MAXDTLS=12,NUMBIN=25)
C
        COMMON/PBLOCK/ACORL1(NUMBIN,7,LMAX,IBLOK*IBLOK,ITYPE)
     &  ,AVACL1(NUMBIN,IBLOK,ITYPE)
     &  ,ACORL2(NUMBIN,7,2*2,LMAX,IBLOK*IBLOK,ITYPE)
     &  ,AVACL2(NUMBIN,2,IBLOK,ITYPE)
     &  ,ACORL3(NUMBIN,7,3*3,LMAX,IBLOK*IBLOK,ITYPE)
     &  ,AVACL3(NUMBIN,3,IBLOK,ITYPE)
        COMMON/ITEM/ITR,IBIN,ITOT
        COMMON/FIT/ACOR(200),SCOR(200),AMM(200),SMM(200),NTMAX
C
        COMPLEX*16 ACORL1,ACORL2,ACORL3
        DIMENSION COR(NUMBIN,LMAX),VAC(NUMBIN),LS(3),AV(NUMBIN)
C
        NBIN=ITOT/IBIN
        LS(1)=LX1
        LS(2)=LX2
        LS(3)=LX3
        EPS=1.0d-10
C
        DO 1 NY=1,1
C
        DO 2 ITY=1,ITYPE
C
        WRITE(6,80)
        WRITE(6,80)
80      FORMAT(' ****************************************************')
        WRITE(6,81)
        WRITE(6,90) NY,ITY
90      FORMAT('      POLYAKOV   LOOPS :  NALITY = ',I4,'   ITY =',I5)
        WRITE(6,81)
        WRITE(6,81)
        WRITE(6,80)
        WRITE(6,80)
        WRITE(6,81)
81      FORMAT(' *')
        WRITE(6,91)
91      FORMAT(' AVERAGE  SPATIAL  LINES ')
        WRITE(6,81)
82      FORMAT(' *************')
C
        IEMAX=NY
        DO 10 IE=1,IEMAX
        DO 10 ID=1,IBLOK
        ACOUNT=3.0*IBIN*LSIZEB/LS(1)
        DO 13 IB=1,NUMBIN
        IF(NY.EQ.1)AV(IB)=AVACL1(IB,ID,ITY)/ACOUNT
        IF(NY.EQ.2)AV(IB)=AVACL2(IB,IE,ID,ITY)/ACOUNT
        IF(NY.EQ.3)AV(IB)=AVACL3(IB,IE,ID,ITY)/ACOUNT
13      CONTINUE
        CALL JACKM(NUMBIN,AV,APLAQ,SPLAQ)
        IF(NY.EQ.1)WRITE(6,110) IE,ID,APLAQ,SPLAQ
110     FORMAT(' IE,  BL =',2I3,'  AV,ER LINE =',2F9.4)
        IF(NY.EQ.2.AND.IE.EQ.1)WRITE(6,111) IE,ID,APLAQ,SPLAQ
111     FORMAT(' IE,  BL =',2I3,'  AV,ER LINE: TrL.TrL =',2F9.4)
        IF(NY.EQ.2.AND.IE.EQ.2)WRITE(6,112) IE,ID,APLAQ,SPLAQ
112     FORMAT(' IE,  BL =',2I3,'  AV,ER LINE: TrL.L =',2F9.4)
        IF(NY.EQ.3.AND.IE.EQ.1)WRITE(6,113) IE,ID,APLAQ,SPLAQ
113     FORMAT(' IE,  BL =',2I3,'  AV,ER LINE: TrL.TrL.TrL =',2F9.4)
        IF(NY.EQ.3.AND.IE.EQ.2)WRITE(6,114) IE,ID,APLAQ,SPLAQ
114     FORMAT(' IE,  BL =',2I3,'  AV,ER LINE: TrL.L.TrL =',2F9.4)
        IF(NY.EQ.3.AND.IE.EQ.3)WRITE(6,115) IE,ID,APLAQ,SPLAQ
115     FORMAT(' IE,  BL =',2I3,'  AV,ER LINE: TrL.L.L =',2F9.4)
10      CONTINUE
C
        DO 20 NPP=1,4
        MOMSQ=NPP-1+NPP/4
        IF(NPP.EQ.1)NPL=1
        IF(NPP.EQ.1)NPU=1
        IF(NPP.EQ.2)NPL=2
        IF(NPP.EQ.2)NPU=3
        IF(NPP.EQ.3)NPL=4
        IF(NPP.EQ.3)NPU=5
        IF(NPP.EQ.4)NPL=6
        IF(NPP.EQ.4)NPU=7
        WRITE(6,81)
        WRITE(6,82)
        WRITE(6,93) MOMSQ
93      FORMAT(' MOM SQ = ',I4)
        WRITE(6,82)
        IEMAX=NY
        IE12=0
        DO 20 IE2=1,IEMAX
        DO 20 IE1=1,IEMAX
        IE12=IE12+1
        IF(IE1.NE.IE2)GOTO20
        ID12=0
        DO 22 ID2=1,IBLOK
        DO 22 ID1=1,IBLOK
        ID12=ID12+1
        IF(ID1.NE.ID2)GOTO22
        WRITE(6,81)
        WRITE(6,82)
        WRITE(6,94) IE1,ID1
94      FORMAT(' IE=',I4,'   BLOCKING LEVEL = ',I4)
        WRITE(6,82)
        DO 24 NT=1,LMAX
        DO 24 IB=1,NBIN
        SUM=0.0
        DO 28 NP=NPL,NPU
        IF(NY.EQ.1)SUM=SUM+REAL(ACORL1(IB,NP,NT,ID12,ITY))
        IF(NY.EQ.2)SUM=SUM+REAL(ACORL2(IB,NP,IE12,NT,ID12,ITY))
        IF(NY.EQ.3)SUM=SUM+REAL(ACORL3(IB,NP,IE12,NT,ID12,ITY))
28      CONTINUE
        COR(IB,NT)=SUM
24      CONTINUE
        IF(ABS(COR(1,1)).LE.EPS) THEN
           write (*,*) 'COR(1,1) =',ABS(COR(1,1)),'< EPS =', EPS 
           GOTO20
        ENDIF
        CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,COR,VAC)
        DO 26 NT=1,MAXDTLS
        WRITE(6,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
102     FORMAT('  DT=',I3,'   AV,ER COR = ',2F8.4,'    E=',2F8.4)
26      CONTINUE
22      CONTINUE
20      CONTINUE
C
2       CONTINUE
C
1       CONTINUE
C
        RETURN
        END
C*********************************************************************
C*********************************************************************
        SUBROUTINE MBPLAN
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER(LX1=64,LX2=24,LX3=24,LX4=24,IBLOK=5,ITYPE=1)
        PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
        PARAMETER(LMAX=LX4/2+1,MAXDTG=12,NUMBIN=25)
        PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
        COMMON/XSTORE/AVAC0(NUMBIN,2,IBLOK,ITYPE)
     &  ,ACOR0R(NUMBIN,2,2,LMAX,IBLOK,IBLOK,ITYPE)
     &  ,ACORER(NUMBIN,6,6,LMAX,IBLOK,IBLOK,ITYPE)
     &  ,ACORTR(NUMBIN,3,3,LMAX,IBLOK,IBLOK,ITYPE)
C
        COMMON/ITEM/ITR,IBIN,ITOT
        COMMON/FIT/ACOR(200),SCOR(200),AMM(200),SMM(200),NTMAX
C
        DIMENSION COR(NUMBIN,LMAX),VAC(NUMBIN),AV(NUMBIN)
C
        NBIN=ITOT/IBIN
        EPS=0.0000000001d0
C
80      FORMAT(' ****************************************************')
81      FORMAT(' *****************************************')
83      FORMAT(' ******')
82      FORMAT(' *')
C
98      FORMAT(' BLOCKING LEVEL = ',I4)
102     FORMAT('  DT=',I3,'   AV,ER COR = ',2F10.5,'  M =',2F8.4)
C
        DO 100 ITY=1,ITYPE
C
        WRITE(6,80)
        WRITE(6,99) ITY
99      FORMAT('    GLUEBALL  CORRELATIONS -  PLANAR:  ITY =',I5)
        WRITE(6,80)
        WRITE(6,82)
C
        WRITE(6,83)
        WRITE(6,91)
91      FORMAT(' AVERAGE  SUPERPLAQUETTES')
        WRITE(6,83)
C
        ACOUNT=IBIN*LSIZE*3*NCOL
        DO 1  ID=1,IBLOK
        DO 2  IB=1,NUMBIN
        AV(IB)=AVAC0(IB,1,ID,ITY)/ACOUNT
2       CONTINUE
        CALL JACKM(NUMBIN,AV,AVER,ERR)
        WRITE(6,92)  ID,AVER,ERR
92      FORMAT(' BLOCKING =',I4,'  AV,ER PLAQ = ',2F12.6)
1       CONTINUE
C
        WRITE(6,81)
        WRITE(6,201)
201     FORMAT('      0++    Glueball     Correlations    ')
        WRITE(6,81)
        WRITE(6,83)
901     FORMAT('   1 x 1 loop  ')
902     FORMAT('   1 x 2 loop  ')
C
        ANORMB=LX4*IBIN
        WRITE(6,82)
        DO 20 IOP=1,2
        WRITE(6,82)
        IF(IOP.EQ.1)WRITE(6,901)
        IF(IOP.EQ.2)WRITE(6,902)
        WRITE(6,82)
        DO 20 ID=1,IBLOK
        WRITE(6,83)
        WRITE(6,98) ID
        WRITE(6,83)
        WRITE(6,82)
        DO 12 NT=1,LMAX
        DO 12 IB=1,NBIN
        COR(IB,NT)=ACOR0R(IB,IOP,IOP,NT,ID,ID,ITY)/ANORMB
        VAC(IB)=AVAC0(IB,IOP,ID,ITY)/ANORMB
12      CONTINUE
        IF(ABS(COR(1,1)).LE.EPS)GOTO15
        CALL JACK(NBIN,0,LX4,NUMBIN,LMAX,COR,VAC)
        DO 14 NT=1,MAXDTG
        WRITE(6,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
14      CONTINUE
15      CONTINUE
20      CONTINUE
10      CONTINUE
C
        WRITE(6,81)
        WRITE(6,202)
202     FORMAT('      2++ (E)   Glueball     Correlations    ')
        WRITE(6,81)
        WRITE(6,83)
C
        DO 30 IOP=1,2
        WRITE(6,82)
        IF(IOP.EQ.1)WRITE(6,901)
        IF(IOP.EQ.2)WRITE(6,902)
        WRITE(6,82)
        DO 30 ID=1,IBLOK
        WRITE(6,83)
        WRITE(6,98) ID
        WRITE(6,83)
        IF(IOP.EQ.1)THEN
        DO 31 NT=1,LMAX
        DO 31 IB=1,NBIN
        COR(IB,NT)=
     &  ACORER(IB,1,1,NT,ID,ID,ITY)+ACORER(IB,2,2,NT,ID,ID,ITY)
31      CONTINUE
        ENDIF
        IF(IOP.EQ.2)THEN
        DO 33 NT=1,LMAX
        DO 33 IB=1,NBIN
        COR(IB,NT)=
     &  ACORER(IB,3,3,NT,ID,ID,ITY)+ACORER(IB,4,4,NT,ID,ID,ITY)
     &  +ACORER(IB,5,5,NT,ID,ID,ITY)+ACORER(IB,6,6,NT,ID,ID,ITY)
33      CONTINUE
        ENDIF
        IF(ABS(COR(1,1)).LE.EPS)GOTO30
        CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,COR,VAC)
        DO 34 NT=1,MAXDTG
        WRITE(6,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
34      CONTINUE
        WRITE(6,83)
30      CONTINUE
C
        WRITE(6,81)
        WRITE(6,203)
203     FORMAT('      2++ (T)   Glueball     Correlations    ')
        WRITE(6,81)
        WRITE(6,83)
C
        WRITE(6,82)
        WRITE(6,902)
        WRITE(6,82)
        DO 40 ID=1,IBLOK
        WRITE(6,83)
        WRITE(6,98) ID
        WRITE(6,83)
        DO 41 NT=1,LMAX
        DO 41 IB=1,NBIN
        COR(IB,NT)=ACORTR(IB,1,1,NT,ID,ID,ITY)
     &  +ACORTR(IB,2,2,NT,ID,ID,ITY)+ACORTR(IB,3,3,NT,ID,ID,ITY)
41      CONTINUE
        IF(ABS(COR(1,1)).LE.EPS)GOTO40
        CALL JACK(NBIN,1,LX4,NUMBIN,LMAX,COR,VAC)
        DO 42 NT=1,MAXDTG
        WRITE(6,102) NT-1,ACOR(NT),SCOR(NT),AMM(NT),SMM(NT)
42      CONTINUE
        WRITE(6,83)
40      CONTINUE
C
100     CONTINUE
C
        RETURN
        END
C*********************************************************************
C create blocked links by multiplying together two smeared 
C links of the previous blocking level -- uses SMEAR1 which smears
C UC11 links once using pointers put in /ASMEAR1/
C*********************************************************************
        SUBROUTINE BLOCK
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER(LX1=64,LX2=24,LX3=24,LX4=24)
        PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
        PARAMETER(IBLOK=5,ITYPE=1)
        PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
        COMMON/ARRAYS/U11(NCOL2,LSIZEB,LX4,4)
        COMMON/ARRAYB/UB11(NCOL2,LSIZEB,3,IBLOK)
        COMMON/ASMEAR1/UC11(NCOL2,LSIZEB,3),
     &  IUP(LSIZEB,3),IDN(LSIZEB,3)
        COMMON/NEXTB/IUPB(LSIZEB,3,IBLOK+1),IDNB(LSIZEB,3,IBLOK+1)
        COMMON/PARAM/PARBS(ITYPE),PARBDS(ITYPE)
C
        DIMENSION A11(NCOL2),B11(NCOL2),C11(NCOL2)
        COMPLEX U11,UB11,A11,B11,C11,UC11
C
        DO 1 I4=1,LX4
        DO 2 ITY=1,ITYPE
        DO 3 IBL=1,IBLOK
C
        IF(IBL.EQ.1)THEN
        DO 8 MU=1,3
        DO 8 NN=1,LSIZEB
        DO 8 IJ=1,NCOL2
        UB11(IJ,NN,MU,1)=U11(IJ,NN,I4,MU)
8       CONTINUE
        GOTO11
        ENDIF
C
        IBLM=IBL-1
        DO 14 KK=1,3
        DO 14 NN=1,LSIZEB
        IUP(NN,KK)=IUPB(NN,KK,IBLM)
        IDN(NN,KK)=IDNB(NN,KK,IBLM)
        DO 14 IJ=1,NCOL2
        UC11(IJ,NN,KK)=UB11(IJ,NN,KK,IBLM)
14      CONTINUE
        PARB=PARBS(ITY)
        PARBD=PARBDS(ITY)
        CALL SMEAR1(PARB,PARBD)
C
        DO 20 MU=1,3
        DO 20 NN=1,LSIZEB
        M1=NN
        M2=IUP(M1,MU)
C
        DO 22 IJ=1,NCOL2
        A11(IJ)=UC11(IJ,M1,MU)
        B11(IJ)=UC11(IJ,M2,MU)
22      CONTINUE
        CALL VMX(1,A11,B11,C11,1)
        DO 24 IJ=1,NCOL2
        UB11(IJ,NN,MU,IBL)=C11(IJ)
24      CONTINUE
C
20      CONTINUE
C
11      CALL PLANARB(I4,IBL,ITY)
        CALL THERMLB(I4,IBL,ITY)
C
3       CONTINUE
2       CONTINUE
1       CONTINUE
C
        DO 40 ITY=1,ITYPE 
        CALL PLNCORB(ITY)
        CALL POTB(ITY)
40      CONTINUE
C
        RETURN
        END
C*********************************************************************
C smear whole time-slice of (possibly blocked and/or smeared)
C links UC11 using given pointers, then output as UC11
C*********************************************************************
        SUBROUTINE SMEAR1(PARA,PARD)
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER(LX1=64,LX2=24,LX3=24,LX4=24)
        PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
        PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
        PARAMETER(IDIAG=1)
C
        COMMON/ASMEAR1/UC11(NCOL2,LSIZEB,3),
     &  IUP(LSIZEB,3),IDN(LSIZEB,3)
        COMMON/DIAGIN/UUC11(NCOL2,LSIZEB,3),
     &  IIUP(LSIZEB,3),IIDN(LSIZEB,3)
        COMMON/DIAGOUT/UDD(NCOL2,LSIZEB,4),IDD(LSIZEB,4)
        DIMENSION UCC11(NCOL2,LSIZEB,3)
        DIMENSION UINT11(NCOL2),UREN11(NCOL2)
        DIMENSION A11(NCOL2),B11(NCOL2),C11(NCOL2),DUM11(NCOL2)
        COMPLEX A11,B11,C11,UINT11,UREN11,UC11,DUM11,UCC11,UDD
        COMPLEX UUC11
C
        DO 11 KK=1,3
        DO 11 NN=1,LSIZEB
        IIUP(NN,KK)=IUP(NN,KK)
        IIDN(NN,KK)=IDN(NN,KK)
11      CONTINUE
        DO 12 MU=1,3
        DO 12 NN=1,LSIZEB
        DO 12 IJ=1,NCOL2
        UUC11(IJ,NN,MU)=UC11(IJ,NN,MU)
12      CONTINUE
C
        DO 1 MU=1,3
        IF(IDIAG.EQ.1) CALL DIAG(MU)
C
        DO 2 NN=1,LSIZEB
        M1=NN
C
        DO 10 IJ=1,NCOL2
        UINT11(IJ)=UC11(IJ,M1,MU)
10      CONTINUE
C
        DO 40 NU=1,3
        IF(MU.EQ.NU)GOTO40
C
        DO 26 IJ=1,NCOL2
        A11(IJ)=UC11(IJ,M1,NU)
26      CONTINUE
        M2=IUP(M1,NU)
        DO 27 IJ=1,NCOL2
        B11(IJ)=UC11(IJ,M2,MU)
27      CONTINUE
        CALL VMX(1,A11,B11,C11,1)
        M3=IUP(M1,MU)
        DO 28 IJ=1,NCOL2
        B11(IJ)=UC11(IJ,M3,NU)
28      CONTINUE
        CALL HERM(1,B11,DUM11,1)
        CALL VMX(1,C11,B11,A11,1)
        DO 21 IJ=1,NCOL2
        UINT11(IJ)=UINT11(IJ)+A11(IJ)*PARA
21      CONTINUE
C
        M2=IDN(M1,NU)
        DO 31 IJ=1,NCOL2
        A11(IJ)=UC11(IJ,M2,NU)
31      CONTINUE
        DO 32 IJ=1,NCOL2
        B11(IJ)=UC11(IJ,M2,MU)
32      CONTINUE
        CALL HERM(1,A11,DUM11,1)
        CALL VMX(1,A11,B11,C11,1)
        M3=IUP(M2,MU)
        DO 33 IJ=1,NCOL2
        B11(IJ)=UC11(IJ,M3,NU)
33      CONTINUE
        CALL VMX(1,C11,B11,A11,1)
        DO 35 IJ=1,NCOL2
        UINT11(IJ)=UINT11(IJ)+A11(IJ)*PARA
35      CONTINUE
C
40      CONTINUE
C
        DO 50 MNU=1,4
        IF(IDIAG.NE.1)GOTO50
C
        DO 42 IJ=1,NCOL2
        A11(IJ)=UDD(IJ,M1,MNU)
42      CONTINUE
        M2=IDD(M1,MNU)
        DO 43 IJ=1,NCOL2
        B11(IJ)=UC11(IJ,M2,MU)
43      CONTINUE
        CALL VMX(1,A11,B11,C11,1)
        M3=IUP(M1,MU)
        DO 44 IJ=1,NCOL2
        B11(IJ)=UDD(IJ,M3,MNU)
44      CONTINUE
        CALL HERM(1,B11,DUM11,1)
        CALL VMX(1,C11,B11,A11,1)
        DO 45 IJ=1,NCOL2
        UINT11(IJ)=UINT11(IJ)+A11(IJ)*PARD
45      CONTINUE
C
50      CONTINUE
C
        NBCOOL=2
        CALL RENORMBS(UINT11,UREN11,NBCOOL)
        DO 60 IJ=1,NCOL2
        UCC11(IJ,M1,MU)=UREN11(IJ)
60      CONTINUE
C
2       CONTINUE
1       CONTINUE
C
        DO 100 MU=1,3
        DO 100 NN=1,LSIZEB
        DO 100 IJ=1,NCOL2
        UC11(IJ,NN,MU)=UCC11(IJ,NN,MU)
100     CONTINUE
C
        RETURN
        END
C*********************************************************************
C form 4 'diagonal' links in plane orth to MU, from each site 
C with end-pt in IDD using links and pointers in /DIAGIN/
C -- output matrices (not suN) and pointers in /DIAGOUT/
C*********************************************************************
        SUBROUTINE DIAG(MU)
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER(LX1=64,LX2=24,LX3=24,LX4=24)
        PARAMETER(LSIZEB=LX1*LX2*LX3)
        PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
        COMMON/DIAGOUT/UDD(NCOL2,LSIZEB,4),IDD(LSIZEB,4)
        COMMON/DIAGIN/UC11(NCOL2,LSIZEB,3),
     &  IUP(LSIZEB,3),IDN(LSIZEB,3)
        COMPLEX C11(NCOL2),F11(NCOL2),DUM11(NCOL2)
        COMPLEX A11(NCOL2),D11(NCOL2),E11(NCOL2)
        COMPLEX AA11(NCOL2),DD11(NCOL2),EE11(NCOL2)
        COMPLEX UDD,UC11
C
        NU=MU+1
        IF(MU.EQ.3)NU=1
        KU=6-MU-NU
C
        DO 1 NN=1,LSIZEB
        M1=NN
C
        DO 2 IJ=1,NCOL2
        A11(IJ)=UC11(IJ,M1,NU)
2       CONTINUE
        M2=IUP(M1,NU)
        DO 3 IJ=1,NCOL2
        AA11(IJ)=UC11(IJ,M2,KU)
3       CONTINUE
        CALL VMX(1,A11,AA11,C11,1)
        DO 4 IJ=1,NCOL2
        D11(IJ)=UC11(IJ,M1,KU)
4       CONTINUE
        M3=IUP(M1,KU)
        DO 5 IJ=1,NCOL2
        DD11(IJ)=UC11(IJ,M3,NU)
5       CONTINUE
        CALL VMX(1,D11,DD11,F11,1)
        DO 6 IJ=1,NCOL2
        UDD(IJ,M1,1)=C11(IJ)+F11(IJ)
6       CONTINUE
        M9=IUP(M3,NU)
        IDD(M1,1)=M9
C
        DO 61 IJ=1,NCOL2
        AA11(IJ)=C11(IJ)+F11(IJ)
61      CONTINUE
        CALL HERM(1,AA11,DUM11,1)
        DO 62 IJ=1,NCOL2
        UDD(IJ,M9,3)=AA11(IJ)
62      CONTINUE
        IDD(M9,3)=M1
C
        M4=IDN(M2,KU)
        DO 7 IJ=1,NCOL2
        AA11(IJ)=UC11(IJ,M4,KU)
7       CONTINUE
        CALL HERM(1,AA11,DUM11,1)
        CALL VMX(1,A11,AA11,C11,1)
        M5=IDN(M1,KU)
        DO 8 IJ=1,NCOL2
        E11(IJ)=UC11(IJ,M5,KU)
8       CONTINUE
        CALL HERM(1,E11,DUM11,1)
        DO 9 IJ=1,NCOL2
        EE11(IJ)=UC11(IJ,M5,NU)
9       CONTINUE
        CALL VMX(1,E11,EE11,F11,1)
        DO 10 IJ=1,NCOL2
        UDD(IJ,M1,2)=C11(IJ)+F11(IJ)
10      CONTINUE
        IDD(M1,2)=M4
C
        DO 63 IJ=1,NCOL2
        AA11(IJ)=C11(IJ)+F11(IJ)
63      CONTINUE
        CALL HERM(1,AA11,DUM11,1)
        DO 64 IJ=1,NCOL2
        UDD(IJ,M4,4)=AA11(IJ)
64      CONTINUE
        IDD(M4,4)=M1
C
1       CONTINUE
C
        RETURN
        END
C**********************************************************
C OK
C  THIS ROUTINE MAKES BLOCKED MATRICES SUN IN CRUDEST
C  POSSIBLE WAY - 1 COOL
C**********************************************************
        SUBROUTINE RENORMBS(UU1,UREN11,NBCOOL)
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
        DIMENSION UU1(NCOL,NCOL),UREN11(NCOL,NCOL)
        COMPLEX UB1,UU1,CSUM,ADUM(NCOL,NCOL),UREN11
        COMPLEX A11(NCOL2),B11(NCOL2),C11(NCOL2)
C
        DO 1 N2=1,NCOL
        DO 1 N1=1,NCOL
        ADUM(N1,N2)=UU1(N1,N2)
1       CONTINUE
C
        DO 20 N2=1,NCOL
C
        DO 10 N3=1,N2-1
C
        CSUM=(0.0,0.0)
        DO 5 N1=1,NCOL
        CSUM=CSUM+ADUM(N1,N2)*CONJG(ADUM(N1,N3))
5       CONTINUE
        DO 6 N1=1,NCOL
        ADUM(N1,N2)=ADUM(N1,N2)-CSUM*ADUM(N1,N3)
6       CONTINUE
C
10      CONTINUE
C
        SUM=0.0
        DO 7 N1=1,NCOL
        SUM=SUM+ADUM(N1,N2)*CONJG(ADUM(N1,N2))
7       CONTINUE
        ANORM=1.0/SQRT(SUM)
        DO 8 N1=1,NCOL
        ADUM(N1,N2)=ADUM(N1,N2)*ANORM
8       CONTINUE
C
20      CONTINUE
C
        DO 21 IK=1,NCOL
        DO 21 IJ=1,NCOL
        IADD=NCOL*(IK-1)
        A11(IJ+IADD)=ADUM(IJ,IK)
21      CONTINUE
C
        DO 22 J=1,NCOL2
        C11(J)=A11(J)
22      CONTINUE
        DO 24 J=1,NCOL
        DO 24 I=1,NCOL
        IJ=I+(J-1)*NCOL
        A11(IJ)=CONJG(UU1(J,I))
24      CONTINUE
C
      DO 30 IBC=1,NBCOOL
C
      CALL VMX(1,A11,C11,B11,1)
C
      IF(NCOL.EQ.2)THEN
      CALL SUBGRB(1,4,2,3,B11,C11)
      ENDIF
      IF(NCOL.EQ.3)THEN
      CALL SUBGRB(1,5,2,4,B11,C11)
      CALL SUBGRB(5,9,6,8,B11,C11)
      CALL SUBGRB(1,9,3,7,B11,C11)
      ENDIF
      IF(NCOL.EQ.4)THEN
      CALL SUBGRB(1,6,2,5,B11,C11)
      CALL SUBGRB(6,11,7,10,B11,C11)
      CALL SUBGRB(11,16,12,15,B11,C11)
      CALL SUBGRB(1,16,4,13,B11,C11)
      CALL SUBGRB(6,16,8,14,B11,C11)
      CALL SUBGRB(1,11,3,9,B11,C11)
      ENDIF
      IF(NCOL.EQ.5)THEN
      CALL SUBGRB(1,7,2,6,B11,C11)
      CALL SUBGRB(7,13,8,12,B11,C11)
      CALL SUBGRB(13,19,14,18,B11,C11)
      CALL SUBGRB(19,25,20,24,B11,C11)
      CALL SUBGRB(7,19,9,17,B11,C11)
      CALL SUBGRB(13,25,15,23,B11,C11)
      CALL SUBGRB(1,13,3,11,B11,C11)
      CALL SUBGRB(1,19,4,16,B11,C11)
      CALL SUBGRB(7,25,10,22,B11,C11)
      CALL SUBGRB(1,25,5,21,B11,C11)
      ENDIF
      IF(NCOL.EQ.6)THEN
      CALL SUBGRB(1,8,2,7,B11,C11)
      CALL SUBGRB(8,15,9,14,B11,C11)
      CALL SUBGRB(15,22,16,21,B11,C11)
      CALL SUBGRB(22,29,23,28,B11,C11)
      CALL SUBGRB(29,36,30,35,B11,C11)
      CALL SUBGRB(1,15,3,13,B11,C11)
      CALL SUBGRB(8,22,10,20,B11,C11)
      CALL SUBGRB(15,29,17,27,B11,C11)
      CALL SUBGRB(22,36,24,34,B11,C11)
      CALL SUBGRB(1,22,4,19,B11,C11)
      CALL SUBGRB(8,29,11,26,B11,C11)
      CALL SUBGRB(15,36,18,33,B11,C11)
      CALL SUBGRB(1,29,5,25,B11,C11)
      CALL SUBGRB(8,36,12,32,B11,C11)
      CALL SUBGRB(1,36,6,31,B11,C11)
      ENDIF
C
30    CONTINUE
C
      DO 500 J=1,NCOL
      DO 500 I=1,NCOL
      IJ=I+(J-1)*NCOL
      UREN11(I,J)=C11(IJ)
500   CONTINUE
C
      RETURN
      END
C*********************************************************************
      SUBROUTINE SUBGRB(I1,I2,I3,I4,B11,C11)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
      COMPLEX B11(NCOL2),C11(NCOL2)
      COMPLEX F11,F12,A11(NCOL2),S11(NCOL2)
C
      F11=(B11(I1)+CONJG(B11(I2)))*0.5
      F12=(B11(I3)-CONJG(B11(I4)))*0.5
      UMAG=SQRT(F11*CONJG(F11)+F12*CONJG(F12))
      UMAG=1./UMAG
      F11=F11*UMAG
      F12=F12*UMAG
      A11(I1)=CONJG(F11)
      A11(I4)=CONJG(F12)
      A11(I3)=-F12
      A11(I2)=F11
      DO 3 IJ2=1,NCOL
      DO 3 IJ1=1,NCOL
      IJ=IJ1+NCOL*(IJ2-1)
      IF(IJ.EQ.I1.OR.IJ.EQ.I2.OR.IJ.EQ.I3.OR.IJ.EQ.I4)GOTO3
      IF(IJ1.EQ.IJ2)THEN
      A11(IJ)=(1.0,0.0)
      ELSE
      A11(IJ)=(0.0,0.0)
      ENDIF
3     CONTINUE
      CALL VMX(1,C11,A11,S11,1)
      CALL VMX(1,B11,A11,C11,1)
      DO 5  IJ=1,NCOL2
      B11(IJ)=C11(IJ)
5     CONTINUE
      DO 6  IJ=1,NCOL2
      C11(IJ)=S11(IJ)
6     CONTINUE
C
      RETURN
      END
C*********************************************************************
C OK
C Here calculate for the simple blocked 1x1 plaquette
C and the simple blocked 1x2 loop
C   PL0R(loop label,time,block level)           : 0++
C*********************************************************************
        SUBROUTINE PLANARB(I4,IBL,ITY)
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER(LX1=64,LX2=24,LX3=24,LX4=24,IBLOK=5,ITYPE=1)
        PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
        PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
        COMMON/PLAN/PL0R(2,LX4,IBLOK,ITYPE),PLER(6,LX4,IBLOK,ITYPE)
     &  ,PLTR(3,LX4,IBLOK,ITYPE)
        COMMON/ARRAYB/UB11(NCOL2,LSIZEB,3,IBLOK)
        COMMON/NEXTB/IUPB(LSIZEB,3,IBLOK+1),IDNB(LSIZEB,3,IBLOK+1)
        COMMON/DUMMY/UC11(NCOL2,LSIZEB,3)
     &  ,IUP(LSIZEB,3),IDN(LSIZEB,3)
        DIMENSION SUM(3),ASUM(3,3),AOP(6)
        COMPLEX A11(NCOL2),B11(NCOL2),C11(NCOL2),DUM11(NCOL2),AKT
        COMPLEX UC11,UB11
        COMPLEX*16 PL0R,PLER,PLTR
C
        ID=IBL
C
        DO 2 KK=1,3
        DO 2 NN=1,LSIZEB
        IUP(NN,KK)=IUPB(NN,KK,ID)
        IDN(NN,KK)=IDNB(NN,KK,ID)
2       CONTINUE
        DO 3 KK=1,3
        DO 3 NN=1,LSIZEB
        DO 3 IJ=1,NCOL2
        UC11(IJ,NN,KK)=UB11(IJ,NN,KK,ID)
3       CONTINUE
C
        DO 4 IP=1,2
        PL0R(IP,I4,ID,ITY)=0.0
4       CONTINUE
        DO 5 IP=1,6
        PLER(IP,I4,ID,ITY)=0.0
5       CONTINUE
        DO 6 IP=1,3
        PLTR(IP,I4,ID,ITY)=0.0
6       CONTINUE
C
        DO 8 IDIR=1,3
        SUM(IDIR)=0.0
8       CONTINUE
C
        DO 10 NN=1,LSIZEB
        DO 10 MU=1,2
        DO 10 NU=MU+1,3
        M1=NN
        IDIR=6-MU-NU
C
        DO 11 IJ=1,NCOL2
        A11(IJ)=UC11(IJ,M1,NU)
        B11(IJ)=UC11(IJ,M1,MU)
11      CONTINUE
        CALL HERM(1,A11,DUM11,1)
        CALL VMX(1,A11,B11,C11,1)
        M2=IUP(M1,MU)
        DO 12 IJ=1,NCOL2
        A11(IJ)=UC11(IJ,M2,NU)
12      CONTINUE
        CALL VMX(1,C11,A11,B11,1)
        M2=IUP(M1,NU)
        DO 13 IJ=1,NCOL2
        A11(IJ)=UC11(IJ,M2,MU)
13      CONTINUE
        CALL HERM(1,A11,DUM11,1)
        CALL TRVMX(1,B11,A11,AKT,1)
C
        SUM(IDIR)=SUM(IDIR)+REAL(AKT)
C
10      CONTINUE
C
        PL0R(1,I4,ID,ITY)=SUM(1)+SUM(2)+SUM(3)
        PLER(1,I4,ID,ITY)=SUM(2)-SUM(3)
        PLER(2,I4,ID,ITY)=SUM(2)+SUM(3)-2.0*SUM(1)
C
        DO 18 NU=1,3
        DO 18 MU=1,3
        ASUM(MU,NU)=0.0
18      CONTINUE
C
        DO 50 NN=1,LSIZEB
        M1=NN
        DO 20 MU=1,3
        DO 20 NU=1,3
        IF(MU.EQ.NU)GOTO20
C
        DO 21 IJ=1,NCOL2
        A11(IJ)=UC11(IJ,M1,NU)
        B11(IJ)=UC11(IJ,M1,MU)
21      CONTINUE
        CALL HERM(1,A11,DUM11,1)
        CALL VMX(1,A11,B11,C11,1)
        M2=IUP(M1,MU)
        DO 22 IJ=1,NCOL2
        A11(IJ)=UC11(IJ,M2,NU)
22      CONTINUE
        CALL VMX(1,C11,A11,B11,1)
        M3=IUP(M2,NU)
        DO 23 IJ=1,NCOL2
        A11(IJ)=UC11(IJ,M3,NU)
23      CONTINUE
        CALL VMX(1,B11,A11,C11,1)
        M3=IUP(M1,NU)
        M2=IUP(M3,NU)
        DO 24 IJ=1,NCOL2
        A11(IJ)=UC11(IJ,M2,MU)
24      CONTINUE
        CALL HERM(1,A11,DUM11,1)
        CALL VMX(1,C11,A11,B11,1)
        DO 25 IJ=1,NCOL2
        A11(IJ)=UC11(IJ,M3,NU)
25      CONTINUE
        CALL HERM(1,A11,DUM11,1)
        CALL TRVMX(1,B11,A11,AKT,1)
C
        ASUM(MU,NU)=ASUM(MU,NU)+REAL(AKT)
        PL0R(2,I4,ID,ITY)=PL0R(2,I4,ID,ITY)+REAL(AKT)
C
20      CONTINUE
50      CONTINUE
        AOP(1)=ASUM(1,2)
        AOP(2)=ASUM(3,1)
        AOP(3)=ASUM(2,3)
        AOP(4)=ASUM(2,1)
        AOP(5)=ASUM(1,3)
        AOP(6)=ASUM(3,2)
C
        PLER(3,I4,ID,ITY)=AOP(2)-AOP(3)+AOP(5)-AOP(6)
        PLER(4,I4,ID,ITY)=-2.0*AOP(1)+AOP(2)+AOP(3)
     &  -2.0*AOP(4)+AOP(5)+AOP(6)
        PLER(5,I4,ID,ITY)=-2.0*AOP(1)+AOP(2)+AOP(3)
     &  +2.0*AOP(4)-AOP(5)-AOP(6)
        PLER(6,I4,ID,ITY)=AOP(3)-AOP(2)+AOP(5)-AOP(6)
        PLTR(1,I4,ID,ITY)=AOP(1)-AOP(4)
        PLTR(2,I4,ID,ITY)=AOP(2)-AOP(5)
        PLTR(3,I4,ID,ITY)=AOP(3)-AOP(6)
C
        RETURN
        END
C*********************************************************************
C OK
C Here we take the wavefunctions from 'PLANAR'
C and calculate correlation functions e.g.
C ACOR2(bin,operator,distance in time,block level oper)
C**
C The ordering of the operators is as follows:
C
C label      0+               label         2+
C  1-2   planar 1-2            1-2     planar 1-2
C
C*********************************************************************
        SUBROUTINE PLNCORB(ITY)
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER(LX1=64,LX2=24,LX3=24,LX4=24,IBLOK=5,NUMBIN=25)
        PARAMETER(ITYPE=1)
        PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4,LMAX=LX4/2+1)
C
        COMMON/XSTORE/AVAC0(NUMBIN,2,IBLOK,ITYPE)
     &  ,ACOR0R(NUMBIN,2,2,LMAX,IBLOK,IBLOK,ITYPE)
     &  ,ACORER(NUMBIN,6,6,LMAX,IBLOK,IBLOK,ITYPE)
     &  ,ACORTR(NUMBIN,3,3,LMAX,IBLOK,IBLOK,ITYPE)
        COMMON/PLAN/PL0R(2,LX4,IBLOK,ITYPE),PLER(6,LX4,IBLOK,ITYPE)
     &  ,PLTR(3,LX4,IBLOK,ITYPE)
        COMMON/ITEM/ITR,IBIN,NTOT
        COMPLEX*16 PL0R,PLER,PLTR
C
        JBIN=(ITR-1)/IBIN+1
C
        DO 25 N4=1,LX4
        DO 25 ID=1,IBLOK
        DO 25 IOP=1,2
        AVAC0(JBIN,IOP,ID,ITY)=AVAC0(JBIN,IOP,ID,ITY)
     &  +REAL(PL0R(IOP,N4,ID,ITY))
25      CONTINUE
C
        DO 30 N4=1,LX4
        DO 30 NT=1,LMAX
        N4X=N4+NT-1
        IF(N4X.GT.LX4)N4X=N4X-LX4
C
        DO 32 ID2=1,IBLOK
        DO 32 ID1=1,IBLOK
C
        DO 34 IOP2=1,2
        DO 34 IOP1=1,2
        ACOR0R(JBIN,IOP1,IOP2,NT,ID1,ID2,ITY)=
     &  ACOR0R(JBIN,IOP1,IOP2,NT,ID1,ID2,ITY)
     &  +REAL(PL0R(IOP1,N4,ID1,ITY)*CONJG(PL0R(IOP2,N4X,ID2,ITY))
     &  +CONJG(PL0R(IOP2,N4,ID2,ITY))*PL0R(IOP1,N4X,ID1,ITY))*0.5
34      CONTINUE
C
        DO 35 IOP2=1,3
        DO 35 IOP1=1,3
        ACORTR(JBIN,IOP1,IOP2,NT,ID1,ID2,ITY)=
     &  ACORTR(JBIN,IOP1,IOP2,NT,ID1,ID2,ITY)
     &  +REAL(PLTR(IOP1,N4,ID1,ITY)*CONJG(PLTR(IOP2,N4X,ID2,ITY))
     &  +CONJG(PLTR(IOP2,N4,ID2,ITY))*PLTR(IOP1,N4X,ID1,ITY))
35      CONTINUE
C
        DO 36 IOP2=1,6
        DO 36 IOP1=1,6
        ACORER(JBIN,IOP1,IOP2,NT,ID1,ID2,ITY)=
     &  ACORER(JBIN,IOP1,IOP2,NT,ID1,ID2,ITY)
     &  +REAL(PLER(IOP1,N4,ID1,ITY)*CONJG(PLER(IOP2,N4X,ID2,ITY))
     &  +CONJG(PLER(IOP2,N4,ID2,ITY))*PLER(IOP1,N4X,ID1,ITY))
36      CONTINUE
C
32      CONTINUE
C
30      CONTINUE
C
c        write (*,*) 'AVAC0:', AVAC0
c        write (*,*) 'ACORTR:', ACORTR

        RETURN
        END
C*********************************************************************
C MOMENTUM PHASE FACTORS
C*********************************************************************
        SUBROUTINE MOMENTM
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER(LX1=64,LX2=24,LX3=24,LX4=24)
C
        COMMON/MOMS/AC1(3,0:50*LX1)
        COMPLEX*16 AC1
C
        PI=4.0*DATAN(1.0d0)
        DO 6 NI=0,LX1*5
        COEFF=(PI*2.0*NI)/LX1
        AC1(1,NI)=CMPLX(DCOS(COEFF),DSIN(COEFF))
6       CONTINUE
        DO 7 NI=0,LX2*5
        COEFF=(PI*2.0*NI)/LX2
        AC1(2,NI)=CMPLX(DCOS(COEFF),DSIN(COEFF))
7       CONTINUE
        DO 8 NI=0,LX3*5
        COEFF=(PI*2.0*NI)/LX3
        AC1(3,NI)=CMPLX(DCOS(COEFF),DSIN(COEFF))
8       CONTINUE
C
        RETURN
        END
C*********************************************************************
C OK
C We calculate all blocked thermal loops in each of the
C directions :    KU=1,2 =x,y
C   ALINE(I,J,KU,***) :
C   I = 1,2  means cosine,sine(p.x) factors included respectively
C   J =   1     2     3
C   p =   0     1     2
C Note:
C  when LXI is not a multiple of the blocked link length, we multiply
C  as many blocked links as possible and then fill the remaining gap
C  with links of the largest possible blocking. The latter are
C  fattened with 'staples' in which the 'horizontal' links are
C  at the original blocking level.
C and:
C  we include an IBLOK+1 level which consists of fattening with
C  staples the largest possible blocking level.
C*********************************************************************
        SUBROUTINE THERMLB(N4,IBLL,ITY)
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER(LX1=64,LX2=24,LX3=24,LX4=24,IBLOK=5,ITYPE=1)
        PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
        PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
        COMMON/LINES/ALINE1(3,7,LX4,IBLOK,ITYPE)
     &  ,ALINE2(3,7,2,LX4,IBLOK,ITYPE)
     &  ,ALINE3(3,7,3,LX4,IBLOK,ITYPE)
        COMMON/ARRAYB/UB11(NCOL2,LSIZEB,3,IBLOK)
        COMMON/NEXTB/IUPB(LSIZEB,3,IBLOK+1),IDNB(LSIZEB,3,IBLOK+1)
        COMMON/DUMMY/UC11(NCOL2,LSIZEB,3)
     &  ,IUP(LSIZEB,3),IDN(LSIZEB,3)
        COMMON/MOMS/AC1(3,0:50*LX1)
        DIMENSION LCNT(IBLOK),LB(IBLOK),IX(3),LS(3)
     &  ,ACT(3),AST(3),A11(NCOL2),B11(NCOL2),C11(NCOL2)
     &  ,D11(NCOL2),E11(NCOL2),UINT11(NCOL2),DUM11(NCOL2)
        COMPLEX UB11,A11,B11,C11,D11,E11,UC11,UINT11,DUM11
        COMPLEX*16 ALINE1,ALINE2,ALINE3,ACT
        COMPLEX*16 CSUM1(7),CSUM2(2,7),CSUM3(3,7),CFACT(7)
        COMPLEX*16 AKT1,AKT2,AKT3,AC1
C
        LS(1)=LX1
        LS(2)=LX2
        LS(3)=LX3
        LB(1)=1
        DO 4 IDD=2,IBLOK
        LB(IDD)=2*LB(IDD-1)
4       CONTINUE
C
        DO 10 KU=1,3
        KUADD=(KU-1)*LSIZEB
        JU=KU+1
        IF(JU.GT.3)JU=JU-3
        IU=KU+2
        IF(IU.GT.3)IU=IU-3
C
C        DO 20 ID=1,IBLOK
        ID=IBLL
C
        DO 11 IB=1,ID
        LCNT(IB)=0
11      CONTINUE
        LREST=LS(KU)
        DO 12 IG=1,ID
        IDG=ID-IG+1
        LCNT(IDG)=LREST/LB(IDG)
        LREST=LREST-LCNT(IDG)*LB(IDG)
12      CONTINUE
        DO 13 IB=1,ID
        IF(LCNT(IB).GE.1)IDS=IB
13      CONTINUE
        IDSM1=IDS-1
        IF(IDS.EQ.1)IDSM1=IDS
C
        DO 14 KK=1,3
        DO 14 NN=1,LSIZEB
        IUP(NN,KK)=IUPB(NN,KK,IDS)
        IDN(NN,KK)=IDNB(NN,KK,IDS)
14      CONTINUE
        DO 15 KK=1,3
        DO 15 NN=1,LSIZEB
        DO 15 IC=1,NCOL2
        UC11(IC,NN,KK)=UB11(IC,NN,KK,IDS)
15      CONTINUE
C
        LSKU=LB(IDS)
        IF(LSKU*(LS(KU)/LSKU).NE.LS(KU)) LSKU=LS(KU)
        DO 16 KK=1,7
        CSUM1(KK)=(0.0,0.0)
        CSUM2(1,KK)=(0.0,0.0)
        CSUM2(2,KK)=(0.0,0.0)
        CSUM3(1,KK)=(0.0,0.0)
        CSUM3(2,KK)=(0.0,0.0)
        CSUM3(3,KK)=(0.0,0.0)
16      CONTINUE
        NN=0
c        DO 30 NK=1,LSKU,2
        DO 30 NK=1,LSKU
        DO 30 NJ=1,LS(JU)
        DO 30 NI=1,LS(IU)
        NN=NN+1
        IX(IU)=NI
        IX(JU)=NJ
        IX(KU)=NK
        M2=IX(1)+LS(1)*(IX(2)-1)+LS(1)*LS(2)*(IX(3)-1)
        DO 17 IJ=1,NCOL2
        A11(IJ)=(0.0,0.0)
17      CONTINUE
        DO 18 N1=1,NCOL
        IJ=N1+NCOL*(N1-1)
        A11(IJ)=(1.0,0.0)
18      CONTINUE
C
        IF(IDS.EQ.ID)THEN
        DO 21 NC=1,LCNT(IDS)
        DO 22 IC=1,NCOL2
        B11(IC)=UC11(IC,M2,KU)
22      CONTINUE
        CALL VMX(1,A11,B11,C11,1)
        M3=IUP(M2,KU)
        M2=M3
        DO 24 IC=1,NCOL2
        A11(IC)=C11(IC)
24      CONTINUE
21      CONTINUE
        ENDIF
C
        DO 40 IG=1,IDS
        IDG=IDS-IG+1
        IF(IDG.EQ.ID)GOTO40
        DO 50 NC=1,LCNT(IDG)
        DO 25 IC=1,NCOL2
        UINT11(IC)=(0.0,0.0)
25      CONTINUE
C
        DO 32 IJ=1,3
        IF(IJ.EQ.KU)GOTO32
        DO 26 IC=1,NCOL2
        B11(IC)=UB11(IC,M2,IJ,IDSM1)
26      CONTINUE
        M3=IUPB(M2,IJ,IDSM1)
        DO 27 IC=1,NCOL2
        C11(IC)=UB11(IC,M3,KU,IDG)
27      CONTINUE
        CALL VMX(1,B11,C11,D11,1)
        M1=IUPB(M2,KU,IDG)
        DO 28 IC=1,NCOL2
        B11(IC)=UB11(IC,M1,IJ,IDSM1)
28      CONTINUE
        CALL HERM(1,B11,DUM11,1)
        CALL VMX(1,D11,B11,C11,1)
        DO 29 IC=1,NCOL2
        UINT11(IC)=UINT11(IC)+C11(IC)
29      CONTINUE
        M3=IDNB(M2,IJ,IDSM1)
        DO 31 IC=1,NCOL2
        B11(IC)=UB11(IC,M3,IJ,IDSM1)
31      CONTINUE
        DO 33 IC=1,NCOL2
        C11(IC)=UB11(IC,M3,KU,IDG)
33      CONTINUE
        CALL HERM(1,B11,DUM11,1)
        CALL VMX(1,B11,C11,D11,1)
        M1=IUPB(M3,KU,IDG)
        DO 34 IC=1,NCOL2
        B11(IC)=UB11(IC,M1,IJ,IDSM1)
34      CONTINUE
        CALL VMX(1,D11,B11,C11,1)
        DO 35 IC=1,NCOL2
        UINT11(IC)=UINT11(IC)+C11(IC)
35      CONTINUE
32      CONTINUE
        DO 36 IC=1,NCOL2
        B11(IC)=UB11(IC,M2,KU,IDG)
36      CONTINUE
        DO 37 IC=1,NCOL2
        UINT11(IC)=UINT11(IC)+B11(IC)
37      CONTINUE
        UMAG=1./6.
        DO 44 IC=1,NCOL2
        B11(IC)=UINT11(IC)*UMAG
44      CONTINUE
C
        CALL VMX(1,A11,B11,C11,1)
        M1=M2
        DO 46 IC=1,NCOL2
        A11(IC)=C11(IC)
46      CONTINUE
        M2=IUPB(M1,KU,IDG)
50      CONTINUE
40      CONTINUE
C
        CALL VMX(1,A11,A11,B11,1)
        CALL VMX(1,B11,A11,E11,1)
        AKT1=(0.0,0.0)
        AKT2=(0.0,0.0)
        AKT3=(0.0,0.0)
        DO 47 N1=1,NCOL
        IJ=N1+NCOL*(N1-1)
        AKT1=AKT1+A11(IJ)
        AKT2=AKT2+B11(IJ)
        AKT3=AKT3+E11(IJ)
47      CONTINUE
C
        CFACT(1)=(1.0,0.0)
        CFACT(2)=AC1(IU,NI)
        CFACT(3)=AC1(JU,NJ)
        CFACT(4)=AC1(JU,NJ)*AC1(IU,NI)
        CFACT(5)=AC1(JU,NJ)*CONJG(AC1(IU,NI))
        CFACT(6)=AC1(IU,NI*2)
        CFACT(7)=AC1(JU,NJ*2)
        DO 48 KK=1,7
        CSUM1(KK)=CSUM1(KK)+AKT1*CFACT(KK)
        CSUM2(1,KK)=CSUM2(1,KK)+AKT1*AKT1*CFACT(KK)
        CSUM2(2,KK)=CSUM2(2,KK)+AKT2*CFACT(KK)
        CSUM3(1,KK)=CSUM3(1,KK)+AKT1*AKT1*AKT1*CFACT(KK)
        CSUM3(2,KK)=CSUM3(2,KK)+AKT2*AKT1*CFACT(KK)
        CSUM3(3,KK)=CSUM3(3,KK)+AKT3*CFACT(KK)
48      CONTINUE
C
30      CONTINUE
C
c        ADIV=1./((LSKU+1)/2)
        ADIV=1./LSKU
        DO 55 KK=1,7
        ALINE1(KU,KK,N4,ID,ITY)=CSUM1(KK)*ADIV
        ALINE2(KU,KK,1,N4,ID,ITY)=CSUM2(1,KK)*ADIV        
        ALINE2(KU,KK,2,N4,ID,ITY)=CSUM2(2,KK)*ADIV
        ALINE3(KU,KK,1,N4,ID,ITY)=CSUM3(1,KK)*ADIV
        ALINE3(KU,KK,2,N4,ID,ITY)=CSUM3(2,KK)*ADIV
        ALINE3(KU,KK,3,N4,ID,ITY)=CSUM3(3,KK)*ADIV
55      CONTINUE
C
20      CONTINUE
10      CONTINUE
C
        RETURN
        END
C*********************************************************************
C OK
C Correlation in t dir of thermal loops in x,y directions
C Here
C ACORL(num bin, dir, dist, block level):  corrln  therm lines
C AVACF(num bin, dir, block level):        averages therm lines
C where
C   dir = 1-3 = x-z
C   dist= correlation distance + 1  in t direction
C*********************************************************************
        SUBROUTINE POTB(ITY)
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER(LX1=64,LX2=24,LX3=24,LX4=24,IBLOK=5,ITYPE=1)
        PARAMETER(LSIZE=LX1*LX2*LX3*LX4)
        PARAMETER(NUMBIN=25,LMAX=LX4/2+1)
        PARAMETER(IIMAX=1)
C
        COMMON/PBLOCK/ACORL1(NUMBIN,7,LMAX,IBLOK*IBLOK,ITYPE)
     &  ,AVACL1(NUMBIN,IBLOK,ITYPE)
     &  ,ACORL2(NUMBIN,7,2*2,LMAX,IBLOK*IBLOK,ITYPE)
     &  ,AVACL2(NUMBIN,2,IBLOK,ITYPE)
     &  ,ACORL3(NUMBIN,7,3*3,LMAX,IBLOK*IBLOK,ITYPE)
     &  ,AVACL3(NUMBIN,3,IBLOK,ITYPE)
        COMMON/LINES/ALINE1(3,7,LX4,IBLOK,ITYPE)
     &  ,ALINE2(3,7,2,LX4,IBLOK,ITYPE)
     &  ,ALINE3(3,7,3,LX4,IBLOK,ITYPE)
        COMMON/ITEM/ITR,IBIN,ITOT
        COMPLEX*16 ALINE1,ALINE2,ALINE3
        COMPLEX*16 ACORL1,ACORL2,ACORL3
        COMPLEX*16 ALLL
C
        write (*,*) 'IIMAX =', IIMAX
C
        JBIN=(ITR-1)/IBIN+1
C
        DO 1 ID=1,IBLOK
        DO 1 N4=1,LX4
        DO 1 I=1,IIMAX
        ALLL=ALINE1(I,1,N4,ID,ITY)
        AVACL1(JBIN,ID,ITY)=AVACL1(JBIN,ID,ITY)+REAL(ALLL)
        ALLL=ALINE2(I,1,1,N4,ID,ITY)
        AVACL2(JBIN,1,ID,ITY)=AVACL2(JBIN,1,ID,ITY)+REAL(ALLL)
        ALLL=ALINE2(I,1,2,N4,ID,ITY)
        AVACL2(JBIN,2,ID,ITY)=AVACL2(JBIN,2,ID,ITY)+REAL(ALLL)
        ALLL=ALINE3(I,1,1,N4,ID,ITY)
        AVACL3(JBIN,1,ID,ITY)=AVACL3(JBIN,1,ID,ITY)+REAL(ALLL)
        ALLL=ALINE3(I,1,2,N4,ID,ITY)
        AVACL3(JBIN,2,ID,ITY)=AVACL3(JBIN,2,ID,ITY)+REAL(ALLL)
        ALLL=ALINE3(I,1,3,N4,ID,ITY)
        AVACL3(JBIN,3,ID,ITY)=AVACL3(JBIN,3,ID,ITY)+REAL(ALLL)
1       CONTINUE
C
        ID12=0
        DO 2 ID2=1,IBLOK
        DO 2 ID1=1,IBLOK
        ID12=ID12+1
C
        DO 3 N4=1,LX4
        DO 3 NT=1,LMAX
        N4X=N4+NT-1
        IF(N4X.GT.LX4) N4X=N4X-LX4
        DO 3 NP=1,7
C
        DO 4 I=1,IIMAX
        ACORL1(JBIN,NP,NT,ID12,ITY)=
     &  ACORL1(JBIN,NP,NT,ID12,ITY)
     &  +(ALINE1(I,NP,N4,ID1,ITY)
     &  *CONJG(ALINE1(I,NP,N4X,ID2,ITY))
     &  +CONJG(ALINE1(I,NP,N4,ID2,ITY))
     &  *ALINE1(I,NP,N4X,ID1,ITY))*0.5
4       CONTINUE
C
3       CONTINUE
2       CONTINUE
C
        ID12=0
        DO 12 ID2=1,IBLOK
        DO 12 ID1=1,IBLOK
        ID12=ID12+1
        IE12=0
        DO 12 IE2=1,2
        DO 12 IE1=1,2
        IE12=IE12+1
C
        DO 13 N4=1,LX4
        DO 13 NT=1,LMAX
        N4X=N4+NT-1
        IF(N4X.GT.LX4) N4X=N4X-LX4
        DO 13 NP=1,7
C
        DO 14 I=1,IIMAX
        ACORL2(JBIN,NP,IE12,NT,ID12,ITY)=
     &  ACORL2(JBIN,NP,IE12,NT,ID12,ITY)
     &  +(ALINE2(I,NP,IE1,N4,ID1,ITY)
     &  *CONJG(ALINE2(I,NP,IE2,N4X,ID2,ITY))
     &  +CONJG(ALINE2(I,NP,IE2,N4,ID2,ITY))
     &  *ALINE2(I,NP,IE1,N4X,ID1,ITY))*0.5
14      CONTINUE
C
13      CONTINUE
12      CONTINUE
C
        ID12=0
        DO 22 ID2=1,IBLOK
        DO 22 ID1=1,IBLOK
        ID12=ID12+1
        IE12=0
        DO 22 IE2=1,3
        DO 22 IE1=1,3
        IE12=IE12+1
C
        DO 23 N4=1,LX4
        DO 23 NT=1,LMAX
        N4X=N4+NT-1
        IF(N4X.GT.LX4) N4X=N4X-LX4
        DO 23 NP=1,7
C
        DO 24 I=1,IIMAX
        ACORL3(JBIN,NP,IE12,NT,ID12,ITY)=
     &  ACORL3(JBIN,NP,IE12,NT,ID12,ITY)
     &  +(ALINE3(I,NP,IE1,N4,ID1,ITY)
     &  *CONJG(ALINE3(I,NP,IE2,N4X,ID2,ITY))
     &  +CONJG(ALINE3(I,NP,IE2,N4,ID2,ITY))
     &  *ALINE3(I,NP,IE1,N4X,ID1,ITY))*0.5
24      CONTINUE
C
23      CONTINUE
22      CONTINUE
C
        RETURN
        END
C*********************************************************************
C OK
C  Here we set up blocked link pointers...
C****************************
C  WARNING: only works if
C     max length blocked link .le. 2*min[lx1,lx2,lx3]
C*********************************************************************
      SUBROUTINE SETUPB
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=64,LX2=24,LX3=24,LX4=24,IBLOK=5)
      PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
C
      COMMON/NEXTB/IUPB(LSIZEB,3,IBLOK+1),IDNB(LSIZEB,3,IBLOK+1)
      DIMENSION LB(IBLOK+1),LB1(IBLOK+1),LB12(IBLOK+1)
C
      LX12=LX1*LX2
      LX123=LX12*LX3
      LB(1)=1
      LB1(1)=LX1
      LB12(1)=LX12
      DO 1 IBL=2,IBLOK+1
      LB(IBL)=2*LB(IBL-1)
      LB1(IBL)=2*LB1(IBL-1)
      LB12(IBL)=2*LB12(IBL-1)
1     CONTINUE
C
      NN=0
      DO 2 L3=1,LX3
      DO 2 L2=1,LX2
      DO 2 L1=1,LX1
      NN=NN+1
      DO 3 ID=1,IBLOK+1
C
      NU1=NN+LB(ID)
      IF((L1+LB(ID)).GT.(LX1*1)) NU1=NU1-LX1
      IF((L1+LB(ID)).GT.(LX1*2)) NU1=NU1-LX1
      IF((L1+LB(ID)).GT.(LX1*3)) NU1=NU1-LX1
      IF((L1+LB(ID)).GT.(LX1*4)) NU1=NU1-LX1
      IUPB(NN,1,ID)=NU1
      ND1=NN-LB(ID)
      IF((L1-LB(ID)).LT.(1-LX1*0)) ND1=ND1+LX1
      IF((L1-LB(ID)).LT.(1-LX1*1)) ND1=ND1+LX1
      IF((L1-LB(ID)).LT.(1-LX1*2)) ND1=ND1+LX1
      IF((L1-LB(ID)).LT.(1-LX1*3)) ND1=ND1+LX1
      IDNB(NN,1,ID)=ND1
C
      NU2=NN+LB1(ID)
      IF((L2+LB(ID)).GT.(LX2*1)) NU2=NU2-LX12
      IF((L2+LB(ID)).GT.(LX2*2)) NU2=NU2-LX12
      IF((L2+LB(ID)).GT.(LX2*3)) NU2=NU2-LX12
      IF((L2+LB(ID)).GT.(LX2*4)) NU2=NU2-LX12
      IUPB(NN,2,ID)=NU2
      ND2=NN-LB1(ID)
      IF((L2-LB(ID)).LT.(1-LX2*0)) ND2=ND2+LX12
      IF((L2-LB(ID)).LT.(1-LX2*1)) ND2=ND2+LX12
      IF((L2-LB(ID)).LT.(1-LX2*2)) ND2=ND2+LX12
      IF((L2-LB(ID)).LT.(1-LX2*3)) ND2=ND2+LX12
      IDNB(NN,2,ID)=ND2
C
      NU3=NN+LB12(ID)
      IF((L3+LB(ID)).GT.(LX3*1)) NU3=NU3-LX123
      IF((L3+LB(ID)).GT.(LX3*2)) NU3=NU3-LX123
      IF((L3+LB(ID)).GT.(LX3*3)) NU3=NU3-LX123
      IF((L3+LB(ID)).GT.(LX3*4)) NU3=NU3-LX123
      IUPB(NN,3,ID)=NU3
      ND3=NN-LB12(ID)
      IF((L3-LB(ID)).LT.(1-LX3*0)) ND3=ND3+LX123
      IF((L3-LB(ID)).LT.(1-LX3*1)) ND3=ND3+LX123
      IF((L3-LB(ID)).LT.(1-LX3*2)) ND3=ND3+LX123
      IF((L3-LB(ID)).LT.(1-LX3*3)) ND3=ND3+LX123
      IDNB(NN,3,ID)=ND3
C
3     CONTINUE
2     CONTINUE
C
      RETURN
      END
C*********************************************************************
C OK
C*********************************************************************
      SUBROUTINE POLYLOOP()
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=64,LX2=24,LX3=24,LX4=24)
      PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
      PARAMETER(NITER=200,NUMBIN=25)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
      COMMON/ARRAYS/U11(NCOL2,LSIZE,4)
      COMMON/NEXT/IUP(LSIZE,4),IDN(LSIZE,4)
      DIMENSION DUM11(NCOL2)
     &,A11(NCOL2),B11(NCOL2),C11(NCOL2),D11(NCOL2)
      COMPLEX U11,A11,B11,C11,D11,DUM11,ACT
C
      DIMENSION AVAC(NUMBIN),AVACSQ(NUMBIN),AVACS(NUMBIN),AVACT(NUMBIN)
      DIMENSION VAL(NUMBIN),AV(NUMBIN)
      dimension icoordvect(4)
      complex cpol1
      complex*16 cpol(4)
C
      icoordvect(1) = LX1-1
      icoordvect(2) = LX2-1
      icoordvect(3) = LX3-1
      icoordvect(4) = LX4-1
C     
      cpol(:) = dcmplx(0.0d0,0.0d0) 
      dnorm = 1.0/dfloat(lsizeb*ncol)
c
      DO ipoint=1, lsizeb
         do idir=1, 4
            a11(:) = u11(:,ipoint,idir)
            ipoint1 = ipoint
            do icoord=2,icoordvect(idir) 
               ipoint2 = iup(ipoint1,idir) 
               b11 = u11(:,ipoint2,idir)
               call vmx(1,a11,b11,c11,1)
               a11(:) = c11(:)
               ipoint2 = ipoint1
            enddo
            ipoint2 = iup(ipoint1,idir) 
            b11 = u11(:,ipoint2,idir)
            call trvmx(1,a11,b11,cpol1,1)
            cpol(idir) = cpol(idir) + cpol1
c            write (*,*) cpol1, cpol(idir)
         enddo
      enddo
C
      cpol(:) = cpol(:)*dnorm
C
 888  format('Poly(' ,i1, ')  =  (',f16.12,','f16.12')')
      do idir=1,4
         write (92,888) idir,dreal(cpol(idir)),dimag(cpol(idir))
      enddo
C
      RETURN
      END
      

C*********************************************************************
C OK
C*********************************************************************
      SUBROUTINE ACTION(IPR,ITER,TOTACT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=64,LX2=24,LX3=24,LX4=24)
      PARAMETER(NITER=200,NUMBIN=25)
      PARAMETER(LSIZE=LX1*LX2*LX3*LX4)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
      COMMON/ASTORE/ACTN(NUMBIN,6)
      COMMON/ARRAYS/U11(NCOL2,LSIZE,4)
      COMMON/NEXT/IUP(LSIZE,4),IDN(LSIZE,4)
      DIMENSION DUM11(NCOL2)
     &,A11(NCOL2),B11(NCOL2),C11(NCOL2),D11(NCOL2)
      COMPLEX U11,A11,B11,C11,D11,DUM11,ACT
C
      DIMENSION AVAC(NUMBIN),AVACSQ(NUMBIN),AVACS(NUMBIN),AVACT(NUMBIN)
      DIMENSION VAL(NUMBIN),AV(NUMBIN)
C
      IBIN=NITER/NUMBIN
C
      IF(IPR.EQ.0)THEN
      JBIN=(ITER-1)/IBIN+1
C
      IF(ITER.EQ.1)THEN
      DO 3 NB=1,NUMBIN
      AVAC(NB)=0.0d0
      AVACSQ(NB)=0.0d0
      AVACS(NB)=0.0d0
      AVACT(NB)=0.0d0
3     CONTINUE
      ENDIF
C
      VACTS=0.0
      VACTT=0.0
C
      DO 10 NN=1,LSIZE
      DO 10 MU=1,3
      M1=NN
C
      DO 6  IJ=1,NCOL2
      A11(IJ)=U11(IJ,M1,MU)
6     CONTINUE
      M2=IUP(M1,MU)
C
      DO 11 NU=MU+1,4
      IPLAQ=6-NU-MU+5*(NU/4)
C
      DO 7  IJ=1,NCOL2
      B11(IJ)=U11(IJ,M2,NU)
7     CONTINUE
      CALL VMX(1,A11,B11,C11,1)
      M3=IUP(M1,NU)
      DO 8  IJ=1,NCOL2
      D11(IJ)=U11(IJ,M3,MU)
8     CONTINUE
      CALL HERM(1,D11,DUM11,1)
      CALL VMX(1,C11,D11,B11,1)
      DO 9  IJ=1,NCOL2
      D11(IJ)=U11(IJ,M1,NU)
9     CONTINUE
      CALL HERM(1,D11,DUM11,1)
      CALL TRVMX(1,B11,D11,ACT,1)
      ANN=1.0/NCOL
      ACT=ANN*REAL(ACT)
      ACTN(JBIN,IPLAQ)=ACTN(JBIN,IPLAQ)+ACT
      IF(NU.NE.4)THEN
      VACTS=VACTS+ACT
      ELSE
      VACTT=VACTT+ACT
      ENDIF
C
11    CONTINUE
10    CONTINUE
C
      VACTS=VACTS/(3.0*LSIZE)
      VACTT=VACTT/(3.0*LSIZE)
      TOTACT = 0.5*(VACTS+VACTT)
      write (91,*) VACTS, VACTT, TOTACT 
C
      AVACS(JBIN)=AVACS(JBIN)+VACTS
      AVACT(JBIN)=AVACT(JBIN)+VACTT
C
      IF(ITER.EQ.NITER)THEN
      DO 17 IP=1,6
      DO 17 JB=1,NUMBIN
      ACTN(JB,IP)=ACTN(JB,IP)/(IBIN*LSIZE)
17    CONTINUE
      ENDIF
C
      ENDIF
C
      IF(IPR.EQ.1)THEN
C
      WRITE(6,80)
80    FORMAT(' *****************************************************')
      WRITE(6,100)
100   FORMAT('                 ACTION            ')
      WRITE(6,80)
      WRITE(6,81)
81    FORMAT(' *  ')
      DO 20 NB=1,NUMBIN
      AV(NB)=AVACS(NB)/IBIN
20    CONTINUE
      CALL JACKM(NUMBIN,AV,AVRS,ERS)
      DO 22 NB=1,NUMBIN
      AV(NB)=AVACT(NB)/IBIN
22    CONTINUE
      CALL JACKM(NUMBIN,AV,AVRT,ERT)
      WRITE(6,110) AVRS,ERS
110   FORMAT('   AVER,ERR   SPACE ACTION =',2F12.8)
      WRITE(6,112) AVRT,ERT
112   FORMAT('   AVER,ERR   TIME  ACTION =',2F12.8)
      WRITE(6,81)
C
      ENDIF
C
      RETURN
      END
C**************************************************************
C OK
C  HERE WE CONSTRUCT A FROZEN GAUGE CONFIGURATION :
C**************************************************************
        SUBROUTINE STARTF
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER(LX1=64,LX2=24,LX3=24,LX4=24)
        PARAMETER(LSIZE=LX1*LX2*LX3*LX4)
        PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
        COMMON/ARRAYS/U11(NCOL2,LSIZE,4)
        COMPLEX U11,A11(NCOL2)
C
        DO 4 IJ=1,NCOL2
        A11(IJ)=(0.0,0.0)
4       CONTINUE
        DO 6 N1=1,NCOL
        NC=N1+NCOL*(N1-1)
        A11(NC)=(1.0,0.0)
6       CONTINUE
C
        DO 1 MU=1,4
        DO 1 NN=1,LSIZE
        DO 2 IJ=1,NCOL2
        U11(IJ,NN,MU)=A11(IJ)
2       CONTINUE
1       CONTINUE
C
        RETURN
        END
C**********************************************************
C  OK
C  THIS ROUTINE REIMPOSES THE UNITARITY CONSTRAINTS ON
C  OUR SU(3) MATRICES
C**********************************************************
        SUBROUTINE RENORM
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER(LX1=64,LX2=24,LX3=24,LX4=24)
        PARAMETER(LSIZE=LX1*LX2*LX3*LX4)
        PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
        COMMON/ARRAYS/U11(NCOL,NCOL,LSIZE,4)
        COMPLEX ADUM(NCOL,NCOL),U11,CSUM
C
        DO 100 MU=1,4
        DO 100 NN=1,LSIZE
C
        DO 1 J=1,NCOL
        DO 1 I=1,NCOL
        ADUM(I,J)=U11(I,J,NN,MU)
1       CONTINUE
C
        DO 20 N2=1,NCOL
C
        DO 10 N3=1,N2-1
C
        CSUM=(0.0,0.0)
        DO 5 N1=1,NCOL
        CSUM=CSUM+ADUM(N1,N2)*CONJG(ADUM(N1,N3))
5       CONTINUE
        DO 6 N1=1,NCOL
        ADUM(N1,N2)=ADUM(N1,N2)-CSUM*ADUM(N1,N3)
6       CONTINUE
C
10      CONTINUE
C
        SUM=0.0
        DO 7 N1=1,NCOL
        SUM=SUM+ADUM(N1,N2)*CONJG(ADUM(N1,N2))
7       CONTINUE
        ANORM=1.0/SQRT(SUM)
        DO 8 N1=1,NCOL
        ADUM(N1,N2)=ADUM(N1,N2)*ANORM
8       CONTINUE
C
20      CONTINUE
C
        DO 2 J=1,NCOL
        DO 2 I=1,NCOL
        U11(I,J,NN,MU)=ADUM(I,J)
2       CONTINUE
C
100     CONTINUE
        RETURN
        END
C***********************************************************
C OK
C  HERE WE SET UP  LINK POINTERS FOR FULL LATTICE
C***********************************************************
      SUBROUTINE SETUP
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=64,LX2=24,LX3=24,LX4=24)
      PARAMETER(LSIZE=LX1*LX2*LX3*LX4)
C
      COMMON/NEXT/IUP(LSIZE,4),IDN(LSIZE,4)
C
      LX12=LX1*LX2
      LX123=LX12*LX3
      LX1234=LX123*LX4
C
      NN=0
      DO 2 L4=1,LX4
      DO 2 L3=1,LX3
      DO 2 L2=1,LX2
      DO 2 L1=1,LX1
      NN=NN+1
C
      NU1=NN+1
      IF((L1+1).GT.LX1) NU1=NU1-LX1
      IUP(NN,1)=NU1
      ND1=NN-1
      IF((L1-1).LT.1) ND1=ND1+LX1
      IDN(NN,1)=ND1
C
      NU2=NN+LX1
      IF((L2+1).GT.LX2) NU2=NU2-LX12
      IUP(NN,2)=NU2
      ND2=NN-LX1
      IF((L2-1).LT.1) ND2=ND2+LX12
      IDN(NN,2)=ND2
C
      NU3=NN+LX12
      IF((L3+1).GT.LX3) NU3=NU3-LX123
      IUP(NN,3)=NU3
      ND3=NN-LX12
      IF((L3-1).LT.1) ND3=ND3+LX123
      IDN(NN,3)=ND3
C
      NU4=NN+LX123
      IF((L4+1).GT.LX4) NU4=NU4-LX1234
      IUP(NN,4)=NU4
      ND4=NN-LX123
      IF((L4-1).LT.1) ND4=ND4+LX1234
      IDN(NN,4)=ND4
C
2     CONTINUE
C
      RETURN
      END
C***********************************************************
C OK
C  HERE WE UPDATE THE LINKS
C***********************************************************
      SUBROUTINE UPDATG(IBIAS,ACTL,BETAG,SRAT,TRAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=64,LX2=24,LX3=24,LX4=24)
      PARAMETER(LSIZE=LX1*LX2*LX3*LX4)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
      COMMON/ARRAYS/U11(NCOL2,LSIZE,4)
      COMMON/NEXT/IUP(LSIZE,4),IDN(LSIZE,4)
      COMMON/HB/BBETAG,IIBIAS
      DIMENSION UINT11(NCOL2),DUM11(NCOL2)
     &,A11(NCOL2),B11(NCOL2),C11(NCOL2)
C
      COMPLEX U11,A11,B11,C11,UINT11,DUM11
C
      BBETAG=BETAG
      IIBIAS=IBIAS
      SUM=0.
      DO 1 NN=1,LSIZE
      M1=NN
C
      DO 10 MU=1,4
C
      DO 11 IJ=1,NCOL2
      UINT11(IJ)=(0.,0.)
11    CONTINUE
      M5=IUP(M1,MU)
C
      DO 20 NU=1,4
      IF(MU.EQ.NU) GO TO 20
      IF(MU.EQ.4.OR.NU.EQ.4)THEN
      CRAT=TRAT
      ELSE
      CRAT=SRAT
      ENDIF
C
      DO 15 IJ=1,NCOL2
      A11(IJ)=U11(IJ,M1,NU)
15    CONTINUE
      M3=IUP(M1,NU)
      DO 17 IJ=1,NCOL2
      B11(IJ)=U11(IJ,M3,MU)
17    CONTINUE
      CALL VMX(1,A11,B11,C11,1)
      DO 19 IJ=1,NCOL2
      A11(IJ)=U11(IJ,M5,NU)
19    CONTINUE
      CALL HERM(1,A11,DUM11,1)
      CALL VMX(1,C11,A11,B11,1)
      DO 23 IJ=1,NCOL2
      UINT11(IJ)=UINT11(IJ)+B11(IJ)*CRAT
23    CONTINUE
      M3=IDN(M1,NU)
      DO 25 IJ=1,NCOL2
      B11(IJ)=U11(IJ,M3,MU)
25    CONTINUE
      DO 27 IJ=1,NCOL2
      A11(IJ)=U11(IJ,M3,NU)
27    CONTINUE
      CALL HERM(1,A11,DUM11,1)
      CALL VMX(1,A11,B11,C11,1)
      M3=IDN(M5,NU)
      DO 31 IJ=1,NCOL2
      A11(IJ)=U11(IJ,M3,NU)
31    CONTINUE
      CALL VMX(1,C11,A11,B11,1)
      DO 32 IJ=1,NCOL2
      UINT11(IJ)=UINT11(IJ)+B11(IJ)*CRAT
32    CONTINUE
20    CONTINUE
C
      CALL HERM(1,UINT11,DUM11,1)
C
      DO 36 IJ=1,NCOL2
      C11(IJ)=U11(IJ,M1,MU)
36    CONTINUE
      CALL VMX(1,UINT11,C11,B11,1)
C
      IF(NCOL.EQ.2)THEN
      CALL SUBGRH(1,4,2,3,B11,C11)
      ENDIF
      IF(NCOL.EQ.3)THEN
      CALL SUBGRH(1,5,2,4,B11,C11)
      CALL SUBGRH(5,9,6,8,B11,C11)
      CALL SUBGRH(1,9,3,7,B11,C11)
      ENDIF
      IF(NCOL.EQ.4)THEN
      CALL SUBGRH(1,6,2,5,B11,C11)
      CALL SUBGRH(6,11,7,10,B11,C11)
      CALL SUBGRH(11,16,12,15,B11,C11)
      CALL SUBGRH(1,16,4,13,B11,C11)
      CALL SUBGRH(6,16,8,14,B11,C11)
      CALL SUBGRH(1,11,3,9,B11,C11)
      ENDIF
      IF(NCOL.EQ.5)THEN
      CALL SUBGRH(1,7,2,6,B11,C11)
      CALL SUBGRH(7,13,8,12,B11,C11)
      CALL SUBGRH(13,19,14,18,B11,C11)
      CALL SUBGRH(19,25,20,24,B11,C11)
      CALL SUBGRH(7,19,9,17,B11,C11)
      CALL SUBGRH(13,25,15,23,B11,C11)
      CALL SUBGRH(1,13,3,11,B11,C11)
      CALL SUBGRH(1,19,4,16,B11,C11)
      CALL SUBGRH(7,25,10,22,B11,C11)
      CALL SUBGRH(1,25,5,21,B11,C11)
      ENDIF
      IF(NCOL.EQ.6)THEN
      CALL SUBGRH(1,8,2,7,B11,C11)
      CALL SUBGRH(8,15,9,14,B11,C11)
      CALL SUBGRH(15,22,16,21,B11,C11)
      CALL SUBGRH(22,29,23,28,B11,C11)
      CALL SUBGRH(29,36,30,35,B11,C11)
      CALL SUBGRH(1,15,3,13,B11,C11)
      CALL SUBGRH(8,22,10,20,B11,C11)
      CALL SUBGRH(15,29,17,27,B11,C11)
      CALL SUBGRH(22,36,24,34,B11,C11)
      CALL SUBGRH(1,22,4,19,B11,C11)
      CALL SUBGRH(8,29,11,26,B11,C11)
      CALL SUBGRH(15,36,18,33,B11,C11)
      CALL SUBGRH(1,29,5,25,B11,C11)
      CALL SUBGRH(8,36,12,32,B11,C11)
      CALL SUBGRH(1,36,6,31,B11,C11)
      ENDIF
C
      DO 40 IJ=1,NCOL2
      U11(IJ,M1,MU)=C11(IJ)
40    CONTINUE
      DO 42 I=1,NCOL
      NC=I+NCOL*(I-1)
      SUM=SUM+REAL(B11(NC))
42    CONTINUE
C
10    CONTINUE
1     CONTINUE
C
      APQ=1.-SUM/(24.*NCOL*LSIZE)
      ACTL=APQ
C
      RETURN
      END
C*****************************************************
C*****************************************************
      SUBROUTINE SUBGRH(I1,I2,I3,I4,B11,C11)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
      COMMON/HB/BETAG,IBIAS,UMAG,A0,A1,A2,A3
C
      COMPLEX B11(NCOL2),C11(NCOL2)
      COMPLEX F11,F12,A11(NCOL2),S11(NCOL2)
      COMPLEX E11,E12,G11,G12
C
      F11=(B11(I1)+CONJG(B11(I2)))*0.5
      F12=(B11(I3)-CONJG(B11(I4)))*0.5
      UMAG=SQRT(F11*CONJG(F11)+F12*CONJG(F12))
      UMAG=1./UMAG
      F11=F11*UMAG
      F12=F12*UMAG
      IF(IBIAS.EQ.1)THEN
      CALL HEATB
      E11=CMPLX(A0,A1)
      E12=CMPLX(A2,A3)
      G11=CONJG(F11)*E11+CONJG(F12)*E12
      G12=-F12*E11+F11*E12
      A11(I1)=G11
      A11(I4)=-CONJG(G12)
      A11(I3)=G12
      A11(I2)=CONJG(G11)
      ENDIF
      IF(IBIAS.EQ.2)THEN
      G11= F11*F11-CONJG(F12)*F12
      G12= F12*F11+CONJG(F11)*F12
      A11(I1)=CONJG(G11)
      A11(I4)=CONJG(G12)
      A11(I3)=-G12
      A11(I2)=G11
      ENDIF
      IF(IBIAS.EQ.3)THEN
      A11(I1)=CONJG(F11)
      A11(I4)=CONJG(F12)
      A11(I3)=-F12
      A11(I2)=F11
      ENDIF
      DO 3 IJ2=1,NCOL
      DO 3 IJ1=1,NCOL
      IJ=IJ1+NCOL*(IJ2-1)
      IF(IJ.EQ.I1.OR.IJ.EQ.I2.OR.IJ.EQ.I3.OR.IJ.EQ.I4)GOTO3
      IF(IJ1.EQ.IJ2)THEN
      A11(IJ)=(1.0,0.0)
      ELSE
      A11(IJ)=(0.0,0.0)
      ENDIF
3     CONTINUE
      CALL VMX(1,C11,A11,S11,1)
      CALL VMX(1,B11,A11,C11,1)
      DO 5  IJ=1,NCOL2
      B11(IJ)=C11(IJ)
5     CONTINUE
      DO 6  IJ=1,NCOL2
      C11(IJ)=S11(IJ)
6     CONTINUE
C
      RETURN
      END
C*****************************************************
C SU2 SUBGROUP HEATBATH
C*****************************************************
      SUBROUTINE HEATB
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LX1=64,LX2=24,LX3=24,LX4=24)
      PARAMETER(LSIZE=LX1*LX2*LX3*LX4)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
      COMMON/HB/BETAG,IBIAS,UMAG,A0,A1,A2,A3
      REAL*8 RNDNUM
C
      PI=4.0*DATAN(1.0d0)
      BTG=BETAG*2.0d0/NCOL
      TEMP=1.d0/BTG
      DD=RNDNUM()
C
      UMAG=UMAG*TEMP
      A1=-DLOG(RNDNUM())*UMAG
      A2=-DLOG(RNDNUM())*UMAG
      A3=DCOS(2.0d0*PI*RNDNUM())
      A3=A3*A3
      A1=A1*A3
      A2=A2+A1
      A0=1.0d0-A2
      A3=RNDNUM()
      A3=A3*A3-1+A2*0.5d0
      IF(A3.GT.0.0d0)THEN
45    X11=-DLOG(RNDNUM())*UMAG
      X22=-DLOG(RNDNUM())*UMAG
      CCC=DCOS(2.0d0*PI*RNDNUM())
      CCC=CCC*CCC
      AA=X11*CCC
      DBB=X22+AA
      XX=RNDNUM()
      XX=XX*XX-1.0d0+DBB*0.5d0
      IF(XX.GT.0.0d0)GOTO45
      A0=1.0d0-DBB
      ENDIF
C
54    CONTINUE
      X1=2.0*RNDNUM()-1.0d0
      X2=2.0*RNDNUM()-1.0d0
      X3=2.0*RNDNUM()-1.0d0
      CC=X1*X1+X2*X2+X3*X3-1.0d0
      IF(CC.LE.0.0)THEN
      A1=X1
      A2=X2
      A3=X3
      ELSE
      GOTO54
      ENDIF
      RAD=(1.d0-A0*A0)
      RAD=(A1*A1+A2*A2+A3*A3)/RAD
      RAD=1.d0/DSQRT(RAD)
      A1=A1*RAD
      A2=A2*RAD
      A3=A3*RAD
C
      RETURN
      END
C***************************************************************
c   here we measure ffdual and action on whole lattice ;
c    - unsymmetrised ffdual for speed
c   we print sites where they are large .
c   vectorised ... except for pieces that are done only
c                  when ipr=2 i.e. usually only last cool.
C**************************************************************
        SUBROUTINE FFDACT
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER(LX1=64,LX2=24,LX3=24,LX4=24)
        PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
        PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
        COMMON/NEXT/IUP(LSIZE,4),IDN(LSIZE,4)
        COMMON/ARRAYS/U11(NCOL2,LSIZE,4)
        DIMENSION PL11(NCOL2,6),DUM11(NCOL2)
     &  ,A11(NCOL2),B11(NCOL2),C11(NCOL2),D11(NCOL2)
        COMPLEX U11,PL11,A11,B11,C11,D11,DUM11,CFF
        COMMON/DUMMY/FFD(LSIZE),ACT(LSIZE)
        DIMENSION IXM(4),IXN(4),IXP(4)
C
        DO 1 NN=1,LSIZE
        M1=NN
C
        ACT(NN)=0.0
        FFD(NN)=0.0
C
      DO 10 MU=1,3
      DO 10 NU=MU+1,4
      IPLAQ=MU+NU-2+NU/4
C
      DO 20 IJ=1,NCOL2
      A11(IJ)=U11(IJ,M1,MU)
20    CONTINUE
      M5=IUP(M1,MU)
      DO 21 IJ=1,NCOL2
      B11(IJ)=U11(IJ,M5,NU)
21    CONTINUE
      CALL VMX(1,A11,B11,C11,1)
      M5=IUP(M1,NU)
      DO 22 IJ=1,NCOL2
      A11(IJ)=U11(IJ,M5,MU)
22    CONTINUE
      CALL HERM(1,A11,DUM11,1)
      CALL VMX(1,C11,A11,B11,1)
      DO 23 IJ=1,NCOL2
      A11(IJ)=U11(IJ,M1,NU)
23    CONTINUE
      CALL HERM(1,A11,DUM11,1)
      CALL VMX(1,B11,A11,C11,1)
      DO 25 IJ=1,NCOL2
      PL11(IJ,IPLAQ)=C11(IJ)
25    CONTINUE
C
10    CONTINUE
C
      DO 30 IJ=1,NCOL2
      C11(IJ)=PL11(IJ,4)
      D11(IJ)=PL11(IJ,3)
30    CONTINUE
      DO 31 IJ=1,NCOL2
      A11(IJ)=C11(IJ)
      B11(IJ)=D11(IJ)
31    CONTINUE
      CALL HERM(1,C11,DUM11,1)
      CALL HERM(1,D11,DUM11,1)
      DO 32 IJ=1,NCOL2
      A11(IJ)=A11(IJ)-C11(IJ)
      B11(IJ)=B11(IJ)-D11(IJ)
32    CONTINUE
      CALL TRVMX(1,A11,B11,CFF,1)
      FFDD=REAL(CFF)
      DO 35 IJ=1,NCOL2
      C11(IJ)=PL11(IJ,5)
      D11(IJ)=PL11(IJ,2)
35    CONTINUE
      DO 36 IJ=1,NCOL2
      A11(IJ)=C11(IJ)
      B11(IJ)=D11(IJ)
36    CONTINUE
      CALL HERM(1,C11,DUM11,1)
      CALL HERM(1,D11,DUM11,1)
      DO 37 IJ=1,NCOL2
      A11(IJ)=A11(IJ)-C11(IJ)
      B11(IJ)=B11(IJ)-D11(IJ)
37    CONTINUE
      CALL TRVMX(1,A11,B11,CFF,1)
      FFDD=FFDD-REAL(CFF)
      DO 40 IJ=1,NCOL2
      C11(IJ)=PL11(IJ,6)
      D11(IJ)=PL11(IJ,1)
40    CONTINUE
      DO 41 IJ=1,NCOL2
      A11(IJ)=C11(IJ)
      B11(IJ)=D11(IJ)
41    CONTINUE
      CALL HERM(1,C11,DUM11,1)
      CALL HERM(1,D11,DUM11,1)
      DO 42 IJ=1,NCOL2
      A11(IJ)=A11(IJ)-C11(IJ)
      B11(IJ)=B11(IJ)-D11(IJ)
42    CONTINUE
      CALL TRVMX(1,A11,B11,CFF,1)
      FFDD=FFDD+REAL(CFF)
C
      FFD(NN)=FFDD*2.0
C
      ACTX=0.0
      DIV3=1.0/NCOL
      DO 46 I=1,6
      DO 46 N1=1,NCOL
      IJ=N1+NCOL*(N1-1)
      ACTX=ACTX+REAL(PL11(IJ,I))
46    CONTINUE
      ACT(NN)=6.0-ACTX*DIV3
      ACTX=0.0
C
1     CONTINUE
C
      ACTTOT=0.0
      FFDTOT=0.0
      FFDMOD=0.0
C
      DO 50 N=1,LSIZE
      FFDC=FFD(N)
      ACTT=ACT(N)
      ACTTOT=ACTTOT+ACTT
      FFDTOT=FFDTOT+FFDC
      FFDMOD=FFDMOD+ABS(FFDC)
50    CONTINUE
C
      ACTMAX=0.0
      FFDMAX=0.0
      FFDMIN=0.0
C
      DO 52 N=1,LSIZE
      FFDC=FFD(N)
      ACTT=ACT(N)
      IF(FFDC.GT.FFDMAX)THEN
      FFDMAX=FFDC
      IMAX=N
      ENDIF
      IF(FFDC.LT.FFDMIN)THEN
      FFDMIN=FFDC
      IMIN=N
      ENDIF
      IF(ACTT.GT.ACTMAX)THEN
      ACTMAX=ACTT
      JMAX=N
      ENDIF
52    CONTINUE
C
        LX123=LX1*LX2*LX3
        JMAX=JMAX-1
        IXP(4)=JMAX/LX123+1
        LX12=LX1*LX2
        JMAX=JMAX-(IXP(4)-1)*LX123
        IXP(3)=JMAX/LX12+1
        JMAX=JMAX-(IXP(3)-1)*LX12
        IXP(2)=JMAX/LX1+1
        JMAX=JMAX-(IXP(2)-1)*LX1
        IXP(1)=JMAX+1
        IMAX=IMAX-1
        IXM(4)=IMAX/LX123+1
        IMAX=IMAX-(IXM(4)-1)*LX123
        IXM(3)=IMAX/LX12+1
        IMAX=IMAX-(IXM(3)-1)*LX12
        IXM(2)=IMAX/LX1+1
        IMAX=IMAX-(IXM(2)-1)*LX1
        IXM(1)=IMAX+1
        IMIN=IMIN-1
        IXN(4)=IMIN/LX123+1
        IMIN=IMIN-(IXN(4)-1)*LX123
        IXN(3)=IMIN/LX12+1
        IMIN=IMIN-(IXN(3)-1)*LX12
        IXN(2)=IMIN/LX1+1
        IMIN=IMIN-(IXN(2)-1)*LX1
        IXN(1)=IMIN+1
C
        WRITE(6,85)ACTTOT,ACTMAX,FFDTOT,FFDMOD,FFDMAX,FFDMIN
85      FORMAT(' TOT/MAX(MOD,MIN) S,Q= ',6F9.2)
        WRITE(6,86)IXP,IXM,IXN
86      FORMAT('SITE SMAX QMAX,MIN = ',4I3,I10,3I3,I10,3I3)
C
        RETURN
        END
C***************************************************************
C   HERE WE SEARCH FOR POINTS WHERE FFD(N,LT) IS A MAXIMUM
C   OR MINIMUM
C**************************************************************
        SUBROUTINE SEARCH
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER(LX1=64,LX2=24,LX3=24,LX4=24)
        PARAMETER(LSIZEB=LX1*LX2*LX3,LSIZE=LSIZEB*LX4)
        PARAMETER(QCUT=0.02,RRCUT=4.0)
        PARAMETER(MAX=999,MAXP=99)
C
        COMMON/DUMMY/FFD(LSIZE)
        COMMON/NEXT/IUP(LSIZE,4),IDN(LSIZE,4)
        DIMENSION DELU(4),DELD(4)
     &  ,IXX(MAX),IYY(MAX),IZZ(MAX),ITT(MAX)
     &  ,IX1(MAXP),IY1(MAXP),IZ1(MAXP),IT1(MAXP)
        DIMENSION Q(MAX),EPS(MAX),IN(MAX),IT(MAX)
        DIMENSION Q1(MAXP),EPS1(MAXP)
C
      NNUM=0
      DO 4 N=1,MAX
      Q(N)=0.0
4     CONTINUE
      DO 6 N=1,MAXP
      Q1(N)=0.0
6     CONTINUE
C
      DO 10 NN=1,LSIZE
      M1=NN
      IF(NNUM.EQ.MAX)GOTO10
C
      FFF=FFD(M1)
      DO 11 MU=1,4
      M2=IUP(M1,MU)
      M3=IDN(M1,MU)
      FFDU=FFD(M2)
      FFDD=FFD(M3)
      DELU(MU)=ABS(FFF)-ABS(FFDU)
      DELD(MU)=ABS(FFF)-ABS(FFDD)
11    CONTINUE
C
      DO 15 MU=1,4
      IF(DELU(MU).LT.0.OR.DELD(MU).LT.0)GOTO10
15    CONTINUE
      IF(ABS(FFF).LT.QCUT)GOTO10
      NNUM=NNUM+1
      Q(NNUM)=FFF
      EPS(NNUM)=DELU(1)+DELD(1)+DELU(2)+DELD(2)
     &+DELU(3)+DELD(3)+DELU(4)+DELD(4)
      IN(NNUM)=NN
C
10    CONTINUE
C
        LX123=LX1*LX2*LX3
        LX12=LX1*LX2
        DO 30 NN=1,MAX
        IF(Q(NN).EQ.0)GOTO30
        JMAX=IN(NN)-1
        IS4=JMAX/(LX123)+1
        ITT(NN)=IS4
        JMAX=JMAX-(IS4-1)*LX123
        IS3=JMAX/(LX12)+1
        IZZ(NN)=IS3
        JMAX=JMAX-(IS3-1)*LX12
        IS2=JMAX/LX1+1
        IYY(NN)=IS2
        JMAX=JMAX-(IS2-1)*LX1
        IS1=JMAX+1
        IXX(NN)=IS1
30      CONTINUE
C
        NUM=0
        DO 32 NN=1,NNUM
        IF(NUM.EQ.MAXP)GOTO32
        DO 33 MM=1,NNUM
        IF(MM.EQ.NN)GOTO33
        IF(ABS(Q(MM)).LT.ABS(Q(NN))) GOTO33
        DELX=ABS(IXX(NN)-IXX(MM))
        IF(DELX.GT.LX1/2) DELX=LX1-DELX
        DELY=ABS(IYY(NN)-IYY(MM))
        IF(DELY.GT.LX2/2) DELY=LX2-DELY
        DELZ=ABS(IZZ(NN)-IZZ(MM))
        IF(DELZ.GT.LX3/2) DELZ=LX3-DELZ
        DELT=ABS(ITT(NN)-ITT(MM))
        IF(DELT.GT.LX4/2) DELT=LX4-DELT
        ISS=DELX*DELX+DELY*DELY+DELZ*DELZ+DELT*DELT
        IF(ISS.GT.RRCUT)GOTO33
        GOTO32
33      CONTINUE
        NUM=NUM+1
        Q1(NUM)=Q(NN)
        IX1(NUM)=IXX(NN)
        IY1(NUM)=IYY(NN)
        IZ1(NUM)=IZZ(NN)
        IT1(NUM)=ITT(NN)
        EPS1(NUM)=EPS(NN)
32      CONTINUE
C
        WRITE(6,93) QCUT,RRCUT
93      FORMAT(' GLOBAL MAXIMA OF ABS(FFDUAL); Q,R.R CUTS =',F10.4,F8.2)
        WRITE(6,94)
94      FORMAT(' FFDUAL  X   Y   Z   T   AV DIFF ')
        DO 35 NN=1,MAXP
        IF(Q1(NN).EQ.0)GOTO35
        QQ=Q1(NN)
        EPSS=EPS1(NN)*0.125
        IS4=IT1(NN)
        IS3=IZ1(NN)
        IS2=IY1(NN)
        IS1=IX1(NN)
        WRITE(6,95) QQ,IS1,IS2,IS3,IS4,EPSS
95      FORMAT(F10.4,4I5,F10.5)
35      CONTINUE
C
        RETURN
        END
C************************************************************
C  jack-knife errors for av(i=1,..,num)
C***********************************************************
      SUBROUTINE JACKM(NUM,AV,AVER,ERR)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION AV(*)
C
      SUM=0.0d0
      DO 10 N=1,NUM
      SUM=SUM+AV(N)
10    CONTINUE
      SUM=SUM/NUM
C
      ESUM=0.0d0
      DO 12 N=1,NUM
      DIF=SUM-AV(N)
      ESUM=ESUM+DIF*DIF
12    CONTINUE
      ESUM=ESUM/NUM
C
      AVER=SUM
      ERR=SQRT(ESUM*NUM)/(NUM-1)
C
      RETURN
      END
C************************************************************
C  jack-knife errors for ratio avu to avd
C***********************************************************
      SUBROUTINE JACKMR(NUM,NMAX,AVU,AVD,AVER,ERR)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION AVU(NMAX),AVD(NMAX)
      DIMENSION DIFU(NMAX),DIFD(NMAX)
C
      AVER=0.0d0
      ERR=0.0d0
C
      SUMU=0.0d0
      SUMD=0.0d0
      DO 10 N=1,NUM
      SUMU=SUMU+AVU(N)
      SUMD=SUMD+AVD(N)
10    CONTINUE
      DO 12 N=1,NUM
      DIFU(N)=SUMU-AVU(N)
      DIFD(N)=SUMD-AVD(N)
12    CONTINUE
      DO 14 N=1,NUM
      IF(DIFD(N).EQ.0.0d0)THEN
      ERR=99.0d0
      GOTO99
      ENDIF
14    CONTINUE
C
      ASUM=0.0d0
      ESUM=0.0d0
      DO 16 N=1,NUM
      DIF=DIFU(N)/DIFD(N)
      ASUM=ASUM+DIF
      ESUM=ESUM+DIF*DIF
16    CONTINUE
      ASUM=ASUM/NUM
      ESUM=ESUM/NUM
      ESUM=ESUM-ASUM*ASUM
C
      IF(SUMD.NE.0.0) AVER=SUMU/SUMD
      IF(ESUM.GT.0.0) ERR=SQRT(ESUM*NUM)
C
99    CONTINUE
      RETURN
      END
C*************************************************************
C fitting masses
C*************************************************************
        SUBROUTINE JACK(NBIN,ISUB,LXI,NUMBIN,LMAX,COR,VAC)
        IMPLICIT REAL*8 (A-H,O-Z)
C
        COMMON/FIT/ACOR(200),SCOR(200),AMM(200),SMM(200),NTMAX
C
        DIMENSION COR(NUMBIN,LMAX),VAC(NUMBIN)
        DIMENSION CORR(200),AW(200)
C
        DO 1 I4=1,LMAX
        ACOR(I4)=0.0d0
        SCOR(I4)=0.0d0
        AMM(I4)=0.0d0
        SMM(I4)=0.0d0
1       CONTINUE
C
        DO 2 IB=1,NBIN
C
        IF(ISUB.EQ.0)THEN
        VACC=0.0d0
        DO 4 JB=1,NBIN
        IF(JB.EQ.IB)GOTO4
        VACC=VACC+VAC(JB)
4       CONTINUE
        ENDIF
        DO 6 NT=1,LMAX
        CORR(NT)=0.0d0
        DO 7 JB=1,NBIN
        IF(JB.EQ.IB)GOTO7
        CORR(NT)=CORR(NT)+COR(JB,NT)
7       CONTINUE
        IF(ISUB.EQ.0) CORR(NT)=CORR(NT)-VACC*VACC/(NBIN-1)
6       CONTINUE
        ANORM=CORR(1)
        IF(ANORM.EQ.0.0)GOTO99
        DO 8 NT=1,LMAX
        CORR(NT)=CORR(NT)/ANORM
8       CONTINUE
        CALL FITT(LXI,CORR,AW,NTMAX)
        DO 9 NT=1,LMAX
        ACOR(NT)=ACOR(NT)+CORR(NT)
        SCOR(NT)=SCOR(NT)+CORR(NT)**2
        AMM(NT)=AMM(NT)+AW(NT)
        SMM(NT)=SMM(NT)+AW(NT)**2
9       CONTINUE
2       CONTINUE
C
        DO 10 NT=1,LMAX
        ACOR(NT)=ACOR(NT)/NBIN
        SCOR(NT)=SCOR(NT)/NBIN
        SCOR(NT)=(SCOR(NT)-ACOR(NT)*ACOR(NT))*NBIN
        IF(SCOR(NT).GT.0.0) SCOR(NT)=SQRT(SCOR(NT))
        AMM(NT)=AMM(NT)/NBIN
        SMM(NT)=SMM(NT)/NBIN
        SMM(NT)=(SMM(NT)-AMM(NT)*AMM(NT))*NBIN
        IF(SMM(NT).GT.0.0) SMM(NT)=SQRT(SMM(NT))
10      CONTINUE
C
        IF(ISUB.EQ.0)THEN
        VACC=0.0d0
        DO 22 JB=1,NBIN
        VACC=VACC+VAC(JB)
22      CONTINUE
        ENDIF
        DO 24 NT=1,LMAX
        CORR(NT)=0.0d0
        DO 26 JB=1,NBIN
        CORR(NT)=CORR(NT)+COR(JB,NT)
26      CONTINUE
        IF(ISUB.EQ.0) CORR(NT)=CORR(NT)-VACC*VACC/NBIN
24      CONTINUE
        ANORM=CORR(1)
        IF(ANORM.EQ.0.0)GOTO99
        DO 28 NT=1,LMAX
        CORR(NT)=CORR(NT)/ANORM
28      CONTINUE
        CALL FITT(LXI,CORR,AW,NTMAX)
        DO 29 NT=1,LMAX
        ACOR(NT)=CORR(NT)
        AMM(NT)=AW(NT)
29      CONTINUE
C
99      RETURN
        END
C*********************************************************************
C effective energies using a local cosh fit - lattice length LT.
C AW(timedif+1) is input corrln function;
C BW(nt) is eff energy from  time differences nt-1 to nt.
C*********************************************************************
        SUBROUTINE FITT(LT,AW,BW,NTMAX)
        IMPLICIT REAL*8 (A-H,O-Z)
C
        DIMENSION AW(200),BW(200)
C
        DO 1 I=1,LT/2+1
        BW(I)=0.0d0
1       CONTINUE
C
        NTMAX=LT/2+1
        DO 2 I4=1,LT/2
C
        IF(AW(I4).LE.0.000001.OR.AW(I4+1).LE.0.000001)THEN
        NTMAX=I4-1
        GOTO99
        ENDIF
        AMSM=0.0d0
        FTM=(AW(I4)/AW(I4+1))
        IF(FTM.GT.1.0)THEN
        AML=DLOG(FTM)
        AMU=DLOG(2.0d0*FTM)
        DO 13 NS=1,20
        AMS=(AML+AMU)/2
        FTS=(EXP(-AMS*(I4-1))+EXP(-(LT-I4+1)*AMS))
     &  /(EXP(-AMS*(I4))+EXP(-(LT-I4)*AMS))
        IF(FTS.LT.FTM)THEN
        AML=AMS
        ELSE
        AMU=AMS
        ENDIF
13      CONTINUE
        AMSM=(AML+AMU)/2
        ELSE
        NTMAX=I4-1
        GOTO99
        ENDIF
        BW(I4)=AMSM
C
2       CONTINUE
C
99      CONTINUE
        RETURN
        END
C**********************************************************
C  VECTOR MATRIX MULTIPLY ... 5*5 COMPLEX
C**********************************************************
        SUBROUTINE VMX(NNN1,A,B,C,NNN2)
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
        COMPLEX A(NCOL,NCOL),B(NCOL,NCOL),C(NCOL,NCOL),CSUM
C
        DO 1 J=1,NCOL
        DO 1 I=1,NCOL
        CSUM=(0.0,0.0)
        DO 2 K=1,NCOL
        CSUM=CSUM+A(I,K)*B(K,J)
2       CONTINUE
        C(I,J)=CSUM
1       CONTINUE
C
        RETURN
        END
C**********************************************************
C  TRACE PRODUCT .. 5*5 COMPLEX
C**********************************************************
        SUBROUTINE TRVMX(NNN1,A,B,CC,NNN2)
        IMPLICIT REAL*8 (A-H,O-Z)
        PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
        COMPLEX A(NCOL,NCOL),B(NCOL,NCOL),CC
C
        CC=(0.0,0.0)
        DO 1 I=1,NCOL
        DO 1 K=1,NCOL
        CC=CC+A(I,K)*B(K,I)
1       CONTINUE
C
        RETURN
        END
C**********************************************************
C  HERMITIAN CONJUGATE
C**********************************************************
      SUBROUTINE HERM(NNN1,A11,DUM11,NNN2)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NCOL=2,NCOL2=NCOL*NCOL)
C
      COMPLEX A11(NCOL,NCOL),DUM11(NCOL,NCOL)
C
      DO 1  I=1,NCOL
      DO 1  J=1,NCOL
      DUM11(I,J)=CONJG(A11(J,I))
2     CONTINUE
1     CONTINUE
      DO 4  I=1,NCOL
      DO 4  J=1,NCOL
      A11(I,J)=DUM11(I,J)
6     CONTINUE
4     CONTINUE
C
      RETURN
      END

