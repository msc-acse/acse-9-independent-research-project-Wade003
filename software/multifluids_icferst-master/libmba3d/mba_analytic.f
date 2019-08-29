      Module mba3d_mba_analytic
C
      use mba3d_ani2
C
      contains
C
C ================================================================
      Subroutine mbaAnalytic(
C ================================================================
c group (M)
     &      nP, MaxP, nF, MaxF, nE, MaxE,
     &      XYP, IPF, IPE, lbF, lbE,
     &      nEStar, 
c group (Dev)
     &      nPv, nFv, nEv, IPV, IFV, IEV, 
     &      flagAuto, status,
c group (Q)
     &      MaxSkipE, MaxQItr,
     &      MetricFunction, Quality, rQuality,
c group (W)
     &      MaxWr, MaxWi, rW, iW,
     &      iPrint, iERR)
C ================================================================
      include 'status.fd'
      include 'magic.fd'
      include 'output.fd'
C ================================================================
C VARIABLES & PARAMETER are decribed in mba_nodal.f except
C
C  MetricFunction - integer function created by the user (see example
C                in file forlibmba.f)
C
C    Integer Function MetricFunction(x, y, z, Metric)
C
C  This routine creates a metric at the given point (x,y, z). The
C  metric is a 3x3 positive definite symmetric tensor:
C                M11   M12   M13
C      Metric =  M12   M22   M23
C                M13   M23   M33
C
C  Only the upper triangular part of array Metric must be defined.
C
C
C *** Authors: K. Lipnikov     (lipnikov@hotmail.com)
C              Yu. Vassilevski (vasilevs@dodo.inm.ras.ru)
C *** Date:   1997 - 2008
C *** Updates: see ChangeLog
C
C ==========================================================
C group (M)
C     Integer MaxP, MaxF, MaxE
      Real*8  XYP(3, *)
      Integer IPE(4, *), IPF(3, *)
      Integer IPV(*), lbF(*), lbE(*)

C group (Dev)
      Integer nFv, nEv, IFV(*), IEV(*) 
      Logical flagAuto
      Integer status

C group (Q)
      Integer  MaxSkipE, MaxQItr
      Real*8   Quality, rQuality
      Logical  flagAnalytic

      Integer  MetricFunction
      EXTERNAL MetricFunction

C group (W)
      Real*8  rW(*)
      Integer iW(*)

C group (Local variables)
      Logical flagFILE
      Integer chanelOUT, statusWork
      Real*8  hStar, tm0, tm1
      Character*80 message
      Character*(LineLenght) output(MaxLines)

      iERR = 0

      iPrintl   = mod(iPrint, 10)
      chanelOUT = iPrint / 10
      flagFILE  = iPrint.GE.10

      If(chanelOUT.GE.100) Call errMes(4201, 'ani_metric',
     &                          'output chanel number is wrong')


      If(flagFILE) Then
         Open(chanelOUT, file='aniMPI.log', status='UNKNOWN')
      End if


      nLines = 0
      nLoop = 1

c ... print header
      If(iPrintl.GE.1) Then
         If(flagFILE) Then
            Write(chanelOUT, 5001) Quality, nEStar, MaxQItr
         Else
            Write(*, 5001) Quality, nEStar, MaxQItr
         End if
      End if

      Call setStatus(flagAuto, status, iPrint)
      statusWork = status


c ... starting clocks
      Call mba3d_seconds(tm0)


c ... memory for sequantial running (is cleaned)
c     iW(1) is overloaded to return colors & nQItr
      ilbP = 1
      iIPF = 1
      iIPE = iIPF + 4 * MaxF
      iIPP = iIPE + 5 * MaxE
      iICP = iIPP + MaxP
      iIHolP = iICP + MaxP
      iIHolF = iIHolP + MaxP
      iIHolE = iIHolF + MaxF 
      iIEP = iIHolE + MaxE
      iIFE = iIEP + MaxP
      iIEE = iIFE + 4 * MaxE
      iL1E = iIEE + 4 * MaxE
      iL2E = iL1E + 2 * MaxE
      iIEPw = iL2E + MaxE
      inEPw = iIEPw + max(4 * MaxE, 3 * MaxF)
      iIPEw = inEPw + MaxP
      iiSE  = iIPEw + 4 * nE

      miLINTRP = MaxWi - iiSE
      If(miLINTRP.LE.MaxP + MaxF) Then
         iERR = 1001
         Write(message,'(A,I10)')
     &        'The approximate size of iW is ', iiSE + MaxP + MaxF
         Call errMes(iERR, 'mbaAnalytic', message)
      End if


c ... memory for sequantial running (is cleaned)
c     rW(1) is overloaded to return total time
      iHesP = 1
      idG = iHesP + 6 * MaxP
      iqE = idG + MaxP
      iHesPw = iqE + MaxE
      iXYPw = iHesPw + 6 * nP
      irSE = iXYPw + 3 * nP

      mrLINTRP = MaxWr - irSE
      If(mrLINTRP.LT.MaxE) Then
         iERR = 1002
         Write(message,'(A,I10)')
     &        'The approximate size of rW is ', irSE + nE
         Call errMes(iERR, 'mbaAnalytic', message)
      End if

      Do n = 1, MaxWr
         rW(n) = 0D0
      End do

      Do n = 1, MaxWi
         iW(n) = 0
      End do


c ... auxiliary structure of the data
      m = iIPF - 1
      Do n = 1, nF
         Do i = 1, 3
            iW(m + i) = IPF(i, n)
         End do
         iW(m + 4) = lbF(n)
         m = m + 4
      End do

      m = iIPE - 1
      Do n = 1, nE
         Do i = 1, 4
            iW(m + i) = IPE(i, n)
         End do
         iW(m + 5) = lbE(n)
         m = m + 5
      End do


c ... metric field is generated by user's formulae in the physical space
      Call iniQ_analytic(nP, XYP, MetricFunction, rW(iHesP))


c ... scale geometry to unit cube
      Call scale2Cube(nP, XYP, .TRUE.)


c ... quality of tetrahedra is computed
      Call makQ(
     &     nLoop,
c group (M)
     &     nP, nE, XYP, IPE, nEv, IEV,
     &     nEStar, hStar,
c group (Q)
     &     rW(iHesP), rW(idG), rW(iqE))


c ... setting the fixed metric for the future loops
      nPw = nP
      nEw = nE
      Call copyMeshData(nP, nE, XYP, rW(iHesP), IPE, 
     &                  rW(iXYPw), rW(iHesPw), iW(iIPEw))


c... runing the basic algorithm for the global grid
      flagAnalytic = .TRUE.
      Call ani2(
c group (M)
     &     nP, MaxP, nF, MaxF, nE, MaxE, 
     &     XYP, iW(iIPF), iW(iIPE), 
     &     nEStar, hStar,
     &     iW(iICP), iW(iIPP), iW(iIEP),
     &     iW(iIFE), iW(iIEE),
     &     iW(iL1E), iW(iL2E),
     &     iW(iIHolP), iW(iIHolF), iW(iIHolE),
     &     iW(iIEPw), iW(inEPw),
     &     miLINTRP, mrLINTRP, iW(iIPEw), iW(iiSE),
     &     flagAuto, statusWork,
c group (Dev)
     &      nPv, nFv, nEv, IPV, IFV, IEV,  
c group (Q)
     &     MaxSkipE, MaxQItr, MaxBasketsGrid,
     &     rW(iHesP), Quality, rQuality, 
     &     rW(idG), rW(iqE), 
     &     nPw, nEw, rW(iXYPw), rW(iHesPw), rW(irSE),
     &     MetricFunction, flagAnalytic,
c group (ERR)
     &     flagFILE, chanelOUT, nLines, output,
     &     iPrintl, nQItr, iERR)

      If(iERR.NE.0 .AND. iERR.NE.1000)
     &   Call errMes(iERR, 'mbaAnalytic', 'memory problems')


      Call mba3d_seconds(tm1)
      tm1 = tm1 - tm0

      If(iPrintl.GE.1) Then
         If(flagFILE) Then
            Do i = 1, nLines
               Write(chanelOUt,'(A)') output(i)
            End do

            Write(chanelOUT, 5003) nQItr, rQuality, nP, nF, nE, tm1
         Else
            Write(*, 5003) nQItr, rQuality, nP, nF, nE, tm1
         End if
      End if


      If(flagFILE) Close(chanelOUT)


c ... rescale geometry back
      Call scale2Cube(nP, XYP, .FALSE.)


c ... original structure of the data
      m = iIPF - 1
      Do n = 1, nF
         Do i = 1, 3
            IPF(i, n) = iW(m + i)
         End do
         lbF(n) = iW(m + 4)
         m = m + 4
      End do

      m = iIPE - 1
      Do n = 1, nE
         Do i = 1, 4
            IPE(i, n) = iW(m + i) 
         End do
         lbE(n) = iW(m + 5)
         m = m + 5
      End do

      Do n = 1, nP
         iW(n) = iW(iICP + n - 1)
      End do
      iW(nP + 1) = nQItr

      rW(1) = tm1
      rW(2) = rQuality
      rW(3) = hStar

      Do n = 1, nE
         rW(3 + n) = rW(iqE + n - 1)
      End do

      Return

 5001 Format(/,
     &    'STONE FLOWER! (1997-2008), version 2.1', /,
     &    'Target: Quality', F5.2, ' with', I9,
     &    ' tetrahedra for at most', I9, ' iterations',/)

 5002 Format(/,'Processor :', I4)

 5003 Format('Total:', I6, ' Q=', E10.4, '  #V,F,E:', I7,I8,I9,
     &       '  tm=', F6.1, 's',/)
      End Subroutine mbaAnalytic



C ================================================================
      Subroutine mbaAnalyticShort(
C ================================================================
c group (M)
     &      nP, MaxP, nF, MaxF, nE, MaxE,
     &      XYP, IPF, IPE, lbF,
     &      nEStar, status,
c group (W)
     &      MaxWr, MaxWi, rW, iW,
     &      iPrint, iERR)
C ================================================================
C group (M)
      Real*8  XYP(3, *)
      Integer IPE(4, *), IPF(3, *)
      Integer lbF(*)
      Integer status

C group (W)
      Real*8  rW(*)
      Integer iW(*)

C group (Local variables)
      Real*8   Quality, rQuality
      Logical  flagAuto
      EXTERNAL MetricFunction_ani
      Character*80 message

C ================================================================
      iERR = 0

      MaxWiAv =  7 * MaxP + nP + 9 * MaxF + 19 * MaxE + 13 * nE
      MaxWrAv = 14 * MaxP + 14 * nP + MaxE

      If(MaxWi.LE.MaxWiAv) Then
         iERR = 1001
         Write(message,'(A,I10)')
     &        'The approximate size of iW is ', MaxWiAv
         Call errMes(iERR, 'mbaAnalyticShort', message)
      End if

      If(MaxWr.LE.MaxWrAv) Then
         iERR = 1002
         Write(message,'(A,I10)')
     &        'The approximate size of rW is ', MaxWrAv
         Call errMes(iERR, 'mbaAnalyticShort', message)
      End if


c ... setup of missing parameter
      nPv = 0
      nFv = 0
      nEv = 0
      flagAuto = .TRUE.

      MaxSkipE = 200
      MaxQItr  = min(50 000, nE) 
      Quality  = 0.3D0

      iIPV = 1
      ilbE = iIPV + MaxF
      iIFV = ilbE + MaxE
      iIEV = iIFV + nFV 
      iiEnd = iIEV + nEV + 1


      Do n = 0, nE - 1
         iW(ilbE + n) = 1
      End do

      Call mbaAnalytic(
c group (M)
     &      nP, MaxP, nF, MaxF, nE, MaxE,
     &      XYP, IPF, IPE, lbF, iW(ilbE),
     &      nEStar, 
c group (Dev)
     &      nPv, nFv, nEv, iW(iIPV), iW(iIFV), iW(iIEV), 
     &      flagAuto, status,
c group (Q)
     &      MaxSkipE, MaxQItr,
     &      MetricFunction_ani, Quality, rQuality, 
c group (W)
     &      MaxWr, MaxWi - iiEnd, rW, iW(iiEnd),
     &      iPrint, iERR)

      Return
      End Subroutine mbaAnalyticShort
C
      End Module mba3d_mba_analytic