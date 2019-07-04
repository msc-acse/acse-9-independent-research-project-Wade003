!*==PROCS127.spg  processed by SPAG 6.72Dc at 17:0!*==PROCS127.spg  processed by SPAG 6.72Dc at 17:0
!
!         module procs127
!           implicit none
!         contains





!
!         module procs127
!           implicit none
!         contains

          PROGRAM MAIN
!          use procs127, only : TEST_R_I, PYTHON_SET_UP_RECBIS 
          IMPLICIT NONE 

!         integer, parameter :: nx=100, ny=100
!         integer, parameter :: nx=1560, ny=1200
         integer, parameter :: nx=100, ny=1
!         integer, parameter :: nonods=nx*ny, nsplt=8
!         integer, parameter :: nonods=nx*ny, nsplt=12
         integer, parameter :: nonods=nx*ny, nsplt=2
         integer, parameter :: mx_ncola = 5*nonods
        
         integer fina(nonods+1),cola(mx_ncola),whichd(nonods),SPLEVS(nsplt)
         integer count,ii,nod,col,ISUB,j,ncola,iexact
         INTEGER ISUB_MAX, ISUB_MIN
         INTEGER havwnod, havmat, na
         real, allocatable :: wnod(:), a(:)
         LOGICAL EXACT

          count=0
          do nod=1,nonods
             fina(nod)=count+1

             if(.false.) then
                do ii=-1,1,1
                   col=nod+ii
                   if((col.ge.1).and.(col.le.nonods)) then
                      count=count+1
                      cola(count)=col
                   endif
                end do
             else ! 5 pt stencil in 2d...
                col=nod-nx
                call set_col(col,count,cola,mx_ncola,nonods)
                col=nod-1
                call set_col(col,count,cola,mx_ncola,nonods)
                col=nod
                call set_col(col,count,cola,mx_ncola,nonods)
                col=nod+1
                call set_col(col,count,cola,mx_ncola,nonods)
                col=nod+nx
                call set_col(col,count,cola,mx_ncola,nonods)
             endif

          end do
          fina(nonods+1)=count+1
          ncola=count

          allocate(wnod(nonods))
          exact=.true. ! exact node balance in each subdomain
          iexact=1
!          exact=.false. ! exact node balance in each subdomain
          havmat=0
          na=ncola*havmat+1
          allocate(A(ncola*havmat+1))
          havwnod=2
          wnod=1.0
          ii=1
          na=0
!          wnod(1:nonods/2)=10.0

!          print *,'fina:',fina
!          print *,'cola:',cola

          SPLEVS(:)=2
!          SPLEVS(1)=5
!          SPLEVS(2)=2
!          SPLEVS(3)=2
!          SPLEVS(4)=2
!          SPLEVS(5)=2
!          SPLEVS(6)=2
!          SPLEVS(7)=2
!          SPLEVS(8)=2

!          SPLEVS(9)=2
!          SPLEVS(10)=2

          call PYTHON_SET_UP_RECBIS(WHICHD, SPLEVS,FINA,COLA, &
              &    WNOD,a, havwnod,havmat,iexact, NSPLT,NCOLA,NONODS,na )
!              &    NSPLT,NCOLA,NONODS, havwnod,WNOD,exact, havmat,a )


          PRINT *,'MAXVAL(WHICHD(:)):',MAXVAL(WHICHD(:))

          ISUB_MAX=0
          ISUB_MIN=100000000
          DO ISUB=1,MAXVAL(WHICHD(:))
             COUNT=0
             DO NOD=1,NONODS
                IF(WHICHD(NOD).EQ.ISUB) COUNT = COUNT+1
             END DO
             PRINT *,'ISUB,COUNT:',ISUB,COUNT
             ISUB_MAX=MAX(ISUB_MAX,COUNT)
             ISUB_MIN=MIN(ISUB_MIN,COUNT)
          END DO
          PRINT *,'ISUB_MAX, ISUB_MIN:',ISUB_MAX, ISUB_MIN

!          PRINT *,'WHICHD:',WHICHD
          print *,'nonods=',nonods
          PRINT *,'WHICHD:'

          do j=1,-ny
             PRINT *,WHICHD( (j-1)*nx +1: j*nx)
          end do


          STOP
          END PROGRAM MAIN



         subroutine one_d_row_stor(ncola, cola,fina, nx,ny,nonods,mx_ncola)
        INTEGER, INTENT(IN) :: nx,ny,nonods,mx_ncola
        INTEGER, INTENT(OUT) :: ncola
        INTEGER, INTENT(OUT) :: cola(mx_ncola),fina(nonods+1)
! local variables...
         integer count,nod,col

          count=0
          do nod=1,nonods
             fina(nod)=count+1

             if(.false.) then
                do ii=-1,1,1
                   col=nod+ii
                   if((col.ge.1).and.(col.le.nonods)) then
                      count=count+1
                      cola(count)=col
                   endif
                end do
             else ! 5 pt stencil in 2d...
                col=nod-nx
                call set_col(col,count,cola,mx_ncola,nonods)
                col=nod-1
                call set_col(col,count,cola,mx_ncola,nonods)
                col=nod
                call set_col(col,count,cola,mx_ncola,nonods)
                col=nod+1
                call set_col(col,count,cola,mx_ncola,nonods)
                col=nod+nx
                call set_col(col,count,cola,mx_ncola,nonods)
             endif

          end do
          fina(nonods+1)=count+1
          ncola=count
          return
          end subroutine 

!
          subroutine set_col(col,count,cola,ncola,nonods)
          integer col,count,ncola,nonods
          integer cola(ncola)
          if((col.ge.1).and.(col.le.nonods)) then
             count=count+1
             cola(count)=col
          endif
          return
          end subroutine set_col
!
!
          SUBROUTINE test_r_i(ORR,OII,  RR,NR,II,NI)
!! SUBROUTINE TO TEST THE REALS AND INTEGERS IN PHTHON
          IMPLICIT NONE
          INTEGER, INTENT(IN) :: NR,NI
!         ! REAL, INTENT(OUT) :: ORR(NR)
!         ! INTEGER, INTENT(OUT) :: OII(NI)
!         ! REAL, INTENT(IN) :: RR(NR)
!         ! INTEGER, INTENT(IN) :: II(NI)
          REAL, INTENT(OUT) :: ORR(NR)
          INTEGER, INTENT(OUT) :: OII(NI)
          REAL, INTENT(IN) :: RR(:)
          INTEGER, INTENT(IN) :: II(:)

          PRINT *,'RR:',RR
          PRINT *,'II:',II
          ORR=RR
          OII=II
          RETURN
          END SUBROUTINE test_r_i


!
! ************RECBIS IS THE KEY SUBROUTINE TO CALL FOR SPLITTING*********** 
! NOTE FOR AUTOENCODERS:
! If we want to use a specified number
! of input neurons (not a bad idea) then we might find the maximum number of
! nodes in each subdomain and just send down zero's into the unused nodes. 
!
!
! *****DOMAIN DECOMPOSITION METHOD***********
! The subdomain splitting method is a neural network - a recurrent mean field theory network. 
! This has the advantage that you can easily modify it (as we have done) to
! as exactly as possible balance the number of nodes in each subdomain. 
!
!
! SET_UP_RECBIS IS A SIMPLIFIED INTERFACE FOR IT **************************




! SET_UP_RECBIS IS A SIMPLIFIED INTERFACE FOR IT **************************

!
!
       SUBROUTINE python_set_up_recbis(WHICHD0, SPLEVS0,FINA0,COLA0, &
                     wnod,a, havwnod,havmat,iexact, NSPLT,NCOLA,NNOD,na) 
! This sub splits the domain into subdomains using recursive bisection. 
! SPLEVS(1..NSPLT) contains the dissection information. 
! NSPLT - THE NUMBER OF RECURSIVE GRAPH CUTS
! SPLEVS(I), i=1,NSPLT: FOR RECURSION =I THE NUMBER OF PARTITIONS is SPLEVS(I)
! FINA,COLA is the matrix sparcity using compressed row storage.  
! NNOD = no of nodes. 
! WHICHD(NOD)= subdomain number of node NOD. *********ONLY THING THAT IS RETURNED BY SUBROUTINE****
! if HAVWNOD.ne.0 then assume we have non-unifrom weight for the nodes 
! in WNOD and we decompose to make these about equal. 
! If HAVWNOD=2 then also use the nodal weight for the graph splitting
! If EXACT balance as best we can the number of nodes in each subdomain. 
! only works for uniform WNOD. 
        IMPLICIT NONE 
! 
  !      INTEGER, INTENT(IN) :: NSPLT,NCOLA,NNOD, havwnod, havmat
!
        INTEGER, INTENT(IN) :: NSPLT,NCOLA,NNOD,na, havwnod, havmat,iexact
!
        INTEGER, INTENT(OUT) :: WHICHD0(NNOD)
        INTEGER, INTENT(IN) :: SPLEVS0(NSPLT)
!        INTEGER, INTENT(IN) :: SPLEVS0(:)
        INTEGER, INTENT(IN) :: COLA0(NCOLA),FINA0(NNOD+1)
        REAL, INTENT(IN) :: WNOD(NNOD),a(na)
  !      REAL, INTENT(IN) :: WNOD(NNOD), A(NCOLA*HAVMAT)
     !   REAL, INTENT(IN) :: WNOD(NNOD)
    !    LOGICAL, INTENT(IN) :: EXACT
!        INTEGER, INTENT(IN) :: COLA0(:),FINA0(:)
! LOCAL VARIABLES...
        LOGICAL IN_PYTHON
        PARAMETER(IN_PYTHON=.false.)
        INTEGER, DIMENSION(:), ALLOCATABLE :: SPLEVS, WHICHD, FINA, COLA
        INTEGER NDOM
    !    real, DIMENSION(:), ALLOCATABLE :: wnod
        logical exact
    !    integer havwnod
     !   integer havmat
     !   real a(1)


     !   print *, "inside fortran subroutine"
     !   print *, "splev", SPLEVS0
     !   print *, "shape fina", size(fina0)
     !   print *, "shape cola", size(cola0)
     !   print *, "nsplit", nsplt
     !   print *, "ncola", ncola
        !print "ii", ii

!        stop


!      
! Convert variables to fortran type variables...
       ALLOCATE(SPLEVS(NSPLT))
       ALLOCATE(WHICHD(NNOD))
       ALLOCATE(FINA(NNOD+1))
       ALLOCATE(COLA(NCOLA))
!       allocate(wnod(nnod))
!       havwnod=2
!       WNOD=1.0
!       exact=.true.
       exact = (iexact==1)
!       havmat=0

     if(IN_PYTHON) then
       COLA(1:NCOLA) = COLA0(1:NCOLA) + 1 
       FINA(1:NNOD+1) = FINA0(1:NNOD+1) + 1
       SPLEVS(1:NSPLT) = SPLEVS0(1:NSPLT)
     else
       COLA(1:NCOLA) = COLA0(1:NCOLA) 
       FINA(1:NNOD+1) = FINA0(1:NNOD+1) 
       SPLEVS(1:NSPLT) = SPLEVS0(1:NSPLT)
     endif
       
       CALL SET_UP_RECBIS(SPLEVS,NSPLT,NDOM, &
                &  FINA,COLA,NCOLA,NNOD,WHICHD, havwnod,WNOD,exact, havmat,a)
     if(IN_PYTHON) then
       WHICHD0(1:NNOD) = WHICHD(1:NNOD) - 1
     else
       WHICHD0(1:NNOD) = WHICHD(1:NNOD)
     endif
       
       END SUBROUTINE python_set_up_recbis
!
!
!       end module




!
! Fortran code... 
       SUBROUTINE SET_UP_RECBIS(SPLEVS,NSPLT,NDOM, &
     &               FINA,COLA,NCOLA,NNOD,WHICHD, havwnod,WNOD,exact, havmat,a) 
! This sub splits the domain into subdomains using recursive bisection. 
! NDOM=no of subdomains. (not really needed but a good check). 
! SPLEVS(1..NSPLT) contains the dissection information. 
! NSPLT - THE NUMBER OF RECURSIVE GRAPH CUTS
! SPLEVS(I), i=1,NSPLT: FOR RECURSION =I THE NUMBER OF PARTITIONS is SPLEVS(I)
! FINA,COLA is the matrix sparcity using compressed row storage.  
! NNOD = no of nodes. 
! WHICHD(NOD)= subdomain number of node NOD. *********ONLY THING THAT IS RETURNED BY SUBROUTINE****
! FOR NESTED BISSECTION THE SUBDOMAINS ARE ORDERED: 
!  1   2   3   4   5   6   7   8   (SUBDOMAIN NUMBERS) 
!    1       2       3       4   (SUBDOMAIN 1 AND 2 MAKES UP SUBDOMAIN 1 ON THIS LEVEL) 
!        1               2
!                1
! In this case NDOM=8, NSPLT=3, SPLEVS(1)=2, SPLEVS(2)=2, SPLEVS(3)=2. 
! if HAVWNOD.ne.0 then assume we have non-unifrom weight for the nodes 
! in WNOD and we decompose to make these about equal. 
! If HAVWNOD=2 then also use the nodal weight for the graph splitting
! If EXACT balance as best we can the number of nodes in each subdomain. 
! only works for uniform WNOD. 
        IMPLICIT NONE 
! 
        INTEGER NDOM,NSPLT,SPLEVS(NSPLT)
        INTEGER NCOLA,NNOD, HAVWNOD, HAVMAT
!
        INTEGER WHICHD(NNOD)
        INTEGER COLA(NCOLA),FINA(NNOD+1)
        REAL WNOD(NNOD), A(NCOLA*HAVMAT) 

! LOCAL VARIABES...
      INTEGER NSUBAL,MAXNLA,MMXTOT
      INTEGER MULLEV
      INTEGER MXNSPL,MXNTRE,NCHAIN,NRMEM
      PARAMETER(MXNSPL=40,MXNTRE=1000)
         INTEGER, DIMENSION(:), ALLOCATABLE :: LWICHD
         INTEGER, DIMENSION(:), ALLOCATABLE :: FITREE,TREE
         REAL, DIMENSION(:), ALLOCATABLE :: SUBALA,LSUBAL
         REAL ALPHA,LODBAL,BETA,TOLER
         REAL, DIMENSION(:), ALLOCATABLE :: X,Y,LX,LY
! MXNTRE=max no of entries in TREE (suggest 1000 
! so that max no of sub-domains is about 500). 
         REAL, DIMENSION(:), ALLOCATABLE :: LA
         INTEGER, DIMENSION(:), ALLOCATABLE :: LCOLA,LFINA
!
        REAL, DIMENSION(:), ALLOCATABLE :: RMEM
        REAL, DIMENSION(:), ALLOCATABLE :: LWNOD
        INTEGER, DIMENSION(:), ALLOCATABLE :: Q,QTEMP
        INTEGER, DIMENSION(:), ALLOCATABLE :: MAP
!
!
        LOGICAL EXACT

        MMXTOT=max(6*NNOD,1000)
        MAXNLA=3*NCOLA

        NRMEM=100*NNOD
        NSUBAL=NNOD

        !print *,'help1' 
        !print *,'nrmem,nNOD:',nrmem,nNOD


        ALLOCATE(X(NNOD),Y(NNOD),Q(NNOD),QTEMP(NNOD),MAP(NNOD), &
     &          LX(MMXTOT),LY(MMXTOT),LWNOD(MMXTOT),LA(NCOLA*3*havmat)) 
        ALLOCATE(RMEM(NRMEM),SUBALA(NSUBAL),LSUBAL(MXNTRE))
        ALLOCATE(LCOLA(MAXNLA),LFINA(MMXTOT+MULLEV)) 
        ALLOCATE(LWICHD(MMXTOT),FITREE(MXNSPL+2),TREE(MXNTRE))
        X=0.0; Y=0.0

      
        SUBALA(:)=1.0
        TOLER=1.E-4
!        NCHAIN=3000
        NCHAIN=300
        ALPHA=0.98
        BETA=0.9
!        LODBAL=0.1
        LODBAL=1.0
!        LODBAL=10.0
!        MULLEV=4
        MULLEV=6
        A=1.0; LA=1.0
        LWNOD(1:NNOD) = WNOD(1:NNOD) 

!
         CALL RECBIS(SPLEVS,NSPLT,MULLEV,  &
     &           RMEM,NRMEM,  &
     &           FINA,COLA,A,NCOLA,NNOD,MMXTOT,NDOM,WHICHD, &
     &           LFINA,LCOLA,LA,MAXNLA, &
     &           NSUBAL, &
     &           LWICHD,MAP  ,WNOD, LWNOD,  Q,QTEMP, &
     &           SUBALA,EXACT,ALPHA,LODBAL,BETA,TOLER,NCHAIN, &
     &           X,Y,LX,LY,HAVMAT, HAVWNOD )

         RETURN
         END
!
!
!
      SUBROUTINE RECBIS(Splevs,Nsplt,Mullev,Rmem,Nrmem,Fina,Cola,A,     &
                      & Ncola,Nnod,Mmxtot,Ndom,Whichd,Lfina,Lcola,La,   &
                      & Maxnla,Nsubal,Lwichd,Map,Wnod,Lwnod,Q,Qtemp,    &
                      & Subala,Exact,Alpha,Lodbal,Beta,Toler,Nchain,X,Y,&
                      & Lx,Ly,Havmat, HAVWNOD)
      IMPLICIT NONE
!*--RECBIS43
!*** Start of declarations inserted by SPAG
      INTEGER i , ii , iii , inlev , isub , maxtot , mxnla2 , Nchain ,  &
            & nlevel , Nrmem , ntree
!*** End of declarations inserted by SPAG
! NB MAP & QTEMP can contain the same storage location.
! This sub splits the domain using recursive bisection.
! SPLEVS(1..NSPLT) contains the dissection information.
! NSPLT - THE NUMBER OF RECURSIVE GRAPH CUTS
! SPLEVS(I), i=1,NSPLT: FOR RECURSION =I THE NUMBER OF PARTITIONS is SPLEVS(I)
! NB LFINA,LCOLA,LMIDPA,LA are defined(working arrays) inside sub.
! WNOD contains the weights of each node.
! MULLEV=no of multi-grid levels.
! MMXTOT is the maximum value of MAXTOT.
! NB it is recommended that MMXTOT=2*NNOD, MAXNLA=2*NCOLA.
! SUBALA() conatains subdomain balancing information.
! HAVMAT=0 if we are not worried about using weights.
! NRMEM=100*NNOD say.
! EXACT WILL EXACTLY BALANCE THE LOAD
      INTEGER Havmat, HAVWNOD
      INTEGER MXNSPL , MXNTRE , Nsubal , Maxnla , Mmxtot
      INTEGER Mullev
      PARAMETER (MXNSPL=100,MXNTRE=100000)
      INTEGER Splevs(Nsplt)
      INTEGER lcount , sub , level , lnnod , micol , Nsplt , lnod
      INTEGER lncola , across , look , lndom , count , count2
      INTEGER accum
      INTEGER Lwichd(Mmxtot) , Map(Nnod)
      INTEGER fitree(MXNSPL+2) , tree(MXNTRE)
      INTEGER Q(Nnod) , Qtemp(Nnod)
      REAL Subala(Nsubal) , lsubal(MXNTRE)
      REAL Alpha , Lodbal , Beta , Toler
      REAL X(Nnod) , Y(Nnod) , Lx(Mmxtot) , Ly(Mmxtot)
! MXNTRE=max no of entries in TREE (suggest 1000
! so that max no of sub-domains is about 500).
      REAL La(Maxnla*Havmat)
!
      INTEGER Lcola(Maxnla) , Lfina(Mmxtot+Mullev)
!
      REAL Rmem(Nrmem)
      REAL Wnod(Nnod) , Lwnod(Mmxtot)
!
      INTEGER Ncola , Nnod , Ndom
!
      INTEGER Whichd(Nnod)
      INTEGER Cola(Ncola) , Fina(Nnod+1)
!
      REAL A(Ncola)
      LOGICAL Exact
!
!
! FORM POINTERS TO TREE.
       !print *,'nrmem:',nrmem
      !WRITE (*,*) 'MXNTRE,NSPLT:' , MXNTRE , Nsplt
      nlevel = Nsplt + 1
      count = 0
      inlev = 1
      DO level = 1 , nlevel
         fitree(level) = count + 1
         IF ( level.EQ.1 ) THEN
            count = count + 1
            inlev = 1
            accum = 1
         ELSE
            inlev = inlev*Splevs(level-1)
            accum = accum + inlev
            count = accum
         ENDIF
      ENDDO
      fitree(nlevel+1) = count + 1
!
      ntree = count
! Now put entries in TREE
      DO level = nlevel , 1 , -1
         IF ( level.EQ.nlevel ) THEN
            ii = 0
            DO count = fitree(level) , fitree(level+1) - 1
               ii = ii + 1
               tree(count) = ii
            ENDDO
            Ndom = ii
         ELSE
! Choose minimum value from future level.
            across = 1
            DO count = fitree(level) , fitree(level+1) - 1
               tree(count) = tree(across+fitree(level+1)-1)
               across = across + Splevs(level)
            ENDDO
         ENDIF
      ENDDO
!
      !WRITE (*,*) 'splevs:' , Splevs
      !WRITE (*,*) 'fitree:' , (fitree(i),i=1,nlevel+1)
      DO level = 1 , -nlevel
         !WRITE (*,*) 'TREE:' ,                                          &
         !          & (tree(count),count=fitree(level),fitree(level+1)-1)
      ENDDO
 
!
!
!
! Split at level 1.
      lndom = Splevs(1)
      !WRITE (*,*) 'LNDOM,SPLEVS(1):' , lndom , Splevs(1)
!
      CALL COPYI(Lfina,Fina,Nnod+1)
      CALL COPYI(Lcola,Cola,Ncola)
      IF ( Havmat.EQ.1 ) CALL COPYR(La,A,Ncola)
      if(havwnod.ne.0) CALL COPYR(Lwnod,Wnod,Nnod)
      CALL COPYR(Lx,X,Nnod)
      CALL COPYR(Ly,Y,Nnod)
      lnnod = Nnod
!      maxtot = max(MIN(4*lnnod,Mmxtot),1000)
!      maxtot = MIN(4*lnnod,Mmxtot)
      maxtot = Mmxtot
      !print *,'1-LNNOD,Mmxtot:',LNNOD,Mmxtot
      CALL SPLIT(Lwichd,maxtot,lnnod,Mullev,lndom,  &
               & Lfina,Lcola,La,Maxnla,Q,Qtemp,Lwnod,lsubal,Exact,Alpha,  &
               & Lodbal,Beta,Toler,Nchain,Lx,Ly,Havmat, HAVWNOD)
!
!
      DO i = 1 , Nnod
         IF ( Lwichd(i).GE.1 ) Whichd(i) = tree(1+Lwichd(i))
         IF ( Lwichd(i).LE.0 ) Whichd(i) = Lwichd(i)
      ENDDO
      !WRITE (*,*) 'after split'
      !WRITE (*,*) 'splevs:' , Splevs
      !WRITE (*,*) 'fitree:' , (fitree(i),i=1,nlevel+1)
      DO level = 1 , nlevel
        ! WRITE (*,*) '1-TREE:' ,                                          &
        !           & (tree(count),count=fitree(level),fitree(level+1)-1)
      ENDDO
!
! count no of nodes in each subdomain
      DO isub = 1 , lndom
         ii = 0
         DO i = 1 , Nnod
            IF ( Lwichd(i).EQ.isub ) ii = ii + 1
         ENDDO
      ENDDO
!
      DO level = 2 , nlevel - 1
         lndom = Splevs(level)
         iii = 0
         DO lcount = fitree(level) , fitree(level+1) - 1
! Look for entries in WHICHD that have values of TREE(LCOUNT),
! form another graph from these and decompose this graph.
            isub = tree(lcount)
            lnnod = 0
            DO i = 1 , Nnod
               Map(i) = 0
            ENDDO
            DO i = 1 , Nnod
               IF ( Whichd(i).EQ.isub ) THEN
                  lnnod = lnnod + 1
!                 MAPBAK(LNNOD)=I
                  Map(i) = lnnod
                  Lwnod(lnnod) = Wnod(i)
                  Lx(lnnod) = X(i)
                  Ly(lnnod) = Y(i)
               ENDIF
            ENDDO
!
            count = 0
            lnod = 0
            DO i = 1 , Nnod
               IF ( Whichd(i).EQ.isub ) THEN
                  lnod = lnod + 1
                  Lfina(lnod) = count + 1
                  DO count2 = Fina(i) , Fina(i+1) - 1
                     micol = Map(Cola(count2))
                     IF ( micol.NE.0 ) THEN
                        count = count + 1
                        Lcola(count) = micol
                        IF ( Havmat.EQ.1 ) La(count) = A(count2)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            Lfina(lnnod+1) = count + 1
            lncola = count
!
            !WRITE (*,*) '*************LNNOD=' , lnnod
!            maxtot = max(4*lnnod,1000)
!            maxtot = max(4*lnnod,1000)
            maxtot = Mmxtot
      !print *,'2-LNNOD,MAxtot:',LNNOD,MAxtot
            !mxnla2 = MIN(2*lncola,Maxnla)
            mxnla2 = Maxnla
            CALL SPLIT(Lwichd,maxtot,lnnod,Mullev,lndom,     &
                     & Lfina,Lcola,La,mxnla2,Q,Qtemp,Lwnod,lsubal,Exact,&
                     & Alpha,Lodbal,Beta,Toler,Nchain,Lx,Ly,Havmat, HAVWNOD)
!
! Use MAP again to obtain MAPBAK
            lnod = 0
            DO i = 1 , Nnod
               IF ( Whichd(i).EQ.isub ) THEN
                  lnod = lnod + 1
!                 MAPBAK(LNOD)=I
                  Map(lnod) = i
               ENDIF
            ENDDO
! Put this decomposition in the next level
            iii = iii + 1
            DO sub = 1 , lndom
               across = (iii-1)*Splevs(level) + sub
               look = fitree(level+1) - 1 + across
               isub = tree(look)
               DO i = 1 , lnnod
!                 IF(LWICHD(I).EQ.SUB) WHICHD(MAPBAK(I))=ISUB
                  IF ( Lwichd(i).LE.0 ) THEN
                     Whichd(Map(i)) = Lwichd(i)
                  ELSE
                     IF ( Lwichd(i).EQ.sub ) Whichd(Map(i)) = isub
                  ENDIF
               ENDDO
            ENDDO
!
         ENDDO ! DO lcount = fitree(level) , fitree(level+1) - 1
!         !print *,'here3'
      ENDDO ! DO level = 2 , nlevel - 1
!         !print *,'here4 exiting RECBIS'
      END
!*==COPYI.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
      SUBROUTINE COPYI(Lcola,Cola,Ncola)
      IMPLICIT NONE
!*--COPYI263
!*** Start of declarations inserted by SPAG
      INTEGER i , Ncola
!*** End of declarations inserted by SPAG
      INTEGER Lcola(Ncola) , Cola(Ncola)
      DO i = 1 , Ncola
         Lcola(i) = Cola(i)
      ENDDO
      END
!*==COPYR.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
      SUBROUTINE COPYR(Lcola,Cola,Ncola)
      IMPLICIT NONE
!*--COPYR277
!*** Start of declarations inserted by SPAG
      INTEGER i , Ncola
!*** End of declarations inserted by SPAG
      REAL Lcola(Ncola) , Cola(Ncola)
      DO i = 1 , Ncola
         Lcola(i) = Cola(i)
      ENDDO
      END
!*==SPLIT_SIMPL.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
      SUBROUTINE SPLIT_SIMPL(Simp_whichd,Nonods,Onsubd,Simp_fina,       &
                           & Simp_cola,Simp_a,Simp_ncola,Simp_wnod)
      IMPLICIT NONE
!*--SPLIT_SIMPL294
!*** Start of declarations inserted by SPAG
      INTEGER maxna , nchain , Nonods
!*** End of declarations inserted by SPAG
      INTEGER havmat
      INTEGER mxnlev
      parameter(mxnlev=4)
! **********NEW INTERFACE*********
! NONODS=number of vertices or nodes in the graph to be decomposed.
! SIMP_WNOD() is the weight of each node - for load balancing.
! ONSUBD=no of subdomains.
! SIMP_WHICHD = domain each node belongs too. THE OUTPUT OF THIS SUB
! SIMP_A contains the edge weights for the graph partitioning and stored
! using compact row storage in SIMP_FINA,SIMP_COLA.
 
! **********THE PREVIOUS INTERFACE*********
! MXNLEV is the maximum number of multi-grid levels. E.G. 4
! NLEVEL is now the number of multi-grid levels e.g. set to 4.
! NONODS=number of vertices or nodes in the graph to be decomposed.
! WNOD() is the weight of each node - for load balancing.
! ONSUBD=no of subdomains.
! SUBAL(1..ONSUBD) contains the load for each processor. e.g. =1.0
! LODBAL gives the importance of balancing the load default is 1
! if EXACT then as near as possible load balance.
! WNOD=weight for each node.
! WHICHD = domain each node belongs too. THE OUTPUT OF THIS SUB
! A contains the edge weights for the graph partitioning and stored
! using compact row storage in FINA,COLA.
! HAVMAT=1 if we strore the matrix.
! ALSO the compact row storage is in FINA, COLA
! RMEM contains a long real working array e.g. NRMEM=100 * NONODS
! Q,QTEMP contains working arrays.
! BETA controls the crytical temp e.g. BETA=0.9
! ALPHA how BETA is annealed donw e.g. ALPHA=0.98
! TOLER tolerence for the iterations e.g. TOLER=1.E-4
! NCHAIN =maximum no of iterations 3000 is suggested.
! MAXTOT=3*nonods is default.
!
      INTEGER maxtot , Onsubd , Simp_ncola
! MAXTOT is the maximum value of TOTNOD allowed.
      INTEGER Simp_whichd(Nonods)
      INTEGER Simp_cola(Simp_ncola) , Simp_fina(Nonods+1)
      REAL Simp_a(Simp_ncola) , Simp_wnod(Nonods)
! Local variables...
      LOGICAL exact
      REAL, allocatable :: rmem(:) , x(:) , y(:) , a(:) , wnod(:)
      INTEGER, allocatable :: q(:) , qtemp(:)
      INTEGER, allocatable :: cola(:) , fina(:) , whichd(:)
      INTEGER nrmem , nlevel , totnod
      REAL alpha , beta , toler , subal(Onsubd) , lodbal
      INTEGER nodlev(mxnlev+1) , ptcola(mxnlev+1)
      integer havwnod
 
      maxtot = 3*Nonods
      exact = .FALSE.
      nlevel = mxnlev
      maxna = 3*Simp_ncola
      nrmem = 100*Nonods
      subal(:) = 1.0
      toler = 1.E-4
      nchain = 3000
      alpha = 0.98
      beta = 0.9
      havmat = 1
      lodbal = 1.0
      havwnod=0
 
 
      ALLOCATE (rmem(nrmem),x(maxtot),y(maxtot),a(maxna),wnod(maxtot))
      ALLOCATE (q(Nonods),qtemp(Nonods))
      ALLOCATE (cola(maxna),fina(maxtot+nlevel))
      ALLOCATE (whichd(maxtot))
      cola(1:Simp_ncola) = Simp_cola(1:Simp_ncola)
      fina(1:Nonods+1) = Simp_fina(1:Nonods+1)
      a(1:Simp_ncola) = Simp_a(1:Simp_ncola)
      wnod(1:Nonods) = Simp_wnod(1:Nonods)
 
      CALL SPLIT(whichd,maxtot,Nonods,nlevel,Onsubd,    &
               & fina,cola,a,maxna,q,qtemp,wnod,subal,exact,alpha, &
               & lodbal,beta,toler,nchain,x,y,havmat, havwnod)
 
      Simp_whichd(1:Nonods) = whichd(1:Nonods)
 
      END
!*==SPLIT.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
 
 
!
!
      SUBROUTINE SPLIT(Whichd,Maxtot,Nonods,Nlevel,Onsubd,   &
                     & Fina,Cola,A,Maxna,Q,Qtemp,Wnod,Subal,Exact,Alpha,&
                     & Lodbal,Beta,Toler,Nchain,X,Y,Havmat, havwnod) 
      IMPLICIT NONE
!*--SPLIT387
!*** Start of declarations inserted by SPAG
      INTEGER i , ii , isub , na , Nchain , need , Nlevel , Nonods 
!*** End of declarations inserted by SPAG
      INTEGER Havmat
      INTEGER MAX_mxnlev
      PARAMETER(MAX_MXNLEV=100) 
! MXNLEV is the maximum number of multi-grid levels. E.G. 4
! NLEVEL is now the number of multi-grid levels e.g. set to 4.
! NONODS=number of vertices or nodes in the graph to be decomposed.
! WNOD() is the weight of each node - for load balancing.
! ONSUBD=no of subdomains.
! SUBAL(1..ONSUBD) contains the load for each processor. e.g. =1.0
! LODBAL gives the importance of balancing the load default is 1
! if EXACT then as near as possible load balance.
! WNOD=weight for each node.
! WHICHD = domain each node belongs too. THE OUTPUT OF THIS SUB
! A contains the edge weights for the graph partitioning and stored
! using compact row storage in FINA,COLA.
! HAVMAT=1 if we strore the matrix.
! ALSO the compact row storage is in FINA, COLA
! RMEM contains a long real working array e.g. NRMEM=100 * NONODS
! Q,QTEMP contains working arrays.
! BETA controls the crytical temp e.g. BETA=0.9
! ALPHA how BETA is annealed donw e.g. ALPHA=0.98
! TOLER tolerence for the iterations e.g. TOLER=1.E-4
! NCHAIN =maximum no of iterations 3000 is suggested.
! MAXTOT=3*nonods is default.
! IF EXACT balce the no of nodes in each subdomain.
      INTEGER nodlev(MAX_mxnlev+1) , ptcola(MAX_mxnlev+1)
!
      INTEGER Maxtot , totnod , Onsubd , Maxna, havwnod
      REAL Lodbal
! MAXTOT is the maximum value of TOTNOD allowed.
      INTEGER mxnlev
      INTEGER Whichd(Maxtot)
      INTEGER Cola(Maxna) , Fina(Maxtot+Nlevel)
      INTEGER Q(Nonods) , Qtemp(Nonods)
      REAL A(Maxna*Havmat) , Alpha , Beta , Toler , Subal(Onsubd)
      REAL WNOD(MAXTOT)
      LOGICAL Exact
!      REAL Rmem(Nrmem)
      REAL X(Maxtot) , Y(Maxtot)

      INTEGER, DIMENSION(:), ALLOCATABLE :: LIST,whichd_new,whichd_long
      REAL, DIMENSION(:), ALLOCATABLE :: V_NEW,B_NEW,TSCALE_NEW

      ALLOCATE(LIST(NONODS))

!
      CALL FORMA(Whichd,Maxtot,totnod,Nonods,nodlev,ptcola,Nlevel,Fina, &
               & Cola,A,Maxna,na,Whichd,Q,Qtemp,LIST,X,Y,Havmat)
! This sub finds the subdomains.
      ALLOCATE(whichd_new(maxtot))
      ALLOCATE(V_NEW(Totnod*Onsubd),B_NEW(Nonods*Onsubd),TSCALE_NEW(NONODS))

      !WRITE (*,*) 'totnod,onsubd=' , totnod , Onsubd
! !print out all the levels
!       DO 82 ILEVEL=NLEVEL,1,-1
!        NCOLA2=PTCOLA(ILEVEL+1)-PTCOLA(ILEVEL)
!        NNOD2=NODLEV(ILEVEL+1)-NODLEV(ILEVEL)
!         CALL DRAXYZ(X(NODLEV(ILEVEL)),Y(NODLEV(ILEVEL)),
!     &     FINA(NODLEV(ILEVEL)+ILEVEL-1),COLA(PTCOLA(ILEVEL)),
!     &     NCOLA2,NNOD2)
!82     CONTINUE
!       STOP

        !print *,'nonods,totnod,onsubd:',nonods,totnod,onsubd
!        stop 221
 
!
!      need = Onsubd*totnod + Onsubd*Nonods + Nonods + 1
!      IF ( Nrmem.LT.need ) THEN
!         WRITE (*,*) 'NOT ENOUGH MEMORY IN SPLIT'
!         WRITE (*,*) 'we need:' , need , ' reals'
!         WRITE (*,*) 'we have ' , Nrmem , ' reals'
!         STOP
!      ENDIF
      !WRITE (*,*) 'INSIDE SPLIT ABOUT TO ENTER NEUR'
      !WRITE (*,*) 'MAXNA,NA,onsubd:' , Maxna , na , Onsubd
     ! WRITE (*,*) 'mxnlev,NLEVEL,ALPHA,LODBAL,BETA,TOLER,NCHAIN:' ,     &
     !           & mxnlev , Nlevel , Alpha , Lodbal , Beta , Toler ,     &
     !           & Nchain
!
! NB V,B,TSCALE are stored in RMEM.
!        CALL NEUR(WHICHD,WHICHD,V,B,TSCALE,
       if(totnod.lt.nonods) stop 3932
       if(maxtot.lt.max(totnod,nonods)) stop 3933
       whichd_new(1:totnod)=whichd(1:totnod)
!      CALL NEUR(Whichd_long,Whichd_new,V_NEW,B_NEW,TSCALE_NEW, nodlev,ptcola,      &
      CALL NEUR(Whichd,Whichd_new,V_NEW,B_NEW,TSCALE_NEW, nodlev,ptcola,      &
              & totnod,Nonods,Onsubd,Fina,Cola,A,na,Nlevel,Alpha,Lodbal,&
              & Beta,Toler,Nchain,Havmat,exact, havwnod,wnod)
! count no of nodes in each subdomain
      !WRITE (*,*) 'inside split onsubd,nonods:' , Onsubd , Nonods
      DO isub = 1 , Onsubd
         ii = 0
         DO i = 1 , Nonods
            IF ( Whichd(i).EQ.isub ) ii = ii + 1
         ENDDO
         !WRITE (*,*) 'no of nodes in sub=' , isub , ' is =' , ii
      ENDDO
!        DO 79 ILEVEL=1,NLEVEL
!          NNOD2=NODLEV(ILEVEL+1)-NODLEV(ILEVEL)
!          MAXNA2=PTCOLA(ILEVEL+1)-PTCOLA(ILEVEL)
!
!            CALL DRAHAL(X(NODLEV(ILEVEL)),Y(NODLEV(ILEVEL)),
!     &        FINA(NODLEV(ILEVEL)+ILEVEL-1),
!     &        COLA(PTCOLA(ILEVEL)),MAXNA2,NNOD2,
!     &        WHICHD(NODLEV(ILEVEL)) )
!79       CONTINUE
      END
!*==COUDOM.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
      SUBROUTINE COUDOM(Whichd,Nonods,Fina,Cola,Ncola)
      IMPLICIT NONE
!*--COUDOM488
!*** Start of declarations inserted by SPAG
      INTEGER i , icol , icount , ii , imax , imin , isub , Ncola ,     &
            & noint , Nonods
!*** End of declarations inserted by SPAG
! this sub counts the number of subdomains.
      INTEGER Whichd(Nonods) , Fina(Nonods+1) , Cola(Ncola)
      imax = -1000
      imin = 1000
      DO i = 1 , Nonods
         imax = MAX(Whichd(i),imax)
         imin = MIN(Whichd(i),imax)
      ENDDO
      DO isub = imax , imin
         ii = 0
         DO i = 1 , Nonods
            IF ( Whichd(i).EQ.isub ) ii = ii + 1
         ENDDO
 
      ENDDO
! count the interface nodes
      noint = 0
      DO i = 1 , Nonods
         DO icount = Fina(i) , Fina(i+1) - 1
            icol = Cola(icount)
            IF ( Whichd(i).NE.Whichd(icol) ) noint = noint + 1
         ENDDO
      ENDDO
      noint = noint/2
      !print * , 'no of edges between subs =' , noint
      END
!*==RECBI2.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
      SUBROUTINE RECBI2(Splevs,Nsplt,Rmem,Nrmem,Fina,Cola,Midpa,A,Ncola,&
                      & Nnod,Ndom,Whichd,Lfina,Lcola,Lmidpa,La,Tempd,Km,&
                      & Kmtemp,Random,Power,Optim,Chopar,Chopa3,Chopa4, &
                      & Halo,Fitree,Tree,Mxntre,Lwichd,Map,Mxldom,Wnod, &
                      & Lwnod)
      IMPLICIT NONE
!*--RECBI2530
!*** Start of declarations inserted by SPAG
      INTEGER i , ii , iii , inlev , isub , Mxntre , nlevel , Nrmem ,   &
            & ntree
!*** End of declarations inserted by SPAG
! This sub splits the domain using recursive bisection.
! SPLEVS(1..NSPLT) contains the dissection information.
! NB LFINA,LCOLA,LMIDPA,LA are defined(working arrays) inside sub.
!         INTEGER MXSPLT
!         PARAMETER(MXSPLT=20)
!         INTEGER SPLEVS(MXSPLT)
      INTEGER Splevs(Nsplt)
      INTEGER lcount , sub , level , lnnod , micol , Nsplt , lnod
      INTEGER lncola , across , look , lndom , count , count2
      INTEGER Mxldom , accum
      INTEGER Lwichd(Nnod) , Map(Nnod)
      INTEGER Fitree(Nsplt+2) , Tree(Mxntre)
      INTEGER Power , Optim
! MXNTRE=max no of entries in TREE (suggest 1000
! so that max no of sub-domains is about 500).
      REAL La(Ncola)
!
      INTEGER Lmidpa(Nnod) , Lcola(Ncola) , Lfina(Nnod+1)
      REAL Lwnod(Nnod)
!
      INTEGER Halo
      REAL Rmem(Nrmem)
      REAL Wnod(Nnod)
!
      INTEGER Ncola , Nnod , Ndom
      REAL Km(Mxldom) , Kmtemp(Mxldom)
      INTEGER Tempd(Nnod)
!
      INTEGER Chopar(Nnod) , Chopa3(Nnod) , Chopa4(Nnod)
      INTEGER Whichd(Nnod)
      INTEGER Midpa(Nnod) , Cola(Ncola) , Fina(Nnod+1)
!
      REAL A(Ncola)
      LOGICAL Random
 
! FORM POINTERS TO TREE.
      nlevel = Nsplt + 1
      count = 0
      inlev = 1
      DO level = 1 , nlevel
         Fitree(level) = count + 1
         IF ( level.EQ.1 ) THEN
            count = count + 1
            inlev = 1
            accum = 1
         ELSE
            inlev = inlev*Splevs(level-1)
            accum = accum + inlev
            count = accum
         ENDIF
      ENDDO
      Fitree(nlevel+1) = count + 1
! NB NTREE might have to be worked out befor hand.
      ntree = count
! Now put entries in TREE
      DO level = nlevel , 1 , -1
         IF ( level.EQ.nlevel ) THEN
            ii = 0
            DO count = Fitree(level) , Fitree(level+1) - 1
               ii = ii + 1
               Tree(count) = ii
            ENDDO
            Ndom = ii
         ELSE
! Choose minimum value from future level.
            across = 1
            DO count = Fitree(level) , Fitree(level+1) - 1
               Tree(count) = Tree(across+Fitree(level+1)-1)
               across = across + Splevs(level)
            ENDDO
         ENDIF
      ENDDO
 
 
! Split at level 1.
      lndom = Splevs(1)
 
      CALL SPLIT2(Rmem,Nrmem,Fina,Cola,Midpa,A,Ncola,Nnod,lndom,Lwichd, &
                & Tempd,Km,Kmtemp,Random,Power,Optim,Chopar,Chopa3,     &
                & Chopa4,Halo,Wnod)
!
!
      DO i = 1 , Nnod
         IF ( Lwichd(i).GE.1 ) Whichd(i) = Tree(1+Lwichd(i))
         IF ( Lwichd(i).LE.0 ) Whichd(i) = Lwichd(i)
      ENDDO
!
!
      DO level = 2 , nlevel - 1
         lndom = Splevs(level)
         iii = 0
         DO lcount = Fitree(level) , Fitree(level+1) - 1
! Look for entries in WHICHD that have values of TREE(LCOUNT),
! form another graph from these and decompose this graph.
            isub = Tree(lcount)
            lnnod = 0
            DO i = 1 , Nnod
               Map(i) = 0
            ENDDO
            DO i = 1 , Nnod
               IF ( Whichd(i).EQ.isub ) THEN
                  lnnod = lnnod + 1
!                 MAPBAK(LNNOD)=I
                  Map(i) = lnnod
                  Lwnod(lnnod) = Wnod(i)
               ENDIF
            ENDDO
!
            count = 0
            lnod = 0
            DO i = 1 , Nnod
               IF ( Whichd(i).EQ.isub ) THEN
                  lnod = lnod + 1
                  Lfina(lnod) = count + 1
                  DO count2 = Fina(i) , Fina(i+1) - 1
                     micol = Map(Cola(count2))
                     IF ( micol.NE.0 ) THEN
                        count = count + 1
                        Lcola(count) = micol
                        La(count) = A(count2)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            Lfina(lnnod+1) = count + 1
            lncola = count
!
            DO i = 1 , lnnod
               DO count2 = Lfina(i) , Lfina(i+1) - 1
                  IF ( Lcola(count2).EQ.i ) Lmidpa(i) = count2
               ENDDO
            ENDDO
 
            CALL SPLIT2(Rmem,Nrmem,Lfina,Lcola,Lmidpa,La,lncola,lnnod,  &
                      & lndom,Lwichd,Tempd,Km,Kmtemp,Random,Power,Optim,&
                      & Chopar,Chopa3,Chopa4,Halo,Lwnod)
!
! Use MAP again to obtain MAPBAK
            lnod = 0
            DO i = 1 , Nnod
               IF ( Whichd(i).EQ.isub ) THEN
                  lnod = lnod + 1
!                 MAPBAK(LNOD)=I
                  Map(lnod) = i
               ENDIF
            ENDDO
! Put this decomposition in the next level
            iii = iii + 1
            DO sub = 1 , lndom
               across = (iii-1)*Splevs(level) + sub
               look = Fitree(level+1) - 1 + across
               isub = Tree(look)
               DO i = 1 , lnnod
!                 IF(LWICHD(I).EQ.SUB) WHICHD(MAPBAK(I))=ISUB
                  IF ( Lwichd(i).LE.0 ) THEN
                     Whichd(Map(i)) = Lwichd(i)
                  ELSE
                     IF ( Lwichd(i).EQ.sub ) Whichd(Map(i)) = isub
                  ENDIF
               ENDDO
            ENDDO
!
         ENDDO
      ENDDO
      END
!*==SPLIT2.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
!
      SUBROUTINE SPLIT2(Rmem,Nrmem,Fina,Cola,Midpa,A,Ncola,Nnod,Ndom,   &
                      & Whichd,Tempd,Km,Kmtemp,Random,Power,Optim,      &
                      & Chopar,Chopa3,Chopa4,Halo,Wnod)
      IMPLICIT NONE
!*--SPLIT2710
!*** Start of declarations inserted by SPAG
      INTEGER ii , iii , iiii , iiiv , iiv , iv , Nrmem
      REAL rbeta , rtr
!*** End of declarations inserted by SPAG
! This sub does the graph splitting.
!
      INTEGER Halo
      REAL Rmem(Nrmem)
!
      INTEGER Ncola , Nnod , Ndom
      REAL Km(Ndom) , Kmtemp(Ndom)
      REAL Wnod(Nnod)
      INTEGER Tempd(Nnod)
      INTEGER nparts , npart3 , npart4 , ndom2 , nnod2
!
      INTEGER Chopar(Nnod) , Chopa3(Nnod) , Chopa4(Nnod)
      INTEGER Whichd(Nnod)
      INTEGER Midpa(Nnod) , Cola(Ncola) , Fina(Nnod+1)
      INTEGER Power , Optim
      INTEGER, ALLOCATABLE :: Noinsu(:)
!
      REAL A(Ncola)
      LOGICAL Random

      ALLOCATE(Noinsu(ndom))
 
      rbeta = 1./(2.**MAX(0,Ndom-2))
      rbeta = SQRT(rbeta)
 
      CALL MESH(rbeta,Wnod,Fina,Cola,Midpa,A,Ncola,Nnod)
 
! This sub finds the subdomains using a replicator approach.
      nparts = INT(REAL(Nnod)**0.75)
      npart3 = INT(REAL(Nnod)**0.5)
      npart4 = INT(REAL(Nnod)**0.25)
 
      ndom2 = Ndom
      nnod2 = Nnod
      ii = Ndom*Nnod + 1
      iii = 2*Ndom*Nnod + 1
      iiii = 2*Ndom*Nnod + Ndom*nparts + 1
      iv = 2*Ndom*Nnod + Ndom*nparts + Ndom*npart3 + 1
      iiv = 2*Ndom*Nnod + Ndom*nparts + Ndom*npart3 + Ndom*npart4 + 1
      iiiv = 2*Ndom*Nnod + Ndom*nparts + Ndom*npart3 + Ndom*npart4 +    &
           & Ndom + 1
      IF ( iiiv+Ndom+1.GT.Nrmem ) THEN
        ! WRITE (0,*) 'ERROR: NOT ENOUGH MEMORY PASSED DOWN FOR '//      &
        !            &'PARTITIONING IN THIS MANNER.'
         GOTO 99999
      ENDIF
 
      CALL ANNEAL(nparts,npart3,npart4,ndom2,nnod2,Rmem(1),Rmem(ii),rtr,&
                & Rmem(iii),Rmem(iiii),Rmem(iv),Rmem(iiv),Rmem(iiiv),   &
                & Fina,Cola,Midpa,A,Ncola,Nnod,Ndom,Whichd,Tempd,Km,    &
                & Kmtemp,Random,Power,Optim,Chopar,Chopa3,Chopa4)
!
! Find a set of interface nodes. -RESCALE MATRIX FOR NEURAL NETWORK.
!
      IF ( Halo.EQ.0 ) THEN
!
         rbeta = 1.
         CALL MESH(rbeta,Wnod,Fina,Cola,Midpa,A,Ncola,Nnod)
!
         CALL NEURI(ndom2,nnod2,Rmem(1),Rmem(nnod2*(ndom2+1)+1),Fina,   &
                  & Cola,Midpa,A,Ncola,Nnod,Ndom,Whichd,                &
!                  & Rmem(nnod2*(ndom2+2)+1),                            &
                  & Noinsu,                            &
                  & Rmem(nnod2*(ndom2+2)+ndom2+2),                      &
                  & Rmem(nnod2*(ndom2+2)+2*ndom2+4),Optim)
 
!
      ENDIF
99999 END
!*==ANNEAL.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
      SUBROUTINE ANNEAL(Nparts,Npart3,Npart4,Ndom2,Nnod2,R,R2,Rtr,      &
                      & Rtrpar,Rtrpa3,Rtrpa4,Rtrsub,Prodm2,Fina,Cola,   &
                      & Midpa,A,Ncola,Nnod,Ndom,Whichd,Tempd,Km,Kmtemp, &
                      & Random,Power,Optim,Chopar,Chopa3,Chopa4)
      IMPLICIT NONE
!*--ANNEAL789
!*** Start of declarations inserted by SPAG
      REAL alpha , FRAC , rrscal
      INTEGER i , ii , isub , itcoun , its , nchain , NOITS
!*** End of declarations inserted by SPAG
!
! FRAC=termination criterion.
!
      INTEGER CHEND
!
      LOGICAL SIMPLE
! If we have no exceptances for CHEND iterations then end.
      PARAMETER (SIMPLE=.FALSE.)
! if KEEPLI then the cooling schedual attempts to keep
! a linear relation ship between the number if excepted
! iteration and ITS.
!        PARAMETER(CHEND=3,FRAC=0.03)
      PARAMETER (CHEND=3,FRAC=0.01)
!
      PARAMETER (NOITS=1000)
!        PARAMETER(NOITS=20)
!
      INTEGER Ncola , Nnod , Ndom
      REAL c , pcost , ncost , initic
      REAL Km(Ndom) , Kmtemp(Ndom)
      INTEGER sub , newsub , oldsub , newnod
      INTEGER Tempd(Nnod) , nexcp , sumex
      INTEGER chain
      INTEGER Nparts , Npart3 , Npart4 , Ndom2 , Nnod2
!
      INTEGER Chopar(Nnod) , Chopa3(Nnod) , Chopa4(Nnod)
      INTEGER Whichd(Nnod)
      INTEGER Midpa(Nnod) , Cola(Ncola) , Fina(Nnod+1)
      INTEGER Power , Optim
!
      REAL R(Ndom2*Nnod2) , R2(Ndom2*Nnod2)
      REAL Prodm2(Ndom2) , Rtrpar(Ndom2*Nparts)
      REAL Rtrsub(Ndom2)
      REAL Rtrpa3(Ndom2*Npart3) , Rtrpa4(Ndom2*Npart4)
      REAL Rtr
      REAL A(Ncola)
      LOGICAL yes , gotc
      LOGICAL Random
 
      alpha = 0
      nchain = 0
      IF ( alpha.LT.0.1 ) alpha = 0.98
 
      rrscal = REAL(Optim)/20.
      nchain = int( rrscal*real((Fina(Nnod+1)-1)*(Ndom-1))/100. )
 
! This sub finds a partition using an annealing approach.
 
!  Initialise, Find a good initial guess.
      Ncola = Fina(Nnod+1) - 1
      CALL INGUES(Fina,Cola,Midpa,A,Ncola,Nnod,Ndom,Whichd,Tempd,Random)
! Find an initial C.
      CALL FININI(nchain,initic,pcost,gotc,Nparts,Npart3,Npart4,Ndom2,  &
                & Nnod2,R,R2,Rtr,Rtrpar,Rtrpa3,Rtrpa4,Rtrsub,Prodm2,    &
                & Fina,Cola,Midpa,A,Ncola,Nnod,Ndom,Whichd,Tempd,Km,    &
                & Kmtemp,Random,Power,Chopar,Chopa3,Chopa4)
!        GOTC=.FALSE.
!        INITIC=100000
      DO i = 1 , Nnod
         Tempd(i) = Whichd(i)
      ENDDO
! work out R & R2 and PCOST
      CALL FIDERC(pcost,Nparts,Npart3,Npart4,Ndom2,Nnod2,R,R2,Rtr,      &
                & Rtrpar,Rtrpa3,Rtrpa4,Rtrsub,Prodm2,Fina,Cola,Midpa,A, &
                & Ncola,Nnod,Ndom,Whichd,Km,Kmtemp,Chopar,Chopa3,Chopa4,&
                & Power)
!
      sumex = 0
!        C=INITIC
      c = pcost*0.2
!        C=PCOST*0.2
      itcoun = 0
      DO its = 1 , NOITS
         itcoun = itcoun + 1
         IF ( (.NOT.gotc) .AND. ((its).EQ.10) ) THEN
            gotc = .TRUE.
            c = pcost*0.2
         ENDIF
!          IF(ITS.EQ.15) C=PCOST*0.2
!          IF(ITS.EQ.50) C=2.0E+11
         nexcp = 0
         DO chain = 1 , nchain
! Find a neighbourhood configuration
            CALL FINEIG(newsub,newnod,oldsub,Nparts,Npart3,Npart4,Ndom2,&
                      & Nnod2,R,R2,Rtr,Rtrpar,Rtrpa3,Rtrpa4,Rtrsub,     &
                      & Prodm2,Fina,Cola,Midpa,A,Ncola,Nnod,Ndom,Whichd,&
                      & Tempd,Power)
            IF ( newsub.NE.oldsub ) THEN
! find cost of neighbourhood config
!              CALL FICOSK(NCOST,PCOST,NEWSUB,OLDSUB,NEWNOD)
               CALL FICOSK(ncost,pcost,newsub,oldsub,newnod,Fina,Cola,  &
                         & Midpa,A,Ncola,Nnod,Ndom,Whichd,Tempd,Km,     &
                         & Kmtemp)
! See if we except this config.
               CALL EXCEPT(pcost,ncost,c,yes)
!               yes=.true.
               IF ( yes ) THEN
                  nexcp = nexcp + 1
                  pcost = ncost
                  Whichd(newnod) = newsub
                  DO sub = 1 , Ndom
                     Km(sub) = Kmtemp(sub)
                  ENDDO
                  IF ( .NOT.SIMPLE )                                    &
                     & CALL UPDERI(ncost,newsub,newnod,oldsub,Nparts,   &
                     & Npart3,Npart4,Ndom2,Nnod2,R,R2,Rtr,Rtrpar,Rtrpa3,&
                     & Rtrpa4,Rtrsub,Prodm2,Fina,Cola,Midpa,A,Ncola,    &
                     & Nnod,Ndom,Whichd,Km,Kmtemp,Chopar,Chopa3,Chopa4, &
                     & Power)
               ELSE
                  Tempd(newnod) = Whichd(newnod)
               ENDIF
            ENDIF
         ENDDO
        ! WRITE (*,*) 'pcost=' , pcost , ' its=' , its , ' NEXCP=' ,     &
        !           & nexcp , ' c=' , c
!
         c = alpha*c
!
         sumex = sumex + nexcp
         IF ( ((its/CHEND)*CHEND.EQ.its) .AND. (its.GT.30) ) THEN
! Check for convergence
            IF ( sumex.LT.CHEND*FRAC*nchain ) GOTO 100
            sumex = 0
         ENDIF
         IF ( itcoun*nchain.GT.500000 ) THEN
            CALL FIDERC(pcost,Nparts,Npart3,Npart4,Ndom2,Nnod2,R,R2,Rtr,&
                      & Rtrpar,Rtrpa3,Rtrpa4,Rtrsub,Prodm2,Fina,Cola,   &
                      & Midpa,A,Ncola,Nnod,Ndom,Whichd,Km,Kmtemp,Chopar,&
                      & Chopa3,Chopa4,Power)
!
            !WRITE (*,*) 'call FIDERC *******************'
            itcoun = 0
         ENDIF
      ENDDO
 100  CALL FIDERC(pcost,Nparts,Npart3,Npart4,Ndom2,Nnod2,R,R2,Rtr,      &
                & Rtrpar,Rtrpa3,Rtrpa4,Rtrsub,Prodm2,Fina,Cola,Midpa,A, &
                & Ncola,Nnod,Ndom,Whichd,Km,Kmtemp,Chopar,Chopa3,Chopa4,&
                & Power)
      !WRITE (*,*) 'pcost=' , pcost
! ********************* OUTPUT THE RESULTS **********************
      DO isub = 1 , Ndom
         ii = 0
         DO i = 1 , Nnod
            IF ( Whichd(i).EQ.isub ) ii = ii + 1
         ENDDO
         !WRITE (*,*) 'no of nodes in sub=' , isub , ' is =' , ii
      ENDDO
!
      END
!*==INGUES.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
      SUBROUTINE INGUES(Fina,Cola,Midpa,A,Ncola,Nnod,Ndom,Whichd,Tempd, &
                      & Random)
      IMPLICIT NONE
!*--INGUES952
!*** Start of declarations inserted by SPAG
      INTEGER i , iarg , icol , idom , indval , iran , istart , its ,   &
            & ival , Ndom , nq , nqtemp
      REAL RAN1 , rran , sub
!*** End of declarations inserted by SPAG
!
      INTEGER Ncola , Nnod
      INTEGER Whichd(Nnod) , Midpa(Nnod) , Cola(Ncola) , Fina(Nnod+1)
      REAL A(Ncola)
!
      INTEGER Tempd(Nnod)
!
!      INTEGER MINVAL,TEMPD(MAXNOD),ROW,COUNT
      INTEGER minval , row , count
      LOGICAL Random
!
      IF ( Random ) THEN
         DO i = 1 , Nnod
!
            iarg = 1
!
            rran = RAN1(iarg)
            IF ( rran.GE.1 ) rran = 0.999999
            iran = INT(rran*REAL(Nnod*Ndom)+1.)
!
            sub = ((iran-1)/Nnod) + 1
            Whichd(i) = sub
            Tempd(i) = sub
         ENDDO
      ELSE
         DO i = 1 , Nnod
            Whichd(i) = 0
            Tempd(i) = 0
         ENDDO
!       write(*,*)'**********************'
!      do 11 i=1,nnod
!        i=286
!        write(*,*)'i=',i,'fina(I)=',fina(I),' fina(i+1)=',fina(i+1)
!11    continue
! Find MIN VALANCY NODE.
         DO idom = 1 , Ndom
!         DO 30 IDOM=1,1
            nq = 0
            istart = 1
 20         nqtemp = nq + 1
            !print * , 'idom :' , idom , ' nq1 :' , nq , ' nqtemp :' ,   &
            !    & nqtemp
!
            minval = 10000
            indval = 0
            DO i = 1 , Nnod
!           write(*,*)'i=',i,'whichd(i)=',whichd(i)
               IF ( Whichd(i).EQ.0 ) THEN
                  ival = 0
!            write(*,*)'fina(i),fina(i+1):',fina(i),fina(i+1)
                  DO count = Fina(i) , Fina(i+1) - 1
!               ival=ival+1
                     IF ( Whichd(Cola(count)).EQ.0 ) ival = ival + 1
                  ENDDO
!             write(*,*)'ival=',ival,' i=',i
                  IF ( ival.LT.minval ) THEN
                     minval = ival
                     indval = i
                  ENDIF
               ENDIF
            ENDDO
!         write(*,*)'minval=',minval
!         write(*,*)'indval=',indval
!         TEMPD(NQ)=INDVAL
!
            Whichd(indval) = idom
            nq = nqtemp
            Tempd(nq) = indval
!         !print*, 'tempd(nq), :', tempd(nq), ' nq :', nq
!
!
!         !print*, 'nq :', nq, ' istart:', istart, ' nqtemp :', nqtemp
!
! GROW DOMAIN FROM A POINT
! find all the nodes that TEMPD is connected to that are not numbered.
            DO its = 1 , Nnod
 
               DO i = istart , nqtemp
                  row = Tempd(i)
!                 write(*,*)'fina(row),fina(row+1):',fina(row),fina(row+1)
                  DO count = Fina(row) , Fina(row+1) - 1
                     icol = Cola(count)
! see if already in TEMPD.
                     IF ( Whichd(icol).EQ.0 ) THEN
!                       write(*,*)'inside'
                        nq = nq + 1
                        Tempd(nq) = icol
!
!                       !print*, 'nq :', nq, ' tempd :', tempd(nq)
!
                        Whichd(icol) = idom
                     ENDIF
                  ENDDO
!                 write(*,*)'nq=',nq
                  IF ( nq.GT.(Nnod/Ndom)+1 ) GOTO 50
               ENDDO
!
               !print * , 'nq :' , nq , ' nqtemp :' , nqtemp
!
               IF ( nq.EQ.nqtemp ) GOTO 20
               istart = nqtemp + 1
               nqtemp = nq
!
               !print * , 'istart :' , istart , ' nqtemp :' , nqtemp
!
               !print * , 'its :' , its , ' nq :' , nq , ' nqtemp :' ,   &
               !    & nqtemp
            ENDDO
!
!
 50      ENDDO
         DO i = 1 , Nnod
            IF ( Whichd(i).EQ.0 ) THEN
               !WRITE (*,*) 'zero at i=' , i
               Whichd(i) = 1
            ENDIF
            Tempd(i) = Whichd(i)
         ENDDO
      ENDIF
      END
!*==NEURI.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
 
!
!
!
!
!
      SUBROUTINE NEURI(Ndom2,Nnod2,V,Tscale,Fina, &
                     & Cola,Midpa,A,Ncola,Nnod,Ndom,Whichd,   &
                     & Noinsu,  &
                     & Stoexp,  &
                     & Sumv,Optim)
      IMPLICIT NONE
!*--NEURI1088
!*** Start of declarations inserted by SPAG
      REAL fcost , ffb , rcost , rr , rscal
      INTEGER icol , iwicol , maxsub , Ncola , Ndom2 , Nnod , Nnod2
!*** End of declarations inserted by SPAG
!
      REAL ANNEAL , FORCE
      INTEGER Ndom , CHAINL , NOITS
      PARAMETER (CHAINL=300)
!
      PARAMETER (ANNEAL=0.7)
      PARAMETER (NOITS=30)
      PARAMETER (FORCE=0.45)
!
      INTEGER Midpa(Nnod) , Cola(Ncola) , Fina(Nnod+1)
      INTEGER Whichd(Nnod)
      INTEGER Optim
      REAL A(Ncola)
!
! THETA=0 means gamma x^T x is used =1 means x^T b is used.
! FORCE is the force with wich the interface nodes repel the other nodes.
! TOLER is the tolerence of the iteration.
! TAKAWA is the force with which the interface nodes are attracted
! to each other.
      REAL kk , maxv , minwei
      REAL Sumv(Ndom+1) , sum
      REAL ttvv , f , ff , Stoexp(Ndom+1)
      REAL V(Nnod2*(Ndom2+1)) , Tscale(Nnod2)
!        REAL V(maxnod*(maxDOM+1)),TSCALE(maxNOD)
      REAL prodm2 , scanu1 , scanu2
!
      REAL alpha , gamma
      INTEGER sub , sub2 , nodsub , i , its , colaco , Noinsu(Ndom+1)
!
      INTEGER count , chain
!
      LOGICAL truble , found
!
! Initialise.
! Work out TSCALE for each node.
      rscal = MIN(0.5,0.05+0.5*REAL(Optim)/200.)
      DO i = 1 , Nnod
         kk = 0.
         DO count = Fina(i) , Fina(i+1) - 1
            kk = kk + A(count)
         ENDDO
!          TSCALE(I)=0.2*KK/REAL(NDOM+1)
!          TSCALE(I)=1.0*KK/REAL(NDOM+1)
!          TSCALE(I)=0.5*KK/REAL(NDOM+1)
!          TSCALE(I)=0.05*KK/REAL(NDOM+1)
         Tscale(i) = rscal*kk/REAL(Ndom+1)
         rr = REAL((Ndom+1)*(Ndom+1)*Nnod)
!        TSCALE(I)=KK*REAL(NDOM-1)/RR
!          TSCALE(I)=KK/RR
      ENDDO
!
      kk = 0.
      minwei = 1.
      DO count = 1 , Fina(Nnod+1) - 1
         kk = kk + A(count)
      ENDDO
!        ALPHA=0.1*REAL(NDOM)*KK/REAL(NNOD**2)
      alpha = 0.5*REAL(Ndom)*kk/REAL(Nnod**2)
      gamma = FORCE*minwei
      scanu1 = alpha*REAL(Ndom-1)/REAL(Ndom)
      scanu2 = -alpha/REAL(Ndom)
!
      DO i = 1 , Nnod*(Ndom+1)
         V(i) = 0.
      ENDDO
      DO i = 1 , Nnod
         V((Whichd(i)-1)*Nnod+i) = 1.
      ENDDO
 
!
      !WRITE (*,*) 'alpha=' , alpha , ' gamma=' , gamma
! The iteration is:
      ffb = 0.
      rcost = 0.
      DO its = 1 , NOITS
         IF ( its.NE.1 ) THEN
            DO i = 1 , Nnod
               Tscale(i) = ANNEAL*Tscale(i)
            ENDDO
         ENDIF
         DO chain = 1 , CHAINL
!
! FF=the functional.
! work out SUMV(SUB)
            DO sub = 1 , Ndom
               sum = 0.
               DO i = 1 , Nnod
                  sum = sum + V((sub-1)*Nnod+i)
               ENDDO
               Sumv(sub) = sum
            ENDDO
!
            fcost = 0.
!
 
            ff = 0.
            DO i = 1 , Nnod
!
               sum = 0.
               DO sub = 1 , Ndom
                  nodsub = (sub-1)*Nnod + i
                  prodm2 = 0.
                  ttvv = 0.
                  DO sub2 = 1 , Ndom
                     IF ( sub.EQ.sub2 ) THEN
                        prodm2 = prodm2 + scanu1*Sumv(sub2)
                     ELSE
                        prodm2 = prodm2 + scanu2*Sumv(sub2)
!
                        DO count = Fina(i) , Fina(i+1) - 1
                           colaco = Cola(count) + (sub2-1)*Nnod
                           ttvv = ttvv + A(count)*V(colaco)
                        ENDDO
                     ENDIF
                  ENDDO
!
                  f = ttvv + prodm2
                  ff = ff + V(nodsub)*f
                  IF ( ABS(f)/Tscale(i).LT.40 ) THEN
                     Stoexp(sub) = EXP(-f/Tscale(i))
                  ELSE
                     Stoexp(sub) = 0.
                  ENDIF
                  sum = sum + Stoexp(sub)
               ENDDO
 
!              F=GAMMA*V(NDOM*NNOD+I)
               f = gamma
               IF ( ABS(f)/Tscale(i).LT.40 ) THEN
                  Stoexp(Ndom+1) = EXP(-f/Tscale(i))
               ELSE
                  Stoexp(Ndom+1) = 0.
               ENDIF
               sum = sum + Stoexp(Ndom+1)
               IF ( sum.LT.0.0000000000000001 ) GOTO 100
!
! Update part that attempts to balance the number of nodes.
               DO sub = 1 , Ndom
                  nodsub = (sub-1)*Nnod + i
                  Sumv(sub) = Sumv(sub) - V(nodsub)
                  V(nodsub) = Stoexp(sub)/sum
                  Sumv(sub) = Sumv(sub) + V(nodsub)
               ENDDO
! Update part for interface nodes
               V(Ndom*Nnod+i) = Stoexp(Ndom+1)/sum
            ENDDO
           ! WRITE (*,*) 'FF=' , ff , ' FFB=' , ffb , ' chain=' , chain ,&
           !            &' ITS=' , its
            IF ( ABS(ff-ffb).LT.0.0001 ) GOTO 100
!        IF(ABS(FF-FFB).LT.0.000001) GOTO 9000
            ffb = ff
!
! Test for convergence.
         ENDDO
 100  ENDDO
!
!
!
!
 
! convert V to WHICHD.
      DO i = 1 , Nnod
         maxv = 0.
         DO sub = 1 , Ndom + 1
            IF ( V((sub-1)*Nnod+i).GT.maxv ) THEN
               maxv = V((sub-1)*Nnod+i)
               maxsub = sub
            ENDIF
         ENDDO
         Whichd(i) = maxsub
      ENDDO
!
! check to see if we have a valide partition.
      truble = .FALSE.
      DO i = 1 , Nnod
         DO count = Fina(i) , Fina(i+1) - 1
            iwicol = Whichd(Cola(count))
            IF ( Whichd(i).LE.Ndom ) THEN
!                 write(*,*)'trouble around i=',i
               IF ( (Whichd(i).NE.iwicol) .AND. (iwicol.LE.Ndom) )      &
                  & truble = .TRUE.
            ENDIF
         ENDDO
      ENDDO
     ! IF ( truble ) WRITE (*,*) 'NOT A VALIDE PARTITION'
     ! IF ( .NOT.truble ) WRITE (*,*) 'A VALIDE PARTITION WAS FOUND'
      DO sub = 1 , Ndom + 1
         Noinsu(sub) = 0
      ENDDO
      DO i = 1 , Nnod
         Noinsu(Whichd(i)) = Noinsu(Whichd(i)) + 1
      ENDDO
      DO sub = 1 , Ndom
         !WRITE (*,*) 'NO OF NODS IN SUB=' , sub , ' IS =' , Noinsu(sub)
      ENDDO
     ! WRITE (*,*) 'NO OF INTERFACE NODES =' , Noinsu(Ndom+1)
! see if interface node is surrounded by other nodes.
      DO i = 1 , Nnod
         IF ( Whichd(i).EQ.Ndom+1 ) THEN
!            write(*,*)'v(i)=',v(i),' i=',i
            found = .FALSE.
            DO count = Fina(i) , Fina(i+1) - 1
               icol = Cola(count)
               IF ( icol.NE.i ) THEN
                  IF ( Whichd(icol).EQ.Ndom+1 ) found = .TRUE.
               ENDIF
            ENDDO
           ! IF ( .NOT.found ) WRITE (*,*) 'I=' , i , 'V(I)=' ,          &
           !                        & V(i+Ndom*Nnod)
         ENDIF
      ENDDO
!
! SET THE INTERFACE NODES TO -1 in WHICHD
      DO i = 1 , Nnod
         IF ( Whichd(i).GT.Ndom ) Whichd(i) = -1
      ENDDO
      END
!*==INTPRO.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
 
      SUBROUTINE INTPRO(Fina,Cola,Midpa,A,Ncola,Nnod,Ndom,Whichd,V,X,   &
                      & Subdon)
      IMPLICIT NONE
!*--INTPRO1319
!*** Start of declarations inserted by SPAG
      REAL A
      INTEGER i , Ncola , Ndom , Nnod
!*** End of declarations inserted by SPAG
! This sub works out which subdomain to attach to the interface nodes.
!
      INTEGER count , sub , sub2
      INTEGER V(Nnod) , X(Nnod)
      INTEGER minint , minsub , nintsu
      INTEGER Whichd(Nnod)
      INTEGER Midpa(Nnod) , Cola(Ncola) , Fina(Nnod+1)
      LOGICAL Subdon(Ndom)
      DO sub = 1 , Ndom
         Subdon(sub) = .FALSE.
      ENDDO
! identify which subdomain has the least in number of
! remaining interface nodes surrounding it.
!
      DO sub2 = 1 , Ndom
         minint = 1000000
         minsub = 0
         DO sub = 1 , Ndom
            IF ( .NOT.Subdon(sub) ) THEN
               nintsu = 0
               DO i = 1 , Nnod
                  V(i) = 0
                  X(i) = 0
                  IF ( Whichd(i).EQ.sub ) V(i) = 1
               ENDDO
! form x=Av
               DO i = 1 , Nnod
                  DO count = Fina(i) , Fina(i+1) - 1
                     X(i) = X(i) + 1*V(Cola(count))
                  ENDDO
               ENDDO
               nintsu = 0
               DO i = 1 , Nnod
!            IF((WHICHD(I).GT.NDOM).AND.(X(I).GT.0)) NINTSU=NINTSU+1
                  IF ( (Whichd(i).LT.0) .AND. (X(i).GT.0) )             &
                     & nintsu = nintsu + 1
               ENDDO
               IF ( nintsu.LT.minint ) THEN
                  minsub = sub
                  minint = nintsu
               ENDIF
            ENDIF
         ENDDO
         Subdon(minsub) = .TRUE.
! Claim all surrounding interface nodes of sub MINSUB
         DO i = 1 , Nnod
            V(i) = 0
            X(i) = 0
            IF ( Whichd(i).EQ.minsub ) V(i) = 1
         ENDDO
! form x=Av
         DO i = 1 , Nnod
            DO count = Fina(i) , Fina(i+1) - 1
               X(i) = X(i) + 1*V(Cola(count))
            ENDDO
         ENDDO
         nintsu = 0
         DO i = 1 , Nnod
            IF ( (Whichd(i).LT.0) .AND. (X(i).GT.0) ) Whichd(i)         &
               & = -minsub
         ENDDO
      ENDDO
      END
!*==FINEIG.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
!
!
!
      SUBROUTINE FINEIG(Newsub,Newnod,Oldsub,Nparts,Npart3,Npart4,Ndom2,&
                      & Nnod2,R,R2,Rtr,Rtrpar,Rtrpa3,Rtrpa4,Rtrsub,     &
                      & Prodm2,Fina,Cola,Midpa,A,Ncola,Nnod,Ndom,Whichd,&
                      & Tempd,Power)
      IMPLICIT NONE
!*--FINEIG1400
!*** Start of declarations inserted by SPAG
      INTEGER iarg , irnod
      REAL RAN1 , rr , rran , rsum
!*** End of declarations inserted by SPAG
! Find a neighbourhood configuration
!
      INTEGER Power
      LOGICAL SIMPLE
      PARAMETER (SIMPLE=.FALSE.)
!        PARAMETER(POWER=2)
!        PARAMETER(NODACR=52,NODOWN=52)
!        PARAMETER(NS=MAXNOD**0.75)
!        PARAMETER(N3=MAXNOD**0.5,N4=MAXNOD**0.25)
!
      INTEGER Ncola , Nnod , Ndom
      INTEGER Whichd(Nnod) , Midpa(Nnod) , Cola(Ncola) , Fina(Nnod+1)
      REAL A(Ncola)
!
      INTEGER Tempd(Nnod)
!
      INTEGER Nparts , Npart3 , Npart4 , Ndom2 , Nnod2
      REAL R(Ndom2*Nnod2) , R2(Ndom2*Nnod2)
      REAL Prodm2(Ndom2) , Rtrpar(Ndom2*Nparts)
      REAL Rtrsub(Ndom2)
      REAL Rtrpa3(Ndom2*Npart3) , Rtrpa4(Ndom2*Npart4)
      REAL Rtr
!
      REAL ranm , prpow
      INTEGER count , sub , iran , Newsub , Oldsub , Newnod
!        INTEGER TEMPD(MAXNOD),SUB2,NODSUB
      INTEGER nodsub
      INTEGER part , part3 , part4
      INTEGER ptpart , ptpar3 , ptpar4
      INTEGER start , finish , isub
      LOGICAL usesim , doloop , missed
!        COMMON /TEMPD/ TEMPD
!        COMMON /R/ R,R2,RTR,RTRPAR,RTRPA3,RTRPA4,RTRSUB,PRODM2
!
!        COMMON /NPART/ NPARTS,NPART3,NPART4
!
!        NPARTS=NNOD**0.75
!        NPART3=NNOD**0.5
!        NPART4=NNOD**0.25
!
      usesim = .FALSE.
      IF ( ABS(Rtr).LT.0.00001 ) usesim = .TRUE.
! If the last nabourgood configuration was not
! excepted then use R2.
      IF ( (.NOT.SIMPLE) .AND. (.NOT.usesim) ) THEN
! Find a neighbour hood configuration which only adds a node
! where it is adjasent to at laest one other
! node in the same sub dom.
!
!            write(*,*)'befor FIDERI'
!            CALL FIDERI(R)
!            write(*,*)'after FIDERI'
!
! Go back to here if missed.
 50      missed = .FALSE.
!
         iarg = 1
         rran = RAN1(iarg)
         IF ( rran.GE.1. ) rran = 0.99999
         ranm = rran*Rtr
!
! Choose which part to look at more closely.
! find out which sub its in
         rsum = 0.
         doloop = .TRUE.
         DO sub = 1 , Ndom
            IF ( doloop ) THEN
               prpow = (2.*Prodm2(sub))**Power
               rr = prpow*Rtrsub(sub)
               IF ( (ranm.GE.rsum) .AND. (ranm.LT.rsum+rr) ) THEN
                  isub = sub
                  doloop = .FALSE.
                  rsum = rsum - rr
               ENDIF
               rsum = rsum + rr
            ENDIF
         ENDDO
         IF ( doloop ) missed = .TRUE.
!
         sub = isub
!            PRPOW=(2.*PRODM2(SUB))**POWER
!
!            rRSUM=RSUM
         doloop = .TRUE.
         start = (sub-1)*Npart4 + 1
         finish = sub*Npart4
         DO part4 = start , finish
            IF ( doloop ) THEN
               rr = prpow*Rtrpa4(part4)
               IF ( (ranm.GE.rsum) .AND. (ranm.LT.rsum+rr) ) THEN
                  ptpar4 = part4 - (sub-1)*Npart4
                  doloop = .FALSE.
                  rsum = rsum - rr
               ENDIF
               rsum = rsum + rr
            ENDIF
         ENDDO
!
         IF ( doloop ) missed = .TRUE.
!
!            RSUM=RRSUM
         doloop = .TRUE.
         start = (sub-1)*Npart3 + ((ptpar4-1)*Npart3)/Npart4 + 1
         finish = (sub-1)*Npart3 + (ptpar4*Npart3)/Npart4
         count = 0
         DO part3 = start , finish
            IF ( doloop ) THEN
               rr = prpow*Rtrpa3(part3)
               IF ( (ranm.GE.rsum) .AND. (ranm.LT.rsum+rr) ) THEN
                  ptpar3 = part3 - (sub-1)*Npart3
                  doloop = .FALSE.
                  rsum = rsum - rr
               ENDIF
               rsum = rsum + rr
            ENDIF
         ENDDO
         IF ( doloop ) missed = .TRUE.
!
!            RSUM=RRSUM
         doloop = .TRUE.
         start = (sub-1)*Nparts + ((ptpar3-1)*Nparts)/Npart3 + 1
         finish = (sub-1)*Nparts + (ptpar3*Nparts)/Npart3
         count = 0
         DO part = start , finish
            IF ( doloop ) THEN
               rr = prpow*Rtrpar(part)
               IF ( (ranm.GE.rsum) .AND. (ranm.LT.rsum+rr) ) THEN
                  ptpart = part - (sub-1)*Nparts
                  doloop = .FALSE.
                  rsum = rsum - rr
               ENDIF
               rsum = rsum + rr
            ENDIF
         ENDDO
         IF ( doloop ) missed = .TRUE.
!
!            RSUM=RRSUM
         doloop = .TRUE.
         start = (sub-1)*Nnod + ((ptpart-1)*Nnod)/Nparts + 1
         finish = (sub-1)*Nnod + (ptpart*Nnod)/Nparts
         DO nodsub = start , finish
!             IF(DOLOOP) THEN
            rr = prpow*R2(nodsub)
            IF ( (ranm.GE.rsum) .AND. (ranm.LT.rsum+rr) ) THEN
               irnod = nodsub - (sub-1)*Nnod
               doloop = .FALSE.
            ENDIF
            rsum = rsum + rr
!             ENDIF
         ENDDO
         IF ( doloop ) missed = .TRUE.
! TRY AGAIN IF MISSED.,
         IF ( missed ) THEN
      !      WRITE (*,*) '******************MISSED'
            GOTO 50
         ENDIF
!
! Have now chosen node IRNOD to be in sub IRSUB
!            write(*,*)'irnod=',irnod
         Newnod = irnod
         Newsub = sub
         Oldsub = Whichd(Newnod)
         Tempd(Newnod) = Newsub
!            stop
      ENDIF
      IF ( SIMPLE .OR. usesim ) THEN
!
         iarg = 1
         rran = RAN1(iarg)
         IF ( rran.GE.1 ) rran = 0.999999
         iran = INT(rran*REAL(Nnod*Ndom)+1.)
!
         Newsub = ((iran-1)/Nnod) + 1
         Newnod = iran - (Newsub-1)*Nnod
         Oldsub = Whichd(Newnod)
         Tempd(Newnod) = Newsub
      ENDIF
      END
!*==FIDERC.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
!
      SUBROUTINE FIDERC(Pcost,Nparts,Npart3,Npart4,Ndom2,Nnod2,R,R2,Rtr,&
                      & Rtrpar,Rtrpa3,Rtrpa4,Rtrsub,Prodm2,Fina,Cola,   &
                      & Midpa,A,Ncola,Nnod,Ndom,Whichd,Km,Kmtemp,Chopar,&
                      & Chopa3,Chopa4,Power)
      IMPLICIT NONE
!*--FIDERC1594
!*** Start of declarations inserted by SPAG
      INTEGER i , j
      REAL rtri , rtrp , rtrp3 , rtrp4
!*** End of declarations inserted by SPAG
! This sub finds the derivative of the cost function=R
! and finds R2 & KM and also finds the cost function PCOST.
!
      INTEGER Power
!        PARAMETER(POWER=2)
!
      INTEGER Ncola , Nnod , Ndom
      INTEGER Whichd(Nnod) , Midpa(Nnod) , Cola(Ncola) , Fina(Nnod+1)
      REAL A(Ncola)
!
      REAL Km(Ndom) , Kmtemp(Ndom)
!
      INTEGER Chopar(Nnod) , Chopa3(Nnod) , Chopa4(Nnod)
!
      INTEGER Nparts , Npart3 , Npart4 , Ndom2 , Nnod2
      REAL R(Ndom2*Nnod2) , R2(Ndom2*Nnod2)
      REAL Prodm2(Ndom2) , Rtrpar(Ndom2*Nparts)
      REAL Rtrsub(Ndom2)
      REAL Rtrpa3(Ndom2*Npart3) , Rtrpa4(Ndom2*Npart4)
      REAL Rtr
!
      REAL prodm1
!        real KM(MAXDOM),KMTEMP(MAXDOM),FF,sp,PRODM1
      REAL Pcost
      REAL prpow
      INTEGER sub , sub2
!        INTEGER CHOPAR(MAXNOD),CHOPA3(MAXNOD),CHOPA4(MAXNOD)
      INTEGER count , row , col
      INTEGER ptpart , ptpar3 , ptpar4
      INTEGER start , finish , st , fi , st2 , fi2 , st3 , fi3 , part , &
            & part3 , part4
!        COMMON /KM/ KM,KMTEMP
!        COMMON /R/ R,R2,RTR,RTRPAR,RTRPA3,RTRPA4,RTRSUB,PRODM2
!        COMMON /CHOPAR/ CHOPAR,CHOPA3,CHOPA4
!        COMMON /NPART/ NPARTS,NPART3,NPART4
! This sub finds the FINDS THE DERIVATIVE OF THE FUNCTIONAL.
!
!        NPARTS=NNOD**0.75
!        NPART3=NNOD**0.5
!        NPART4=NNOD**0.25
!
      DO sub = 1 , Ndom
         DO i = 1 , Nnod
            R(i+(sub-1)*Nnod) = 0.
         ENDDO
!
         DO col = 1 , Nnod
            IF ( Whichd(col).EQ.sub ) THEN
               DO count = Fina(col) , Fina(col+1) - 1
                  row = Cola(count)
                  R(row+(sub-1)*Nnod) = R(row+(sub-1)*Nnod) + A(count)
               ENDDO
            ENDIF
         ENDDO
!
         Km(sub) = 0.
         DO i = 1 , Nnod
            R2((sub-1)*Nnod+i) = R((sub-1)*Nnod+i)**Power
            IF ( Whichd(i).EQ.sub ) THEN
               R2((sub-1)*Nnod+i) = 0.
               Km(sub) = Km(sub) + R(i+(sub-1)*Nnod)
            ENDIF
         ENDDO
      ENDDO
!
      DO sub = 1 , Ndom
         prodm1 = 1.
         DO sub2 = 1 , Ndom
            IF ( sub2.NE.sub ) prodm1 = prodm1*Km(sub2)
         ENDDO
         Prodm2(sub) = prodm1
      ENDDO
!
      Pcost = 1.
      DO sub = 1 , Ndom
         Kmtemp(sub) = Km(sub)
         Pcost = Pcost*Km(sub)
      ENDDO
 
!
! WORK OUT RTRPAR & RTR.
! Work out CHOPA3() & CHOPA4()
      Rtr = 0.
      DO sub = 1 , Ndom
         prpow = (2.*Prodm2(sub))**Power
         start = (sub-1)*Npart4 + 1
         finish = sub*Npart4
         rtrp4 = 0.
         DO part4 = start , finish
            ptpar4 = part4 - (sub-1)*Npart4
            st3 = (sub-1)*Npart3 + ((ptpar4-1)*Npart3)/Npart4 + 1
            fi3 = (sub-1)*Npart3 + (ptpar4*Npart3)/Npart4
            rtrp3 = 0.
            DO part3 = st3 , fi3
               ptpar3 = part3 - (sub-1)*Npart3
               st2 = (sub-1)*Nparts + ((ptpar3-1)*Nparts)/Npart3 + 1
               fi2 = (sub-1)*Nparts + (ptpar3*Nparts)/Npart3
               rtrp = 0.
               DO part = st2 , fi2
                  ptpart = part - (sub-1)*Nparts
                  st = (sub-1)*Nnod + ((ptpart-1)*Nnod)/Nparts + 1
                  fi = (sub-1)*Nnod + (ptpart*Nnod)/Nparts
                  rtri = 0.
                  DO i = st , fi
                     rtri = rtri + R2(i)
!
                     IF ( sub.EQ.1 ) THEN
                        j = i - (sub-1)*Nnod
                        Chopar(j) = part - (sub-1)*Nparts
                        Chopa3(j) = part3 - (sub-1)*Npart3
                        Chopa4(j) = part4 - (sub-1)*Npart4
                     ENDIF
                  ENDDO
                  Rtrpar(part) = rtri
                  rtrp = rtrp + rtri
               ENDDO
               Rtrpa3(part3) = rtrp
               rtrp3 = rtrp3 + rtrp
            ENDDO
            Rtrpa4(part4) = rtrp3
            rtrp4 = rtrp4 + rtrp3
         ENDDO
         Rtrsub(sub) = rtrp4
         Rtr = Rtr + prpow*rtrp4
      ENDDO
!
      !WRITE (*,*) 'rtr=' , Rtr
      END
!*==UPDERI.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
!
 
      SUBROUTINE UPDERI(Ncost,Newsub,Newnod,Oldsub,Nparts,Npart3,Npart4,&
                      & Ndom2,Nnod2,R,R2,Rtr,Rtrpar,Rtrpa3,Rtrpa4,      &
                      & Rtrsub,Prodm2,Fina,Cola,Midpa,A,Ncola,Nnod,Ndom,&
                      & Whichd,Km,Kmtemp,Chopar,Chopa3,Chopa4,Power)
      IMPLICIT NONE
!*--UPDERI1739
! This sub updates THE NEW R & finds Rtr & RTRPAR.
!
      REAL SMALL
      INTEGER Power
!        PARAMETER(POWER=2)
      PARAMETER (SMALL=1000.)
!
      INTEGER Ncola , Nnod , Ndom
      INTEGER Whichd(Nnod) , Midpa(Nnod) , Cola(Ncola) , Fina(Nnod+1)
      REAL A(Ncola)
!
      REAL Km(Ndom) , Kmtemp(Ndom)
!
      INTEGER Chopar(Nnod) , Chopa3(Nnod) , Chopa4(Nnod)
!
      INTEGER Nparts , Npart3 , Npart4 , Ndom2 , Nnod2
      REAL R(Ndom2*Nnod2) , R2(Ndom2*Nnod2)
      REAL Prodm2(Ndom2) , Rtrpar(Ndom2*Nparts)
      REAL Rtrsub(Ndom2)
      REAL Rtrpa3(Ndom2*Npart3) , Rtrpa4(Ndom2*Npart4)
      REAL Rtr
!
!        real KM(MAXDOM),KMTEMP(MAXDOM),FF,PRODM1
      REAL prodm1
      REAL prpow , r1pr2p
      INTEGER sub , sub2 , nodsub , colaco
!        INTEGER CHOPAR(MAXNOD),CHOPA3(MAXNOD),CHOPA4(MAXNOD)
      REAL Ncost
      INTEGER Newsub , Newnod , Oldsub
      INTEGER count
      INTEGER ptpart , ptpar3 , ptpar4
!        COMMON /KM/ KM,KMTEMP
!        COMMON /R/ R,R2,RTR,RTRPAR,RTRPA3,RTRPA4,RTRSUB,PRODM2
!        COMMON /CHOPAR/ CHOPAR,CHOPA3,CHOPA4
!        COMMON /NPART/ NPARTS,NPART3,NPART4
! This sub finds the FINDS THE DERIVATIVE OF THE FUNCTIONAL.
!
!        NPARTS=NNOD**0.75
!        NPART3=NNOD**0.5
!        NPART4=NNOD**0.25
!
! WORK OUT R.
      IF ( Ncost.LT.SMALL ) THEN
         DO sub = 1 , Ndom
            prodm1 = 1.
            DO sub2 = 1 , Ndom
               IF ( sub2.NE.sub ) prodm1 = prodm1*Km(sub2)
            ENDDO
            Prodm2(sub) = prodm1
         ENDDO
      ELSE
         DO sub = 1 , Ndom
            Prodm2(sub) = Ncost/Km(sub)
         ENDDO
      ENDIF
!
!
!               prob=.FALSE.
      DO count = Fina(Newnod) , Fina(Newnod+1) - 1
         colaco = Cola(count)
!
         nodsub = (Newsub-1)*Nnod + colaco
         ptpart = Chopar(colaco) + (Newsub-1)*Nparts
         ptpar3 = Chopa3(colaco) + (Newsub-1)*Npart3
         ptpar4 = Chopa4(colaco) + (Newsub-1)*Npart4
! Must update RTRPAR as well
! Work out new R(NOSUB)
         R(nodsub) = R(nodsub) + A(count)
         r1pr2p = -R2(nodsub)
         IF ( Whichd(colaco).NE.Newsub ) THEN
            R2(nodsub) = R(nodsub)**Power
         ELSE
            R2(nodsub) = 0.
         ENDIF
         r1pr2p = r1pr2p + R2(nodsub)
         Rtrsub(Newsub) = Rtrsub(Newsub) + r1pr2p
         Rtrpar(ptpart) = Rtrpar(ptpart) + r1pr2p
         Rtrpa3(ptpar3) = Rtrpa3(ptpar3) + r1pr2p
         Rtrpa4(ptpar4) = Rtrpa4(ptpar4) + r1pr2p
!
! takeaway a 1.
         nodsub = (Oldsub-1)*Nnod + colaco
         ptpart = Chopar(colaco) + (Oldsub-1)*Nparts
         ptpar3 = Chopa3(colaco) + (Oldsub-1)*Npart3
         ptpar4 = Chopa4(colaco) + (Oldsub-1)*Npart4
! Must update RTRPAR as well
! Work out new R(NOSUB)
         R(nodsub) = R(nodsub) - A(count)
         r1pr2p = -R2(nodsub)
         IF ( Whichd(colaco).NE.Oldsub ) THEN
            R2(nodsub) = R(nodsub)**Power
         ELSE
            R2(nodsub) = 0.
         ENDIF
         r1pr2p = r1pr2p + R2(nodsub)
         Rtrsub(Oldsub) = Rtrsub(Oldsub) + r1pr2p
         Rtrpar(ptpart) = Rtrpar(ptpart) + r1pr2p
         Rtrpa3(ptpar3) = Rtrpa3(ptpar3) + r1pr2p
         Rtrpa4(ptpar4) = Rtrpa4(ptpar4) + r1pr2p
      ENDDO
!
!
! Must update RTR as well
      Rtr = 0.
      DO sub = 1 , Ndom
         prpow = (2.*Prodm2(sub))**Power
         Rtr = Rtr + prpow*Rtrsub(sub)
!                 RTR=RTR+ RTRSUB(SUB)*(PRODM2(SUB)**POWER)
      ENDDO
!              RTR=4.*RTR
!
      END
!*==FININI.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
 
      SUBROUTINE FININI(Nchain,Initic,Pcost,Gotc,Nparts,Npart3,Npart4,  &
                      & Ndom2,Nnod2,R,R2,Rtr,Rtrpar,Rtrpa3,Rtrpa4,      &
                      & Rtrsub,Prodm2,Fina,Cola,Midpa,A,Ncola,Nnod,Ndom,&
                      & Whichd,Tempd,Km,Kmtemp,Random,Power,Chopar,     &
                      & Chopa3,Chopa4)
      IMPLICIT NONE
!*--FININI1864
!*** Start of declarations inserted by SPAG
      INTEGER its
!*** End of declarations inserted by SPAG
!
! This sub finds the initial C parameter INITIC.
!
      REAL EXCEPI
      LOGICAL SIMPLE
      PARAMETER (SIMPLE=.FALSE.)
      PARAMETER (EXCEPI=0.9)
!
      INTEGER Ncola , Nnod , Ndom
      INTEGER Whichd(Nnod) , Midpa(Nnod) , Cola(Ncola) , Fina(Nnod+1)
      REAL A(Ncola)
!
      REAL Km(Ndom) , Kmtemp(Ndom)
      INTEGER Tempd(Nnod)
      INTEGER Power
!
      INTEGER Chopar(Nnod) , Chopa3(Nnod) , Chopa4(Nnod)
!
      INTEGER Nparts , Npart3 , Npart4 , Ndom2 , Nnod2
      REAL R(Ndom2*Nnod2) , R2(Ndom2*Nnod2)
      REAL Prodm2(Ndom2) , Rtrpar(Ndom2*Nparts)
      REAL Rtrsub(Ndom2)
      REAL Rtrpa3(Ndom2*Npart3) , Rtrpa4(Ndom2*Npart4)
      REAL Rtr
!
! EXCEPI=initial probability of exceptance with decreased cost fu.
! if MAXBAS then find INITIC based on max difference in cost.
      REAL Initic , Pcost , dif , nloacp
      REAL mxdifc , ncost
!        real KM(MAXDOM),KMTEMP(MAXDOM)
!        INTEGER TEMPD(MAXNOD),I,IRAN,SAM
      INTEGER Nchain , chain , newsub , newnod , oldsub , sub
      LOGICAL Gotc
      LOGICAL Random
!        COMMON /TEMPD/ TEMPD
!        COMMON /KM/ KM,KMTEMP
      mxdifc = 0.
! work out R & R2.
      Ncola = Fina(Nnod+1) - 1
      CALL FIDERC(Pcost,Nparts,Npart3,Npart4,Ndom2,Nnod2,R,R2,Rtr,      &
                & Rtrpar,Rtrpa3,Rtrpa4,Rtrsub,Prodm2,Fina,Cola,Midpa,A, &
                & Ncola,Nnod,Ndom,Whichd,Km,Kmtemp,Chopar,Chopa3,Chopa4,&
                & Power)
      !WRITE (*,*) 'PCOST=' , Pcost
! Find an initial INITIC s.t. [1/(e^( {ncost-pcost}/INITIC))] =0.98
!        DO 20 ITS=1,70
      DO its = 1 , 7
         DO chain = 1 , Nchain
! Find a neighbourhood configuration
            CALL FINEIG(newsub,newnod,oldsub,Nparts,Npart3,Npart4,Ndom2,&
                      & Nnod2,R,R2,Rtr,Rtrpar,Rtrpa3,Rtrpa4,Rtrsub,     &
                      & Prodm2,Fina,Cola,Midpa,A,Ncola,Nnod,Ndom,Whichd,&
                      & Tempd,Power)
            IF ( newsub.NE.oldsub ) THEN
! find cost of neighbourhood config
!              CALL FICOSK(NCOST,PCOST,NEWSUB,OLDSUB,NEWNOD)
               CALL FICOSK(ncost,Pcost,newsub,oldsub,newnod,Fina,Cola,  &
                         & Midpa,A,Ncola,Nnod,Ndom,Whichd,Tempd,Km,     &
                         & Kmtemp)
               dif = 0.
               IF ( ncost.LT.Pcost ) dif = ABS(Pcost-ncost)
               IF ( dif.GT.mxdifc ) mxdifc = dif
               Pcost = ncost
               Whichd(newnod) = newsub
               DO sub = 1 , Ndom
                  Km(sub) = Kmtemp(sub)
               ENDDO
               IF ( .NOT.SIMPLE )                                       &
                  & CALL UPDERI(ncost,newsub,newnod,oldsub,Nparts,      &
                  & Npart3,Npart4,Ndom2,Nnod2,R,R2,Rtr,Rtrpar,Rtrpa3,   &
                  & Rtrpa4,Rtrsub,Prodm2,Fina,Cola,Midpa,A,Ncola,Nnod,  &
                  & Ndom,Whichd,Km,Kmtemp,Chopar,Chopa3,Chopa4,Power)
            ENDIF
!            write(*,*)'pcost=',pcost
         ENDDO
      ENDDO
!
      nloacp = LOG(EXCEPI)
!        INITIC= - MXDIFC/NLOACP
      Initic = Pcost
      !WRITE (*,*) 'INITIC=' , Initic
      Gotc = .TRUE.
      IF ( ABS(Initic).LT.0.00001 ) Gotc = .FALSE.
!         STOP
      END
!*==FICOSK.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
 
!
!
!
!
!
      SUBROUTINE FICOSK(Ncost,Pcost,Newsub,Oldsub,Newnod,Fina,Cola,     &
                      & Midpa,A,Ncola,Nnod,Ndom,Whichd,Tempd,Km,Kmtemp)
      IMPLICIT NONE
!*--FICOSK1963
! find the cost function associated with configuration TEMPD.
! using KM from previous configuration.
!
      INTEGER Ncola , Nnod , Ndom
      INTEGER Whichd(Nnod) , Midpa(Nnod) , Cola(Ncola) , Fina(Nnod+1)
      REAL A(Ncola)
!
      REAL Km(Ndom) , Kmtemp(Ndom)
      INTEGER Tempd(Nnod)
!
      REAL Ncost , Pcost
!           real KM(MAXDOM),KMTEMP(MAXDOM)
      REAL moxro1 , moxro2
!           INTEGER TEMPD(MAXNOD),SUB
      INTEGER sub
      INTEGER Newsub , Oldsub , Newnod
!
      INTEGER count
!           COMMON /TEMPD/ TEMPD
!
!           COMMON /KM/ KM,KMTEMP
!
! Find NCOST from KM and PCOST and NEWNOD,NEWSUB
!            IF(PCOST.LT.1.0E+08) THEN
      DO sub = 1 , Ndom
         Kmtemp(sub) = Km(sub)
      ENDDO
!            ENDIF
!
 
      moxro1 = 0.
      moxro2 = 0.
      DO count = Fina(Newnod) , Fina(Newnod+1) - 1
         IF ( Whichd(Cola(count)).EQ.Newsub ) moxro1 = moxro1 + A(count)
         IF ( Whichd(Cola(count)).EQ.Oldsub ) moxro2 = moxro2 + A(count)
      ENDDO
! NB A(MIDPA(COUNT))=0.
!
      Kmtemp(Newsub) = Km(Newsub) + 2.*moxro1
      Kmtemp(Oldsub) = Km(Oldsub) - 2.*moxro2
 
!
      Ncost = Kmtemp(1)
      DO sub = 2 , Ndom
         Ncost = Ncost*Kmtemp(sub)
      ENDDO
      END
!*==MESH.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
!
!
!
      SUBROUTINE MESH(Beta,Wnod,Fina,Cola,Midpa,A,Ncola,Nnod)
      IMPLICIT NONE
!*--MESH2021
!*** Start of declarations inserted by SPAG
      INTEGER i , Ncola , Nnod
!*** End of declarations inserted by SPAG
!
      REAL Beta , Wnod(Nnod)
      INTEGER count
      INTEGER Midpa(Nnod) , Cola(Ncola) , Fina(Nnod+1)
      REAL A(Ncola)
!
      CALL CLRRL(Nnod,Wnod,1.0)
!
      DO i = 1 , Nnod
         DO count = Fina(i) , Fina(i+1) - 1
            A(count) = Beta
         ENDDO
         A(Midpa(i)) = 0.
      ENDDO
! Now scale the matrix.
      DO i = 1 , Nnod
         DO count = Fina(i) , Fina(i+1) - 1
            A(count) = A(count)/(Wnod(Cola(count))*Wnod(i))
         ENDDO
      ENDDO
      END SUBROUTINE MESH
!*==USEDIV.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!      idum = 1
!      do 10 i = 1, 30
!        z = ran1(idum)
!10    write(*,*) idum, z, z*1000000
!      end
!
!
!
!
! ******************************************************
! ******************************************************
! ****** DECOMPOSE USEING MULTI-LEVEL NEURAL NET *******
! ******************************************************
! ******************************************************
! ******************************************************
 
 
 
!
!
      SUBROUTINE USEDIV(Rmem,Nrmem,Fina,Cola,A,Ncola,Nonods,Whichd,     &
                      & Renum,Imem,Nimem,X,Y,Ndom,Havmat,Splevs,Nsplt,  &
                      & Mxnspl,Optim, havwnod)
      IMPLICIT NONE
!*--USEDIV2074
!*** Start of declarations inserted by SPAG
      INTEGER i , ichoic , Ncola , Ndom , Nimem , Nonods , Nrmem ,      &
            & Nsplt , nsubal, havwnod
!*** End of declarations inserted by SPAG
! This sub allows the user to input information about
      INTEGER Havmat , Optim
      REAL Rmem(Nrmem) , A(Ncola*Havmat)
      REAL X(Nonods) , Y(Nonods)
      INTEGER Fina(Nonods+1) , Cola(Ncola) , Whichd(Nonods)
      INTEGER Renum(Nonods) , Imem(Nimem)
      INTEGER Mxnspl , MXNTRE
!      PARAMETER(MXNSPL=40,MXNTRE=1000)
      PARAMETER (MXNTRE=1000)
! the partition wanted.
      INTEGER Splevs(Mxnspl) , mullev
      REAL subala(MXNTRE)
      REAL alpha , lodbal , beta , toler
      INTEGER nchain
      LOGICAL defaul , topcon , sidcon
      LOGICAL exact
      INTEGER hypdim , proacr , prodow
! MULLEV=no of multi-grid levels.
! SPLEVS(1..NSPLT) contains the dissection information.
      !WRITE (*,*) 'INSIDE USEDIV'
      ichoic = 1
      mullev = 4
!        MULLEV=1
!
      GOTO 100
!
      !WRITE (*,*) '**********************************************'
      !WRITE (*,*) '+ THIS PROGRAM WILL NOW SPLIT THE MESH'
      !WRITE (*,*) '+ RORDER THE EQUATIONS'
      !WRITE (*,*) '+ MAP THE SUBDOMAINS TO A PARALLEL COMPUTER'
      !WRITE (*,*) '  ARCHITECTURE'
      !WRITE (*,*) '*** BUT IT NEEDS SOME INFORMATION TO DO THIS *'
      !WRITE (*,*) '**********************************************'
      !WRITE (*,*) ' '
      !WRITE (*,*) 'ENTER THE NUMBER OF RECURSIVE GRAPH CUTS:'
    !  READ * , Nsplt
      DO i = 1 , Nsplt
         !WRITE (*,*) 'FOR RECURSION =' , i , ' ENTER THE '
         !WRITE (*,*) 'NUMBER OF PARTITIONS:'
        ! READ * , Splevs(i)
         IF ( i.EQ.1 ) THEN
            Ndom = Splevs(i)
         ELSE
            Ndom = Ndom*Splevs(i)
         ENDIF
      ENDDO
!
      !WRITE (*,*) ' '
      !WRITE (*,*) 'FOR THE GRAPH CUT OPTIMIZATION :'
      !WRITE (*,*) 'PARTITIONING THE GRID ON SMALLER GRIDS'
      !WRITE (*,*) 'CAN HELP THE OPTIMIZATION'
      !WRITE (*,*) 'ENTER NUMBER OF MULTI-GRID LEVELS 1,..,5'
     ! READ * , mullev
     ! WRITE (*,*) 'USE DEFAULT  oprimization PARAMETERS(1=YES'
     ! WRITE (*,*) 'OR 0 OTHERWISE'
     ! READ * , ichoic
!
!
 100  DO i = 1 , Nsplt
!
         IF ( i.EQ.1 ) THEN
            Ndom = Splevs(i)
         ELSE
            Ndom = Ndom*Splevs(i)
         ENDIF
!
      ENDDO
!
      IF ( ichoic.EQ.1 ) THEN
         lodbal = 1.
         beta = 0.9
!          BETA=0.8
         toler = 0.0001
!          NCHAIN=1000
         nchain = 700
      ELSE
         !WRITE (*,*) 'LODBAL=1. is the default .gt.1 then more'
         !WRITE (*,*) 'importance is placed on load balancing.'
         !WRITE (*,*) 'BETA=0.9 is the default but this controles'
         !WRITE (*,*) 'how close to the critical temp the optimization'
         !WRITE (*,*) 'is performed WHEN BETA=0.5 - HALF CRYTICAL TEMP'
         !WRITE (*,*) 'IS USED.'
         !WRITE (*,*) 'TOLER =the solution tolerence 0.0001is suggested'
         !WRITE (*,*) 'NCHAIN =maximum no of its 3000 is suggested'
         !WRITE (*,*) 'ENTER LODBAL:'
        ! READ * , lodbal
         !WRITE (*,*) 'ENTER BETA:'
       !  READ * , beta
         !WRITE (*,*) 'ENTER TOLER:'
        ! READ * , toler
         !WRITE (*,*) 'ENTER NCHAIN:'
         !READ * , nchain
      ENDIF
! ALPHA is not eventually used.
      alpha = 1.
!
      !WRITE (*,*) 'ONCE WE HAVE A GRAPH PARTITION WE MUST'
      !WRITE (*,*) 'MAP IT ONTO PROCESSORS CONNECTED'
      !WRITE (*,*) 'WITH A PARTICULAR TAPOLOGY'
      !WRITE (*,*) ' '
      !WRITE (*,*) 'THE DEFAULT MAPPING WORKS WELL IN SOME CASES'
      !WRITE (*,*) 'DO YOU WANT A DEFAULT MAPPING: 1, O-OTHERWISE'
!        READ *,ICHOIC
      ichoic = 1
      IF ( ichoic.EQ.1 ) THEN
         defaul = .TRUE.
      ELSE
         defaul = .FALSE.
         !WRITE (*,*) '2 ARCHITECTURES ARE AVAILABLE MESH CONNECTED'
         !WRITE (*,*) 'SIMILAR TO A N X M F.D grid or a hypercube.'
!
         !WRITE (*,*) ' '
         !WRITE (*,*) 'HYPERCUBE DIMENSION S.T. =2^(no of domains)'
         !WRITE (*,*) 'ENTER THE HYPERCUBE DIMENSIONS: 0 FOR MESH'
       !  READ * , hypdim
         IF ( hypdim.EQ.0 ) THEN
            !WRITE (*,*) ' '
            !WRITE (*,*) 'ENTER THE PARTICULARS OF THE MESH ARCHITECTURE'
            !WRITE (*,*) 'SUPPOSING THE MESH HAS DIM: PRODOW X PROACR'
            !WRITE (*,*) 'THE PROCESSORS WILL BE NUMBERED'
            !WRITE (*,*) '1,2,...,PROACR'
            !WRITE (*,*) '.......,2*PROACR'
            !WRITE (*,*) '.................'
            !WRITE (*,*) '.......,PRODOW*PROACR'
            !WRITE (*,*) 'ENTER PROACR'
           ! READ * , proacr
            !WRITE (*,*) 'ENTER PRODOW'
            !READ * , prodow
            !WRITE (*,*) ' '
            !WRITE (*,*) 'IS THE MESH WRAPPED ROUND FROM THE TOP(1=YES'
       !     READ * , ichoic
            topcon = .FALSE.
            IF ( ichoic.EQ.1 ) topcon = .TRUE.
            !WRITE (*,*) 'MESH WRAPPED ROUND FROM THE SIDES ?(1=YES'
       !     READ * , ichoic
            sidcon = .FALSE.
            IF ( ichoic.EQ.1 ) sidcon = .TRUE.
         ENDIF
      ENDIF
!
      nsubal = MXNTRE
! For reordering equations....
! For the subdomain to proc mapping.......
      CALL PREPAR(Splevs,Nsplt,mullev,Rmem,Nrmem,X,Y,Fina,Cola,A,Ncola, &
                & Nonods,Ndom,Whichd,subala,nsubal,exact,alpha,lodbal,  &
                & beta,toler,nchain,Renum,defaul,hypdim,topcon,sidcon,  &
                & proacr,prodow,Imem,Nimem,Havmat,havwnod)
 
      END
!*==PREPAR.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
! For reordering equations....
! For the subdomain to proc mapping.......
      SUBROUTINE PREPAR(Splevs,Nsplt,Mullev,Rmem,Nrmem,X,Y,Fina,Cola,A, &
                      & Ncola,Nonods,Ndom,Whichd,Subala,Nsubal,Exact,   &
                      & Alpha,Lodbal,Beta,Toler,Nchain,Renum,Defaul,    &
                      & Hypdim,Topcon,Sidcon,Proacr,Prodow,Imem,Nimem,  &
                      & Havmat,havwnod)
      IMPLICIT NONE
!*--PREPAR2241
!*** Start of declarations inserted by SPAG
      INTEGER i , ii , imem2 , inod1 , inod2 , isub , Nchain , Nimem ,  &
            & nimem2 , Nrmem , nrmem2 , Nsubal
!*** End of declarations inserted by SPAG
!
! THIS SUB DIVIDS UP THE DOMAIN INTO SUBDOMAINS
! & ALLOCATES SUBDOMAINS TO PROCESSORS
! & MINIMIZES THE BANDWIDTH ALSO INCREASES CONNECTIVITY
! FOR AN ITERATIVE SOLUTION.
! AFTER CALLING THIS ROUTINE SUBDOMAIN J SHOULD BE
! ATTACHED TO PROC J
! %%%%%%%%%% WHICHD(OLD NOD)=SUBDOMAIN OF %%%%%%%%%%%%%
! %%%%%%%%%% RENUM(OLD NODE)=NEW NODE NO  %%%%%%%%%%%%%
      INTEGER Havmat,havwnod
      INTEGER maxnla , mmxtot , Mullev
      INTEGER Splevs(Nsplt)
      INTEGER Nsplt
      REAL Subala(Nsubal)
      REAL Alpha , Lodbal , Beta , Toler
      REAL X(Nonods) , Y(Nonods)
! MXNTRE=max no of entries in TREE (suggest 1000
! so that max no of sub-domains is about 500).
!
      INTEGER lcola , lfina
!
      REAL Rmem(Nrmem)
!
      REAL A(Ncola*Havmat)
      LOGICAL Exact
! This is for REORGL.....
      INTEGER Ndom , Nonods , nnod , Ncola , mxcols
      INTEGER Fina(Nonods+1) , Cola(Ncola)
      INTEGER Whichd(Nonods) , Renum(Nonods)
! If DEFAUL then map using rordering routine.
! COLSUB contains the subdomain to subdomain communication matrix.
! MXCOLS= maximum value of NCOLSU
!
      INTEGER ptimem
! For the subdomain to proc mapping.......
      LOGICAL Defaul , Topcon , Sidcon
      INTEGER Hypdim , Proacr , Prodow , Imem(Nimem)
      INTEGER lwichd , map , finsub , colsub , blasub
      INTEGER q , qtemp
      INTEGER ptrmem , la , wnod , lwnod , rmem2
      INTEGER lx , ly , maxnq
!
      nnod = Nonods
!
! SPLIT THE DOMAIN INTO SUBDOMAINS
!         MMXTOT=2*NNOD
!         MAXNLA=2*NCOLA
      mmxtot = MAX(1000,INT(1.5*nnod))
      maxnla = MAX(1000,INT(1.5*Ncola))
      mxcols = Ndom*Ndom
! INTEGERS....
      ptimem = 1
      lwichd = ptimem
      ptimem = ptimem + mmxtot
      map = ptimem
      ptimem = ptimem + Nonods
      lfina = ptimem
      ptimem = ptimem + mmxtot + Mullev
      lcola = ptimem
      ptimem = ptimem + maxnla
      q = ptimem
      ptimem = ptimem + Nonods
! give MAP & QTEMP same storage space.
      qtemp = map
 
! REALS....
      ptrmem = 1
      la = ptrmem
      ptrmem = ptrmem + maxnla*Havmat
      wnod = ptrmem
      ptrmem = ptrmem + Nonods
      lwnod = ptrmem
      ptrmem = ptrmem + mmxtot
      lx = ptrmem
      ptrmem = ptrmem + mmxtot
      ly = ptrmem
      ptrmem = ptrmem + mmxtot
      rmem2 = ptrmem
      nrmem2 = Nrmem - rmem2 + 1
!
      CALL RECBIS(Splevs,Nsplt,Mullev,Rmem(rmem2),nrmem2,Fina,Cola,A,   &
                & Ncola,Nonods,mmxtot,Ndom,Whichd,Imem(lfina),          &
                & Imem(lcola),Rmem(la),maxnla,Nsubal,Imem(lwichd),      &
                & Imem(map),Rmem(wnod),Rmem(lwnod),Imem(q),Imem(qtemp), &
                & Subala,Exact,Alpha,Lodbal,Beta,Toler,Nchain,X,Y,      &
                & Rmem(lx),Rmem(ly),Havmat, havwnod)
! count no of nodes in each subdomain
      !WRITE (*,*) 'after RECBIS ndom,nonods:' , Ndom , Nonods
      DO isub = 1 , Ndom
         ii = 0
         DO i = 1 , Nonods
            IF ( Whichd(i).EQ.isub ) ii = ii + 1
         ENDDO
         !WRITE (*,*) 'no of nodes in sub=' , isub , ' is =' , ii
      ENDDO
!
!
      ptimem = 1
      inod1 = ptimem
      maxnq = Nonods + 2000
      ptimem = ptimem + maxnq
      inod2 = ptimem
      ptimem = ptimem + Nonods
!
      finsub = ptimem
      ptimem = ptimem + Ndom + 1
      colsub = ptimem
      ptimem = ptimem + mxcols
      blasub = ptimem
      ptimem = ptimem + Ndom
      imem2 = ptimem
      nimem2 = Nimem - imem2 + 1
! REORDER EQUATIONS -INCLUDING THE SUBDOMAIN TO PROCESSOR MAPPING.
! RENUM(OLD NODE)=NEW NODE NO.
!     &           Q,BLANK,BLASUB,
! For the subdomain to proc mapping.......
      CALL REORG2(Fina,Cola,Ncola,Nonods,Ndom,Whichd,Renum,Imem(finsub),&
                & Imem(colsub),mxcols,Imem(inod1),Imem(inod2),          &
                & Imem(blasub),maxnq,Defaul,Hypdim,Topcon,Sidcon,Proacr,&
                & Prodow,Imem(imem2),nimem2)
      !WRITE (*,*) 'OLDNOD,RENUM(OLDNOD):' , 1 , Renum(1)
      !WRITE (*,*) 'NDOM,DEFAUL' , Ndom , Defaul
      END
!*==FORMA.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
! **************END OF MAIN ********************************
!
!
!
!
      SUBROUTINE FORMA(Cgrid,Maxtot,Totnod,Nonods,Nodlev,Ptcola,Nlevel, &
                     & Fina,Cola,A,Maxna,Na,Color,Q,Qtemp,List,X,Y,     &
                     & Havmat)
      IMPLICIT NONE
!*--FORMA2379
!*** Start of declarations inserted by SPAG
      INTEGER inca , Maxna , maxna2 , ncolo , ncolor , nnod , nod
!*** End of declarations inserted by SPAG
! COLOR & CGRID can contain the same storage space.
! Q & LIST can contain the same storage space.
! This subroutine forms the matrix A from the
! finnest grid information contained in the matrix A.
! CGRID(Fine grid node)=0 if it is not the cause grid node
! = corresponding cause grid node number, otherwise,
! this number varies between 1 and NSHORT.
! NONODS= number of nodes on the finest grid.
      INTEGER Havmat
      INTEGER ONE
      PARAMETER (ONE=1)
      INTEGER Nonods
      INTEGER Totnod , Na , Nlevel , Maxtot
      INTEGER Cgrid(Maxtot)
      INTEGER Nodlev(Nlevel+1) , Ptcola(Nlevel+1)
      INTEGER Fina(Maxtot+Nlevel) , Cola(Maxna)
      INTEGER Color(Maxtot) , Q(Nonods) , Qtemp(Nonods) , List(Nonods)
 
      REAL A(Maxna*Havmat)
      REAL X(Maxtot) , Y(Maxtot)
!
      INTEGER level
      INTEGER fnods , cnods , nca , mxnca , fncola
      LOGICAL zero
      LOGICAL quick
!
      !WRITE (*,*) 'maxtot=' , Maxtot
      zero = .TRUE.
      quick = .TRUE.
!
      Nodlev(1) = 1
      Nodlev(2) = Nonods + 1
      Ptcola(1) = 1
      Ptcola(2) = Fina(Nonods+1)
      Na = Ptcola(2) - 1
      DO level = 2 , Nlevel
! Colour nodes on level LEVEL-1
         nnod = Nodlev(level) - Nodlev(level-1)
         maxna2 = Ptcola(level) - Ptcola(level-1)
!
!          IF(LEVEL.EQ.-level) THEN
!            CALL DRAXYZ(X(NODLEV(LEVEL-1)),Y(NODLEV(LEVEL-1)),
!     &              FINA(NODLEV(LEVEL-1)+LEVEL-2),
!     &                 COLA(PTCOLA(LEVEL-1)),MAXNA2,NNOD)
!           ENDIF
!
         CALL COLOR2(Color(Nodlev(level-1)),ncolor,Q,Qtemp,nnod,        &
                   & Fina(Nodlev(level-1)+level-2),Cola(Ptcola(level-1))&
                   & ,maxna2,zero,quick)
!
!
! Work out CGRID & NCOLO & NODLEV(LEVEL+1).
         ncolo = 0
         DO nod = Nodlev(level-1) , Nodlev(level) - 1
            IF ( Color(nod).EQ.ONE ) THEN
               ncolo = ncolo + 1
               Cgrid(nod) = ncolo
            ELSE
               Cgrid(nod) = 0
            ENDIF
         ENDDO
         Nodlev(level+1) = Nodlev(level) + ncolo
!
! Work out A for level LEVEL.
         fnods = Nodlev(level) - Nodlev(level-1)
         cnods = Nodlev(level+1) - Nodlev(level)
         fncola = Ptcola(level) - Ptcola(level-1)
         mxnca = Maxna - Ptcola(level) + 1
         !WRITE (*,*) '999777 LEVEL=' , level
         IF ( Havmat.EQ.1 ) THEN
            CALL LOCALA(A(Ptcola(level)),A(Ptcola(level-1)),fnods,cnods,&
                      & Cgrid(Nodlev(level-1)),Cola(Ptcola(level)),     &
                      & Cola(Ptcola(level-1)),                          &
                      & Fina(Nodlev(level)+level-1),                    &
                      & Fina(Nodlev(level-1)+level-2),nca,mxnca,fncola, &
                      & List,Nonods,X(Nodlev(level-1)),                 &
                      & Y(Nodlev(level-1)),X(Nodlev(level)),            &
                      & Y(Nodlev(level)),Havmat)
         ELSE
            CALL LOCALA(A,A,fnods,cnods,Cgrid(Nodlev(level-1)),         &
                      & Cola(Ptcola(level)),Cola(Ptcola(level-1)),      &
                      & Fina(Nodlev(level)+level-1),                    &
                      & Fina(Nodlev(level-1)+level-2),nca,mxnca,fncola, &
                      & List,Nonods,X(Nodlev(level-1)),                 &
                      & Y(Nodlev(level-1)),X(Nodlev(level)),            &
                      & Y(Nodlev(level)),Havmat)
         ENDIF
!
!           write(*,*)'JUST AFTER LOCALA call checka'
!       CALL CHECKA(A(PTCOLA(LEVEL)),
!     &            FINA(NODLEV(LEVEL)+LEVEL-1),
!     &            COLA(PTCOLA(LEVEL)),NCA,cnods )
!
         Na = Na + nca
         Ptcola(level+1) = Ptcola(level) + nca
         inca = Fina(Nodlev(level)+cnods+level-1)                       &
              & - Fina(Nodlev(level)+level-1)
         !WRITE (*,*) 'NCA,INCA:' , nca , inca
!
!
      ENDDO
!
      Totnod = Nodlev(Nlevel+1) - 1
      IF ( Totnod.GT.Maxtot ) THEN
         !WRITE (*,*) 'MAXTOT.LT.TOTNOD; Totnod,Maxtot:',Totnod,Maxtot
         STOP
      ENDIF
!
      END
!*==LOCALA.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
      SUBROUTINE LOCALA(Ca,Fa,Fnods,Cnods,Cgrid,Ccola,Fcola,Cfina,Ffina,&
                      & Nca,Mxnca,Fncola,List,Mxlist,Fx,Fy,Cx,Cy,Havmat)
      IMPLICIT NONE
!*--LOCALA2500
!*** Start of declarations inserted by SPAG
      REAL c2nod
      INTEGER i , iinod , inod , nlist , nod
!*** End of declarations inserted by SPAG
      INTEGER Havmat
      INTEGER Fnods , Cnods , Nca , Mxnca , Fncola
      INTEGER Cgrid(Fnods) , Ccola(Mxnca) , Fcola(Fncola)
      INTEGER Cfina(Cnods+1) , Ffina(Fnods+1)
      INTEGER ccount , count , count2 , cnod
      INTEGER Mxlist , List(Mxlist)
      REAL Fx(Fnods) , Fy(Fnods) , Cx(Cnods) , Cy(Cnods)
      REAL Fa(Fncola*Havmat) , Ca(Mxnca*Havmat)
      LOGICAL yes
! This sub finds the matrix on the causer mesh
! from the matrix on the finner mesh.
! Work out A for level LEVEL.
      ccount = 0
      DO nod = 1 , Fnods
         cnod = Cgrid(nod)
         IF ( cnod.NE.0 ) THEN
            Cx(cnod) = Fx(nod)
            Cy(cnod) = Fy(nod)
            nlist = 0
            Cfina(cnod) = ccount + 1
            DO count = Ffina(nod) , Ffina(nod+1) - 1
               inod = Fcola(count)
               DO count2 = Ffina(inod) , Ffina(inod+1) - 1
                  iinod = Fcola(count2)
                  c2nod = Cgrid(iinod)
                  IF ( c2nod.NE.0 ) THEN
! Is C2NOD in list of cause grid nodes surrounding cause node CNOD.
                     yes = .TRUE.
                     DO i = 1 , nlist
                        IF ( c2nod.EQ.List(i) ) yes = .FALSE.
                     ENDDO
                     IF ( yes ) THEN
                        nlist = nlist + 1
                        List(nlist) = c2nod
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
! Put list in assending order
! - now called GEM_IBUBLE to avoid clash with fluidity's version
! - (in case the two become different somehow...)
            CALL GEM_IBUBLE(List,nlist)
            DO i = 1 , nlist
               ccount = ccount + 1
               Ccola(ccount) = List(i)
            ENDDO
            Cfina(cnod+1) = ccount + 1
         ENDIF
      ENDDO
      Nca = ccount
!
! Now find the cause A (-CA) from the fine A (-FA).
      IF ( Havmat.EQ.1 ) THEN
         DO count = 1 , Nca
            Ca(count) = 1.
         ENDDO
         DO cnod = 1 , Cnods
            DO count = Cfina(cnod) , Cfina(cnod+1) - 1
               IF ( Ccola(count).EQ.cnod ) Ca(count) = 0.
            ENDDO
         ENDDO
         IF ( Nca.GT.Mxnca ) THEN
            !WRITE (*,*) 'NOT ENOUGH MEMORY SENT DOWN'
            !WRITE (*,*) 'INSIDE LOCALA'
            STOP
         ENDIF
      ENDIF
      END
!*==GEM_IBUBLE.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
! - now called GEM_IBUBLE, to avoid clash with fluidity's version
! - (will the two ever become different? -well, never mind...)
!
      SUBROUTINE GEM_IBUBLE(List,Nlist)
      IMPLICIT NONE
!*--GEM_IBUBLE2582
!*** Start of declarations inserted by SPAG
      INTEGER i , ii , j , Nlist
!*** End of declarations inserted by SPAG
      INTEGER List(Nlist)
      DO i = 1 , Nlist
         DO j = 2 , Nlist
            IF ( List(j-1).GT.List(j) ) THEN
! SWOP
               ii = List(j-1)
               List(j-1) = List(j)
               List(j) = ii
            ENDIF
         ENDDO
      ENDDO
      END
!*==NEUR.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
      SUBROUTINE NEUR(Whichd,Cgrid,V,B,Tscale,Nodlev,Ptcola,Totnod,     &
                    & Nonods,Onsubd,Fina,Cola,A,Na,Nlevel,Alpha,Lodbal, &
                    & Beta,Toler,Nchain,Havmat,EXACT, havwnod, wnod)
      IMPLICIT NONE
!*--NEUR2606
!*** Start of declarations inserted by SPAG
      INTEGER icstar , ifini2 , ifleve , ifstar , iileve , ilevel ,     &
            & istar2 , istart , Na , ncola , Nlevel , nonod2
      REAL rmean
      LOGICAL EXACT
!*** End of declarations inserted by SPAG
      INTEGER Havmat, havwnod
      REAL Beta , Toler
      INTEGER Onsubd , Nonods , Totnod
      INTEGER Nchain
      REAL V(Totnod*Onsubd)
      REAL B(Nonods*Onsubd)
      REAL Lodbal , beta2 , Alpha
      REAL A(Na*Havmat) , Tscale(Nonods)
      real wnod(totnod) 
      INTEGER i , fnods , cnods, ifstart, icstart
      INTEGER Cola(Na) , Fina(Totnod+Nlevel)
      INTEGER Whichd(Totnod) , Cgrid(Totnod)
      INTEGER Nodlev(Nlevel+1) , Ptcola(Nlevel+1)
      LOGICAL satur , mixup
! WHICHD(NOD) can overwrite CGRID to save space.
! WHICHD(NOD) contains the subdomain that node NOD is in.
! NONODS=number of nodes on the finnest grid level.
! TOTNOD=total number of nodes on all grid levels.
! if exact balance the no of nodes in each subdomain. 
! only works for uniform WNOD 
    !  WRITE (*,*) 'INSIDE NEUR ONSUBD,totnod,nlevel:' , Onsubd ,        &
    !            & Totnod , Nlevel
!      print *,'--wnod(1:nonods):',wnod(1:nonods)

      IF(TOTNOD.EQ.0) RETURN ! NOTHING TO DO AS EVERYTHING HAS ZERO LENGTH
!

      IF(HAVWNOD.NE.0) then
         DO ilevel = 1 , nlevel-1
            iileve = ilevel + 1
            fnods = Nodlev(iileve) - Nodlev(iileve-1)
            cnods = Nodlev(iileve+1) - Nodlev(iileve)
            istart = Nodlev(iileve-1)
            ncola = Ptcola(iileve) - Ptcola(iileve-1)
            ifstart = Nodlev(iileve-1)
            icstart = Nodlev(iileve)
            ifleve = iileve - 1
            CALL MAP_SIMPL(WNOD(ifstart),WNOD(icstart),fnods,cnods,Cgrid(istart), &
                           Fina(Nodlev(ifleve)+ifleve-1),Cola(Ptcola(ifleve)),ncola)
!            if(ilevel==nlevel-1) then
!               print *,'nonods,icstart,fnods,cnods:',nonods,icstart,fnods,cnods
!               PRINT *,'WNOD(icstart:icstart+cnods-1):',WNOD(icstart:icstart+cnods-1)
!               stop 292
!            endif
         END DO
      ENDIF ! IF(HAVWNOD.NE.0) then
!
! NB TOTNOD is the total number of cause and fine grid nodes.
!
!
      DO ilevel = Nlevel , 1 , -1
! LEVEL 1 is the fine mesh.
!
         beta2 = Beta
!
!
 50      IF ( ilevel.NE.Nlevel ) THEN
            mixup = .FALSE.
            iileve = ilevel + 1
            fnods = Nodlev(iileve) - Nodlev(iileve-1)
            cnods = Nodlev(iileve+1) - Nodlev(iileve)
            ifstar = (Nodlev(iileve-1)-1)*Onsubd + 1
            icstar = (Nodlev(iileve)-1)*Onsubd + 1
            istart = Nodlev(iileve-1)
            ifleve = iileve - 1
            ncola = Ptcola(iileve) - Ptcola(iileve-1)
!
            CALL MAP(V(ifstar),V(icstar),fnods,cnods,Onsubd,            &
                   & Fina(Nodlev(ifleve)+ifleve-1),Cola(Ptcola(ifleve)),&
                   & Cgrid(istart),ncola)
         ELSE
            mixup = .TRUE.
         ENDIF
!
         nonod2 = Nodlev(ilevel+1) - Nodlev(ilevel)
         istart = Nodlev(ilevel)
         istar2 = (Nodlev(ilevel)-1)*Onsubd + 1
         ifini2 = (Nodlev(ilevel+1)-1)*Onsubd
         ncola = Ptcola(ilevel+1) - Ptcola(ilevel)
!
         !WRITE (*,*) '**** ILEVEL=' , ilevel
!
         !WRITE (*,*) 'HAVMAT=' , Havmat
         IF ( Havmat.EQ.1 ) THEN
            CALL NEURAL(V(istar2),B,Tscale,mixup,nonod2,Onsubd,         &
                      & Fina(Nodlev(ilevel)+ilevel-1),                  &
                      & Cola(Ptcola(ilevel)),A(Ptcola(ilevel)),ncola,   &
                      & Lodbal,beta2,Toler,Nchain,Havmat, havwnod,wnod(istart))
         ELSE
          !print *,'istar2,Totnod,Onsubd,NONODS=',istar2,Totnod,Onsubd,NONODS
            CALL NEURAL(V(istar2),B,Tscale,mixup,nonod2,Onsubd,         &
                      & Fina(Nodlev(ilevel)+ilevel-1),                  &
                      & Cola(Ptcola(ilevel)),A,ncola,Lodbal,beta2,Toler,&
                      & Nchain,Havmat, havwnod,wnod(istart))
         ENDIF

!
! Test the validity of solution V ****************
         rmean = 0.
         DO i = istar2 , ifini2
            rmean = rmean + MAX(V(i),1.-V(i))
         ENDDO
         rmean = rmean/REAL(ifini2-istar2)
         satur = .TRUE.
         IF ( rmean.LT.0.7 ) satur = .FALSE.
         !WRITE (*,*) 'THE SOLUTION SATURATION :RMEAN=' , rmean
         !!print *,'v:',v
         !WRITE (*,*) 'SATUR=' , satur
!
         IF ( .NOT.satur ) THEN
            beta2 = beta2*0.8
            GOTO 50
         ENDIF
!
! *************************************************
!
         mixup = .FALSE.
!
      ENDDO
!
!
      !WRITE (*,*) 'INSIDE NEUR HERE 4'
!
! MAKE SURE ALL NODES ARE SATURATED ***************.
! WORK OUT WHICHD *********************************.
      DO ilevel = Nlevel , 1 , -1
         cnods = Nodlev(ilevel+1) - Nodlev(ilevel)
         icstar = (Nodlev(ilevel)-1)*Onsubd + 1
         istart = Nodlev(ilevel)
         CALL FIWICD(V(icstar),Whichd(istart),cnods,Onsubd,exact,havwnod,wnod(istart))
      ENDDO
! *************************************************.
      !WRITE (*,*) 'GOING OUT OF NEUR ONSUBD:' , Onsubd
      RETURN
      END
!*==CHECKA.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
      SUBROUTINE CHECKA(A,Fina,Cola,Ncola,Nonods)
      IMPLICIT NONE
!*--CHECKA2724
!*** Start of declarations inserted by SPAG
      INTEGER i , icol , j
!*** End of declarations inserted by SPAG
! This sub checks to see if A is O.K.
      INTEGER Ncola , Nonods , Fina(Nonods+1) , Cola(Ncola)
      INTEGER count
      REAL A(Ncola)
      IF ( Ncola.LT.Fina(Nonods+1)-1 ) THEN
         !WRITE (*,*) 'INSIDE CHECKA'
         !WRITE (*,*) 'PROBLEM WITH BOUNDS'
         STOP
      ENDIF
      DO i = 1 , Nonods
         DO count = Fina(i) , Fina(i+1) - 1
            icol = Cola(count)
            IF ( icol.NE.i ) THEN
               IF ( A(count).LT.0.00001 ) THEN
                  !WRITE (*,*) 'PROBLEM WITH A:'
                  !WRITE (*,*) 'COUNT,I,ICOL:' , count , i , icol
                  !WRITE (*,*) 'A(COUNT):' , A(count)
                  !WRITE (*,*) 'FINA(I+1)-FINA(I):' , Fina(i+1) - Fina(i)
                  !WRITE (*,*) 'A(J),J=FINA(I),FINA(I+1)-1' ,            &
                  !          & (A(j),j=Fina(i),Fina(i+1)-1)
                  STOP
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      END
!*==FIWICD.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
      SUBROUTINE FIWICD(V,Whichd,Nonods,Onsubd,exact,havwnod,wnod)
      IMPLICIT NONE
!*--FIWICD2761
!*** Start of declarations inserted by SPAG
!*** End of declarations inserted by SPAG
! This sub finds whichd from the neuron values V.
! if exact calculate which subdomain every node belong too to exactly balance the load as far as possible. 
      INTEGER Nonods , Onsubd , Whichd(Nonods) , sub , maxsub, havwnod
      REAL V(Nonods*Onsubd) , wnod(nonods), maxv
      logical exact
! local variables
      INTEGER i , ii , isub, n_nods_bal
      DO i = 1 , Nonods
         maxv = 0.
         DO sub = 1 , Onsubd
            IF ( V((sub-1)*Nonods+i).GT.maxv ) THEN
               maxv = V((sub-1)*Nonods+i)
               maxsub = sub
            ENDIF
!             V((SUB-1)*NONODS+I)=0.
         ENDDO
!          V((MAXSUB-1)*NONODS+I)=1.
         Whichd(i) = maxsub
      ENDDO

      if(exact) then ! try to make balance the number of nodes. 
         whichd=0
         DO isub = 1 , Onsubd
            call critical_bal_size(isub,v((isub-1)*nonods+1),onsubd,nonods,whichd,havwnod,wnod) ! only consider subdomains that have not been filled in. 
         end do
      endif

! count no of nodes in each subdomain
      !WRITE (*,*) 'onsubd,nonods:' , Onsubd , Nonods
      DO isub = 1 , Onsubd
         ii = 0
         DO i = 1 , Nonods
            IF ( Whichd(i).EQ.isub ) ii = ii + 1
         ENDDO
         !WRITE (*,*) 'no of nodes in sub=' , isub , ' is =' , ii
      ENDDO
      END
!*==MAP.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
      SUBROUTINE critical_bal_size(isub,v,onsubd,nonods,whichd,havwnod,wnod) ! only consider subdomains that have not been filled in.
! calculate which subdomain every node belong too to exactly balance the load as far as possible. 
      implicit none
      integer isub,onsubd,nonods,havwnod
      real wnod(nonods) 
      real v(nonods)
      integer whichd(nonods)
! local variables...
      integer nits
      parameter(nits=40) 
      integer its,nod
      real n_nods_bal
      real low_crit_bal,high_crit_bal
      real crit_bal,rwsum,rcount,tol_near_one

      if(isub==onsubd) then
         do nod=1,nonods
            if(whichd(nod)==0) whichd(nod)=isub
         end do
         return ! simple end to make sure we have all the nodes attached to subdomains. 
      endif

      if(havwnod.ne.0) then
         rwsum = sum(wnod) 
      else
         rwsum = real(nonods)
      endif

      n_nods_bal=rwsum/real(onsubd)  

      tol_near_one = 0.99* rwsum/real(nonods) 

      low_crit_bal=0.0
      high_crit_bal=1.0    

      crit_bal=1.0/real(onsubd)

      do its=1,nits

         rcount=0.0
         do nod=1,nonods
            if(whichd(nod)==0) then
               if(havwnod.ne.0) then
                  if(v(nod).gt.crit_bal) rcount=rcount+wnod(nod) 
               else
                  if(v(nod).gt.crit_bal) rcount=rcount+1.0
               endif
            endif
         end do
         if(rcount.gt.n_nods_bal) then
            low_crit_bal=crit_bal
         else
            high_crit_bal=crit_bal
         endif
         !print *,'its,low_crit_bal,high_crit_bal,isub,count,n_nods_bal:', &
         !         its,low_crit_bal,high_crit_bal,isub,count,n_nods_bal

!         if( (count > n_nods_bal-1.000001).and.(count < n_nods_bal+1.000001) ) exit 
         if( (rcount > n_nods_bal-tol_near_one).and.(rcount < n_nods_bal+tol_near_one) ) exit 
!         if( (count > n_nods_bal-0.99).and.(count < n_nods_bal+0.99) ) exit 
!         if( (count > n_nods_bal-0.1).and.(count < n_nods_bal+0.1) ) exit 
         crit_bal=0.5*(low_crit_bal + high_crit_bal) 
      end do ! do its=1,nits

      rcount=0
      do nod=1,nonods
         if(whichd(nod)==0) then
            if(v(nod).gt.crit_bal) then
               rcount=rcount+1
               whichd(nod)=isub
            end if
         endif
      end do

      !print *,'final isub,nod of nodes,n_nods_bal:',isub,count,n_nods_bal
!      !print *,'-v:',v
      
      return
      end 
!
!
!
!
      SUBROUTINE MAP_SIMPL(Vl,Vs,Nonods,Nshort,Cgrid, &
                           FINAF,COLAF,NCOLAF)
      IMPLICIT NONE
!*--MAP2797
!*** Start of declarations inserted by SPAG
      INTEGER Nonods, Nshort, NCOLAF
      LOGICAL SIMPLE_INTERP
      PARAMETER(SIMPLE_INTERP=.FALSE.)
!*** End of declarations inserted by SPAG
      INTEGER Onsubd
      INTEGER Cgrid(Nonods) 
      REAL Vl(Nonods) , Vs(Nshort)
      integer FINAF(NONODS+1),COLAF(NCOLAF)
! CGRID(Fine grid node)=0 if it is not the cause grid node
! = corresponding cause grid node number, otherwise,
! this number varies between 1 and NSHORT.
! Lcal variables...
      INTEGER COUNT,ICOUNT,J, ndcgri , nod 
      REAL RSUM

      DO nod = 1 , Nonods
!
         ndcgri = Cgrid(nod)
!         IF ( ndcgri.NE.0 ) Vl(NOD) = Vs(ndcgri)
         IF ( ndcgri.NE.0 )  THEN
            IF(SIMPLE_INTERP) THEN
               Vs(ndcgri) = Vl(NOD)
            ELSE
               RSUM=0.0
               ICOUNT=0
               DO COUNT=FINAF(NOD),FINAF(NOD+1)-1
                  J=COLAF(COUNT)
                  RSUM=RSUM + Vl(J)
                  ICOUNT=ICOUNT+1
               END DO
               ICOUNT=ICOUNT+1
               RSUM=RSUM + Vl(NOD)
               Vs(ndcgri) = RSUM/REAL(ICOUNT)
            ENDIF ! ENDOF IF THEN ELSE IF(SIMPLE_INTERP) THEN
         ENDIF
      ENDDO
      END SUBROUTINE MAP_SIMPL
!
!
!      
!
      SUBROUTINE MAP(Vl,Vs,Nonods,Nshort,Onsubd,Fina,Cola,Cgrid,Ncola)
      IMPLICIT NONE
!*--MAP2797
!*** Start of declarations inserted by SPAG
      INTEGER icol , isub , Ncola , ndcgri , nod , Nonods , Nshort
      REAL rsum
!*** End of declarations inserted by SPAG
      INTEGER Onsubd
      INTEGER Cgrid(Nonods) , Fina(Nonods+1) , Cola(Ncola)
      REAL Vl(Nonods*Onsubd) , Vs(Nshort*Onsubd)
! CGRID(Fine grid node)=0 if it is not the cause grid node
! = corresponding cause grid node number, otherwise,
! this number varies between 1 and NSHORT.
      logical adjust_small
      parameter(adjust_small=.true.) 
      INTEGER count

      DO isub = 1 , Onsubd
         DO nod = 1 , Nonods
            Vl(nod+(isub-1)*Nonods) = 0.
         ENDDO
      ENDDO
!
      DO nod = 1 , Nonods
!
         ndcgri = Cgrid(nod)
         IF ( ndcgri.NE.0 ) THEN
            DO count = Fina(nod) , Fina(nod+1) - 1
               icol = Cola(count)
               DO isub = 1 , Onsubd
                  Vl(icol+(isub-1)*Nonods) = Vl(icol+(isub-1)*Nonods)   &
                   & + Vs(ndcgri+(isub-1)*Nshort)
               ENDDO
            ENDDO
         ENDIF
!
      ENDDO
!
!
! Now normalise VL
      DO nod = 1 , Nonods
         rsum = 0.
         DO isub = 1 , Onsubd
            rsum = rsum + Vl(nod+(isub-1)*Nonods)
         ENDDO
!          IF(rSUM.LT.0.000001) THEN
!            write(*,*)'rSUM,NOD:',rSUM,NOD
!            write(*,*)'ABOUT TO DIVIDE BY ZERO'
!            STOP
!          ENDIF
         DO isub = 1 , Onsubd
            Vl(nod+(isub-1)*Nonods) = Vl(nod+(isub-1)*Nonods)/rsum
         ENDDO
      ENDDO
      
      if(adjust_small) then
      if(nonods.le.onsubd) then
         call random_number(vl)
!         !print *,'vl:',vl
         vl = 1./real(onsubd)  + 0.01*(0.5-vl/real(onsubd))
!         !print *,'-vl:',vl
!         stop 2011
         DO nod = 1 , Nonods
            rsum = 0.
            DO isub = 1 , Onsubd
               
               rsum = rsum + Vl(nod+(isub-1)*Nonods)
            ENDDO
            DO isub = 1 , Onsubd
               Vl(nod+(isub-1)*Nonods) = Vl(nod+(isub-1)*Nonods)/rsum
            ENDDO
         ENDDO
      endif
      endif

      !print *,'*********nonods,vl:',nonods,vl
      END
!*==NEURAL.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
      SUBROUTINE NEURAL(V,B,Tscale,Mixup,Nonods,Onsubd,Fina,Cola,A,     &
                      & Ncola,Lodbal,Beta,Toler,Nchain,Havmat, havwnod,wnod)
      IMPLICIT NONE
!*--NEURAL2855
!*** Start of declarations inserted by SPAG
      INTEGER ichain , iichain , inits , iseed , kk , Ncola , NINITS, havwnod
      REAL r , rninle
!*** End of declarations inserted by SPAG
      INTEGER Havmat
      REAL Beta , SUGEST , TOLINT , Toler
      INTEGER Onsubd , Nonods
      INTEGER Nchain , MXNDOM
      real toler_small
      PARAMETER (MXNDOM=32)
      PARAMETER (SUGEST=0.1, NINITS=1, toler_small=1.e-15)
!        PARAMETER(TOLINT=0.0001,TOLER=0.001)
!        PARAMETER(TOLINT=0.0001,TOLER=0.0001)
      PARAMETER (TOLINT=0.0001)
      REAL V(Nonods*Onsubd)
      REAL B(Nonods*Onsubd) , Tscale(Nonods)
      real wnod(nonods)
      REAL A(Ncola*Havmat)
      REAL sumv(MXNDOM) , stoexp(MXNDOM) , rsum
      REAL Lodbal
      REAL ttvv , f
      REAL alpha
      REAL RAN1 , rran , maxdif, maxgl, wsum
      INTEGER sub , nodsub , i , j, colaco
      INTEGER Cola(Ncola) , Fina(Nonods+1)
      INTEGER count
      LOGICAL Mixup
      REAL, DIMENSION(:), ALLOCATABLE :: wwnod ! the node weights defined internally. 
! LODBAL=1. is the default .gt.1 then more importance
! is placed on load balancing.
! BETA=0.9 is the default but this number controles
! how close to the critical temp the optimization is performed
! WHEN BETA=0.5 - HALF CRYTICAL TEMP IS USED.
! TOLER =the solution tolerence 0.0001 is suggested.
! NCHAIN =maximum no of iterations 3000 is suggested.
!
! If MIXUP then start with a random solution else use one already
! in V.
! WHICHD(FNOD) initially contains the number of cause grid
! nodes that fine grid node FNOD takes its value from (IS SURROUNDED BY).
! But eventually contains the subdomain that node FNOD is in.
!
! Initialise.
! Randomly initialise only on causest grid *******************
      !WRITE (*,*) 'inside neural'
      IF ( Mixup ) THEN
         iseed = 1
         DO sub = 1 , Onsubd
            DO i = 1 , Nonods
               V((sub-1)*Nonods+i) = 1./REAL(Onsubd)
            ENDDO
            DO i = 1 , Nonods
               rran = RAN1(iseed)
               IF ( rran.GE.1 ) rran = 0.999999
!
               V((sub-1)*Nonods+i) = rran*SUGEST + 1./REAL(Onsubd)      &
                                   & - 0.5*SUGEST
!
            ENDDO
         ENDDO
      ENDIF
      !WRITE (*,*) 'inside neural HERE1'
!
      !WRITE (*,*) 'fina(1):' , Fina(1)
      !WRITE (*,*) 'fina(NONODS+1):' , Fina(Nonods+1)
      !WRITE (*,*) 'ncola:' , Ncola
      !WRITE (*,*) 'NONODS:' , Nonods
      !WRITE (*,*) 'BETA,LODBAL,NCHAIN:' , Beta , Lodbal , Nchain

      allocate(wwnod(nonods))
! 
      if(havwnod.ne.0) then
         wwnod=wnod
      else
         wwnod=1.0
      endif
!
      wsum=sum(wwnod)
!
! Work out ALPHA.
      ttvv = 0.
      IF ( Havwnod.EQ.2 ) THEN
         do i=1,nonods
            DO count = Fina(i) , Fina(i+1) - 1
               j = cola(count)
               ttvv = ttvv + 0.5*(wwnod(i) + wwnod(j))
            ENDDO
!            ttvv = ttvv - wwnod(i)  ! take away diagonal
         end do
      ELSE IF ( Havmat.EQ.1 ) THEN
         do i=1,nonods
            DO count = Fina(i) , Fina(i+1) - 1
            ttvv = ttvv + A(count)
         ENDDO
         end do
      ELSE
         do i=1,nonods
            DO count = Fina(i) , Fina(i+1) - 1
            ttvv = ttvv + 1.
         ENDDO
!            ttvv = ttvv - 1.  ! take away diagonal
         end do
      ENDIF
!      alpha = Lodbal*REAL(Onsubd)*ttvv/(REAL(Nonods)**2)
      alpha = Lodbal*REAL(Onsubd)*ttvv/ (wsum**2) 
!      WRITE (*,*) 'havwnod,alpha:' , havwnod,alpha
!      WRITE (*,*) 'wwnod:',wwnod
!      stop 2921
!        ALPHA=100.
!
!
! Work out the source B() for all the neurons.
      DO sub = 1 , Onsubd
         rninle = wsum/REAL(Onsubd)
         DO i = 1 , Nonods
            ttvv = 0.
            IF ( Havwnod.EQ.2 ) THEN
               DO count = Fina(i) , Fina(i+1) - 1
                  j = cola(count)
                  ttvv = ttvv + 0.5*(wwnod(i)+wwnod(j))
               ENDDO
!               ttvv = ttvv - wwnod(i) ! take away diagonal
            ELSE IF ( Havmat.EQ.1 ) THEN
               DO count = Fina(i) , Fina(i+1) - 1
                  ttvv = ttvv + A(count)
               ENDDO
            ELSE
               DO count = Fina(i) , Fina(i+1) - 1
                  ttvv = ttvv + 1.
               ENDDO
!               ttvv = ttvv - 1. ! take away diagonal
            ENDIF
            B(i+(sub-1)*Nonods) = 0.5*(ttvv-alpha*wwnod(i)*rninle)
         ENDDO
      ENDDO
!
!
! Work out TSCALE for all the nodes.
      DO i = 1 , Nonods
         ttvv = 0.
         
         IF ( Havwnod.EQ.2 ) THEN
            DO count = Fina(i) , Fina(i+1) - 1
              j = cola(count)
              ttvv = ttvv + 0.5*(wwnod(i)+wwnod(j))
            ENDDO
!            ttvv = ttvv - wwnod(i) ! take away diagonal
         ELSE IF ( Havmat.EQ.1 ) THEN
            DO count = Fina(i) , Fina(i+1) - 1
               ttvv = ttvv + A(count)
            ENDDO
         ELSE
            DO count = Fina(i) , Fina(i+1) - 1
               ttvv = ttvv + 1.
            ENDDO
!            ttvv = ttvv - 1. ! subtract out the digonal value
         ENDIF
         Tscale(i) = Beta*ttvv/REAL(Onsubd) ! does not depend on load balancing
!          if(ttvv.lt.0.00001) then
!            write(*,*)'problem with tscale(i),i:',tscale(i),i
!            write(*,*)'FINA(I+1)-FINA(I):',FINA(I+1)-FINA(I)
!            DO 1865 COUNT=FINA(I),FINA(I+1)-1
!              write(*,*)'count,A(COUNT):',count,A(COUNT)
!              write(*,*)'COLA(COUNT):',COLA(COUNT)
!1865         CONTINUE
!            stop
!          endif
      ENDDO
!
!
!
!
! NOW START ITERATION
      iichain = Nchain
!       IF(ILEVEL.EQ.NLEVEL) IICHAIN=500
      DO ichain = 1 , Nchain
!        DO 210 ICHAIN=1,IICHAIN
! work out SUMV(SUB)
         DO sub = 1 , Onsubd
            rsum = 0.
            DO i = 1 , Nonods
               rsum = rsum + V((sub-1)*Nonods+i)*wwnod(i)
            ENDDO
            sumv(sub) = rsum
         ENDDO
!        write(*,*)'inside neural HERE5'
!
         maxgl = 0.
         DO i = 1 , Nonods
! TOTNOD is the total number of cause and fine grid nodes.
!
            DO inits = 1 , NINITS
               rsum = 0.
                  DO sub = 1 , Onsubd
                     nodsub = (sub-1)*Nonods + i
                     ttvv = 0.
                     if( Havwnod.EQ.2 ) THEN
                        DO count = Fina(i) , Fina(i+1) - 1
                           j=Cola(count)
                           colaco = j + (sub-1)*Nonods
                           ttvv = ttvv + 0.5*(wwnod(i)+wwnod(j))*V(colaco)
                        ENDDO
                        colaco = i + (sub-1)*Nonods
!                        ttvv = ttvv - wwnod(i)*V(colaco) ! dont include the diagonal
                     else if ( Havmat.EQ.1 ) THEN
                     DO count = Fina(i) , Fina(i+1) - 1
                        colaco = Cola(count) + (sub-1)*Nonods
                        ttvv = ttvv + A(count)*V(colaco)
                     ENDDO
                     else
                     DO count = Fina(i) , Fina(i+1) - 1
                        colaco = Cola(count) + (sub-1)*Nonods
                        ttvv = ttvv + V(colaco)
                     ENDDO
                        colaco = i + (sub-1)*Nonods
!                       ttvv = ttvv - V(colaco) ! dont include the diagonal
                     endif
!
! C is the source from the surrounding other layers
! in multi-grid approach.
! B is the source for each neuron.
                     f = -ttvv + alpha*wwnod(i)*sumv(sub) + B(nodsub)
                     stoexp(sub) = EXP(-f/Tscale(i))
                     rsum = rsum + stoexp(sub)
                  ENDDO
!
! Update part that attempts to balance the number of nodes.
               maxdif = 0.
               DO sub = 1 , Onsubd
                  nodsub = (sub-1)*Nonods + i
                  sumv(sub) = sumv(sub) - V(nodsub)*wwnod(i)
                  !if(rsum==0.0) 
                  !!print *,'rsum=',rsum
                  r = stoexp(sub)/max(toler_small, rsum) 
                  maxdif = MAX(maxdif,ABS(r-V(nodsub)))
                  maxgl = MAX(maxdif,maxgl)
                  V(nodsub) = r
                  sumv(sub) = sumv(sub) + V(nodsub)*wwnod(i)
               ENDDO
!
               IF ( maxdif.LT.TOLINT ) GOTO 50
            ENDDO
!         write(*,*)'ichain,ilevel,I,INITS:',ichain,ilevel,I,INITS
!
 50      ENDDO
!
 
!         write(*,*)'ichain,maxgl,TOLER,sumv(1),sumv(2):',
!     &            ichain,maxgl,TOLER,sumv(1),sumv(2)
         kk = ichain
         IF ( maxgl.LT.Toler ) GOTO 100
!
      ENDDO
!
 100  continue! WRITE (*,*) 'NO OF ITERATIONS=' , kk
      END
!*==MULPV.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
      SUBROUTINE MULPV(Nlong,Nshort,Xcount,X,Vec,Fina,Cola,Ncola)
      IMPLICIT NONE
!*--MULPV3074
!*** Start of declarations inserted by SPAG
      INTEGER i , icount
!*** End of declarations inserted by SPAG
      INTEGER Nlong , Nshort , Ncola , cnod , colj
      INTEGER Fina(Nlong+1) , Cola(Ncola)
      INTEGER Xcount(Nlong)
!         integer COLOR(NLONG)
      REAL X(Nlong) , Vec(Nshort)
! This subroutine performs X=P*vec multiplication where
! P is the prolongation operator.
! P is not explicitly formed.
! COLOR contains the colouring on the finner grid level.
! XCOUNT(FNOD) contains the number of cause grid
! nodes that fine grid node FNOD takes its value from (IS SURROUNDED BY).
      DO i = 1 , Nlong
         X(i) = 0.
!           XCOUNT(I)=0
      ENDDO
!
      cnod = 0
      DO i = 1 , Nlong
!           IF(COLOR(I).EQ.1) THEN
         IF ( Xcount(i).EQ.1 ) THEN
            cnod = cnod + 1
            DO icount = Fina(i) , Fina(i+1) - 1
               colj = Cola(icount)
!               XCOUNT(COLJ)=XCOUNT(COLJ)+1
               X(colj) = X(colj) + Vec(cnod)
            ENDDO
         ENDIF
      ENDDO
!
      DO i = 1 , Nlong
         X(i) = X(i)/REAL(Xcount(i))
      ENDDO
      END
!*==MULPTV.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
      SUBROUTINE MULPTV(Nlong,Nshort,Xcount,X,Vec,Fina,Cola,Ncola)
      IMPLICIT NONE
!*--MULPTV3118
!*** Start of declarations inserted by SPAG
      INTEGER i , icount
      REAL r
!*** End of declarations inserted by SPAG
! This subroutine performs X=P^T*vec multiplication where
! P^T is the RESTRICTION operator.
! P is not explicitly formed.
! COLOR contains the colouring on the finner grid level.
! XCOUNT(FNOD) contains the number of cause grid
! nodes that fine grid node FNOD takes its value from (IS SURROUNDED BY).
      INTEGER Nlong , Nshort , Ncola , cnod , colj
      INTEGER Fina(Nlong+1) , Cola(Ncola)
      INTEGER Xcount(Nlong)
      REAL X(Nshort) , Vec(Nlong)
!
      cnod = 0
      DO i = 1 , Nlong
         IF ( Xcount(i).EQ.1 ) THEN
            cnod = cnod + 1
            r = 0.
            DO icount = Fina(i) , Fina(i+1) - 1
               colj = Cola(icount)
               r = r + Vec(colj)/REAL(Xcount(colj))
            ENDDO
            X(cnod) = r
         ENDIF
      ENDDO
      END
!*==GETMES.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
!
!
      SUBROUTINE GETMES(Cola,Fina,Aa,Nonods,Nodacr,Nodown,Maxna)
      IMPLICIT NONE
!*--GETMES3156
!*** Start of declarations inserted by SPAG
      INTEGER i , iacr , jdow , Maxna
      REAL SCALE
!*** End of declarations inserted by SPAG
      INTEGER Nodacr , Nodown , Nonods
      PARAMETER (SCALE=1.)
!        PARAMETER(ONSUBD=2)
!        PARAMETER(NODACR=16,NODOWN=16,NONODS=NODACR*NODOWN)
      REAL Aa(Maxna) , rp
!        REAL XCOR(NONODS),YCOR(NONODS),DSQRT(NONODS)
      INTEGER Cola(Maxna) , Fina(Nonods+1)
      INTEGER row , count
!        COMMON /MIDPA/ MIDPA,COLA,FINA,AA
!        COMMON /XYCOR/ XCOR,YCOR,MXX,MXY
      row = 0
      count = 0
!        MXX=REAL(NODACR-1)
!        MXY=REAL(NODOWN-1)
      DO i = 1 , Maxna
         Aa(i) = 0.
      ENDDO
!
!
      DO jdow = 1 , Nodown
         DO iacr = 1 , Nodacr
            row = row + 1
!
!            XCOR(ROW)=REAL(IACR)-1.
!            YCOR(ROW)=REAL(JDOW)-1.
!
            Fina(row) = count + 1
            rp = 0
!
            IF ( jdow.NE.1 ) THEN
               count = count + 1
               Cola(count) = row - Nodacr
               Aa(count) = -SCALE
               rp = rp + 1.
            ENDIF
!
            IF ( iacr.NE.1 ) THEN
               count = count + 1
               Cola(count) = row - 1
               Aa(count) = -SCALE
               rp = rp + 1.
            ENDIF
!
            count = count + 1
            Cola(count) = row
!                MIDPA(ROW)=COUNT
!
            IF ( iacr.NE.Nodacr ) THEN
               count = count + 1
               Cola(count) = row + 1
               Aa(count) = -SCALE
               rp = rp + 1.
            ENDIF
!
            IF ( jdow.NE.Nodown ) THEN
               count = count + 1
               Cola(count) = row + Nodacr
               Aa(count) = -SCALE
               rp = rp + 1.
            ENDIF
!
!                AA(MIDPA(ROW))=SP*RP
!
         ENDDO
      ENDDO
      Fina(Nonods+1) = count + 1
! apply boundary condition.
!            AA(FINA(1))=INFINY
      DO count = 1 , Fina(Nonods+1) - 1
         Aa(count) = -Aa(count)
      ENDDO
!
!         DO 35 I=1,NONODS
!           XCOR(I)=XCOR(I)/MXX
!           YCOR(I)=YCOR(I)/MXY
!35       CONTINUE
!          MXX=1.
!          MXY=1.
      END
!*==OUTPUT.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
!
      SUBROUTINE OUTPUT(Whichd,Nodacr,Nodown,Nonods,Totnod,Onsubd,V,    &
                      & Mored,Fgrid,Cgrid,Nodlev,Nlevel)
      IMPLICIT NONE
!*--OUTPUT3249
!*** Start of declarations inserted by SPAG
      INTEGER i , ilevel , ilevel2 , Nlevel
!*** End of declarations inserted by SPAG
      INTEGER Onsubd , Nodacr , Nodown , Nonods , Totnod
      INTEGER jdow , iacr , Whichd(Totnod) , fnod
      INTEGER Mored(Totnod) , Cgrid(Totnod) , Fgrid(Totnod)
      INTEGER Nodlev(Nlevel+1)
      REAL V(Onsubd*Totnod)
!
      !WRITE (*,*) 'nodacr,nodown,nonods:' , Nodacr , Nodown , Nonods
      !WRITE (*,*) 'onsubd,totnod' , Onsubd , Totnod
      !WRITE (*,*) ' '
      DO jdow = 1 , Nodown
!            write(6,'(40(1X,I3))')
!            write(6,'(60(0X,I2))')
!            write(*,*)
     !    WRITE (6,'(60(I2))') (Whichd((jdow-1)*Nodacr+iacr),iacr=1,     &
      !                      & Nodacr)
      ENDDO
!          write(*,*)(V(I),I=1,ONSUBD*TOTNOD)
!
! We will put the decomposition in MORED.
      DO i = 1 , Nonods
         Mored(i) = Whichd(i)
      ENDDO
!
! !print the other levels
      DO ilevel = 2 , Nlevel
         DO i = 1 , Nonods
            Mored(i) = 0
         ENDDO
         DO i = Nodlev(ilevel) , Nodlev(ilevel+1) - 1
            fnod = Fgrid(i)
            DO ilevel2 = 3 , ilevel
               fnod = Fgrid(fnod)
            ENDDO
! FNOD is the finest grid node.
            Mored(fnod) = Whichd(i)
         ENDDO
!
      !   WRITE (*,*) 'PARTITION FOR LEVEL:' , ilevel
         DO jdow = 1 , Nodown
       !     WRITE (6,'(40(1X,I2))') (Mored((jdow-1)*Nodacr+iacr),iacr=1,&
        !                          & Nodacr)
         ENDDO
      ENDDO
      END
!*==RAN1.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
!
!
!      idum = 1
!      do 10 i = 1, 30
!        z = ran1(idum)
!10    write(*,*) idum, z, z*1000000
!      end
!
!
      FUNCTION RAN1(Idum)
      IMPLICIT NONE
!*--RAN13313
!*** Start of declarations inserted by SPAG
      INTEGER IA1 , IA2 , IA3 , IC1 , IC2 , IC3 , Idum , iff , ix1 ,    &
            & ix2 , ix3 , j , M1 , M2 , M3
      REAL r , RAN1 , RM1 , RM2
!*** End of declarations inserted by SPAG
      DIMENSION r(97)
      PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      PARAMETER (M3=243000,IA3=4561,IC3=51349)
      logical new_ran
      parameter(new_ran=.true.)
      DATA iff/0/
      real rh1

      if(new_ran) then
         call random_number(rh1)
         ran1=rh1
!        !print *,'ran1:',ran1
         return
      endif

      IF ( Idum.LT.0 .OR. iff.EQ.0 ) THEN
         iff = 1
         ix1 = MOD(IC1-Idum,M1)
         ix1 = MOD(IA1*ix1+IC1,M1)
         ix2 = MOD(ix1,M2)
         ix1 = MOD(IA1*ix1+IC1,M1)
         ix3 = MOD(ix1,M3)
         DO j = 1 , 97
            ix1 = MOD(IA1*ix1+IC1,M1)
            ix2 = MOD(IA2*ix2+IC2,M2)
            r(j) = (FLOAT(ix1)+FLOAT(ix2)*RM2)*RM1
         ENDDO
         Idum = 1
      ENDIF
      ix1 = MOD(IA1*ix1+IC1,M1)
      ix2 = MOD(IA2*ix2+IC2,M2)
      ix3 = MOD(IA3*ix3+IC3,M3)
      j = 1 + (97*ix3)/M3
!      IF ( j.GT.97 .OR. j.LT.1 ) PAUSE
      IF ( j.GT.97 .OR. j.LT.1 ) then
         !print *,'j,idum=',j,idum
         STOP 27
      endif
      RAN1 = r(j)
      r(j) = (FLOAT(ix1)+FLOAT(ix2)*RM2)*RM1
      END
!*==REORG2.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
! **********************************************************
! **** THIS PART RE-ORDERS THE EQUATIONS *******************
! **********************************************************
!
!
! For the subdomain to proc mapping.......
      SUBROUTINE REORG2(Fina,Cola,Ncola,Nonods,Ndom,Whichd,Renum,Finsub,&
                      & Colsub,Mxcols,Q,Blank,Blasub,Maxnq,Defaul,      &
                      & Hypdim,Topcon,Sidcon,Proacr,Prodow,Imem,Nimem)
      IMPLICIT NONE
!*--REORG23361
!*** Start of declarations inserted by SPAG
      INTEGER i , iedist , iiqn , imem2 , insub , Nimem , nimem2 ,      &
            & nintnd
!*** End of declarations inserted by SPAG
!
      INTEGER Maxnq , Ndom , Nonods , Ncola , Mxcols
      INTEGER Fina(Nonods+1) , Cola(Ncola)
      INTEGER Whichd(Nonods) , Renum(Nonods)
      INTEGER Colsub(Mxcols) , Finsub(Ndom+1)
! If DEFAUL then map using rordering routine.
! COLSUB contains the subdomain to subdomain communication matrix.
! MXCOLS= maximum value of NCOLSU
      INTEGER Q(Maxnq) , Blank(Nonods) , Blasub(Ndom)
!
      INTEGER sub , sub2 , ncolsu , newsub , stasub
      INTEGER count , count2
      LOGICAL found , zero
! For the subdomain to proc mapping.......
      LOGICAL Defaul , Topcon , Sidcon
      INTEGER Hypdim , Proacr , Prodow , Imem(Nimem)
! This sub reorders the equations for entire global system of
! equations.
! RENUM(OLDNOD no)=new nod no.
 
      IF ( Ndom.LE.1 ) THEN
! for entire system.
         zero = .TRUE.
         DO i = 1 , Nonods
            Blank(i) = 1
         ENDDO
         iiqn = MAX(500,Nonods/2)
         CALL REORD2(Blank,Fina,Cola,Nonods,Ncola,Renum,Q(500),         &
                   & Q(Nonods/2+501),iiqn,zero)
      ENDIF
!
      IF ( Ndom.GT.1 ) THEN
! REORDER THE SUBDOMAINS ****************
         nintnd = 0
         IF ( nintnd.EQ.0 ) THEN
! This is for BFBGS
! find pointer for subdomains, then reorder subdomains.
            count2 = 0
            DO sub = 1 , Ndom
               Blasub(sub) = 1
               Finsub(sub) = count2 + 1
               DO i = 1 , Nonods
! IX will be BLANK   IX(I)=0
                  Blank(i) = 0
                  IF ( Whichd(i).EQ.sub ) Blank(i) = 1
               ENDDO
               DO i = 1 , Nonods
! IL will be Q                 IL(I)=0
                  Q(i) = 0
                  DO count = Fina(i) , Fina(i+1) - 1
                     Q(i) = Q(i) + Blank(Cola(count))
                  ENDDO
               ENDDO
! see which subdomains to connect it with.
               DO i = 1 , Nonods
                  IF ( Q(i).GE.1 ) THEN
                     sub2 = Whichd(i)
! Look to see if we have sub2 connected to subdomain SUB already.
                     found = .FALSE.
                     DO count = Finsub(sub) , count2
                        IF ( Colsub(count).EQ.sub2 ) found = .TRUE.
                     ENDDO
! Put into pointers?
                     IF ( .NOT.found ) THEN
                        count2 = count2 + 1
                        Colsub(count2) = sub2
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
            Finsub(Ndom+1) = count2 + 1
            ncolsu = count2
!             CALL SAVSUB(NDOM,NCOLSU,FINSUB,COLSUB)
! Re-order subdomains.
! RESUB(OLDSUB)=NEWSUB
            IF ( Defaul ) THEN
               zero = .TRUE.
               iiqn = MAX(500,Nonods/2)
!     &      Q,Q(NONODS/2 +1),IIQN,ZERO)
               CALL REORD2(Blasub,Finsub,Colsub,Ndom,ncolsu,Renum,Q(500)&
                         & ,Q(Nonods/2+501),iiqn,zero)
            ELSE
               imem2 = Ndom*Ndom + 1
               iedist = 1
               nimem2 = Nimem - imem2 + 1
               CALL MAPROC(Hypdim,Topcon,Sidcon,Proacr,Prodow,          &
                         & Imem(imem2),nimem2,Finsub,Colsub,Imem(iedist)&
                         & ,Ndom,Ndom,ncolsu,Blasub,Renum)
!     &                 NCOLDO,WICNOD,WICPRO)
            ENDIF
            DO i = 1 , Nonods
               newsub = Renum(Whichd(i))
               Whichd(i) = newsub
            ENDDO
         ENDIF
! FINISHED REORDERING SUBDOMAINS ********
! ***************************************
!
! reorder equations within each subdomain.
         zero = .TRUE.
         insub = 0
         DO sub = 1 , Ndom
            stasub = insub
            DO i = 1 , Nonods
               Blank(i) = 0
               IF ( Whichd(i).EQ.sub ) THEN
                  Blank(i) = 1
                  insub = insub + 1
               ENDIF
            ENDDO
            iiqn = MAX(500,Nonods/2)
            CALL REORD2(Blank,Fina,Cola,Nonods,Ncola,Renum,Q(500),      &
                      & Q(Nonods/2+501),iiqn,zero)
            zero = .FALSE.
         ENDDO
!
         IF ( nintnd.GT.1 ) THEN
! now reorder interface nodes.
            insub = 0
            DO sub = -Ndom , -1 , 1
               stasub = insub
               found = .FALSE.
               DO i = 1 , Nonods
                  Blank(i) = 0
                  IF ( Whichd(i).EQ.sub ) THEN
                     Blank(i) = 1
                     insub = insub + 1
                     found = .TRUE.
                  ENDIF
               ENDDO
               IF ( found ) THEN
                  iiqn = MAX(500,Nonods/2)
                  CALL REORD2(Blank,Fina,Cola,Nonods,Ncola,Renum,Q(500),&
                            & Q(Nonods/2+501),iiqn,zero)
               ENDIF
               zero = .FALSE.
            ENDDO
         ENDIF
      ENDIF
      END
!*==SAVSUB.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
      SUBROUTINE SAVSUB(Ndom,Ncolsu,Finsub,Colsub)
      IMPLICIT NONE
!*--SAVSUB3510
!*** Start of declarations inserted by SPAG
      INTEGER i
!*** End of declarations inserted by SPAG
      INTEGER Ndom , Ncolsu
      INTEGER Finsub(Ndom+1) , Colsub(Ncolsu)
 
     ! OPEN (2,FILE='rubish',STATUS='UNKNOWN')
      !WRITE (2,*) Ndom
      !WRITE (2,*) Ncolsu
      !WRITE (2,*) (Finsub(i),i=1,Ndom+1)
      !WRITE (2,*) (Colsub(i),i=1,Ncolsu)
     ! CLOSE (2)
      END
!*==REORD2.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
      SUBROUTINE REORD2(Blank,Fina,Cola,Nonods,Ncola,Renum2,Q,Qtemp,    &
                      & Iiqn,Zero)
      IMPLICIT NONE
!*--REORD23532
!*** Start of declarations inserted by SPAG
      INTEGER i , iav , icol , ifront , ii , iii , Iiqn , imxbwd ,      &
            & ival , j , nod , nodq , nq , nqtemp
!*** End of declarations inserted by SPAG
! This subroutine mimimizes the bandwidth for a direct solution
! and maximizes the connectivity for an iterative solution.
! IT works by forming fronts, the nodes contained in which
! are renumbered.
! The nodes with the minimum mean of the surrounding
! numbered nodes, the front nods can be taken as
! negative or egnored.
! RENUM2(OLDNOD no)=new nod no.
! It will begin renumbering from the maximum entry in RENUM2.
      LOGICAL QUICK
      PARAMETER (QUICK=.TRUE.)
!
      INTEGER Nonods , Ncola
      INTEGER Q(Iiqn) , Qtemp(Iiqn)
      INTEGER Blank(Nonods) , Fina(Nonods+1)
      INTEGER Cola(Ncola)
      INTEGER Renum2(Nonods)
!
      REAL minodn , mean , froval
      INTEGER mxval , nval , mxnod , mxrenu
      INTEGER minval , minnod
      INTEGER nodgot , count , row
! Q,QTEMP contains the work spase.
      INTEGER LARGE
      PARAMETER (LARGE=10000000)
      LOGICAL found , Zero
! This sub reorders the equations so an iterative OR DIRECT solver
! can solve the equations efficiently.
! When BLANK(I)=0 node I is egnored (or renumbered).
! When BLANK(I)=1 node I is not egnored(or renumbered).
! if ZERO then
      IF ( Zero ) THEN
         DO i = 1 , Nonods
            Renum2(i) = 0
         ENDDO
      ENDIF
!
! Find MXRENU = max node number in RENUM2
 100  mxrenu = 0
      DO i = 1 , Nonods
         mxrenu = MAX(Renum2(i),mxrenu)
      ENDDO
!
      nod = mxrenu + 1
!
      minval = LARGE
      DO i = 1 , Nonods
         IF ( Blank(i).NE.0 ) THEN
            IF ( Renum2(i).EQ.0 ) THEN
               nval = 0
               DO count = Fina(i) , Fina(i+1) - 1
                  icol = Cola(count)
                  IF ( Blank(icol).NE.0 ) THEN
                     IF ( Renum2(icol).EQ.0 ) nval = nval + 1
                  ENDIF
               ENDDO
               IF ( nval.LT.minval ) THEN
                  minval = nval
                  minnod = i
               ENDIF
            ENDIF
         ENDIF
      ENDDO
!
!         RENUM(NOD)=MXNOD
      Renum2(minnod) = nod
!
 
      nq = 1
      Q(1) = minnod
      DO ifront = 1 , 100000
! Find another front Q.
! find all the nodes(put in Q) that QTEMP is
! connected to that are not numbered.
! and not already in Q.
         DO i = 1 , nq
            Qtemp(i) = Q(i)
         ENDDO
         nqtemp = nq
         nq = 0
         DO i = 1 , nqtemp
            row = Qtemp(i)
            DO count = Fina(row) , Fina(row+1) - 1
               icol = Cola(count)
! see if already in Q.
               IF ( Blank(icol).NE.0 ) THEN
                  IF ( Renum2(icol).EQ.0 ) THEN
                     found = .FALSE.
                     DO ii = 1 , nq
                        IF ( icol.EQ.Q(ii) ) found = .TRUE.
                     ENDDO
                     IF ( .NOT.found ) THEN
                        nq = nq + 1
                        Q(nq) = icol
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
! GIVE NOD NOS TO Q
         IF ( (ifront.EQ.1) .AND. (nq.GT.2) ) THEN
! Creat a spacial ordering. *************
! Find node in Q with max valancy
            mxval = -1
            DO j = 1 , nq
               i = Q(j)
               IF ( Blank(i).NE.0 ) THEN
                  IF ( Renum2(i).EQ.0 ) THEN
                     nval = 0
                     DO count = Fina(i) , Fina(i+1) - 1
                        icol = Cola(count)
                        IF ( Blank(icol).NE.0 ) THEN
                           IF ( Renum2(icol).EQ.0 ) nval = nval + 1
                        ENDIF
                     ENDDO
                     IF ( nval.GT.mxval ) THEN
                        mxval = nval
                        mxnod = i
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
!
            iav = INT(REAL(nq)/2.+0.6)
! Put to -1 for loops 910,920.
            Renum2(mxnod) = -1
            iii = 0
            DO i = 1 , nq
               IF ( iii.LT.iav-1 ) THEN
                  IF ( Renum2(Q(i)).EQ.0 ) THEN
                     iii = iii + 1
                     nod = nod + 1
                     Renum2(Q(i)) = -nod
                  ENDIF
               ENDIF
            ENDDO
            nod = nod + 1
            Renum2(mxnod) = -nod
            DO i = 1 , nq
               IF ( Renum2(Q(i)).EQ.0 ) THEN
                  nod = nod + 1
                  Renum2(Q(i)) = -nod
               ENDIF
            ENDDO
! ***************************************
         ELSE
            DO j = 1 , nq
               IF ( .NOT.QUICK ) THEN
                  minodn = LARGE
                  DO i = 1 , nq
! find node not NUMBERED connected to the minimum of the
! MEAN of the surrounding nodes.
! ALSO TAKE INTO ACCOUNT THE VALANCY.
                     nodq = Q(i)
                     IF ( Blank(nodq).NE.0 ) THEN
                        IF ( Renum2(nodq).EQ.0 ) THEN
                           ival = 0
                           mean = 0.
                           froval = 0.
                           DO count = Fina(nodq) , Fina(nodq+1) - 1
                              icol = Cola(count)
                              IF ( Blank(icol).EQ.1 ) THEN
!                     IF(RENUM2(ICOL).NE.0) THEN
                                 IF ( Renum2(icol).GT.0 ) THEN
                                    ival = ival + 1
                                    mean = mean + Renum2(icol)
! Front nodes count as a zero.
                                 ENDIF
                                 IF ( Renum2(icol).LT.0 )               &
                                    & froval = froval + 1.
                              ENDIF
                           ENDDO
                           IF ( ival.NE.0 ) THEN
                              mean = mean/REAL(ival) - 0.01*froval
                              IF ( mean.LT.minodn ) THEN
                                 minodn = mean
                                 nodgot = nodq
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDDO
               ELSE
                  nodgot = Q(j)
               ENDIF
!
               nod = nod + 1
               Renum2(nodgot) = -nod
            ENDDO
!
         ENDIF
!
         DO i = 1 , nq
            Renum2(Q(i)) = -Renum2(Q(i))
         ENDDO
         IF ( nq.EQ.0 ) THEN
! See if we have numbered all the nodes, if not
! then goto 9770 else goto 9700
            found = .FALSE.
            DO i = 1 , Nonods
               IF ( Blank(i).EQ.1 ) THEN
                  IF ( Renum2(i).EQ.0 ) found = .TRUE.
               ENDIF
            ENDDO
            IF ( found ) GOTO 100
            IF ( .NOT.found ) GOTO 200
         ENDIF
      ENDDO
!
! find maximum bandwidth.
 200  imxbwd = 0
      DO i = 1 , Nonods
         IF ( Blank(i).NE.0 ) THEN
            ii = Renum2(i)
            DO count = Fina(i) , Fina(i+1) - 1
               icol = Cola(count)
               IF ( Blank(icol).NE.0 )                                  &
                  & imxbwd = MAX(imxbwd,ABS(ii-Renum2(icol)))
            ENDDO
         ENDIF
      ENDDO
 
      END
!*==READOM.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
 
!
!
! **********************************************************
! **** THIS PART MAPS THE SUBDOMAINS TO PROCESSORS *********
! **********************************************************
!
!
!
      SUBROUTINE READOM(Ndom,Ncolsu,Finsub,Colsub)
      IMPLICIT NONE
!*--READOM3772
!*** Start of declarations inserted by SPAG
      INTEGER i , iincol , indom
!*** End of declarations inserted by SPAG
      INTEGER Ndom , Ncolsu
      INTEGER Finsub(Ndom+1) , Colsub(Ncolsu)
 
    !  OPEN (2,FILE='rubish',STATUS='UNKNOWN')
 !     READ (2,*) indom
      IF ( indom.NE.Ndom ) THEN
  !       WRITE (*,*) 'THERE IS AN ERROR'
         STOP
      ENDIF
   !   READ (2,*) iincol
    !  READ (2,*) (Finsub(i),i=1,indom+1)
     ! READ (2,*) (Colsub(i),i=1,iincol)
    !  CLOSE (2)
     ! WRITE (*,*) 'no of edges of domain to domain connectivity=' ,     &
     !           & iincol/2
      END
!*==MAPROC.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
      SUBROUTINE MAPROC(Hypdim,Topcon,Sidcon,Proacr,Prodow,Imem,Nimem,  &
                      & Findom,Coldom,Edist,Ndom,Nprocs,Ncoldo,Wicnod,  &
                      & Wicpro)
      IMPLICIT NONE
!*--MAPROC3801
!*** Start of declarations inserted by SPAG
      INTEGER Ncoldo , Ndom , Nimem , Nprocs
!*** End of declarations inserted by SPAG
! This sub finds the PROCESSOR TO SUBDOMAIN MAPPING.
! NB nonods=nprocs.
! If HYPDIM=0 assume not a hypercube but a net.
! TOPCON & SIDCON are for the net.
! PROACR,PRODOW are the dimensions of the net.
! If TOPCON then the top row of nodes are connected to the bottom.
! If SIDCON then the r.h.s column of nodes are connected to the l.h.s
! IF(hypdim.gt.0) NCOLA=no of edges of hypercube * 2. = nprocs*log(base 2)(nprocs).
! log(base 2)(nprocs)=HYPDIM.
! WICPRO(SUBDOM NO)=PROCESSOR NUMBER
! WICNOD(PROC NO)=SUBDOMAIN NUMBER.
      INTEGER Hypdim , Proacr , Prodow
      LOGICAL Topcon , Sidcon
      INTEGER Findom(Ndom+1) , Coldom(Ncoldo)
      INTEGER Edist(Nprocs*Nprocs)
      INTEGER Imem(Nimem)
      INTEGER Wicnod(Nprocs) , Wicpro(Ndom)
! Define the processor connectivity, put in EDIST matrix.
      CALL GLBPRO(Hypdim,Topcon,Sidcon,Proacr,Prodow,Nprocs,Edist,Imem, &
                & Nimem)
!
      IF ( Nimem.LT.Ndom*Ndom ) THEN
      !   WRITE (*,*) 'NOT ENOUGH MEM IN MAPROC'
         STOP
      ENDIF
! solve mapping problem.
      CALL MAPPRO(Imem,Findom,Coldom,Edist,Ndom,Nprocs,Ncoldo,Wicnod,   &
                & Wicpro)
 
      END
!*==GLBPRO.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
      SUBROUTINE GLBPRO(Hypdim,Topcon,Sidcon,Proacr,Prodow,Nprocs,      &
                      & Fulmat,Imem,Nimem)
      IMPLICIT NONE
!*--GLBPRO3843
!*** Start of declarations inserted by SPAG
      INTEGER need , Nimem , nz
!*** End of declarations inserted by SPAG
! This sub finds the global inter processor communication matrix.
! NB nonods=nprocs.
! If HYPDIM=0 assume not a hypercube but a net.
! TOPCON & SIDCON are for the net.
! PROACR,PRODOW are the dimensions of the net.
! If TOPCON then the top row of nodes are connected to the bottom.
! If SIDCON then the r.h.s column of nodes are
! connected to the l.h.s
      INTEGER Hypdim , Proacr , Prodow , Nprocs
      INTEGER Fulmat(Nprocs*Nprocs) , Imem(Nimem)
! NIMEM=MAX(NPROCS*HYPDIM,4*NPROCS)*2 +(NPROCS+1)*2
      INTEGER ncola , fina , cola , fina2 , cola2 , front , lfront
      LOGICAL Topcon , Sidcon
!
! IF(hypdim.gt.0) NCOLA=no of edges of hypercube * 2. = nprocs*log(base 2)(nprocs).
! log(base 2)(nprocs)=HYPDIM.
!           !WRITE(*,*)'inside glbpro'
      ncola = MAX(Nprocs*Hypdim,4*Nprocs)
      IF ( Hypdim.GT.0 ) ncola = Nprocs*Hypdim
! ALlocate space
      fina = 1
      cola = fina + Nprocs + 1
      fina2 = cola + ncola
      cola2 = fina2 + Nprocs + 1
      need = fina2 + Nprocs + ncola
      IF ( need.GT.Nimem ) THEN
         !WRITE (*,*) 'NOT ENOUGH MEMORY STOPED IN GLBPRO'
         STOP
      ENDIF
!
      IF ( Hypdim.GT.0 ) THEN
         CALL HYPPTR(Hypdim,Nprocs,Imem(fina),Imem(cola),Imem(fina2),   &
                   & Imem(cola2),ncola,nz)
      ELSE
!            write(*,*)NCOLA,NPROCS,PRODOW,PROACR,
!     &                     TOPCON,SIDCON,fina,cola,nimem
         CALL DEFDOM(Imem(fina),Imem(cola),ncola,Nprocs,Prodow,Proacr,  &
                   & Topcon,Sidcon,nz)
      ENDIF
!
      front = fina2
      lfront = front + Nprocs
!           write(*,*)'befor getful'
      CALL GETFUL(Fulmat,Imem(fina),Imem(cola),Nprocs,ncola,Imem(front),&
                & Imem(lfront),Nprocs)
!           write(*,*)'after getful'
      END
!*==GETFUL.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
      SUBROUTINE GETFUL(Fulmat,Fina,Cola,Nonods,Ncola,Front,Lfront,     &
                      & Mxnfro)
      IMPLICIT NONE
!*--GETFUL3902
!*** Start of declarations inserted by SPAG
      INTEGER ifr , inod , j , jnod , Ncola , nod , Nonods
!*** End of declarations inserted by SPAG
      INTEGER MXNFRT
      PARAMETER (MXNFRT=100000)
      INTEGER Mxnfro
      INTEGER count , Front(Mxnfro) , Lfront(Mxnfro)
      INTEGER nfront , lnfron , idist
      INTEGER Fina(Nonods+1) , Cola(Ncola)
      INTEGER Fulmat(Nonods*Nonods)
! MXNFRO=max number of nodes in a front, it can not exceed the
! maximum distance between two nodes.
! MXNFRT=max number of fronts.
! This sub finds a full matrix from the matrix in FINA,COLA.
! It works using fronts.  It assums that a node does not point to its self.
      DO inod = 1 , Nonods
! Work out the INOD'th row of the dense matrix.
         DO j = 1 , Nonods
            Fulmat((inod-1)*Nonods+j) = 0
         ENDDO
         Fulmat((inod-1)*Nonods+inod) = -1
! Initialize front
         nfront = 1
         Front(nfront) = inod
!
         DO idist = 1 , MXNFRT
            lnfron = nfront
            DO ifr = 1 , nfront
               Lfront(ifr) = Front(ifr)
            ENDDO
!
            nfront = 0
            DO ifr = 1 , lnfron
! LFRONT contains the previous front.
               jnod = Lfront(ifr)
               DO count = Fina(jnod) , Fina(jnod+1) - 1
                  nod = Cola(count)
                  IF ( Fulmat((inod-1)*Nonods+nod).EQ.0 ) THEN
                     nfront = nfront + 1
                     Front(nfront) = nod
                     Fulmat((inod-1)*Nonods+nod) = idist
                  ENDIF
               ENDDO
            ENDDO
            IF ( nfront.EQ.0 ) GOTO 50
         ENDDO
!
 50      Fulmat((inod-1)*Nonods+inod) = 0
      ENDDO
!
!         write(*,*)'FULMAT:'
!         DO 19 JNOD=1,NONODS
!           jnod=1
!          do 19 i=1,9
!           write(*,*)(FULMAT((JNOD-1)*NONODS+ (i-1)*9+j),j=1,9)
!19       CONTINUE
!         STOP
      END
!*==DEFDOM.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
      SUBROUTINE DEFDOM(Findom,Coldom,Ncoldo,Nonods,Nodown,Nodacr,      &
                      & Topcon,Sidcon,Nz)
      IMPLICIT NONE
!*--DEFDOM3968
!*** Start of declarations inserted by SPAG
      INTEGER iacr , jdow , Ncoldo , Nodacr , Nodown , Nonods , Nz
!*** End of declarations inserted by SPAG
      INTEGER Findom(Nonods+1) , Coldom(Ncoldo) , row , count
      LOGICAL Topcon , Sidcon
! This subroutine finds the pointers for cray T3D or domain connectivity.
! If TOPCON then the top row of nodes are connected to the bottom.
! If SIDCON then the r.h.s column of nodes are connected to the l.h.s
! column of nodes.
!
!         write(*,*)'inside defdom'
!         write(*,*)'NCOLDO,NONODS,NODOWN,NODACR,TOPCON,SIDCON:'
!         write(*,*)NCOLDO,NONODS,NODOWN,NODACR,TOPCON,SIDCON
      count = 0
      row = 0
      DO jdow = 1 , Nodown
         DO iacr = 1 , Nodacr
            row = row + 1
            Findom(row) = count + 1
!
            IF ( jdow.NE.1 ) THEN
               count = count + 1
               Coldom(count) = row - Nodacr
            ENDIF
            IF ( (jdow.EQ.1) .AND. Topcon ) THEN
               count = count + 1
               Coldom(count) = Nonods - Nodacr + iacr
            ENDIF
!
            IF ( iacr.NE.1 ) THEN
               count = count + 1
               Coldom(count) = row - 1
            ENDIF
            IF ( (iacr.EQ.1) .AND. Sidcon ) THEN
               count = count + 1
               Coldom(count) = row + Nodacr - 1
            ENDIF
! miss out the middle node.
!                COUNT=COUNT+1
!                COLDOM(COUNT)=ROW
!                MIDPA(ROW)=COUNT
!
            IF ( iacr.NE.Nodacr ) THEN
               count = count + 1
               Coldom(count) = row + 1
            ENDIF
            IF ( (iacr.EQ.Nodacr) .AND. Sidcon ) THEN
               count = count + 1
               Coldom(count) = row - Nodacr + 1
            ENDIF
!
            IF ( jdow.NE.Nodown ) THEN
               count = count + 1
               Coldom(count) = row + Nodacr
            ENDIF
            IF ( (jdow.EQ.Nodown) .AND. Topcon ) THEN
               count = count + 1
               Coldom(count) = row - (Nonods-Nodacr)
            ENDIF
!
         ENDDO
      ENDDO
      Findom(Nonods+1) = count + 1
      Nz = count
!         write(*,*)'going out of defdom'
      END
!*==DEFPRO.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
      SUBROUTINE DEFPRO(Edist,Nprocs,Nacr,Ndow)
      IMPLICIT NONE
!*--DEFPRO4041
!*** Start of declarations inserted by SPAG
      INTEGER ii , ij , iproc , ji , jj , jproc , Nacr , Ndow , Nprocs
!*** End of declarations inserted by SPAG
      INTEGER Edist(Nprocs*Nprocs)
      DO ji = 1 , Ndow
         DO ii = 1 , Nacr
            iproc = (ji-1)*Nacr + ii
            DO jj = 1 , Ndow
               DO ij = 1 , Nacr
                  jproc = (jj-1)*Nacr + ij
                  Edist((iproc-1)*Nprocs+jproc) = ABS(ji-jj)            &
                   & + ABS(ii-ij)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      END
!*==HYPPTR.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
      SUBROUTINE HYPPTR(Hypdim,Nprocs,Finpro,Colpro,Finloc,Colloc,      &
                      & Ncolpr,Nz)
      IMPLICIT NONE
!*--HYPPTR4066
!*** Start of declarations inserted by SPAG
      INTEGER i , ilev , Ncolpr , nploc , npro , Nprocs , Nz
!*** End of declarations inserted by SPAG
      INTEGER count , Hypdim
      INTEGER Finpro(Nprocs+1) , Colpro(Ncolpr)
      INTEGER Finloc(Nprocs+1) , Colloc(Ncolpr)
! This sub finds the pointers for the HYPERCUBE
! HYPDIM=dimensions of the hypercube NB NPROCS=2^HYPDIM.
! The connectivity is contained in the pointers FINPRO,COLPRO,NZ
! NZ=length of COLPRO.
!          write(*,*)'inside hypptr'
      npro = 2
      Finpro(1) = 1
      Finpro(2) = 2
      Finpro(3) = 3
      Colpro(1) = 2
      Colpro(2) = 1
      Ncolpr = 2
!
      DO ilev = 2 , Hypdim
         nploc = npro
         DO i = 1 , npro
            Finloc(i) = Finpro(i)
            DO count = Finpro(i) , Finpro(i+1) - 1
               Colloc(count) = Colpro(count)
            ENDDO
         ENDDO
         Finloc(npro+1) = Finpro(npro+1)
 
         npro = 2*nploc
         CALL NEXTLE(Finpro,Colpro,Ncolpr,npro,Finloc,Colloc,Ncolpr,    &
                   & nploc,Nz)
      ENDDO
!
      IF ( npro.NE.Nprocs ) THEN
         !WRITE (*,*) 'PROBLEM WITH THE HYPERCUBE'
         STOP
      ENDIF
!          write(*,*)'going out of hypptr'
!          stop
      END
!*==NEXTLE.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
      SUBROUTINE NEXTLE(Finpro,Colpro,Ncolpr,Nprocs,Finloc,Colloc,      &
                      & Ncollo,Nploc,Nz)
      IMPLICIT NONE
!*--NEXTLE4115
!*** Start of declarations inserted by SPAG
      INTEGER i , ii , Ncollo , Ncolpr , Nploc , Nprocs , Nz
!*** End of declarations inserted by SPAG
      INTEGER count , count2
      INTEGER Finpro(Nprocs+1) , Colpro(Ncolpr)
      INTEGER Finloc(Nploc+1) , Colloc(Ncollo)
! This sub is called recursively to build up the links on the HYPERCUBE.
      count2 = 0
      DO i = 1 , Nploc
         Finpro(i) = count2 + 1
         DO count = Finloc(i) , Finloc(i+1) - 1
            count2 = count2 + 1
            Colpro(count2) = Colloc(count)
         ENDDO
         count2 = count2 + 1
         Colpro(count2) = i + Nploc
      ENDDO
!
      DO ii = 1 , Nploc
         i = ii + Nploc
         Finpro(i) = count2 + 1
!
         count2 = count2 + 1
         Colpro(count2) = ii
         DO count = Finloc(i) , Finloc(i+1) - 1
            count2 = count2 + 1
            Colpro(count2) = Colloc(count) + Nploc
         ENDDO
      ENDDO
      Nz = count2
      Finpro(Nprocs+1) = count2 + 1
      END
!*==MAPPRO.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
      SUBROUTINE MAPPRO(Ax,Findom,Coldom,Edist,Nonods,Nprocs,Ncoldo,    &
                      & Wicnod,Wicpro)
      IMPLICIT NONE
!*--MAPPRO4155
!*** Start of declarations inserted by SPAG
      INTEGER i , icurf , its , j , nchain , Ncoldo , nexcp , NOITS ,   &
            & Nonods , Nprocs
!*** End of declarations inserted by SPAG
      REAL ALPHA , BETA
      INTEGER CHEND
!         PARAMETER(NOITS=10000,ALPHA=0.995,BETA=0.01)
      PARAMETER (NOITS=10000,ALPHA=0.98,BETA=0.01)
      PARAMETER (CHEND=5)
      REAL c , pcost , ncost
      INTEGER curf , newf , lastf(CHEND) , ikeep(CHEND)
!         INTEGER X(NONODS*NPROCS)
      INTEGER Ax(Nonods*Nprocs)
      INTEGER Findom(Nonods+1) , Coldom(Ncoldo)
      INTEGER Edist(Nprocs*Nprocs)
      INTEGER length , rnod , rproc
      INTEGER Wicnod(Nprocs) , Wicpro(Nonods)
      INTEGER proc , chain , mxexce
      LOGICAL yes , same
      length = Nonods*Nprocs
! Initialize X & WICNOD,WICPRO
!         DO 10 I=1,LENGTH
!           X(I)=0
!10       CONTINUE
      DO proc = 1 , Nprocs
!           X(PROC+(PROC-1)*NONONDS)=1
         Wicnod(proc) = proc
         Wicpro(proc) = proc
      ENDDO
!
      CALL FINFAX(curf,Ax,Findom,Coldom,Ncoldo,Nonods,Nprocs,Edist,     &
                & Wicnod)
      !WRITE (*,*) 'befor we randomise CURF=' , curf
!
      nchain = 10 + 1*Nonods*Nprocs
! NCHAIN=maximum length of a markolf chain.
      mxexce = Nonods + INT(BETA*REAL(nchain))
!
! This sub finds the initial C.
      CALL GETC(c,curf,nchain,mxexce,Ax,Findom,Coldom,Edist,Nonods,     &
              & Nprocs,Ncoldo,Wicnod,Wicpro)
      !WRITE (*,*) 'initial curf=' , curf
      !WRITE (*,*) 'initia C=' , c
!
      DO i = 1 , CHEND
         lastf(i) = 0
      ENDDO
!
      DO its = 1 , NOITS
         nexcp = 0
         DO chain = 1 , nchain
!             write(*,*)'CURF=',CURF
! Generate a nabourhood configuration ** make sure it is different from old one.
            CALL GENNAB(rproc,rnod,Nprocs,Nonods,Wicnod)
! Find new cost of this nabourhood configuration.
            CALL FINEWF(newf,curf,rnod,rproc,Wicnod,Wicpro,Nprocs,      &
                      & Nonods,Edist,Findom,Coldom,Ncoldo,Ax)
! See if we except this config.
            pcost = curf
            ncost = newf
            CALL EXCEPT(ncost,pcost,c,yes)
            IF ( yes ) THEN
               nexcp = nexcp + 1
! Update everything Change X,AX.
               CALL UPXAX(rnod,rproc,Wicnod,Wicpro,Nonods,Nprocs,curf,  &
                        & newf,Ax,Edist,Findom,Coldom,Ncoldo)
               IF ( nexcp.GT.mxexce ) GOTO 50
            ENDIF
         ENDDO
 50      continue!    WRITE (*,*) 'ITS,CURF,NEXCP,CHAIN:' , its , curf , nexcp ,     &
      !             & chain
!
         c = ALPHA*c
!
         IF ( its.GT.30 ) THEN
! Put CURF into LASTF and shift LASTF down.
            DO i = 1 , CHEND
               ikeep(i) = lastf(i)
            ENDDO
            lastf(CHEND) = curf
            DO i = 1 , CHEND - 1
               lastf(i) = ikeep(i+1)
            ENDDO
! See if the SAME
            same = .TRUE.
            DO i = 1 , CHEND
               DO j = 1 , CHEND
                  IF ( lastf(i).NE.lastf(j) ) same = .FALSE.
               ENDDO
            ENDDO
!
! Check for convergence Ie) see if all the entries in LASTF are the SAME.
            IF ( same ) GOTO 100
         ENDIF
!
      ENDDO
!
! ********************* OUTPUT THE RESULTS **********************
!
 100  CALL FINFAX(icurf,Ax,Findom,Coldom,Ncoldo,Nonods,Nprocs,Edist,    &
                & Wicnod)
      IF ( icurf.NE.curf ) THEN
       !  WRITE (*,*) '********* THERE WAS A PROBLEM: ICURF,CURF' ,      &
        !           & icurf , curf
         curf = icurf
      ENDIF
!
      !WRITE (*,*) 'IN THE END CURF=' , curf
      DO i = 1 , -Nonods
         !WRITE (*,*) 'NOD, WICPRO(NOD):' , i , Wicpro(i)
      ENDDO
      DO i = 1 , -Nprocs
         !WRITE (*,*) 'PROC, WICNOD(PROC):' , i , Wicnod(i)
      ENDDO
      END
!*==GETC.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
      SUBROUTINE GETC(C,Curf,Nchain,Mxexce,Ax,Findom,Coldom,Edist,      &
                    & Nonods,Nprocs,Ncoldo,Wicnod,Wicpro)
      IMPLICIT NONE
!*--GETC4278
!*** Start of declarations inserted by SPAG
      INTEGER Nchain , Ncoldo , ncost , Nonods , Nprocs , ntry
      REAL pcost , rrlog
!*** End of declarations inserted by SPAG
      REAL C , mean
      INTEGER Curf , newf
!         INTEGER X(NONODS*NPROCS)
      INTEGER Ax(Nonods*Nprocs)
      INTEGER Findom(Nonods+1) , Coldom(Ncoldo)
      INTEGER Edist(Nprocs*Nprocs)
      INTEGER rnod , rproc
      INTEGER Wicnod(Nprocs) , Wicpro(Nonods)
      INTEGER chain , Mxexce , count
! This subroutine finds the initial value of the temperature like parameter C
!         NTRY=MXEXCE
      ntry = Nchain
      mean = 0.
      count = 0
      DO chain = 1 , ntry
! Generate a nabourhood configuration ** make sure it is different from old one.
         CALL GENNAB(rproc,rnod,Nprocs,Nonods,Wicnod)
! Find new cost of this nabourhood configuration.
         CALL FINEWF(newf,Curf,rnod,rproc,Wicnod,Wicpro,Nprocs,Nonods,  &
                   & Edist,Findom,Coldom,Ncoldo,Ax)
! See if we except this config.
         pcost = Curf
         ncost = newf
         IF ( newf.GT.Curf ) THEN
            mean = mean + REAL(newf-Curf)
            count = count + 1
         ENDIF
! Update everything Change X,AX.
         CALL UPXAX(rnod,rproc,Wicnod,Wicpro,Nonods,Nprocs,Curf,newf,Ax,&
                  & Edist,Findom,Coldom,Ncoldo)
      ENDDO
!           write(*,*)'got to here'
      mean = mean/REAL(count)
      rrlog = LOG(1./0.8)
      C = mean/rrlog
!           write(*,*)'C=',C,' RRLOG=',RRLOG
!           STOP
      END
!*==FINFAX.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
 
!
!
!
      SUBROUTINE FINFAX(Curf,Ax,Findom,Coldom,Ncoldo,Nonods,Nprocs,     &
                      & Edist,Wicnod)
      IMPLICIT NONE
!*--FINFAX4329
!*** Start of declarations inserted by SPAG
      INTEGER i , iproc , jproc , Ncoldo , Nonods , Nprocs
!*** End of declarations inserted by SPAG
      INTEGER Curf , Ax(Nonods*Nprocs)
      INTEGER Findom(Nonods+1) , Coldom(Ncoldo)
      INTEGER Edist(Nprocs*Nprocs)
      INTEGER edis , count , col
      INTEGER Wicnod(Nprocs)
! This sub finds the functional CURF and the vector AX
! NB FINDOM & COLDOM do not include the middle value.
      DO i = 1 , Nonods*Nprocs
         Ax(i) = 0
      ENDDO
      DO iproc = 1 , Nprocs
         DO jproc = 1 , Nprocs
            edis = Edist((iproc-1)*Nprocs+jproc)
!
            i = Wicnod(jproc)
            DO count = Findom(i) , Findom(i+1) - 1
               col = Coldom(count)
               Ax(col+(iproc-1)*Nonods) = Ax(col+(iproc-1)*Nonods)      &
                & + edis
            ENDDO
!
         ENDDO
      ENDDO
!
! Work out CURF
      Curf = 0
      DO iproc = 1 , Nprocs
         i = Wicnod(iproc)
         Curf = Curf + Ax(i+(iproc-1)*Nonods)
      ENDDO
      END
!*==GENNAB.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
      SUBROUTINE GENNAB(Rproc,Rnod,Nprocs,Nonods,Wicnod)
      IMPLICIT NONE
!*--GENNAB4370
!*** Start of declarations inserted by SPAG
      INTEGER iarg
      REAL RAN1
!*** End of declarations inserted by SPAG
      INTEGER Rproc , Rnod , Nprocs , Nonods , Wicnod(Nprocs)
      INTEGER iran
      REAL rran
!
 100  iarg = 3
      rran = RAN1(iarg)
      IF ( rran.GE.1 ) rran = 0.999999
      iran = INT(rran*REAL(Nprocs*Nonods)+1.)
!
      Rproc = INT((iran-1)/Nonods) + 1
      Rnod = iran - (Rproc-1)*Nonods
      IF ( Wicnod(Rproc).EQ.Rnod ) GOTO 100
      END
!*==FINEWF.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
      SUBROUTINE FINEWF(Newf,Curf,Rnod,Rproc,Wicnod,Wicpro,Nprocs,      &
                      & Nonods,Edist,Findom,Coldom,Ncoldo,Ax)
      IMPLICIT NONE
!*--FINEWF4395
!*** Start of declarations inserted by SPAG
      INTEGER i , Ncoldo , Nonods , Nprocs
!*** End of declarations inserted by SPAG
! Find new cost of this nabourhood configuration.
! WICNOD(PROC) returns the domain or 'nod' that processor proc contains.
! WICPRO(NOD) returns the processor that domain NOD is on.
      INTEGER Newf , Curf , Rnod , Rproc , Wicnod(Nprocs) ,             &
            & Wicpro(Nonods) , e1
      INTEGER count , col , knod , kproc , exf
      INTEGER Edist(Nprocs*Nprocs)
      INTEGER Findom(Nonods+1) , Coldom(Ncoldo)
      INTEGER Ax(Nonods*Nprocs)
      LOGICAL rr , rk , kr , kk
!
      knod = Wicnod(Rproc)
      kproc = Wicpro(Rnod)
!
      exf = 0
      exf = exf + Ax(Rnod+(Rproc-1)*Nonods)
      exf = exf - Ax(knod+(Rproc-1)*Nonods)
      exf = exf - Ax(Rnod+(kproc-1)*Nonods)
      exf = exf + Ax(knod+(kproc-1)*Nonods)
!
      exf = 2*exf
!
! Now for r^T A r.
      rr = .FALSE.
      rk = .FALSE.
      kr = .FALSE.
      kk = .FALSE.
      i = Rnod
      DO count = Findom(i) , Findom(i+1) - 1
         col = Coldom(count)
         IF ( col.EQ.Rnod ) rr = .TRUE.
         IF ( col.EQ.knod ) rk = .TRUE.
      ENDDO
      i = knod
      DO count = Findom(i) , Findom(i+1) - 1
         col = Coldom(count)
         IF ( col.EQ.Rnod ) kr = .TRUE.
         IF ( col.EQ.knod ) kk = .TRUE.
      ENDDO
!
      e1 = Edist((Rproc-1)*Nprocs+kproc)
! The first entry
      IF ( rr ) exf = exf - e1
      IF ( rk ) exf = exf + e1
! The 2ND entry
      IF ( kr ) exf = exf - (-e1)
      IF ( kk ) exf = exf - (e1)
! The 3RD entry
      IF ( rr ) exf = exf - (e1)
      IF ( rk ) exf = exf - (-e1)
! The 4TH entry
      IF ( kr ) exf = exf + e1
      IF ( kk ) exf = exf - e1
!
      Newf = Curf + exf
      END
!*==UPXAX.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
      SUBROUTINE UPXAX(Rnod,Rproc,Wicnod,Wicpro,Nonods,Nprocs,Curf,Newf,&
                     & Ax,Edist,Findom,Coldom,Ncoldo)
      IMPLICIT NONE
!*--UPXAX4462
!*** Start of declarations inserted by SPAG
      INTEGER i , iproc , Ncoldo , Nonods , Nprocs
!*** End of declarations inserted by SPAG
      INTEGER Rnod , Rproc , knod , kproc
      INTEGER Wicnod(Nprocs) , Wicpro(Nonods)
      INTEGER Ax(Nprocs*Nonods)
      INTEGER Edist(Nprocs*Nprocs)
      INTEGER count , col , Curf , Newf , eir , eik
      INTEGER Findom(Nonods+1) , Coldom(Ncoldo)
! Update AX as well
! WICNOD(PROC) returns the domain or 'nod' that processor proc contains.
! WICPRO(NOD) returns the processor that domain NOD is on.
!
      knod = Wicnod(Rproc)
      kproc = Wicpro(Rnod)
!
      Curf = Newf
!
      i = Rnod
      DO count = Findom(i) , Findom(i+1) - 1
         DO iproc = 1 , Nprocs
            eir = Edist((Rproc-1)*Nprocs+iproc)
            eik = Edist((kproc-1)*Nprocs+iproc)
            col = Coldom(count) + (iproc-1)*Nonods
            Ax(col) = Ax(col) + eir - eik
         ENDDO
      ENDDO
!
      i = knod
      DO count = Findom(i) , Findom(i+1) - 1
         DO iproc = 1 , Nprocs
            eir = Edist((Rproc-1)*Nprocs+iproc)
            eik = Edist((kproc-1)*Nprocs+iproc)
            col = Coldom(count) + (iproc-1)*Nonods
            Ax(col) = Ax(col) + eik - eir
         ENDDO
      ENDDO
!
!         I=RNOD+(RPROC-1)*NONODS
!         X(I)=X(I)+1
!         I=KNOD+(RPROC-1)*NONODS
!         X(I)=X(I)-1
!         I=RNOD+(KPROC-1)*NONODS
!         X(I)=X(I)-1
!         I=KNOD+(KPROC-1)*NONODS
!         X(I)=X(I)+1
!
      Wicnod(Rproc) = Rnod
      Wicpro(Rnod) = Rproc
      Wicnod(kproc) = knod
      Wicpro(knod) = kproc
      END
!*==UPXAX2.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
!
      SUBROUTINE UPXAX2(Rnod,Rproc,Wicnod,Wicpro,Nonods,Nprocs,Curf,    &
                      & Newf,Ax,Edist,Findom,Coldom,Ncoldo)
      IMPLICIT NONE
!*--UPXAX24523
!*** Start of declarations inserted by SPAG
      INTEGER i , idisp , iproc , Ncoldo , Nonods , Nprocs
!*** End of declarations inserted by SPAG
      INTEGER Rnod , Rproc , knod , kproc
      INTEGER Wicnod(Nprocs) , Wicpro(Nonods)
      INTEGER Ax(Nprocs*Nonods)
      INTEGER Edist(Nprocs*Nprocs)
      INTEGER count , col , Curf , Newf , eir , eik
      INTEGER Findom(Nonods+1) , Coldom(Ncoldo)
! Update AX as well
! WICNOD(PROC) returns the domain or 'nod' that processor proc contains.
! WICPRO(NOD) returns the processor that domain NOD is on.
!
      knod = Wicnod(Rproc)
      kproc = Wicpro(Rnod)
!
      Curf = Newf
!
      DO iproc = 1 , Nprocs
         eir = Edist((iproc-1)*Nprocs+Rproc)
         eik = Edist((iproc-1)*Nprocs+kproc)
         idisp = (iproc-1)*Nonods
         i = Rnod
         DO count = Findom(i) , Findom(i+1) - 1
            col = Coldom(count) + idisp
            Ax(col) = Ax(col) + eir - eik
         ENDDO
         i = knod
         DO count = Findom(i) , Findom(i+1) - 1
            col = Coldom(count) + idisp
            Ax(col) = Ax(col) + eik - eir
         ENDDO
      ENDDO
!
!         I=RNOD+(RPROC-1)*NONODS
!         X(I)=X(I)+1
!         I=KNOD+(RPROC-1)*NONODS
!         X(I)=X(I)-1
!         I=RNOD+(KPROC-1)*NONODS
!         X(I)=X(I)-1
!         I=KNOD+(KPROC-1)*NONODS
!         X(I)=X(I)+1
!
      Wicnod(Rproc) = Rnod
      Wicpro(Rnod) = Rproc
      Wicnod(kproc) = knod
      Wicpro(knod) = kproc
      END
!*==EXCEPT.spg  processed by SPAG 6.72Dc at 17:05 on 27 Apr 2019
!
!
!
      SUBROUTINE EXCEPT(Pcost,Ncost,C,Yes)
      IMPLICIT NONE
!*--EXCEPT4578
!*** Start of declarations inserted by SPAG
      REAL RAN1
!*** End of declarations inserted by SPAG
!
      REAL C , Pcost , Ncost , aprob , ran
      LOGICAL Yes
      IF ( Ncost.GE.Pcost ) THEN
         Yes = .TRUE.
         aprob = 9999.
      ELSEIF ( ABS(C).LT.0.00000000000000001 ) THEN
         Yes = .FALSE.
         aprob = -9991.
      ELSE
         IF ( (Pcost-Ncost)/C.LT.40 ) THEN
            aprob = 1/EXP((Pcost-Ncost)/C)
         ELSE
            aprob = 0.
         ENDIF
         IF ( aprob.GE.1 ) Yes = .TRUE.
         IF ( aprob.LT.1 ) THEN
            ran = RAN1(3)
!                write(*,*)'RAN=',RAN,' APROB=',APROB
!     &               ,' pcost=',PCOST,' NCOST=',NCOST
! RAN is between 0 and 1.
            Yes = .FALSE.
            IF ( ran.LT.aprob ) Yes = .TRUE.
         ENDIF
      ENDIF
!          write(*,*)' APROB=',APROB
      END
!
!


      SUBROUTINE COLOR2( WHICHC,NCOLOR,Q,QTEMP,NNOD, &
     &     FINA,COLA,NA,ZERO,QUICK)
!      use FLDebug
      IMPLICIT NONE
!     
!     - This subroutine colours the equations and puts the results in WHICHC()
!     - The pointers are FINA,COLA
!     - The number of nodes is NNOD.
!     - NCOLOR = no of colours needed. 
!     - WHICHC  contains the colouring. 
!     - Q,QTEMP are working arrays.
!     - if not ZERO then some of the nodes have already been coloured. 
      LOGICAL QUICK
!     PARAMETER(QUICK=.TRUE.)
!     
      INTEGER NCOLOR,NNOD,NA
      INTEGER Q(NNOD), QTEMP(NNOD), WHICHC(NNOD)
      INTEGER FINA(NNOD+1),COLA(NA)
      LOGICAL ZERO
!
      INTEGER VAL, COUNT, COL, QCOLOR, ROW
      INTEGER MINVAL,MAXCON,MAXCOL
      LOGICAL DONE,FOUND
      INTEGER I,INCOLO,NODE,ICOUNT,NQ,ITS,J,NODQ,NUM,NQTEMP,ICOL,II
!
!
      IF(ZERO) THEN
        DO I=1,NNOD
          WHICHC(I)=0
        END DO
       ENDIF
       INCOLO=0
!
! - find node with min valency that is not already coloured.
!
!2100  CONTINUE
      NQ = 0
      DO WHILE ( NQ .EQ. 0 ) 
!
      MINVAL = 100000
      NODE=0
      ICOUNT=0
      DO I = 1, NNOD
       IF(WHICHC(I).EQ.0) THEN
         ICOUNT=ICOUNT+1
         VAL = FINA(I+1) - FINA(I)
         IF( VAL .LT. MINVAL ) THEN
            MINVAL = VAL
            NODE   = I
         END IF
        END IF
      END DO  ! DO I = 1, NNOD
!     ewrite(3,*)'node,minval,nnod:',node,minval,nnod, 
!    & ' NDS NOT COLED=',ICOUNT
!
      IF(NODE.EQ.0) THEN
!       ewrite(3,*)'HAVE COLOURED ALL THE NODES'
         NQ=-1
      ELSE ! IF(NODE.EQ.0) THEN
!        GOTO 9700
!      ENDIF
         NQ   = 1
!
         Q(1) = NODE
!
         DO ITS = 1, 100000
!        ewrite(3,*)'front no, nq,NDS COLD:',its,nq,INCOLO
!        ewrite(3,*)'q:',(q(iijj),iijj=1,nq)
!
! - colour Q
            DO J = 1, NQ
              MAXCON = -1
              IF(.NOT.QUICK) THEN
               DO I = 1, NQ
! - find node not coloured connected to max number of coloured nodes. 
                  NODQ = Q(I)
                  IF( WHICHC(NODQ) .EQ. 0 ) THEN
                     NUM  = 0
                     DO COUNT = FINA(NODQ),FINA(NODQ+1)-1
                        IF( WHICHC(COLA(COUNT)) .NE. 0 ) NUM = NUM+1
                     END DO ! DO COUNT = FINA(NODQ),FINA(NODQ+1)-1
!
                     IF( NUM .GT. MAXCON ) THEN
                        MAXCON = NUM
                        QCOLOR = NODQ
                     END IF
                  END IF
               END DO ! DO I = 1, NQ
             ELSE
               QCOLOR=Q(J)
             ENDIF
!
! - Colour node QCOLOR
! - find MAXCOL = max color number surrounding node.  
             MAXCOL = 0
             DO COUNT = FINA(QCOLOR),FINA(QCOLOR+1)-1
                MAXCOL = MAX( WHICHC( COLA(COUNT) ), MAXCOL )
             END DO ! DO COUNT = FINA(QCOLOR),FINA(QCOLOR+1)-1
!
             DONE = .FALSE.
             DO COL = 1, MAXCOL
!
! - look for color 
               IF( .NOT. DONE ) THEN
                  FOUND = .FALSE.
                  DO COUNT = FINA(QCOLOR),FINA(QCOLOR+1)-1
                   IF( COL .EQ. WHICHC(COLA(COUNT)) ) FOUND=.TRUE.
                  END DO ! DO COUNT = FINA(QCOLOR),FINA(QCOLOR+1)-1
!
                  IF( .NOT. FOUND ) THEN
                     DONE           = .TRUE.
                     WHICHC(QCOLOR) = COL
                  END IF
               END IF
             END DO ! DO COL = 1, MAXCOL
!
             IF( .NOT. DONE ) WHICHC(QCOLOR) = MAXCOL + 1
           END DO ! DO J = 1, NQ
           INCOLO=INCOLO+NQ
!
! - find Q
           DO I = 1, NQ
              QTEMP(I) = Q(I)
           END DO ! DO I = 1, NQ
!
           NQTEMP = NQ
!
! - find all the nodes that QTEMP is connected 
! to that are not coloured.
! - and not already in Q. 
!
           NQ = 0
           DO I = 1, NQTEMP
              ROW = QTEMP(I)
              DO COUNT = FINA(ROW), FINA(ROW+1)-1
                 ICOL = COLA(COUNT)
!
! - see if already in Q.
                 IF( WHICHC(ICOL) .EQ. 0 ) THEN
                    FOUND = .FALSE.
                    DO II = 1, NQ
                       IF( ICOL .EQ. Q(II) ) FOUND = .TRUE.
                    END DO ! DO II = 1, NQ
                    IF( .NOT. FOUND ) THEN
                       NQ    = NQ + 1
                       Q(NQ) = ICOL
                    END IF
                 END IF
              END DO ! DO COUNT = FINA(ROW), FINA(ROW+1)-1
           END DO ! DO I = 1, NQTEMP
!
!         IF( NQ .EQ. 0 ) GO TO 2100
           IF( NQ .EQ. 0 ) EXIT ! exit its loop 
!
        END DO ! DO ITS = 1, 100000
!
      ENDIF ! IF(NODE.EQ.0) THEN ELSE

      
      END DO ! DO WHILE ( NQ .EQ. 0 )
!
! 9700  CONTINUE
!
! - Now find NCOLOR
!
        NCOLOR = 0
        DO I = 1, NNOD
           NCOLOR = MAX( NCOLOR, WHICHC(I) )
        END DO ! DO I = 1, NNOD
!        ewrite(3,*) '   MAXIMUM NO. OF NODAL COLOURS IS : ', NCOLOR
        RETURN
        END SUBROUTINE COLOR2
!
!
!
!
    SUBROUTINE CLRRL( N, A, CONST )
!     ------------------------------
! - this subroutine clears a real array
!
!   -------------------------------
! - date last modified : 03/05/1998
!   -------------------------------
!
      IMPLICIT NONE
!
      INTEGER N, I
!
      REAL    A(N), CONST
!
      DO I = 1, N
!
         A(I) = CONST
!
      END DO
!
      RETURN
      END SUBROUTINE CLRRL

!
!
!
!



