	PROGRAMM NUM_SV

	INCLUDE 'parameters' 

	COMMON/BASE/ Q(3,LINK,Ngr),Y(3,LINK,Ngr),W(3,LINK,Ngr),
     *	Z(3,LINK,Ngr),
     *  BW(3,LINK,Ngr),H(3,LINK,Ngr),YB(LINK,Ngr),FK(2,Ngr),FKD(2,Ngr)
     
       Common/K/ C(2,LINK,Ngr),CS(2,LINK,Ngr),CB(2,LINK,Ngr),
     *            S(2,LINK,Ngr),Seqv(LINK,Ngr)
c S(*,*,*) - suspended concentration, Seqv(*,*) - equlibrium conc.
        Common/pch/ flam,a12,a21,a13,a31,fKs,fKb,ros,ro,epsb
        Common/pchm/ d_part,dmu,dnu,z_obm,ER,grav,Zobm,delvv
	COMMON/GRANC/C_GRAN(LINK),CS_GRAN(LINK),S_GRAN(LINK)
	  COMMON/CGRN/ C_GR(LINK,500),TGR_C(LINK,500),IT_C(LINK)
	  COMMON/SGRN/ S_GR(LINK,500),TGR_S(LINK,500),IT_S(LINK)
	  COMMON/CSGRN/ CS_GR(LINK,500),TGR_CS(LINK,500),IT_CS(LINK)
	
        COMMON/GRN/ QGR(LINK,500),TGR(LINK,500),IT(LINK)
	COMMON/GRNY/  YGR(LINK,500),TGRY(LINK,500),ITY(LINK)
	COMMON/GRAN/QGRAN(LINK),QGRNO(node)
	COMMON/GRANY/YGRAN(LINK),YGRNO(node)
c  QGRAN,YGRAN - gr.yslovija po brancham, QGRNO,YGRNO - gr.ysl.v yzlax
C   Q - discharge, Y - level,W - crossection, BW - width, H - depth
c  (time,link,point) 
c  time: 1 - n, 2 - n+1, 3 - n-1 
C   YB(link,point) - bottom level
C   FK - rashodnaja funkcija, FKD - proizvodnaja rash. f-ii
c   IL(*) - number of points on corresponding link
	COMMON/COEF/ COE(LINK,Ngr,10),RCOE(LINK,Ngr,6),IL(LINK)
C COE(L,i) - Array for linearised equatons coefficiants A(i),B(i),...,
C i=1,IL(L)-1, L - branch number, A(1),B(1),...,G1(1);A(2),B(2),... 
C RCOE(L,i) - Array for coefficients E(i),F(i),H(i),E1(i),F1(i),H1(i)
	COMMON/SWEEP/ LIN(LINK,2),DNY(node)
C LIN(L,1) - initial NODE of branch L (point IL(L)), LIN(L,2) - end
c NODE of branch L (point 1)
c DNY(node) - array for level elevation in nodes
	Common/SCAL/phi,teta,bet,eps,delt,delx,fi,fn
        COMMON/CB/ TA(LINK,Ngr),ALF(LINK,Ngr)
	COMMON/dY/  dY(2,LINK,Ngr),dQ(2,LINK,Ngr)
	COMMON/MATR/ A(NODE,NODE),B(NODE),MM(node,2),
     *               M_IN(node,LINK),MOUT(node,LINK)
	COMMON/SOLUT/MGRUP(node),MDAM(node,2),NGRUP(node,2),NDAM
c MM(i,*), i-node number
c MM(i,1)-kol-vo vtokov
c MM(i,2)-kol-vo vytokov
c M_IN(i,*)-vtok numbers
c MOUT(i,*)-vytok numbers
c   MGRUP(*) - numbers of dams nodes and last node
        COMMON/CROSS/ WPR(LINK,Ngr,n_cross)
	COMMON /YBS/YBASE(LINK)
C dY(*,*,*)- delY, dQ(*,*,*) - delQ at every grid point on all branches
c  LINKS - real number of branches 
c  M - real number of nodes
	DIMENSION VVV(LINK,ngr),CH(28,1000),TCH(1000)
	DIMENSION QTR1(4,5000),QTR2(4,5000),QTR3(4,5000),QTR4(4,5000),
     *	QTR5(4,5000), QTR6(4,5000),QTR0(4,5000)
	DIMENSION Smean(4,100),Qmean(4,100),CBB(2500),
     %	Q_mean(100),Q_m3(100,50)

	iqtr=1
	idud=1
	jchek=0
c---------READING HYDROLOGICAL DATA--------------
 7	continue
	CALL HYD_INP(LINKS,M,KGR)  ! zadanie nach. uslov. Q,W,H,BW
        write(*,*)'QGR=',(QGR(i,1),i=1,LINKS)
        write(*,*)'IT=',(IT(i),i=1,LINKS)
	write(*,*)'links=',links
c	pause
c        write(*,*)'YGR=',(YGR(i,1),i=1,LINKS)
c        write(*,*)'ITY=',(ITY(i),i=1,LINKS)
c	write(*,*)'IT_C=',(IT_C(i),i=1,LINKS)
c	write(*,*)'C_GR=',(C_GR(2,i),i=1,31)
c	write(*,*)'TGR_C=',(TGR_C(1,i),i=1,31)
c	pause
	open(1,file='config.riv')
	read(1,*)tm
	read(1,*)delt
	read(1,*)ts
c tm - time of calculation        
c ts - interval zapisi dannyh
	close(1)
        	write(*,*)'tm=',tm,' ts=',ts
c---------READING COEFFICIENT DATA--------------
              OPEN(3,file='nusl') 
                     read(3,*) ktrnsprt
                     read(3,*)
              read(3,*) delx,phi,teta,bet,eps,fi,fn 
                     read(3,*)
              read(3,*) d_part,dmu,dnu,z_obm,ER,grav,Zobm,delvv
                    read(3,*)
	      read(3,*) flam,a12,a21,a13,a31,fKs,fKb,ros,ro,epsb 	
              read(3,*)
              read(3,*) (IL(i),i=1,LINKS)
		read(3,*)
		read(3,*)(LIN(j,1),j=1,LINKS)
		read(3,*)(LIN(j,2),j=1,LINKS)
		read(3,*)
		read(3,*)NDAM ! kol-vo grupp
		read(3,*)
		read(3,*)(MDAM(i,1),i=1,NDAM-1) ! verhn. uzly plotin
		read(3,*)
		read(3,*)(MDAM(i,2),i=1,NDAM-1)  ! nizhnie uzly plotin
		read(3,*)
		read(3,*)(NGRUP(i,1),i=1,NDAM) ! verhn. uzly grupp
		read(3,*)
		read(3,*)(NGRUP(i,2),i=1,NDAM)  ! nizhnie uzly grupp
		read(3,*)
		read(3,*)(MGRUP(i),i=1,NDAM)
              CLOSE(3) 
	vtr=24.*3600
          flam=flam/vtr !(1/sek)
           a12=a12/vtr !(1/sek)
	   a21=a21/vtr
           a13=a13/vtr !(1/sek)
           a31=a31/vtr !(1/sek)
c	write(*,*)'a12=',a12,' flam=',flam,' a13=',a13,' a31=',a31
c	pause              
	write(*,*)'MGRUP=',(MGRUP(i),i=1,NDAM+1)
c	pause
	              OPEN(3,file='node')
              read(3,*)
	write(*,*)'M=',M
	DO i=1,M
	read(3,*) k,MM(k,1),MM(k,2),(M_IN(k,j),j=1,MM(k,1)),
     +            (MOUT(k,j),j=1,MM(k,2))
	ENDDO
              CLOSE(3)
	WRITE(*,*)'MM='
	vsc=sqrt(grav*2.7)
	tsc=330000/vsc/86400
	
        		DO L=1,LINKS
        		do i=1,IL(L)
        		H(1,L,i)=Y(1,L,i)-YB(L,i)
        	enddo
        	ENDDO



c------------- CROSSECTION READING (WPR)-------------------------------
c	write(*,*)'CROSSECTION READING'
	call WPR_RE(delY,LINKS)
c	write(*,*)(YBASE(i),i=1,links)
c        PAUSE'WPR is written'
c---------------------------------------------------------------------
	J=1
c	WRITE(*,*)'j=',J

        CALL NEWCROSS(delY,LINKS,J)

!         definition of max point number on the links

	CALL INITIAL(LINKS,LMAX)

        KTS=int(tm*86400/delt)
        kint=int(ts*86400/delt)
        write(*,*)'KTS=',kts,' kint=',kint
	pause

                 DO ISTEP=1,KTS ! kol-vo shagov po vremeni
        
!!!!!!!! ZADANIE VERHNIH GRANICHNYH USLOVIJ Q !!!!!!!!!!!!!!!!!!!!       
	call BOUND_Q(delt,LINKS,istep,M)
c	if(istep.ge.23) then
c	write(*,*)'TGR(!,1-3)=',TGR(1,1),TGR(1,2),TGR(1,3)
c	write(*,*)'QGR(!,1-3)=',QGR(1,1),QGR(1,2),QGR(1,3)
c	write(*,*)'istep=',istep
c	write(*,*)'QGRAN(1)=',QGRAN(1)
c	pause
c	endif
	call bound_dam(LINKS,istep,M)
	Q_balans_gr=Q_balans_gr+(Q(2,3,IL(3))+QGRAN(4))*delt
!!!!!!! ZADANIE  GRANICHNYH USLOVIJ FOR CONCENTRATION(C)       
        DO L=1,LINKS
	if(IT_C(L).ne.0) then ! 1
         tim=istep*delt/86400.
       if(tim.le.TGR_C(L,1)) then ! 2
        C_GRAN(L)=C_GR(L,1)
        goto 310
        else
        do ik=1,IT_C(L)-1
        if(tim.gt.TGR_C(L,ik).and.tim.le.TGR_C(L,ik+1)) then ! 3
        C_GRAN(L)=C_GR(L,ik)
        goto 310
        else
          C_GRAN(L)=C_GR(L,IT_C(L))
        end if ! 3
        enddo
        end if ! 2
        endif ! 1
 310    continue        
        ENDDO
!!!!!!! ZADANIE  GRANICHNYH USLOVIJ FOR SEDIMENT(S)       
        DO L=1,LINKS
	if(IT_S(L).ne.0) then ! 1
         tim=istep*delt/86400.
       if(tim.le.TGR_S(L,1)) then ! 2
        S_GRAN(L)=S_GR(L,1)
        goto 320
        else
        do ik=1,IT_S(L)-1
        if(tim.gt.TGR_S(L,ik).and.tim.le.TGR_S(L,ik+1)) then ! 3
        S_GRAN(L)=S_GR(L,ik)
        goto 320
        else
          S_GRAN(L)=S_GR(L,IT_S(L))
        end if ! 3
        enddo
        end if ! 2
        endif ! 1
 320    continue        
        ENDDO
!!!!!!! ZADANIE  GRANICHNYH USLOVIJ FOR CONCENR. IN SEDIMENT(CS)       
        DO L=1,LINKS
	if(IT_CS(L).ne.0) then ! 1
         tim=istep*delt/86400.
       if(tim.le.TGR_CS(L,1)) then ! 2
        CS_GRAN(L)=CS_GR(L,1)
        goto 330
        else
        do ik=1,IT_CS(L)-1
        if(tim.gt.TGR_CS(L,ik).and.tim.le.TGR_CS(L,ik+1)) then ! 3
        CS_GRAN(L)=CS_GR(L,ik)
        goto 330
        else
          CS_GRAN(L)=CS_GR(L,IT_CS(L))
        end if ! 3
        enddo
        end if ! 2
        endif ! 1
 330    continue        
        ENDDO
	call bound_y(delt,LINKS,istep,M)
c	write(*,*)(YGRNO(i),i=1,M)
c	pause

                    k=1 ! 1 iteration
	CALL ITER_1 (LINKS)

 1            CONTINUE    ! new iteration
c	write(*,*)'k=',k
c	pause
             CALL COEFN(LINKS,LMAX,istep)
             CALL MATRIX(LINKS,M,istep) 
             CALL SOLSYST(M)  ! solution Ax=B
c		write(*,*)'LEVEL ELEVATION AT NODES'
c                 write(*,*)'dNY=',(dNY(i),i=1,M)
              CALL BACKSW(LINKS) ! obratnaja progonka
		eps1=eps
               CALL SNORM(LINKS,eps,eps1,JJ)  ! sravnenie 
               if(JJ.EQ.1) then
                      go to 10
               else
                      go to 20
               end if

 10               k=k+1
	CALL ITER_2 (LINKS)
           write(*,*)'k=',k
c    Vychisl. novyh W,BW,FK,FKD
                   J=2
        CALL NEWCROSS(delY,LINKS,J)
                   GO TO 1       ! new iteration !

C                           NEXT TIME STEP
 20     CONTINUE
	CALL STEP_DEL (LINKS)
        CALL NEWCROSS(delY,LINKS,2)
	CALL NEW_STEP (LINKS)
c 	if(ktrnsprt.ne.0) then
c	CALL TRANSP(delx,delt,LINKS,istep,M) !for transp.
c	endif
	itv=int(istep*delt/(ts*86400)) ! output intermediate data
	itvp= int((istep-1)*delt/(ts*86400))

	        if(itv.GT.itvp) then
 		th=delt*istep/86400
        	to=1
c	QTR2(1,iqtr)=Q(2,2,5) 	! p.Kiev base
c	QTR2(2,iqtr)=Y(2,2,5)   ! 20 km
!	QTR2(3,iqtr)=Y(2,2,4)   ! 18 km ot Kiev dam

	QTR2(1,iqtr)=Q(2,2,27) 	! p.Kiev delx=500m
!	QTR2(2,iqtr)=Y(2,2,27)   ! 19.5 km ot Kiev dam
	QTR2(2,iqtr)=Y(2,2,28)   ! 20 km ot Kiev dam

	QTR3(1,iqtr)=Q(2,2,114) 	! Kijliv, 63 km
	QTR3(2,iqtr)=Y(2,2,114)   

	iqtr=iqtr+1
	endif        

	write(*,*)'new time step',' ISTEP=',istep+1

               END DO ! (ISTEP, po vremeni) 
 		th=delt*KTS/86400
 50       continue        
        to=10E+10

	open(1,file='post Kiev')
		do i=1,iqtr-1
	write(1,*)i,QTR2(1,i),QTR2(2,i) 
		enddo
	close(1)
	open(1,file='Kijliv')
		do i=1,iqtr-1
	write(1,*)i,QTR3(1,i),QTR3(2,i)
		enddo
	close(1)

 800	CALL HYD_OUT(LINKS,th,M,to,istep)  ! M - node number
 850	format((f6.1),x,5(f12.4))
 400     format((i5),6(F12.3))        
c 400     format((i5),4(E12.3))        

         STOP
         END
