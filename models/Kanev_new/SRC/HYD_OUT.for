	SUBROUTINE HYD_OUT(LINKS,th,M,to,istep)
	INCLUDE 'parameters' 
	COMMON/BASE/ Q(3,LINK,Ngr),Y(3,LINK,Ngr),W(3,LINK,Ngr),
     *	Z(3,LINK,Ngr),
     *  BW(3,LINK,Ngr),H(3,LINK,Ngr),YB(LINK,Ngr),FK(2,Ngr),FKD(2,Ngr)
     
       Common/K/ C(2,LINK,Ngr),CS(2,LINK,Ngr),CB(2,LINK,Ngr),
     *            S(2,LINK,Ngr),Seqv(LINK,Ngr)

	COMMON/COEF/ COE(LINK,Ngr,10),RCOE(LINK,Ngr,6),IL(LINK)
          WRITE(*,*)' WRITING HYD FILE. '
	lu=1
        OPEN(UNIT=lu,STATUS='UNKNOWN',FILE='hyd1',err=13)
c-----------------
	lu1=2
        OPEN(UNIT=lu1,ACCESS='APPEND',STATUS='UNKNOWN',
     *   FILE='1.hyp',err=13)
c         write(*,'(A)') STR
        j=1
	Np=0
	do L=1,LINKS
	Np=Np+IL(L)
	enddo
	 WRITE (lu,*) 'To= ',th
	WRITE(lu,*) LINKS,M
         WRITE (lu,*) Np
         WRITE (lu,*) 'Branch Grid Discharge Level Bottom level'
c
c	 WRITE (lu1,*) 'To= ',th
	 WRITE (lu1,*) th
c	WRITE(lu1,*) LINKS,M
c         WRITE (lu1,*) Np
         WRITE (lu1,*) 'Branch Grid C	Cvol  CS  CB   S'

         DO 11 L=1,LINKS
           DO 11 K=1,IL(L)
	vv=Q(2,L,k)/W(2,L,k)
           WRITE (lu,*)
     +    L,K, Q(2,L,k),Y(2,L,k),YB(L,k)
 11	continue
		do L=1,LINKS
		do k=1,IL(L)
c	if(istep.le.720) then
c	if(L.le.4) then
c            WRITE (lu1,*)
c    +    L,K, Q(2,L,k),W(2,L,k),vv
c     +    L,K,C(2,L,k)*27.02/1000,CB(2,L,k)
c	endif
c	else
		vv=Q(2,L,k)/W(2,L,k)
	Cvol=(C(2,L,k)+CS(2,L,k)*S(2,L,k))/1000.
		WRITE(lu1,300)
c     +   L,k,Q(2,L,k),C(2,L,k)/1000.,CS(2,L,k),CB(2,L,k),S(2,L,k)
     +   L,k,C(2,L,k)/1000.,Cvol,CS(2,L,k),CB(2,L,k),S(2,L,k)
! конц. в pCi/l
c     +   L,k,Q(2,L,k),vv,Seqv(L,k),S(2,L,k)
c     +   L,k,Q(2,L,k),Y(2,L,k),vv
c	endif
		enddo
		enddo
c---------  end loop  ------------    
	ii=0
        write(lu,*)'Boundary Condition:',ii  !!!
        if(to.eq.10E+10) then
	 WRITE (lu1,*) 'To= ',to
        endif
	GOTO 131      
13      write(*,*) ' Can not open file'
        goto 14 
131        write(*,*)' NORMAL END HYD FILE    '
 	close(lu)
 	close(lu1) 
c 		open(3,file='cntr1')
c	write(3,*)'conc. solution        conc. vzvesi      conc.bottom'
c	write(3,*)'       Branch       Grid         Conc. solution
c     *     Conc. susp.	  Conc. bottom'
c           DO  L=1,LINKS
c           DO  K=1,IL(L)
c	           WRITE (3,*)
c     +    L,K, C(2,L,k),CS(2,L,k),CB(2,L,k)
c     +    L,K, C(2,L,k)*27.02/1000,CB(2,L,k)
c		ENDDO
c		ENDDO
c	close(3)
c 		open(3,file='bot_layer')
c	write(3,*)'       Branch       Grid         Bottom layer'
c           DO  L=1,LINKS
c           DO  K=1,IL(L)
c	           WRITE (3,*) L,K, Z(2,L,k)
c		ENDDO
c		ENDDO
c	close(3)
c 		open(3,file='sedim1')
c	write(3,*)'       Branch       Grid     Concentration 
c     *  of sediment'
c           DO  L=1,LINKS
c           DO  K=1,IL(L)
c	           WRITE (3,*) L,K, S(2,L,k)
c		ENDDO
c		ENDDO
c	close(3)
 200     format(21a,i5)        
 300     format(2(i5),5(E12.3))        
14        return
	 end
