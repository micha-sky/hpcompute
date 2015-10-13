	SUBROUTINE MATRIX(LINKS,M,istep)
c M - number of nodes
	INCLUDE 'parameters' 
	COMMON/BASE/ Q(3,LINK,Ngr),Y(3,LINK,Ngr),W(3,LINK,Ngr),
     *	Z(3,LINK,Ngr),
     *  BW(3,LINK,Ngr),H(3,LINK,Ngr),YB(LINK,Ngr),FK(2,Ngr),FKD(2,Ngr)
	COMMON/COEF/ COE(LINK,Ngr,10),RCOE(LINK,Ngr,6),IL(LINK)
	COMMON/SWEEP/ LIN(LINK,2),DNY(node)
	COMMON/GRAN/QGRAN(LINK),QGRNO(node)
	COMMON/GRANY/YGRAN(LINK),YGRNO(node)
	COMMON/MATR/ A(NODE,NODE),B(NODE),MM(node,2),
     *               M_IN(node,LINK),MOUT(node,LINK)
c		write(*,*)'matrix works'
                     DO i=1,M
                     DO j=1,M
               A(i,j)=0.
                      enddo
                       enddo
		DO I=1,M
	IF(MM(i,2).NE.0) THEN  ! kol-vo vytokov ne 0

	if(MM(i,1).ne.0) then  ! kol-vo vtokov ne 0
	DO iv=1,MM(i,1)
	j=M_IN(i,iv)  ! j - nomer vtoka
	jin=LIN(j,1)  ! nomer yzla, iz kotorogo vytekaet j 
	A(i,jin)=RCOE(j,1,4)
	ENDDO		
	endif		! Opredelenie koef. Aij dla vtokov
c---------
	DO ivyt=1,MM(i,2)
	j=MOUT(i,ivyt)	 ! j - nomer vytoka
	jend=LIN(j,2) 	 ! nomer yzla, v kotoryj vpadaet j 
	A(i,jend)=-RCOE(j,1,3)
	ENDDO		! Opredelenie koef. Aij dla vytokov
c----------
	sum_vt=0
		DO iv=1,MM(i,1)
		j=M_IN(i,iv)  ! j - nomer vtoka
		sum_vt=sum_vt+RCOE(j,1,6)
		ENDDO	
	sum_vyt=0
	DO ivyt=1,MM(i,2)
	j=MOUT(i,ivyt)	 ! j - nomer vytoka
	sum_vyt=sum_vyt+RCOE(j,1,1)
	ENDDO
	A(i,i)=sum_vt-sum_vyt
c------------ opredelenie B(i)
		sumbvt=0
		DO iv=1,MM(i,1)
		j=M_IN(i,iv)  ! j - nomer vtoka
		sumbvt=sumbvt+Q(2,j,IL(j))
		ENDDO
		sumbvyt=0
	DO ivyt=1,MM(i,2)
	j=MOUT(i,ivyt)	 ! j - nomer vytoka
	sumbvyt=sumbvyt+Q(2,j,1)
	ENDDO
		sumfvt=0
		DO iv=1,MM(i,1)
		j=M_IN(i,iv)  ! j - nomer vtoka
		sumfvt=sumfvt+RCOE(j,1,5)
		ENDDO
	sumfvyt=0
	DO ivyt=1,MM(i,2)
	j=MOUT(i,ivyt)	 ! j - nomer vytoka
	sumfvyt=sumfvyt+RCOE(j,1,2)
	ENDDO
	B(i)=-QGRNO(i)-sumbvt+sumbvyt-sumfvt+sumfvyt
	if(istep.eq.1.and.i.eq.5.or.istep.eq.2399.and.i.eq.5) then
c	write(*,*)'istep=',istep,' i=',i
c	write(*,*)'sumbvt=',sumbvt,' sumbvyt=',sumbvyt,
c     *            ' sumfvt=',sumfvt,' sumfvyt=',sumfvyt
c	pause
	endif
		ELSE
	A(i,i)=1
	j=M_IN(i,1) 	! nomer odnogo iz branch.,vtek.v yzel i
	B(i)=YGRNO(i)-Y(2,j,IL(j))
	ENDIF
		ENDDO ! po I
	if(istep.eq.1.or.istep.eq.2399) then
c	write(*,*)'istep=',istep
c	write(*,*)'A(i,j)'
c	do i=1,M
c	write(*,*)(A(i,j),j=1,M)
c	enddo
c	write(*,*)'QGRNO(5)=',QGRNO(5)
c	write(*,*)'B(i)'
c	write(*,*)(B(i),i=1,M)
c	pause
	endif
         RETURN
         END