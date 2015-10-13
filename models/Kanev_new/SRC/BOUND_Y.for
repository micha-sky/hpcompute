             SUBROUTINE BOUND_Y(delt,LINKS,istep,M)

	INCLUDE 'parameters' 
                    	
	COMMON/BASE/ Q(3,LINK,Ngr),Y(3,LINK,Ngr),W(3,LINK,Ngr),
     *	Z(3,LINK,Ngr),
     *  BW(3,LINK,Ngr),H(3,LINK,Ngr),YB(LINK,Ngr),FK(2,Ngr),FKD(2,Ngr)
	COMMON/SWEEP/ LIN(LINK,2),DNY(node)

	COMMON/GRNY/  YGR(LINK,500),TGRY(LINK,500),ITY(LINK)
	COMMON/GRANY/YGRAN(LINK),YGRNO(node)
	COMMON/MATR/ A(NODE,NODE),B(NODE),MM(node,2),
     *               M_IN(node,LINK),MOUT(node,LINK)

	DO L=1,LINKS
	if(ITY(L).ne.0) then  !  1
	tim=istep*delt/86400.
	if(tim.le.TGRY(L,1)) then ! 2
	YGRAN(L)=YGR(L,1)
	goto 300
	else
	do ik=1,ITY(L)-1
	if(tim.gt.TGRY(L,ik).and.tim.le.TGRY(L,ik+1)) then ! 3
	YGRAN(L)=YGR(L,ik)
	goto 300
	else
	YGRAN(L)=YGR(L,ITY(L))
	endif ! 3
	enddo
	endif  ! 2
	endif  ! 1
 300	continue
	ENDDO
		DO i=1,M
	IF(MM(i,2).EQ.0) THEN  ! kol-vo vytokov = 0
	YGRNO(i)=YGRAN(M_IN(i,1))
	endif
		ENDDO	
        RETURN
        END
