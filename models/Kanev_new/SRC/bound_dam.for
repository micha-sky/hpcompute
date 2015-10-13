             SUBROUTINE BOUND_dam(LINKS,istep,M)

	INCLUDE 'parameters' 
                    	
	COMMON/BASE/ Q(3,LINK,Ngr),Y(3,LINK,Ngr),W(3,LINK,Ngr),
     *	Z(3,LINK,Ngr),
     *  BW(3,LINK,Ngr),H(3,LINK,Ngr),YB(LINK,Ngr),FK(2,Ngr),FKD(2,Ngr)
	COMMON/SWEEP/ LIN(LINK,2),DNY(node)

       Common/K/ C(2,LINK,Ngr),CS(2,LINK,Ngr),CB(2,LINK,Ngr),
     *            S(2,LINK,Ngr),Seqv(LINK,Ngr)
	COMMON/COEF/ COE(LINK,Ngr,10),RCOE(LINK,Ngr,6),IL(LINK)
        COMMON/GRN/ QGR(LINK,500),TGR(LINK,500),IT(LINK)
	COMMON/GRAN/QGRAN(LINK),QGRNO(node)
	COMMON/GRANC/C_GRAN(LINK),CS_GRAN(LINK),S_GRAN(LINK)
	COMMON/MATR/ A(NODE,NODE),B(NODE),MM(node,2),
     *               M_IN(node,LINK),MOUT(node,LINK)
	COMMON/SOLUT/MGRUP(node),MDAM(node,2),NGRUP(node,2),NDAM

	if(NDAM.NE.1) then
		do i=1,NDAM-1
	jup=MDAM(i,1)    ! nomer verhnego uzla plotiny
	ivt=M_IN(jup,1)   ! nomer br., vtek. v uzel jup
	jlow=MDAM(i,2)   ! nomer nizhnego uzla plotiny
c (predpolagaetsa kol-vo vtek. bran. = 1 )
	IF(Q(2,ivt,IL(ivt)).LE.0) THEN
	QGRNO(jlow)=0
	C_GRAN(jlow)=0
	CS_GRAN(jlow)=0
	S_GRAN(jlow)=0
	ELSE
	QGRNO(jlow)=Q(2,ivt,IL(ivt))
	C_GRAN(jlow)=C(2,ivt,IL(ivt))
	CS_GRAN(jlow)=CS(2,ivt,IL(ivt))
	S_GRAN(jlow)=S(2,ivt,IL(ivt))
	ENDIF
		enddo ! (po i, i=1,NDAM)	
	endif
        RETURN
        END
