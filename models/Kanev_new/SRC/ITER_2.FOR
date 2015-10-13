             SUBROUTINE iter_2(LINKS)

	INCLUDE 'parameters' 
                    	
	COMMON/BASE/ Q(3,LINK,Ngr),Y(3,LINK,Ngr),W(3,LINK,Ngr),
     *	Z(3,LINK,Ngr),
     *  BW(3,LINK,Ngr),H(3,LINK,Ngr),YB(LINK,Ngr),FK(2,Ngr),FKD(2,Ngr)
	COMMON/SWEEP/ LIN(LINK,2),DNY(node)

	COMMON/dY/  dY(2,LINK,Ngr),dQ(2,LINK,Ngr)
	COMMON/COEF/ COE(LINK,Ngr,10),RCOE(LINK,Ngr,6),IL(LINK)

             DO L=1,LINKS  ! pereprisvoenie dY,dQ 
             DO i=1,IL(L)  ! i compute Y(2,*,*),Q(2,*,*)
                        dY(2,L,i)=dY(1,L,i)
                        dQ(2,L,i)=dQ(1,L,i)
                        Y(2,L,i)=Y(2,L,i)+dY(2,L,i) 
                        Q(2,L,i)=Q(2,L,i)+dQ(2,L,i) 
             END DO
             END DO

	RETURN
	END

