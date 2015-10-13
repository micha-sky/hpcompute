             SUBROUTINE iter_1(LINKS)

	INCLUDE 'parameters' 
                    	
	COMMON/BASE/ Q(3,LINK,Ngr),Y(3,LINK,Ngr),W(3,LINK,Ngr),
     *	Z(3,LINK,Ngr),
     *  BW(3,LINK,Ngr),H(3,LINK,Ngr),YB(LINK,Ngr),FK(2,Ngr),FKD(2,Ngr)
	COMMON/SWEEP/ LIN(LINK,2),DNY(node)

	COMMON/dY/  dY(2,LINK,Ngr),dQ(2,LINK,Ngr)
	COMMON/COEF/ COE(LINK,Ngr,10),RCOE(LINK,Ngr,6),IL(LINK)

                       DO L=1,LINKS  ! zadanie dY,dQ na 1-oj iteracii
                       DO i=1,IL(L)
                         do j=1,2
                                dY(j,L,i)=0.
                                dQ(j,L,i)=0.
			end do
                       END DO
                       END DO

	RETURN
	END

