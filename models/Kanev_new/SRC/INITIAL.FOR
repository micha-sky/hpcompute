	SUBROUTINE INITIAL(LINKS,LMAX)

	INCLUDE 'parameters' 
	
	COMMON/BASE/ Q(3,LINK,Ngr),Y(3,LINK,Ngr),W(3,LINK,Ngr),
     *	Z(3,LINK,Ngr),
     *  BW(3,LINK,Ngr),H(3,LINK,Ngr),YB(LINK,Ngr),FK(2,Ngr),FKD(2,Ngr)
	COMMON/COEF/ COE(LINK,Ngr,10),RCOE(LINK,Ngr,6),IL(LINK)

!         definition of max point number on the links
                   LMAX=IL(1)  
                   DO i=2,LINKS 
                      if(IL(i).GT.LMAX) then
                      LMAX=IL(i)
                      end if
                   ENDDO
                      DO L=1,LINKS                     
                            DO i=1,IL(L)
                              do j=2,3
                   Y(J,L,I)=Y(1,L,i)
                   Q(J,L,i)=Q(1,L,i)
                   W(J,L,i)=W(1,L,i)
                   BW(J,L,i)=BW(1,L,i)
                   H(J,L,i)=H(1,L,i)
                                 end do ! po j
                              END DO ! po i
                            END DO ! po L
		RETURN
		END
