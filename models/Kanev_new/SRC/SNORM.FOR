	SUBROUTINE SNORM(LINKS,eps,eps1,JJ)
	INCLUDE 'parameters' 
	COMMON/dY/  dY(2,LINK,Ngr),dQ(2,LINK,Ngr)
	COMMON/COEF/ COE(LINK,Ngr,10),RCOE(LINK,Ngr,6),IL(LINK)
c  sravnenie norm vektorov dY and dQ
c  JJ - priznak povtornoj iteracii
c  JJ=1 trebuetca povtornaja iteracia, JJ=0 - ne trebuetca
c		write(*,*)'subr. snorm'
c		do L=1,links
c		write(*,*)'L=',l
c                write(*,*)'dY1=',(dY(1,L,i),i=1,IL(L))
c                write(*,*)'dQ1=',(dQ(1,L,i),i=1,IL(L))
c                write(*,*)'dY2=',(dY(2,L,i),i=1,IL(L))
c                write(*,*)'dQ2=',(dQ(2,L,i),i=1,IL(L))
c		end do
                  fNY1=abs(dY(1,1,1)) ! Norma dY na k+1 iteracii
                  fNQ1=abs(dQ(1,1,1)) ! Norma dQ na k+1 iteracii
                  fNY2=abs(dY(2,1,1)) ! Norma dY na k iteracii
                  fNQ2=abs(dQ(2,1,1)) ! Norma dQ na k iteracii
                   do L=1,LINKS
                   DO i=1,IL(L)
                     if(abs(dY(1,L,i)).GT.fNY1) then
                     fNY1=abs(dY(1,L,i))
	iiy=i
                     end if
                     if(abs(dQ(1,L,i)).GT.fNQ1) then
                     fNQ1=abs(dQ(1,L,i))
        	iiq=i
                     end if
                     if(abs(dY(2,L,i)).GT.fNY2) then
                     fNY2=abs(dY(2,L,i))
                     end if
                     if(abs(dQ(2,L,i)).GT.fNQ2) then
        	ii2=i
                     fNQ2=abs(dQ(2,L,i))
                     end if
                 END DO
                 END DO
c		write(*,*)'fNY1=',fNY1,' iiy=',iiy
c		write(*,*)'fNQ1=',fNQ1,' iiq=',iiq
c		write(*,*)'fNY2=',fNY2
c		write(*,*)'fNQ2=',fNQ2,' ii2=',ii2
c		write(*,*)'eps=',eps,' eps1=',eps1
                   PY=eps1*fNY2+eps
                   PQ=eps1*fNQ2+eps
c		write(*,*)'PY=',PY
c		write(*,*)'PQ=',PQ
c		pause
        if(fNY1.gt.PY.OR.fNQ1.gt.PQ) then
                    JJ=1
        else
                    JJ=0
        end if
                    RETURN
                    END