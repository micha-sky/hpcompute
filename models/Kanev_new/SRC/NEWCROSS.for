	SUBROUTINE NEWCROSS(delY,LINKS,J)
	INCLUDE 'parameters' 

	COMMON/BASE/ Q(3,LINK,Ngr),Y(3,LINK,Ngr),W(3,LINK,Ngr),
     *	Z(3,LINK,Ngr),
     *  BW(3,LINK,Ngr),H(3,LINK,Ngr),YB(LINK,Ngr),FK(2,Ngr),FKD(2,Ngr)
	COMMON/COEF/ COE(LINK,Ngr,10),RCOE(LINK,Ngr,6),IL(LINK)
        COMMON/CROSS/ WPR(LINK,Ngr,n_cross)
	COMMON /dY/  dY(2,LINK,Ngr),dQ(2,LINK,Ngr)
	COMMON /YBS/ YBASE(LINK)

C SUBROUTINE CALCULATES CROSSECTIONS AND WIDTH                
c WPR - crossection table 
c W - crosssection, BW - width, H - depth, DH - depth step,
c IL - number of points on the branches, LINKS -  number of branches
c J -  time step number
		DO L=1,LINKS
		DO i=1,IL(L)
           H(J,L,i)=Y(J,L,i)-YB(L,i)
	DEL= Y(J,L,i)-YBASE(L)
	if(del.lt.0) then
	write(*,*)'from newcross'
	write(*,*)'L=',L,' i=',i, ' DEL=',DEL,' Y=',Y(J,L,i),
     *             ' YBASE=',YBASE(L),' Q=',Q(j,L,i)
	pause
	endif
	if(L.eq.3) then
c	write(*,*)'L=',L,' DEL=',DEL,' Y=',Y(J,L,i),' YBASE=',YBASE(L)
c		pause
	endif
			ih=int(DEL/delY)+1	
			fint=(DEL-(ih-1)*delY)/delY
		W(J,L,i)=WPR(L,i,ih)+fint*(WPR(L,i,ih+1)-WPR(L,i,ih))
		BW(J,L,i)=(WPR(L,i,ih+1)-WPR(L,i,ih))/delY  !!!
!		BW(J,L,i)=(W(J,L,i)-WPR(L,i,ih))/(DEL-(ih-1)*delY)  
c		BW(J,L,i)=W(J,L,i)/H(J,L,i)
        if  (BW(J,L,i).le.0) then
        write(*,*)'from NEWCROSS '
        write(*,*)'ERROR negative wids'
        write(*,*)'L=',l,' i=',i,' ih=',ih
	write(*,*)'Y=',Y(j,L,i),' YBASE=',YBASE(L),' Q=',Q(j,L,i)
	write(*,*)'WPR(ih)=',WPR(L,i,ih),' WPR(ih+1)=',WPR(L,i,ih+1)
        write(*,*)'L-1=',l-1,' i=',IL(L-1)
	write(*,*)'Y=',Y(j,L-1,IL(L-1)),' Q=',Q(j,L-1,IL(L-1))
	pause
       write(*,*)'WPR(ih)=',WPR(L,i,ih),' WPR(ih-1)=',
     +        WPR(L,i,ih-1),' WPR(ih+1)=',wpr(L,i,ih+1)
	pause
         endif
		END DO
		END DO

         RETURN
         END
