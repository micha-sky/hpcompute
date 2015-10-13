             SUBROUTINE BOUND_Q(delt,LINKS,istep,M)

	INCLUDE 'parameters' 
                    	
	COMMON/BASE/ Q(3,LINK,Ngr),Y(3,LINK,Ngr),W(3,LINK,Ngr),
     *	Z(3,LINK,Ngr),
     *  BW(3,LINK,Ngr),H(3,LINK,Ngr),YB(LINK,Ngr),FK(2,Ngr),FKD(2,Ngr)
	COMMON/SWEEP/ LIN(LINK,2),DNY(node)

	COMMON/GRN/   QGR(LINK,500),TGR(LINK,500),IT(LINK)
	COMMON/GRAN/QGRAN(LINK),QGRNO(node)
	COMMON/MATR/ A(NODE,NODE),B(NODE),MM(node,2),
     *               M_IN(node,LINK),MOUT(node,LINK)

        tim=istep*delt/86400.

	DO L=1,LINKS
c	 write(*,*)'IT=',IT(L)
	 if(IT(L).ne.0) then ! 1
        if(tim.le.TGR(L,1)) then ! 2
         QGRAN(L)=QGR(L,1)
c         write(*,*)'tim<TGR(L,2):L,ik,tim,TGR2',L,ik,tim,TGR(L,2)
         goto 200
        else
         do ik=1,IT(L)-1

          if(tim.gt.TGR(L,ik).and.tim.lt.TGR(L,ik+1)) then ! 3

		  del_TT=(TGR(L,ik+1)-TGR(L,ik))*0.6*86400 
c           write(*,*)'del_TT=',del_TT                     
	      nd=int(del_TT/delt)                          
c           write(*,*)'nd=',nd                             
		  if(nd.ne.0) then                         
	       delQ=(QGR(L,ik+1)-QGR(L,ik))/nd             
		  else                                     
		   delQ=0                                  
		  endif                                    
     	if((tim+delt/86400).le.(TGR(L,ik)+del_TT/86400)) then          
            QGRAN(L)=QGRAN(L)+delQ                         
	endif
	if((tim.le.(TGR(L,ik)+del_TT/86400)).and.
     %	((tim+delt/86400).gt.(TGR(L,ik)+del_TT/86400))) then          
            QGRAN(L)=QGR(L,ik+1)
	endif                                              
            goto 200
c		 endif                                   
          endif ! 3
c          if(ik.eq.(IT(L)-1)) then
          if(tim.ge.TGR(L,IT(L))) then
           QGRAN(L)=QGR(L,IT(L))
           goto 200
          endif
         enddo
        end if ! 2
       endif ! 1
 200   continue        
      ENDDO		! po L
		DO i=1,M
	IF(MM(i,2).NE.0) THEN  ! kol-vo vytokov ne 0
	QGRNO(i)=0
	do j=1,MM(i,2)
	QGRNO(i)=QGRNO(i)+QGRAN(MOUT(i,j))
	enddo
	endif
		ENDDO ! (po i)
c        write(*,*)'QGRNO=',(QGRNO(i),i=1,M)
c        pause
c	write(*,*)'bound finished'	
        RETURN
        END
