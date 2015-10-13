	SUBROUTINE BACKSW(LINKS)
C                BACKWARD SWEEP
C  compute dY(1,L,i),dQ(1,L,i) via DNY(node)
	INCLUDE 'parameters' 
	COMMON/COEF/ COE(LINK,Ngr,10),RCOE(LINK,Ngr,6),IL(LINK)
	COMMON/SWEEP/ LIN(LINK,2),DNY(node)
	COMMON/dY/  dY(2,LINK,Ngr),dQ(2,LINK,Ngr)
c		write(*,*)'backsw works'
             DO L=1,LINKS
c	write(*,*)'L=',l
c	write(*,*)'LIN(L,1)=',LIN(L,1),'LIN(L,2)',LIN(L,2)
c	pause
               MINT=LIN(L,1)  ! node for point 1 on L (nachalo)
               MEND=LIN(L,2)   ! node for point IL(L) on L (konec)
               dY(1,L,1)=DNY(MINT)  ! dely(1) for branch L (nachalo)
               dY(1,L,IL(L))=DNY(MEND) ! dely(il) (konec)
c	write(*,*)'dY(1,L,1)',dY(1,L,1),' dY(1,L,IL(L))=',dY(1,L,IL(L))
c	pause
                dQ(1,L,1)=RCOE(L,1,1)*dY(1,L,1)+RCOE(L,1,2)+
     *          RCOE(L,1,3)*dY(1,L,IL(L)) 
c	write(*,*)'E=',RCOE(L,1,1),' F=',RCOE(L,1,2),' H=',RCOE(L,1,3)
c	write(*,*)'E1=',RCOE(L,1,4),' F1=',RCOE(L,1,5),
c     *      ' H1=',RCOE(L,1,6)
                dQ(1,L,IL(L))=RCOE(L,1,4)*dY(1,L,1)+RCOE(L,1,5)+
     *          RCOE(L,1,6)*dY(1,L,IL(L)) 
c      write(*,*)'dQ(1,L,1)=',dQ(1,L,1),' dQ(1,L,IL(L))=',dQ(1,L,IL(L))
c	pause
                    DO i=1,IL(L)-2
                 A=COE(L,i,1)
                 B=COE(L,i,2)
                 C=COE(L,i,3)
                 D=COE(L,i,4)
                 G=COE(L,i,5)
                 A1=COE(L,i,6)
                 B1=COE(L,i,7)
                 C1=COE(L,i,8)
                 D1=COE(L,i,9)
                 G1=COE(L,i,10)
            DENOM=A*B1-A1*B
            RL=(C*B1-C1*B)/DENOM
            RMM=(D*B1-D1*B)/DENOM
            RN=(G*B1-G1*B)/DENOM
          dY(1,L,i+1)=dY(1,L,i)*RL+dQ(1,L,i)*RMM+RN
          dQ(1,L,i+1)=dY(1,L,i+1)*RCOE(L,i+1,1)+RCOE(L,i+1,2)+
     *       dY(1,L,IL(L))*RCOE(L,i+1,3)
                    END DO ! po i
              END DO ! po L
c                    DO ik=1,IL(1)-1 !!!experiment
c	jk=IL(1)-ik
c                 A=COE(1,jk,1)
c                 B=COE(1,jk,2)
c                 C=COE(1,jk,3)
c                 D=COE(1,jk,4)
c                 G=COE(1,jk,5)
c                 A1=COE(1,jk,6)
c                 B1=COE(1,jk,7)
c                 C1=COE(1,jk,8)
c                 D1=COE(1,jk,9)
c                 G1=COE(1,jk,10)
c            DENOM=C*D1-C1*D
c            RL1=(A*D1-A1*D)/DENOM
c            RMM1=(B*D1-B1*D)/DENOM
c            RN1=(G1*D-G*D1)/DENOM
c          dY(2,1,jk)=dY(1,1,jk+1)*RL1+dQ(1,1,jk+1)*RMM1+RN1
c          dQ(2,1,jk)=dY(1,1,jk)*RCOE(1,jk,1)+RCOE(1,jk,2)+
c     *       dY(1,1,IL(1))*RCOE(1,jk,3)
c                    END DO ! po ik
c	dY(2,1,IL(1))=dY(1,1,IL(1))
c	dQ(2,1,IL(1))=dQ(1,1,IL(1))
c	goto 10
        L=1
        goto 10
               write(*,*)'dY1=',(dY(1,L,i),i=181,270)
c        pause
               write(*,*)'dY1=',(dY(1,L,i),i=271,360)
c        pause
               write(*,*)'dY1=',(dY(1,L,i),i=361,450)
c        pause
               write(*,*)'dY1=',(dY(1,L,i),i=451,540)
c        pause
               write(*,*)'dY1=',(dY(1,L,i),i=541,650)
c        pause
               write(*,*)'dY1=',(dY(1,L,i),i=651,740)
c        pause
               write(*,*)'dY1=',(dY(1,L,i),i=741,830)
        pause
               write(*,*)'dY1=',(dY(1,L,i),i=831,930)
        pause
               write(*,*)'dY1=',(dY(1,L,i),i=931,956)
        pause
 10	continue
        goto 20
        	write(*,*)'i=',181,270
	pause 
                write(*,*)'dQ1=',(dQ(1,L,i),i=181,270)
        pause
c	goto 20
        write(*,*)'i=',271,300
                write(*,*)'dQ1=',(dQ(1,L,i),i=271,300)
        pause
        write(*,*)'i=',301,350
                write(*,*)'dQ1=',(dQ(1,L,i),i=301,350)
        pause
        write(*,*)'i=',350,400
                write(*,*)'dQ1=',(dQ(1,L,i),i=350,400)
        pause
        write(*,*)'i=',401,450
                write(*,*)'dQ1=',(dQ(1,L,i),i=401,450)
        pause
        goto 20
                write(*,*)'dQ1=',(dQ(1,L,i),i=451,640)
        pause
                write(*,*)'dQ1=',(dQ(1,L,i),i=641,730)
        pause
                write(*,*)'dQ1=',(dQ(1,L,i),i=731,830)
        pause
                write(*,*)'dQ1=',(dQ(1,L,i),i=831,920)
        pause
                write(*,*)'dQ1=',(dQ(1,L,i),i=931,956)
        pause
c                write(*,*)'dY2=',(dY(2,L,i),i=741,830)
                write(*,*)'dQ2=',(dQ(2,L,i),i=831,956)
c        pause
c		write(*,*)'backsw finished'
c	pause
 20	continue
         RETURN
         END