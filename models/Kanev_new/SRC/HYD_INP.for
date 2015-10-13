	SUBROUTINE HYD_INP(LINKS,M,kgr)
	INCLUDE 'parameters' 
	COMMON/BASE/ Q(3,LINK,Ngr),Y(3,LINK,Ngr),W(3,LINK,Ngr),
     *	Z(3,LINK,Ngr),
     *  BW(3,LINK,Ngr),H(3,LINK,Ngr),YB(LINK,Ngr),FK(2,Ngr),FKD(2,Ngr)
	  COMMON/CGRN/ C_GR(LINK,500),TGR_C(LINK,500),IT_C(LINK)
	  COMMON/SGRN/ S_GR(LINK,500),TGR_S(LINK,500),IT_S(LINK)
	  COMMON/CSGRN/ CS_GR(LINK,500),TGR_CS(LINK,500),IT_CS(LINK)
	
	COMMON/GRN/   QGR(LINK,500),TGR(LINK,500),IT(LINK)
	COMMON/GRNY/  YGR(LINK,500),TGRY(LINK,500),ITY(LINK)
          WRITE(*,*)' READING HYD FILE. '
	lu=1
        OPEN(UNIT=lu,STATUS='OLD',FILE='hyd',err=13)
         write(*,'(A)') STR
        j=1
	 READ (lu,'(4x,g15.0)',ERR=12,end=13) th
         read (lu,*) LINKS,M  ! number of branches, number of nodes
         write(*,*) LINKS, M
         GOTO 14
	  
 12       write(*,*)' format violation'
          WRITE(*,*)' ERR in HYD FILE. '
          BACKSPACE(lu)
          READ (lu,'(A)')str
C          WRITE (*,'(A)')str
	    STOP
 13       write(*,*)' ABNORMAL END FILE    '
          STOP
	
 14     CONTINUE
   	  WRITE(*,*)' end first stage of HYD'
c --------   loop  --------------
 15	continue
	  READ (lu,*,ERR=12,END=23) IPO
        READ (lu,'(A)',ERR=12,END=13)str
               write(*,*)'IPO=',ipo
        DO 11 i=1,IPO
          READ (lu,fmt=*,ERR=12,END=13)
     +    K ,L, Q(1,k,l) ,Y(1,k,l),YB(k,l)
c          READ (lu,fmt=*,ERR=12,END=13)
c     +    K ,L, q ,ar, h
c          write(*,*) k,l, q(1,k,l),H(1,k,l),BW(1,k,l)
 11     continue
c---------  end loop  ------------          
 23 	write(*,*)'NORMAL END HYD FILE'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ÇÖêïçàÖ ÉêÄçàóçõÖ ìëãéÇàü (êÄëïéÑ)
        OPEN(UNIT=lu,STATUS='OLD',FILE='q_bc.dat')
	read(lu,*)kgr   ! kol-vo gr. uslov.
c	write(*,*)kgr   
	if(kgr.eq.0) goto 50
	DO j=1,kgr
	read(lu,*)num  ! branch number
	read(lu,*)npoint
	read(lu,*)IT(num) ! kol-vo vrem. tochek
c	write(*,*)'from hyd_inp',' IT=',IT(num),' num=',num
c	pause
	DO i=1,IT(num)
	read(lu,*)TGR(num,i),QGR(num,i)
c	write(*,*)'TGR=',TGR(num,i),' QGR=',QGR(num,i)
c	pause
	enddo
	ENDDO
	close(lu)
  50   	continue

!   çàÜçàÖ à ÇçìíêÖççàÖ ÉêÄçàóçõÖ  ìëãéÇàü (ìêéÇÖçú)
        OPEN(UNIT=lu,STATUS='OLD',FILE='y_bc.dat')
	read(lu,*)kgry
c	write(*,*)'kol-vo gr.usl. for Y =',kgry
c	pause
	if(kgry.eq.0) goto 60

	DO j=1,kgry
	read(lu,*)num
	read(lu,*)npoint
	read(lu,*)ITY(num)
c	write(*,*)'num=',num,' ITY(num)=',ITY(num)
c	pause
	DO i=1,ITY(num)
	read(lu,*)TGRY(num,i),YGR(num,i)
c	write(*,*)'T=',TGRY(num,i),' Y=',YGR(num,i)
c	write(*,*)'i=',i,' num=',num
c	pause
	enddo
	ENDDO
	close(lu)
 60 		continue
!    ÉêÄçàóçõÖ ìëãéÇàü (CONCENTRATION)
c        OPEN(UNIT=lu,STATUS='OLD',FILE='c_bc.dat')
c	read(lu,*)kgr_c   ! kol-vo gr. uslov. for concentr.
c	DO j=1,kgr_c
c	read(lu,*)num  ! branch number
c	read(lu,*)npoint
c	read(lu,*)IT_C(num) ! kol-vo vrem. tochek
c	DO i=1,IT_C(num)
c	read(lu,*)TGR_C(num,i),C_GR(num,i) ! conc. v Bq/m3
! 1 Bq=27.02 pCi
c	write(*,*)'do',C_GR(num,i)
c	C_GR(num,i)=C_GR(num,i)*1000/27.02 ! conc. v Bq/m3
c	write(*,*)'posle',C_GR(num,i)
c	pause
c	enddo
c	ENDDO
c	close(lu)
!    ÉêÄçàóçõÖ ìëãéÇàü (SEDIMENT)
c        OPEN(UNIT=lu,STATUS='OLD',FILE='s_bc.dat')
c	read(lu,*)kgr_s   ! kol-vo gr. uslov. for sediment
c	DO j=1,kgr_s
c	read(lu,*)num  ! branch number
c	read(lu,*)npoint
c	read(lu,*)IT_S(num) ! kol-vo vrem. tochek
c	DO i=1,IT_S(num)
c	read(lu,*)TGR_S(num,i),S_GR(num,i) 
c	enddo
c	ENDDO
c	close(lu)
!    ÉêÄçàóçõÖ ìëãéÇàü (CONC. IN SEDIMENT)
c        OPEN(UNIT=lu,STATUS='OLD',FILE='cs_bc.dat')
c	read(lu,*)kgr_cs   ! kol-vo gr. uslov. for sediment
c	DO j=1,kgr_cs
c	read(lu,*)num  ! branch number
c	read(lu,*)npoint
c	read(lu,*)IT_CS(num) ! kol-vo vrem. tochek
c	DO i=1,IT_CS(num)
c	read(lu,*)TGR_CS(num,i),CS_GR(num,i) 
c	enddo
c	ENDDO
c	close(lu)
c	write(*,*)'from HYD_INP'
c		write(*,*)'C_GR1=',(C_GR(1,i),i=1,IT_C(1))
c		write(*,*)'C_GR2=',(C_GR(2,i),i=1,IT_C(2))
c		write(*,*)'C_GR4=',(C_GR(4,i),i=1,IT_C(4))
c	pause
        return
	end
