	SUBROUTINE COEFN(LINKS,LMAX,istep)

	INCLUDE 'parameters' 

c SUBROUTINE COMPUTE COEFF. A(i),B(i),...,G1(i);
c E(i),F(i),H(i),..., H1(i) FOR ALL BRABCHES AND GRIDS
C LINKS - number of branches
C    LMAX=maxIL(L) - max number of points on the  link
C COE(L,i) - Array for linearised equatons coefficiants A(i),B(i),...,
C i=1,IL(L)-1, L - branch number, A(1),B(1),...,G1(1);A(2),B(2),... 
c RCOE(L,i) - Array for coefficients E(i),F(i),H(i),E1(i),F1(i),H1(i)

	COMMON/BASE/ Q(3,LINK,Ngr),Y(3,LINK,Ngr),W(3,LINK,Ngr),
     *	Z(3,LINK,Ngr),
     *  BW(3,LINK,Ngr),H(3,LINK,Ngr),YB(LINK,Ngr),FK(2,Ngr),FKD(2,Ngr)
	COMMON/COEF/ COE(LINK,Ngr,10),RCOE(LINK,Ngr,6),IL(LINK)
        Common/pchm/ d_part,dmu,dnu,z_obm,ER,grav,Zobm,delvv
	Common/SCAL/phi,teta,bet,eps,delt,delx,fi,fn

        REAL CF(10,LMAX),CQ(2,LMAX),CW(2,LMAX),CBW(2,LMAX),
     *   CY(2,LMAX),CE(6,LMAX) 

c CF(1,i) - linearised equatons coefficiants A(i),
c CF(2,i) - B(i),...,CF(10,i) - G1(i), i=1,IL(L)-1
c CE(1,i) - E(i), CE(2,i) - F(i), CE(3,i) - H(i)
c CE(4,i) - E1(i), CE(5,i) - F1(i), CE(6,i) - H1(i)
c rashody, urovni i t.d. dla brancha L  
c CFK and CFKD - rashodnye harakteristiky rusla

c 		write(*,*)'coefn works'
c	write(*,*)'LMAX=',LMAX
c	pause
	do i=1,LMAX
	do j=1,10
			CF(j,i)=0
	enddo
	enddo
	do i=1,LMAX
	do j=1,6
			CE(j,i)=0
	enddo
	enddo
	do i=1,LMAX
	do j=1,2
			CQ(j,i)=0
			CW(j,i)=0
			CBW(j,i)=0
			CY(j,i)=0
	enddo
	enddo
                DO L=1,LINKS !! nachalo cykla po L
c	write(*,*)'L=',L	
             N=IL(L) ! number of points on branch L
             DO J=1,2
c	write(*,*)'j=',j
       CALL PRFKN(L,N,J) ! compute rash. harakt. rusla i proizv.
             DO i=1,IL(L)
! vychisl. massivov dla dannogo brancha
              CQ(j,i)=Q(j,L,i)  ! rashody
              CY(j,i)=Y(j,L,i)  ! urovni
              CW(j,i)=W(j,L,i)  ! crossections
              CBW(j,i)=BW(j,L,i) ! width
            END DO ! po i
            END DO ! po J
	 j=2
C    Compute coefficients A(i),B(i),...; i=1,IL(L)-1
                     alfa=1
                 DO i=1,IL(L)-1  ! cykl po i dla vychisl. coef.
c	write(*,*)'i=',i
      CF(1,i)=phi*CBW(2,i+1)/delt    ! A(i)
      CF(2,i)=teta/delx              ! B(i)
      CF(3,i)=(1.-phi)*CBW(2,i)/delt !C(i)
      CF(4,i)=-CF(2,i)               !D(i)
      CF(5,i)=phi/delt*(CW(2,i+1)-CW(1,i+1))+(1-phi)/delt*
     * (CW(2,i)-CW(1,i))+            !G(i)  
     * teta/delx*(CQ(2,i+1)-CQ(2,i))+(1-teta)/delx*(CQ(1,i+1)-CQ(1,i)) 
c	write(*,*)'phi=',phi,' delt=',delt
      CF(6,i)=-alfa*teta*CBW(2,i+1)/(CW(2,i+1)**2)*CQ(2,i+1)*(teta/delx*
     1 (CQ(2,i+1)-CQ(2,i))+(1-teta)/delx*(CQ(1,i+1)-CQ(1,i)))
     2 -alfa*teta*CBW(2,i+1)/delx*     !!         A1(i)
     2 (teta/4*(CQ(2,i)/CW(2,i)+CQ(2,i+1)/CW(2,i))**2+(1-teta)/4*
     2 (CQ(1,i)/CW(1,i)+CQ(1,i+1)/CW(1,i+1))**2) !
     2 +alfa*teta*CBW(2,i+1)/2/(CW(2,i+1)**2)*CQ(2,i+1)*(CQ(2,i)/
     2 CW(2,i)+
     2 CQ(2,i+1)/CW(2,i+1))*(teta/delx*(CW(2,i+1)-CW(2,i))+(1-teta)/
     2 delx*(CW(1,i+1)-CW(1,i)))
     3 +teta*grav/delx*(teta/2*(CW(2,i)+CW(2,i+1))+(1-teta)/2*
     3 (CW(1,i)+CW(1,i+1)))+teta*grav*CBW(2,i+1)/2*
     3 (teta/delx*(CY(2,i+1)-CY(2,i))+(1-teta)/delx*(CY(1,i+1)-CY(1,i)))
     4 -teta*grav*(1-bet)*FKD(2,i+1)/(FK(2,i+1)**3)*CQ(2,i+1)*
     4 ABS(CQ(2,i+1))
     4 *(teta*(CW(2,i)+CW(2,i+1))+(1-teta)*(CW(1,i)+CW(1,i+1)))
      CF(6,i)=CF(6,i)+teta*grav*CBW(2,i+1)/2*(teta*(bet*CQ(2,i)*
     4 ABS(CQ(2,i))/
     4 (FK(2,i)**2)+(1-bet)*CQ(2,i+1)*ABS(CQ(2,i+1))/(FK(2,i+1)**2))
     4 +(1-teta)*(bet*CQ(1,i)*ABS(CQ(1,i))/(FK(1,i)**2)+
     4 (1-bet)*CQ(1,i+1)*ABS(CQ(1,i+1))/(FK(1,i+1)**2))) !!  A1(i)
      CF(7,i)=phi/delt+teta/delx*alfa*(teta*(CQ(2,i)/CW(2,i)+CQ(2,i+1)/
     2 CW(2,i+1))+(1-teta)*(CQ(1,i)/CW(1,i)+CQ(1,i+1)/CW(1,i+1)))
     2 +alfa*teta/CW(2,i+1)*(teta/delx*(CQ(2,i+1)-CQ(2,i))+(1-teta)/
     2 delx*(CQ(1,i+1)-CQ(1,i)))
     3 -alfa*teta/2/CW(2,i+1)*(CQ(2,i)/CW(2,i)+CQ(2,i+1)/CW(2,i+1))*
     3 (teta/delx*(CW(2,i+1)-CW(2,i))+(1-teta)/delx*(CW(1,i+1)-CW(1,i)))
     5 +grav*teta*(1-bet)*ABS(CQ(2,i+1))/(FK(2,i+1)**2)
     5 *(teta*(CW(2,i)+CW(2,i+1))+(1-teta)*(CW(1,i)+CW(1,i+1))) !! B1(i)
       CF(8,i)=-alfa*teta*CBW(2,i)*CQ(2,i)/(CW(2,i)**2)*(teta/delx*
     2 (CQ(2,i+1)-CQ(2,i))+(1-teta)/delx*(CQ(1,i+1)-CQ(1,i)))
     3 +alfa*teta*CBW(2,i)/delx*(teta/4*(CQ(2,i)/CW(2,i)+   !
     3 CQ(2,i+1)/CW(2,i+1))**2+(1-teta)/4*(CQ(1,i)/CW(1,i)+ !
     3 CQ(1,i+1)/CW(1,i+1))**2)+alfa*teta/2* !
     3 (CBW(2,i)*CQ(2,i)**2/CW(2,i)**3+CBW(2,i)*CQ(2,i)*CQ(2,i+1)/ !
     3 (CW(2,i+1)*CW(2,i)**2))*  !
     3 (teta/delx*(CW(2,i+1)-CW(2,i))+  !
     3 (1-teta)/delx*(CW(1,i+1)-CW(1,i)))  !
       CF(8,i)=CF(8,i)-teta*grav/delx*(teta/2*(CW(2,i+1)+CW(2,i+1))  !
     4 +(1-teta)/2*(CW(1,i)+CW(1,i+1)))+  !
     4 teta*grav*CBW(2,i)/2/delx*(teta*(CY(2,i+1)-CY(2,i))+(1-teta)*
     4 (CY(1,i+1)-CY(1,i)))
     5 -grav*teta*bet*FKD(2,i)*CQ(2,i)*ABS(CQ(2,i))/
     5 (FK(2,i)**3)*(teta*(CW(2,i)+CW(2,i+1))+     !
     5 (1-teta)*(CW(1,i)+CW(1,i+1)))    !
      CF(8,i)=CF(8,i)+teta*grav*CBW(2,i)/2*(teta*(bet*CQ(2,i)*
     5 ABS(CQ(2,i))/(FK(2,i)**2)+ !
     5 (1-bet)*CQ(2,i+1)*ABS(CQ(2,i+1))/(FK(2,i+1)**2))+
     5 (1-teta)*(bet*CQ(1,i)*ABS(CQ(1,i))/(FK(1,i)**2)
     5 +(1-bet)*CQ(1,i+1)*ABS(CQ(1,i+1))/(FK(1,i+1)**2))) !! C1(i)
       CF(9,i)=(1-phi)/delt+(-teta/delx*alfa*(teta*(CQ(2,i)/CW(2,i)+
     1 CQ(2,i+1)/CW(2,i+1))+(1-teta)*(CQ(1,i)/CW(1,i)+CQ(1,i+1)/
     1 CW(1,i+1)))
     1 +alfa*teta/CW(2,i)*(teta/delx*(CQ(2,i+1)-CQ(2,i))+(1-teta)/
     2 delx*(CQ(1,i+1)-CQ(1,i)))-alfa*teta/2/CW(2,i)*(CQ(2,i)/CW(2,i)+
     2 CQ(2,i+1)/CW(2,i+1))*
     3 (teta/delx*(CW(2,i+1)-CW(2,i))+
     3 (1-teta)/delx*(CW(1,i+1)-CW(1,i))))+
     4 grav*teta*bet*ABS(CQ(2,i))/(FK(2,i)**2)*(teta*(CW(2,i)+
     4 CW(2,i+1))+
     5 (1-teta)*(CW(1,i)+CW(1,i+1))) !!       D1(i)
      CF(10,i)=phi/delt*(CQ(2,i+1)-CQ(1,i+1))+(1-phi)/delt*(CQ(2,i)-
     1 CQ(1,i))+
     1 alfa*(teta*(CQ(2,i)/CW(2,i)+CQ(2,i+1)/CW(2,i+1))+(1-teta)*
     2 (CQ(1,i)/CW(1,i)+
     2 CQ(1,i+1)/CW(1,i+1)))*(teta/delx*(CQ(2,i+1)-CQ(2,i))+(1-teta)/
     3 delx*
     3 (CQ(1,i+1)-CQ(1,i)))-alfa*(teta/4*(CQ(2,i)/CW(2,i)+CQ(2,i+1)
     3 /CW(1,i+1))**2+
     4 (1-teta)/4*(CQ(1,i)/CW(1,i)+CQ(1,i+1)/CW(1,i+1))**2)*(teta/delx*
     5 (CW(2,i+1)-CW(2,i))+(1-teta)/delx*
     5 (CW(1,i+1)-CW(1,i)))+grav*(teta/2*(CW(2,i)+CW(2,i+1))
     6 +(1-teta)/2*(CW(1,i)+CW(1,i+1)))*
     7 (teta/delx*(CY(2,i+1)-CY(2,i))+(1-teta)/delx*(CY(1,i+1)-CY(1,i)))
      CF(10,i)=CF(10,i)+ grav*(teta/2*(CW(2,i)+CW(2,i+1))+(1-teta)/2*
     1 (CW(1,i)+CW(1,i+1)))*
     1 (teta*(bet*CQ(2,i)*ABS(CQ(2,i))/(FK(2,i)**2)+(1-bet)*CQ(2,i+1)*
     1 ABS(CQ(2,i+1))/(FK(2,i+1)**2))+
     2 (1-teta)*(bet*CQ(1,i)*ABS(CQ(1,i))/(FK(1,i)**2)+(1-bet)*
     3 CQ(1,i+1)*ABS(CQ(1,i+1))/(FK(1,i+1)**2)))  !! G1(i)
		do iii=3,5
		CF(iii,i)=-CF(iii,i)
		enddo
		do iii=8,10
		CF(iii,i)=-CF(iii,i)
		enddo
                           END DO ! konec cykla po i
c        pause'coefn. A,B,... calc.'
c----------------------------------------------------------------------
c             COMPUTE RECURSION RELATION COEFF E,F,--- ETC.
               DO ij=1,IL(L)-1
                i=IL(L)-ij
            IF(IJ.GT.1) GO TO 50
! zadanie nachalnyh znachenij pri i=IL(L)-1
            DENOM=-CF(2,i)*CF(9,i)+CF(7,i)*CF(4,i)
            CE(1,i)=(CF(8,i)*CF(2,i)-CF(3,i)*CF(7,i))/DENOM ! E(IL-1)
            CE(2,i)=(CF(10,i)*CF(2,i)-CF(5,i)*CF(7,i))/DENOM ! F(IL-1)
            CE(3,i)=(CF(1,i)*CF(7,i)-CF(6,i)*CF(2,i))/DENOM ! H(IL-1) 
            DENOM=-CF(4,i)*CF(7,i)+CF(9,i)*CF(2,i) 
            CE(4,i)=(CF(3,i)*CF(9,i)-CF(8,i)*CF(4,i))/DENOM ! E1(IL-1)
            CE(5,i)=(CF(5,i)*CF(9,i)-CF(10,i)*CF(4,i))/DENOM ! F1(IL-1)
            CE(6,i)=(CF(6,i)*CF(4,i)-CF(1,i)*CF(9,i))/DENOM ! H1(IL-1)
             GO TO 90
 50          DENOM=CF(1,i)*CF(7,i)-CF(6,i)*CF(2,i)
            RL=(CF(3,i)*CF(7,i)-CF(8,i)*CF(2,i))/DENOM
            RMM=(CF(4,i)*CF(7,i)-CF(9,i)*CF(2,i))/DENOM
            RN=(CF(5,i)*CF(7,i)-CF(10,i)*CF(2,i))/DENOM
            DENOM=(CF(1,i)+CF(2,i)*CE(1,i+1))*RMM-CF(4,i)
            CE(1,i)=(CF(3,i)-RL*(CF(1,i)+CF(2,i)*CE(1,i+1)))/DENOM
            CE(2,i)=(CF(5,i)-(CF(1,i)+CF(2,i)*CE(1,i+1))*RN-
     *      CF(2,i)*CE(2,i+1))/DENOM
            CE(3,i)=-CF(2,i)*CE(3,i+1)/DENOM
            CE(4,i)=CE(4,i+1)*(RL+RMM*CE(1,i))
            CE(5,i)=CE(5,i+1)+CE(4,i+1)*(RMM*CE(2,i)+RN)
            CE(6,i)=CE(6,i+1)+CE(3,i)*CE(4,i+1)*RMM
 90       CONTINUE
                    END DO    ! po ij
C        WRITING CALCULATED COEFFICIENTS CF(J,i) (A(i),B(i),...) 
C        TO ARRAY COE(L,i,*) AND COEF. CE(J,i) (E(i),F(i),H(i),...)
C        TO ARRAY RCOE(L,i,*)

		DO i=1,IL(L)
			do j=1,10
			COE(L,i,j)=CF(j,i)
			enddo
		do k=1,6
		RCOE(L,i,k)=CE(k,i)
		enddo
 		  ENDDO ! po i
		  ENDDO ! po L

         RETURN
         END
