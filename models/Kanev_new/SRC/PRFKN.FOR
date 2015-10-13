	SUBROUTINE PRFKN(L,N,J)
	INCLUDE 'parameters' 
	COMMON/BASE/ Q(3,LINK,Ngr),Y(3,LINK,Ngr),W(3,LINK,Ngr),
     *	Z(3,LINK,Ngr),
     *  BW(3,LINK,Ngr),H(3,LINK,Ngr),YB(LINK,Ngr),FK(2,Ngr),FKD(2,Ngr)
c   Rashodnaja harakteristika rusla: K=cW(h**1/2); c - coeff. Shezy

        Common/pchm/ d_part,dmu,dnu,z_obm,ER,grav,Zobm,delvv
	Common/SCAL/phi,teta,bet,eps,delt,delx,fi,fn
c  c=1/n(R**1/6), n=0.02 - 0.04; R =~ h - glubina
                 DO i=1,N
c         write(*,*)'prfkn works'
	hreal=H(J,L,i)
	hmean=W(J,L,i)/BW(J,L,i)
c	write(*,*)'hmean=',hmean,' h=',hreal
c	write(*,*)'L=',l,' i=',i
c           FK(J,i)=1./fn*W(J,L,i)*(hreal**(2./3.))  
c           FKD(J,i)=(5./3.)*BW(J,L,i)/0.02*(hreal**(2./3.))
c	write(*,*)'h=',hmean,' W=',W(J,L,i),' BW=',BW(J,L,i)
	if(hmean.le.0) then
	write(*,*)'from prfkn'
	write(*,*)'h=',hmean,' i=',i
	pause
	write(*,*)'Y=',(Y(2,2,k),k=1,15)
	pause
	write(*,*)'Q=',(Q(2,2,k),k=1,15)
	pause
	endif
           FK(J,i)=1./fn*W(J,L,i)*(hmean**(2./3.))  
           FKD(J,i)=(5./3.)*BW(J,L,i)/fn*(hmean**(2./3.))
c	write(*,*)'i=',i
c	write(*,*)'FK=',FK(J,i)
c        pause
                  END DO
c        write(*,*)'prfkn finished'
         RETURN
         END