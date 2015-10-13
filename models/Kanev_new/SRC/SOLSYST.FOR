	SUBROUTINE SOLSYST(M)
	INCLUDE 'parameters' 
	COMMON/SWEEP/ LIN(LINK,2),DNY(node)
	COMMON/MATR/ A(NODE,NODE),B(NODE),MM(node,2),
     *               M_IN(node,LINK),MOUT(node,LINK)
	COMMON/SOLUT/MGRUP(node),MDAM(node,2),NGRUP(node,2),NDAM
c   M - number of nodes
c   LINKS - number of branches
c   NDAM - kol-vo plotin
c   MGRUP(*) - numbers of dams nodes and last node
C   subr. computes array DNY, level elevation at nodes
c	write(*,*)'sols works'
c	do i=1,M
c	write(*,*)'i=',i
c	write(*,*)(A(i,j),j=1,M),' B(i)=',B(i)
c	pause
c	enddo
			DO igrup=1,NDAM
			NY1=NGRUP(igrup,1) 
			NY2=NGRUP(igrup,2) 
		DO k=NY1,NY2
	if(A(k,k).ne.0) then
		DO i=k+1,NY2
		if(A(i,k).NE.0) then
		alfa=A(i,k)/A(k,k)
			do j=k,NY2
			A(i,j)=A(i,j)-alfa*A(k,j)
			enddo ! po j
		B(i)=B(i)-alfa*B(k) 
		endif
		ENDDO   ! po i
		endif
		ENDDO ! po k
c-------------------------------------------
		DNY(NY2)=B(NY2)/A(NY2,NY2)
		DO i=1,NY2-1
			j=NY2-i
		sum=0
			do k=j+1,NY2
			sum=sum+A(j,k)*DNY(k)
			enddo   ! po k   
		DNY(j)=(B(j)-sum)/A(j,j)
		ENDDO ! po i
			ENDDO   ! po igrup
c-------------------------------------------
c	write(*,*)'DNY=',(DNY(i),i=1,M)
c	pause
c	write(*,*)'sols finished'
         RETURN
         END