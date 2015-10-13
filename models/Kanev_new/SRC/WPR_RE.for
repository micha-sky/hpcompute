         SUBROUTINE WPR_RE(delY,LINKS)

	INCLUDE 'parameters' 

        COMMON/CROSS/ WPR(LINK,Ngr,n_cross)
	COMMON/COEF/ COE(LINK,Ngr,10),RCOE(LINK,Ngr,6),IL(LINK)
	COMMON /YBS/YBASE(LINK)

c           CROSSECTION READING (WPR)
c---------------- 1 branch --------------------------------------
        open(1,file='wpr1.txt')  !! L=1
        read(1,*)L,npoint,ny,YBASE(L),delY
c L - номер бранча, npoint - кол. точек в бранче, 
c ny - кол. слоев по глубине, YBASE(L) - базовый уровень бранча
c        write(*,*)'L=',L,' npoint=',npoint,' ny=',ny,
c     *   ' YBASE(L)=',YBASE(L),' delY=',delY
c        pause
        	do i=1,npoint
        	do j=1,ny
	read(1,*)ii,jj,WPR(L,i,j)
        	enddo
        	enddo
         close(1)
			do i=1,npoint
			do ik=1,30
		WPR(L,i,ny+ik)=2*WPR(L,i,ny+ik-1)-WPR(L,i,ny+ik-2)
			enddo
			enddo
c---------------- 2 branch -------------------------------------------        
        open(1,file='wpr2.txt')  !! L=2
        read(1,*)L,npoint,ny,YBASE(L),delY
c L - номер бранча, npoint - кол. точек в бранче, 
c ny - кол. слоев по глубине, YBASE(L) - базовый уровень бранча
c        write(*,*)'L=',L,' npoint=',npoint,' ny=',ny,
c     *   ' YBASE(L)=',YBASE(L),' delY=',delY
c        pause
        	do i=1,npoint
        	do j=1,ny
	read(1,*)ii,jj,WPR(L,i,j)
        	enddo
        	enddo
         close(1)
			do i=1,npoint
			do ik=1,30
		WPR(L,i,ny+ik)=2*WPR(L,i,ny+ik-1)-WPR(L,i,ny+ik-2)
			enddo
			enddo
c---------------- 3 branch -------------------------------------------        
c        open(1,file='wpr3.txt')  !! L=3
c        read(1,*)L,npoint,ny,YBASE(L),delY
c L - номер бранча, npoint - кол. точек в бранче, 
c ny - кол. слоев по глубине, YBASE(L) - базовыйауровень бранча
c        write(*,*)'L=',L,' npoint=',npoint,' ny=',ny,
c     *   ' YBASE(L)=',YBASE(L),' delY=',delY
c        pause
c        	do i=1,npoint
c        	do j=1,ny
c	read(1,*)ii,jj,WPR(L,i,j)
c        	enddo
c        	enddo
c         close(1)
c			do i=1,npoint
c			do ik=1,30
c		WPR(L,i,ny+ik)=2*WPR(L,i,ny+ik-1)-WPR(L,i,ny+ik-2)
c			enddo
c			enddo
c---------------- 4 branch -------------------------------------------        
c        open(1,file='wpr4.txt')  !! L=4
c        read(1,*)L,npoint,ny,YBASE(L),delY
c L - номер бранча, npoint - кол. точек в бранче, 
c ny - кол. слоев по глубине, YBASE(L) - базовый уровень бранча
c        write(*,*)'L=',L,' npoint=',npoint,' ny=',ny,
c     *   ' YBASE(L)=',YBASE(L),' delY=',delY
c        pause
c        	do i=1,npoint
c        	do j=1,ny
c	read(1,*)ii,jj,WPR(L,i,j)
c        	enddo
c        	enddo
c         close(1)
c			do i=1,npoint
c			do ik=1,30
c		WPR(L,i,ny+ik)=2*WPR(L,i,ny+ik-1)-WPR(L,i,ny+ik-2)
c			enddo
c			enddo
c---------------- 5 branch -------------------------------------------        
c        open(1,file='wpr5.txt')  !! L=5
c        read(1,*)L,npoint,ny,YBASE(L),delY
c L - номер бранча, npoint - кол. точек в бранче, 
c ny - кол. слоев по глубине, YBASE(L) - базовый уровень бранча
c        write(*,*)'L=',L,' npoint=',npoint,' ny=',ny,
c     *   ' YBASE(L)=',YBASE(L),' delY=',delY
c        pause
c        	do i=1,npoint
c        	do j=1,ny
c	read(1,*)ii,jj,WPR(L,i,j)
c        	enddo
c        	enddo
c         close(1)
c			do i=1,npoint
c			do ik=1,30
c		WPR(L,i,ny+ik)=2*WPR(L,i,ny+ik-1)-WPR(L,i,ny+ik-2)
c			enddo
c			enddo
c---------------- 6 branch -------------------------------------------        
c        open(1,file='wpr6.txt')  !! L=6
c        read(1,*)L,npoint,ny,YBASE(L),delY
c L - номер бранча, npoint - кол. точек в бранче, 
c ny - кол. слоев по глубине, YBASE(L) - базовый уровень бранча
c        write(*,*)'L=',L,' npoint=',npoint,' ny=',ny,
c     *   ' YBASE(L)=',YBASE(L),' delY=',delY
c        pause
c        	do i=1,npoint
c        	do j=1,ny
c	read(1,*)ii,jj,WPR(L,i,j)
c        	enddo
c        	enddo
c         close(1)
c			do i=1,npoint
c			do ik=1,30
c		WPR(L,i,ny+ik)=2*WPR(L,i,ny+ik-1)-WPR(L,i,ny+ik-2)
c			enddo
c			enddo
c---------------- 7 branch -------------------------------------------        
c        open(1,file='wpr7.txt')  !! L=7
c        read(1,*)L,npoint,ny,YBASE(L),delY
c L - номер бранча, npoint - кол. точек в бранче, 
c ny - кол. слоев по глубине, YBASE(L) - базовый уровень бранча
c        write(*,*)'L=',L,' npoint=',npoint,' ny=',ny,
c     *   ' YBASE(L)=',YBASE(L),' delY=',delY
c        	do i=1,npoint
c        	do j=1,ny
c	read(1,*)ii,jj,WPR(L,i,j)
c        	enddo
c        	enddo
c         close(1)
c			do i=1,npoint
c			do ik=1,30
c		WPR(L,i,ny+ik)=2*WPR(L,i,ny+ik-1)-WPR(L,i,ny+ik-2)
c			enddo
c			enddo
c---------------- 8 branch -------------------------------------------        
c        open(1,file='wpr8.txt')  !! L=8
c        read(1,*)L,npoint,ny,YBASE(L),delY
c L - номер бранча, npoint - кол. точек в бранче, 
c ny - кол. слоев по глубине, YBASE(L) - базовый уровень бранча
c        write(*,*)'L=',L,' npoint=',npoint,' ny=',ny,
c     *   ' YBASE(L)=',YBASE(L),' delY=',delY
c        	do i=1,npoint
c        	do j=1,ny
c	read(1,*)ii,jj,WPR(L,i,j)
c        	enddo
c        	enddo
c         close(1)
c			do i=1,npoint
c			do ik=1,30
c		WPR(L,i,ny+ik)=2*WPR(L,i,ny+ik-1)-WPR(L,i,ny+ik-2)
c			enddo
c			enddo
           RETURN
           END
