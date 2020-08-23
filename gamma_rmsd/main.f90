!***********************************************
!    A program for calculate error
!    author: Baihua Wu (wubaihua@pku.edu.cn)   
!***********************************************
   
    
    
program gamma_ad
    implicit real*8(a-h,o-z)
    real*8 rmsd
    real*8,allocatable :: data_test(:,:),data_benchmark(:,:),&
                            x_test(:),x_benchmark(:),&
                            itp_test(:,:),itp_benchmark(:,:),&
                            x_itp(:),error(:,:) 
    character*500 c1,c2,cname
    
    write(*,"(a)")"***********************************************"
    write(*,"(a)")"A program for calculate error"
    write(*,"(a)")"author: Baihua Wu (wubaihua@pku.edu.cn)"   
    write(*,"(a)")"***********************************************"
    
    write(*,*) "input the job name:"
    read(*,*) cname
    open(99,file="result_"//trim(cname)//".dat")
    write(*,*) "input the number of states:"
    read(*,*) n_site
    
    write(*,*) "input the benchmark file:"
    read(*,*) c2
    open(200,file=trim(c2))
    rewind(200)
    nx_bm=-1
    do while(.true.)
        read(200,*,iostat=ierr)
        nx_bm=nx_bm+1
        if(ierr<0)then
            exit
        end if
    end do
    rewind(200)
    
    allocate(x_benchmark(nx_bm))
    allocate(data_benchmark(nx_bm,n_site))
    do i=1,nx_bm
        read(200,*) x_benchmark(i),data_benchmark(i,:)
    end do
    
    write(*,*) "input the x_start, x_end, N_grid"
    read(*,*) x_s,x_e,n_grid
    allocate(x_itp(n_grid))
    allocate(itp_benchmark(n_grid,n_site))
    do i=1,n_grid
        x_itp(i)=x_s+(i-1)*(x_e-x_s)/(n_grid-1)
    end do
    do i=1,n_site
        call ITP(n_grid,x_itp,itp_benchmark(:,i),nx_bm,x_benchmark,data_benchmark(:,i))
    end do
    
    !open(121,file="itp_benchmark.dat")
    !do i=1,n_grid
    !    write(121,"(<n_site+1>E18.8E3)") x_itp(i),itp_benchmark(i,:)
    !end do
    !stop
    
171 write(*,*) "input the gamma:"
    read(*,*) gamma
    
    write(*,*) "input the test file:"
    read(*,*) c1
    open(100,file=trim(c1))
    rewind(100)
    nx_te=-1
    do while(.true.)
        read(100,*,iostat=ierr)
        nx_te=nx_te+1
        if(ierr<0)then
            exit
        end if
    end do
    rewind(100)
    
    allocate(x_test(nx_te))
    allocate(data_test(nx_te,n_site))
    do i=1,nx_te
        read(100,*) x_test(i),data_test(i,:)
    end do
    allocate(itp_test(n_grid,n_site))
    allocate(error(n_grid,n_site))
    do i=1,n_site
        call ITP(n_grid,x_itp,itp_test(:,i),nx_te,x_test,data_test(:,i))
    end do
    
    !error=abs(itp_test-itp_benchmark)
    rmsd=0
    do i=1,n_grid
        do j=1,n_site
            rmsd=rmsd+(itp_test(i,j)-itp_benchmark(i,j))**2
        end do
    end do
    rmsd=sqrt(rmsd/(n_grid*n_site))
        
    
    

    
    write(99,*) gamma,rmsd
!    
    write(*,*) "continue? 1 for yes, 0 for no:"
    read(*,*) iif
    select case(iif)
    case(1)
        close(100)
        deallocate(x_test)
        deallocate(data_test)
        deallocate(itp_test)
        deallocate(error)
        goto 171
    case(0)
        close(100)
        close(200)
        close(99)
    end select
    
        
    
        
    
    
    
    
    
    
    
    
    
    
    
end program

!program test
!    implicit real*8(a-h,o-z)
!    real*8 x0(5),y0(5),x(6),y(6)
!    
!    x0=(/0.0_8, 3.0_8, 5.0_8, 8.5_8, 10.0_8/)
!    y0=(/1.0_8, 5.4_8, 8.3_8, 11.5_8, 20.0_8/)
!    x=(/0.0_8, 2.0_8, 4.0_8, 6.0_8, 8.0_8, 10.0_8/)
!    
!    call ITP(6,x,y,5,x0,y0)
!    open(11,file="x0y0.dat")
!    open(12,file="xy.dat")
!
!    do i=1,5
!        write(11,*) x0(i),y0(i)
!    end do
!    
!    do i=1,6
!        write(12,*) x(i),y(i)
!    end do
!
!end program




subroutine ITP(n,x,y,n0,x0,y0)
    implicit real*8(a-h,o-z)
    integer n,n0
    real*8 x(n),y(n),x0(n0),y0(n0)
    
    
    
        
    
    
    do i=1,n
        do j=1,n0-1
            if(x0(j)<=x(i) .and. x0(j+1)>=x(i))then
                if(x0(j)==x(i))then
                    y(i)=y0(j)
                elseif(x0(j+1)==x(i))then
                    y(i)=y0(j+1)
                else
                    y(i)=y0(j)+(y0(j+1)-y0(j))*(x(i)-x0(j))/(x0(j+1)-x0(j))
                end if
                exit
            end if
        end do
    end do
    
    
  
end subroutine
