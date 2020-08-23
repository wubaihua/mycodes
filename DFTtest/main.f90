program dfttest
    implicit real*8(a-h,o-z)
    real*8,allocatable :: x(:),y(:),z(:),wa(:),wr(:),R(:)
    real*8,parameter :: pi=3.1415926535897932
    
    NA=5810
    NR=1000
    allocate(x(NA))
    allocate(y(NA))
    allocate(z(NA))
    allocate(wa(NA))
    allocate(R(NR))
    allocate(wr(NR))
    !
    call LD5810(X,Y,Z,Wa,N) 
    !
    !!write(*,*) "x=",x
    !!write(*,*) "y=",y
    !!write(*,*) "z=",z
    !!write(*,*) "w=",w
    !!write(*,*) sum(w)
    !
    !c=0
    !do i=1,N
    !    c=c+w(i)*(x(i)**2+y(i)**2)
    !end do
    !c=c*4*pi
    !e=8*pi/3
    !
    !write(*,*) c,e
    
    P=1
    call chebyshev(NR,P,wr,R)
    c=0
    do i=1,NA
        do j=1,NR
            c=c+4*pi*wr(j)*wa(i)*(R(j)*x(i))**2*exp(-(R(j)*x(i))**2-(R(j)*y(i))**2-(R(j)*z(i))**2)
        end do
    end do
    
    e=0.5*pi**1.5
    
    write(*,*) c,e
    
    
    
pause   
end program
