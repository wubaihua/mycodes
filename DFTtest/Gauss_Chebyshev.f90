subroutine chebyshev(N,P,w,R)
    implicit real*8(a-h,o-z)
    integer N !number of grids
    real*8 P ! parameter for tansform coordinate
    real*8 w(N),R(N),x(N)
    real*8,parameter :: pi=3.1415926535897932
    
    do i=1,N
        x(i)=cos(i*pi/(N+1))
        w(i)=(2*pi*P**3*(1+x(i))**2.5)/((N+1)*(1-x(i))**3.5)
        R(i)=((1+x(i))*P)/(1-x(i))
    end do
    
    
    
    
end subroutine
