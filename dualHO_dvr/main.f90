program dvr_dual_ho
    use LAPACK95
    implicit real*8(a-h,o-z)
    real*8,allocatable :: T(:,:),V(:,:),H(:,:),E(:),R(:),work(:),lambda(:,:),U(:,:)
    complex*8,allocatable :: pro(:,:),c(:),c0(:)
    complex*8,parameter :: im=(0.0_8,1.0_8)
    real*8,parameter :: pi=3.14159265358979
    real*8 mass
    
    
    open(21,file="result.dat")
    
    Nstate=2
    dx=0.01
    x_min=-10
    x_max=10
    ngrid=int((x_max-x_min)/dx)
    
    allocate(T(Nstate*Ngrid,Nstate*Ngrid))
    allocate(V(Nstate*Ngrid,Nstate*Ngrid))
    allocate(H(Nstate*Ngrid,Nstate*Ngrid))
    allocate(lambda(Nstate*Ngrid,Nstate*Ngrid))
    allocate(U(Nstate*Ngrid,Nstate*Ngrid))
    allocate(c(Nstate*Ngrid))
    allocate(c0(Nstate*Ngrid))
    allocate(pro(Nstate*Ngrid,Nstate*Ngrid))
    allocate(E(Nstate*Ngrid))
    allocate(R(Ngrid))
    allocate(work(3*Nstate*Ngrid))
    
    do i=1,Ngrid
        R(i)=x_min+(i-1)*dx
    end do
    
    
    omega1=1
    omega2=2
    eps=1
    delta=1
    mass=1
    q=1
    p0=0
    
    ttot=10
    dt=0.01
    
    T=0
    V=0
    
    c=0
    
    
    do i=1,Ngrid
        V(i,i)=0.5*omega1**2*(R(i)+q)**2+eps
        !V(i,i)=0.5*R(i)**2
        V(i+Ngrid,i+Ngrid)=0.5*omega2**2*R(i)**2-eps
        V(i+Ngrid,i)=delta
        V(i,Ngrid+i)=delta
        !write(21,"(5E18.8)") R(i), V(i,i), V(i+Ngrid,i+Ngrid), V(i+Ngrid,i), V(i,Ngrid+i)
        c(i)=(omega1/pi)**0.25*exp(-omega1*(R(i)+q)**2/2)*exp(-im*p0*(R(i)+q))
        !write(21,*) R(i),real(c(i)),imag(c(i))
        T(i,i)=-2
        T(i+Ngrid,i+Ngrid)=-2
        if(i>1)then
            T(i,i-1)=1
            T(i-1,i)=1
            T(i+Ngrid,i-1+Ngrid)=1
            T(i-1+Ngrid,i+Ngrid)=1
        end if
    end do
    !do i=1,Ngrid
    !    write(*,*) V(i,:)
    !end do
    
    H=-T/(2*mass*dx**2)+V
    U=H
    call dsyev('V','L',Nstate*Ngrid,U,Nstate*Ngrid,E,work,3*Nstate*Ngrid,info)
    lambda=0
    !do i=1,Nstate*Ngrid
    !    lambda(i,i)=E(i)
    !end do
    
    c0=c
    
    
    !write(*,*) E(1:12)
    !write(*,*) E(1+Ngrid:6+Ngrid)
    !
    !do i=1,Ngrid
    !    
    !    write(21,"(11E18.8)") R(i),H(i,1:10)
    !end do
    !!do i=1,10
    !    write(*,*) E(:)
    !!end do
    time=0
    rho1=sum(real(c(1:Ngrid))**2+imag(c(1:Ngrid))**2)*dx
    rho2=sum(real(c(1+Ngrid:2*Ngrid))**2+imag(c(1+Ngrid:2*Ngrid))**2)*dx
    write(21,*) time,rho1,rho2
    
    
    
    do i=1,int(ttot/dt)
        !c=c-im*matmul(H,c)*dt
        !cd=cd-im*matmul(lambda,cd)*dt
        time=time+dt
        pro=0
        do k=1,Nstate*Ngrid
            pro(k,k)=cos(E(k)*time)-im*sin(E(k)*time)
        end do
        pro=matmul(U,matmul(pro,transpose(U)))
        c=matmul(pro,c0)
        !c=matmul(transpose(U),cd)
        rho1=sum(real(c(1:Ngrid))**2+imag(c(1:Ngrid))**2)*dx
        rho2=sum(real(c(1+Ngrid:2*Ngrid))**2+imag(c(1+Ngrid:2*Ngrid))**2)*dx
        write(21,*) time,rho1,rho2
    
        
        
    end do
    
        
    
    
    
    
    
    
    
    
    
    
    
 pause   
end program
