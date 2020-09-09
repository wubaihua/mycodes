subroutine cal_thermmass(itype,x,beta,mass,thermmass)
    implicit none
    real*8,external :: hessian
    integer itype
    real*8 beta,Q,x,thermmass,mass,u,g
    
    
    g=hessian(itype,x)
    if(g<0)then
        u=beta*sqrt(-hessian(itype,x)/mass)
        Q=tanh(0.5*u)/(0.5*u)
    else
        u=beta*sqrt(hessian(itype,x)/mass)
        Q=(0.5*u)/tanh(0.5*u)
    end if
    
   thermmass=mass*Q
    
end subroutine






subroutine pild_force(u,f,Nbeads,itype)
    implicit real*8(a-h,o-z)
    integer Nbeads,itype
    real*8 u(Nbeads),f(Nbeads),x(Nbeads)
    real*8,external :: Vpot,dVpot
    
    call staging(x,u,Nbeads)
    
    f(1)=0
    do i=1,Nbeads
        f(1)=f(1)+dVpot(itype,x(i))/Nbeads
    end do
    
    do i=2,Nbeads
        f(i)=dVpot(itype,x(i))/Nbeads+f(i-1)*real(i-2)/real(i-1)
    end do
    

end subroutine



subroutine pild_A(dt,mass,omega_ad,m4p,u,p,Nbeads)
    implicit real*8(a-h,o-z)
    real*8 dt,omega_ad,mass
    real*8,intent(in) :: m4p(Nbeads)
    real*8,intent(inout) :: u(Nbeads),p(Nbeads)
    real*8 u0(Nbeads),p0(Nbeads)
    
    u0=u
    p0=p
    
    
    
    u(1)=u(1)+p(1)*dt/mass
    
    do i=2,Nbeads
        u(i)=u0(i)*cos(omega_ad*dt)+p0(i)*sin(omega_ad*dt)/(omega_ad*m4p(i))
        p(i)=p0(i)*cos(omega_ad*dt)-u0(i)*sin(omega_ad*dt)*omega_ad*m4p(i)
    end do
    !write(*,*) u(:)
end subroutine


subroutine pild_B(dt,beta,mass,u,p,itype,Nbeads)
    implicit real*8(a-h,o-z)
    integer itype,Nbeads
    real*8 dt,thermmass,mass,beta,x(Nbeads)
    real*8,intent(in) :: u(Nbeads)
    real*8,intent(inout) :: p(Nbeads)
    real*8,allocatable :: f(:)
    
    
    call staging(x,u,Nbeads)
    allocate(f(Nbeads))
    
    call cal_thermmass(itype,x(1),beta,mass,thermmass)
    
    call pild_force(u,f,Nbeads,itype)
    
    p(1)=p(1)-thermmass*f(1)*dt/mass
    p(2:Nbeads)=p(2:Nbeads)-f(2:Nbeads)*dt
    deallocate(f)
    !write(*,*) p(:)
end subroutine


subroutine pild_O(dt,p,omega_ad,m4p,beta,Nbeads)
    implicit real*8(a-h,o-z)
    integer itype,Nbeads
    real*8 dt,omega_ad,beta
    real*8 gamma_lang,c1,c2,eta(Nbeads)
    real*8,intent(inout) :: p(Nbeads)
    real*8,intent(in) :: m4p(Nbeads)
    
    
    
    !do i=1,Nbeads
    !    gamma_lang=omega
    !    !gamma_lang=1
    !    c1=exp(-gamma_lang*dt)
    !    c2=sqrt(1-c1**2)
    !    call box_muller(eta,x2, 1.0_8, 0.0_8)
    !    p(i)=p(i)*c1+c2*sqrt(m4p(i))*eta/sqrt(beta)
    !    
    !end do
    
    do i=1,Nbeads
        call box_muller(eta(i),x2, 1.0_8, 0.0_8)
    end do
        
    gamma_lang=2*omega_ad
    c1=exp(-gamma_lang*dt)
    c2=sqrt(1-c1**2)
    
    p(2:Nbeads)=p(2:Nbeads)*c1+c2*sqrt(m4p(2:Nbeads))*eta(2:Nbeads)/sqrt(beta)
        
    
    
  
end subroutine


subroutine pild_corr_fun(Nbeads,beta,mass,m4pild,itype,omega_ad,u,p,corr_fun,x0,p0)
    implicit real*8(a-h,o-z)
    integer Nbeads,itype
    real*8 beta,mass,m4pild(Nbeads),omega_ad,u(Nbeads),x(Nbeads),p(Nbeads),corr_fun(2),x0,p0
    real*8,external :: Vpot
    real*8,parameter :: pi=3.1415926535897932
     
    call staging(x,u,Nbeads)
    
    
    corr_fun(1)=x0*x(1)
    corr_fun(2)=p0*p(1)





end subroutine
    


