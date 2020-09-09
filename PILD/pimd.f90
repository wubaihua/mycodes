subroutine staging(x,u,Nbeads)
    implicit real*8(a-h,o-z)
    integer Nbeads
    real*8 u(Nbeads)
    real*8 x(Nbeads)
    
    !N=size(u)
    
    do i=1,Nbeads
        if(i==1)then
            x(i)=u(1)
        else
            x(i)=u(1)
            do j=i,Nbeads
                x(i)=x(i)+real(i-1)*u(j)/real(j-1)
            end do
        end if
    end do
    
    
    
end subroutine

subroutine box_muller(x1,x2,sigma,miu)
    implicit none
    real*8, parameter :: pi = 3.14159265358979323846
    !real*8, parameter :: hbar    = 1.0
    real*8 x1,x2,u1,u2,sigma,miu
    call RANDOM_NUMBER(u1)
    call RANDOM_NUMBER(u2)
    x1=sqrt(-2*log(u1))*cos(2*pi*u2)
    x2=sqrt(-2*log(u1))*sin(2*pi*u2)
    
    x1=x1*sigma+miu
    x2=x2*sigma+miu
end subroutine

subroutine pimd_mass(mass,m4u,m4p,Nbeads)
    implicit real*8(a-h,o-z)
    real*8 mass
    real*8,intent(inout) :: m4u(Nbeads),m4p(Nbeads)
    
    
    
    do i=1,Nbeads
        if(i==1)then
            m4u(i)=0
            m4p(i)=mass
        else
            m4u(i)=mass*real(i)/real(i-1)
            m4p(i)=m4u(i)
        end if
    end do
    
end subroutine








! dynamics part

subroutine pimd_force(u,f,Nbeads,itype)
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



subroutine pimd_A(dt,omega,m4p,u,p,Nbeads)
    implicit real*8(a-h,o-z)
    real*8 dt,omega
    real*8,intent(in) :: m4p(Nbeads)
    real*8,intent(inout) :: u(Nbeads),p(Nbeads)
    real*8 u0(Nbeads),p0(Nbeads)
    
    u0=u
    p0=p
    
    
    
    u(1)=u(1)+p(1)*dt/m4p(1)
    
    do i=2,Nbeads
        u(i)=u0(i)*cos(omega*dt)+p0(i)*sin(omega*dt)/(omega*m4p(i))
        p(i)=p0(i)*cos(omega*dt)-u0(i)*sin(omega*dt)*omega*m4p(i)
    end do
    
end subroutine


subroutine pimd_B(dt,omega,u,p,itype,Nbeads)
    implicit real*8(a-h,o-z)
    integer itype,Nbeads
    real*8 dt
    real*8,intent(in) :: u(Nbeads)
    real*8,intent(inout) :: p(Nbeads)
    real*8,allocatable :: f(:)
    
    
    
    allocate(f(Nbeads))
    
    
    
    call pimd_force(u,f,Nbeads,itype)
    
    p=p-f*dt
    deallocate(f)
    
end subroutine


subroutine pimd_O(dt,p,omega,m4p,beta,Nbeads)
    implicit real*8(a-h,o-z)
    integer itype,Nbeads
    real*8 dt,omega,beta
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
        
    gamma_lang=omega
    c1=exp(-gamma_lang*dt)
    c2=sqrt(1-c1**2)
    
    p(:)=p(:)*c1+c2*sqrt(m4p(:))*eta(:)/sqrt(beta)
        
    
    
  
end subroutine



subroutine pimd_estimators(est,u,p,nbeads,m4p,beta,itype)
    implicit real*8(a-h,o-z)
    integer Nbeads,itype
    real*8 est(6),u(Nbeads),p(Nbeads),m4p(Nbeads),beta
    real*8 E_tot,K_prim,K_vir
    real*8 x(Nbeads)
    real*8,external :: Vpot,dVpot
    
    call staging(x,u,Nbeads)
    
    xc=sum(x)/Nbeads
    
    E_tot=0
    
    V=0
    K_vir=1/(2*beta)
    !K_vir=0
    K_prim=Nbeads/(2*beta)
    
    omega_P2=real(Nbeads)/(beta**2)
    
    do i=1,Nbeads
        
        K_vir=K_vir+0.5*(x(i)-xc)*dVpot(itype,x(i))/real(Nbeads)
        !K_vir=K_vir+0.5*x(i)*dVpot(itype,x(i))/real(Nbeads)
        !
        V=V+Vpot(itype,x(i))/Nbeads
        !E_tot=
    end do
    do i=2,Nbeads
        K_prim=K_prim-0.5*Nbeads*(x(i)-x(i-1))**2/beta**2
    end do
    K_prim=K_prim-0.5*Nbeads*(x(1)-x(Nbeads))**2/beta**2
    
    est(1)=V
    est(2)=K_prim
    est(3)=K_vir
    est(4)=V+K_prim
    est(5)=V+K_vir  
    est(6)=1./(2*K_vir)
    !write(*,*) "test1=",est(3)


end subroutine