program pild
    implicit real*8(a-h,o-z)
    real*8,allocatable :: p(:),x(:),u(:) ,m4p(:),m4u(:),m4pild(:)
    real*8 mass,est(6)
    real*8,allocatable :: pot_est(:),k_prim_est(:),k_vir_est(:),E_prim_est(:),E_vir_est(:)
    integer number(5)
    real*8,allocatable :: corre_fun(:,:,:)
    
    
    
    
    open(23,file="result_traj.dat")
    
    
    
    Nbeads=30
    beta=8
    itype=2
    
    Ntraj=100
    x0=0
    mass=1
    
    
    
    
    omega=sqrt(real(Nbeads))/beta
    
    gamma_ad=1E-5
    omega_ad=omega/sqrt(gamma_ad)
    
    igrid=100
    
    ttot_pimd=1000
    dt_pimd=0.01
    
   
    
    ttot_pild=10
    dt_pild=0.0001
    
    
    
    allocate(p(Nbeads))
    allocate(x(Nbeads))
    allocate(u(Nbeads))
    allocate(m4p(Nbeads))
    allocate(m4u(Nbeads))
    allocate(m4pild(Nbeads))
    allocate(corre_fun(2,int(ttot_pild/dt_pild),Ntraj))
    do itraj=1,Ntraj
        u=x0
    
        do i=1,nBeads
            call box_muller(p(i),x22,sqrt(mass/beta),0.0_8)
        end do
        call pimd_mass(mass,m4u,m4p,Nbeads)
        m4pild=m4u*gamma_ad
        t=0
        nstep=int(ttot_pimd/dt_pimd)
    
        do itime=1,nstep
        
            call pimd_B(dt_pimd/2,omega,u,p,itype,Nbeads)
            call pimd_A(dt_pimd/2,omega,m4p,u,p,Nbeads)
            call pimd_O(dt_pimd,p,omega,m4p,beta,Nbeads)
            call pimd_A(dt_pimd/2,omega,m4p,u,p,Nbeads)
            call pimd_B(dt_pimd/2,omega,u,p,itype,Nbeads)
        
            
            
            t=t+dt
        end do
        
        x0=u(1)
        call cal_thermmass(itype,x0,beta,mass,thermmass)
        call box_muller(p(1),x2,sqrt(thermmass/beta),0.0_8)
        do j=2,Nbeads
            call box_muller(p(j),x2,sqrt(m4pild(j)/beta),0.0_8)
        end do
        p0=p(1)
        
        t=0
        nstep=int(ttot_pild/dt_pild)
        
        
        
        do itime=1,nstep
        
            call pild_B(dt_pild/2,beta,mass,u,p,itype,Nbeads)
            call pild_A(dt_pild/2,mass,omega_ad,m4pild,u,p,Nbeads)
            call pild_O(dt_pild,p,omega,m4pild,beta,Nbeads)
            call pild_A(dt_pild/2,mass,omega_ad,m4pild,u,p,Nbeads)
            call pild_B(dt_pild/2,beta,mass,u,p,itype,Nbeads)
            t=t+dt
            !write(*,*) u(:)
            !write(*,*) p(:)
            call staging(x,u,Nbeads)
            corre_fun(1,itime,itraj)=(x0**2+0.25*beta/thermmass-beta**2*p0**2/(4*thermmass**2))*x(1)**2
            corre_fun(2,itime,itraj)=mass*p0*p(1)/thermmass
            
            !write(*,*) t
        end do
        
        
        
    end do
    
    t=0
    do i=1,nstep-1
        t=t+dt_pild
        write(23,*) t,sum(corre_fun(1,i,:))/Ntraj,sum(corre_fun(2,i,:))/Ntraj
    end do
    
    
   
    

    
    
    
    
    
    
pause   
end program



