program pild
    implicit real*8(a-h,o-z)
    real*8,allocatable :: p(:),x(:),u(:) ,m4p(:),m4u(:),m4pild(:)
    real*8 mass,est(6)
    real*8,allocatable ::thm_traj(:),u0(:,:)
    integer number(5)
    real*8,allocatable :: corre_fun(:,:,:)
    
    
    
    
    open(23,file="result_traj.dat",status='replace')
    
    
    
    Nbeads=50
    beta=8
    itype=3
    
    Ntraj=100
    x0=0
    mass=1
    
    
    
    
    omega=sqrt(real(Nbeads))/beta
    
    gamma_ad=1E-6
    omega_ad=omega/sqrt(gamma_ad)
    
    igrid=100
    
    ttot_pimd=500
    dt_pimd=0.001
    
   
    
    ttot_pild=14
    dt_pild=0.00001
    
    allocate(thm_traj(Ntraj))
    allocate(u0(Ntarj,Nbeads))
    allocate(p(Nbeads))
    allocate(x(Nbeads))
    allocate(u(Nbeads))
    allocate(m4p(Nbeads))
    allocate(m4u(Nbeads))
    allocate(m4pild(Nbeads))
    allocate(corre_fun(2,int(ttot_pild/dt_pild),Ntraj))
    call pimd_mass(mass,m4u,m4p,Nbeads)
    m4pild=m4u*gamma_ad
    do itraj=1,Ntraj
        u=0
        
        do i=1,nBeads
            call box_muller(p(i),x22,sqrt(mass/beta),0.0_8)
        end do
       
        t=0
        nstep=int(ttot_pimd/dt_pimd)
        ntm=0
        thermmass=0
        do itime=1,nstep
        
            call pimd_B(dt_pimd/2,omega,u,p,itype,Nbeads)
            call pimd_A(dt_pimd/2,omega,m4p,u,p,Nbeads)
            call pimd_O(dt_pimd,p,omega,m4p,beta,Nbeads)
            call pimd_A(dt_pimd/2,omega,m4p,u,p,Nbeads)
            call pimd_B(dt_pimd/2,omega,u,p,itype,Nbeads)
        
            
            
            t=t+dt_pimd

            if(t>300)then
                ! thermmass0=0
                ! do ib=1,Nbeads
                !     call cal_thermmass(itype,u(ib),beta,mass,tm)
                !     thermmass0=thermmass0+tm
                ! end do
                call cal_thermmass(itype,u(1),beta,mass,tm)
                ! write(*,*) tm
                thermmass=thermmass+tm
                ntm=ntm+1
            end if




        end do

       thm_traj(itraj)=thermmass/ntm
    !    write(*,*) thm_traj(it)
       u0(itraj,:)=u(:)

    end do
    thermmass=sum(thm_traj(:))/Ntraj
    write(*,*) thermmass
    ! stop
    do itraj=1,Ntraj    
        x0=u0(itraj,1)
        u=u0(itraj,:)
        ! call cal_thermmass(itype,x0,beta,mass,thermmass)
        call box_muller(p(1),x2,sqrt(thermmass/beta),0.0_8)
        do j=2,Nbeads
            call box_muller(p(j),x2,sqrt(m4pild(j)/beta),0.0_8)
            ! call box_muller(u(j),x2,sqrt(1.0_8/(beta*m4pild(j)))/omega_ad,0.0_8)
        end do
        ! p(2:Nbeads)=p(2:Nbeads)*sqrt(gamma_ad)
        p0=p(1)
        
        t=0
        nstep=int(ttot_pild/dt_pild)
        
        
        
        do itime=1,nstep
            ! call cal_thermmass(itype,u(1),beta,mass,thermmass)
            call pild_B(dt_pild/2,beta,mass,u,p,itype,Nbeads,thermmass)
            call pild_A(dt_pild/2,mass,omega_ad,m4pild,u,p,Nbeads)
            call pild_O(dt_pild,p,omega_ad,m4pild,beta,Nbeads)
            call pild_A(dt_pild/2,mass,omega_ad,m4pild,u,p,Nbeads)
            ! call cal_thermmass(itype,u(1),beta,mass,thermmass)
            call pild_B(dt_pild/2,beta,mass,u,p,itype,Nbeads,thermmass)
            t=t+dt_pild
            !write(*,*) u(:)
            !write(*,*) p(:)
            ! call cal_thermmass(itype,u(1),beta,mass,thermmass)
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
    
    
   
    

    
    
    
    
    
    
! pause   
end program



