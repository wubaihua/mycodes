program pimd
    implicit real*8(a-h,o-z)
    real*8,allocatable :: p(:),x(:),u(:) ,m4p(:),m4u(:)
    real*8 mass,est(6)
    real*8,allocatable :: pot_est(:),k_prim_est(:),k_vir_est(:),E_prim_est(:),E_vir_est(:)
    integer number(5)
    
    number=(/2,10,50,100,200/)
    
    inum=1
    n_cont=5
    
    open(11,file="potential.dat")
    open(12,file="kinetic_prim.dat")
    open(13,file="kinetic_vir.dat")
    open(14,file="E_prim.dat")
    open(15,file="E_vir.dat")
    
    !open(23,file="result_traj.dat")
    
    
    
60  Nbeads=number(inum)
    beta=8
    itype=2
    
    Ntraj=20
    x0=0
    mass=1
    
    
    omega=sqrt(real(Nbeads))/beta
    
    igrid=100
    
    ttot=1000
    dt=0.01
    
    t_count=200
    
    
    allocate( pot_est(Ntraj))
    allocate( K_prim_est(Ntraj))
    allocate( K_vir_est(Ntraj))
    allocate( E_prim_est(Ntraj))
    allocate( E_vir_est(Ntraj))
    pot_est(:)=0
    K_prim_est=0
    K_vir_est=0
    E_prim_est=0
    E_vir_est=0
    
    allocate(p(Nbeads))
    allocate(x(Nbeads))
    allocate(u(Nbeads))
    allocate(m4p(Nbeads))
    allocate(m4u(Nbeads))
    do itraj=1,Ntraj
        u=x0
    
        do i=1,nBeads
            call box_muller(p(i),x22,sqrt(mass/beta),0.0_8)
        end do
        call pimd_mass(mass,m4u,m4p,Nbeads)
        t=0
        nstep=int(ttot/dt)
    
        do itime=1,nstep
        
            call pimd_B(dt/2,omega,u,p,itype,Nbeads)
            call pimd_A(dt/2,omega,m4p,u,p,Nbeads)
            call pimd_O(dt,p,omega,m4p,beta,Nbeads)
            call pimd_A(dt/2,omega,m4p,u,p,Nbeads)
            call pimd_B(dt/2,omega,u,p,itype,Nbeads)
        
            !if(mod(itime,igrid)==0)then
            if(t>t_count)then
                call pimd_estimators(est,u,p,nbeads,m4p,beta,itype)
                !write(23,"(7E18.8E3)") t,est(:)
                !write(*,*) "test2=",est(3)
            
                pot_est(itraj)=pot_est(itraj)+est(1)
                K_prim_est(itraj)=K_prim_est(itraj)+est(2)
                K_vir_est(itraj)=K_vir_est(itraj)+est(3)
                !write(*,*) K_vir_est(itraj)
                E_prim_est(itraj)=E_prim_est(itraj)+est(4)
                E_vir_est(itraj)=E_vir_est(itraj)+est(5)
            end if
            
            t=t+dt
        end do
    
    end do
    
110 format(I4,<Ntraj>E18.8E3)    
    
    write(11,110) Nbeads,pot_est(:)/int((ttot-t_count)/dt)
    write(12,110) Nbeads,K_prim_est(:)/int((ttot-t_count)/dt)
    write(13,110) Nbeads,K_vir_est(:)/int((ttot-t_count)/dt)
    write(14,110) Nbeads,E_prim_est(:)/int((ttot-t_count)/dt)
    write(15,110) Nbeads,E_vir_est(:)/int((ttot-t_count)/dt)
    
    deallocate(p)
    deallocate(x)
    deallocate(u)
    deallocate(m4p)
    deallocate(m4u)
    deallocate( pot_est)
    deallocate( K_prim_est)
    deallocate( K_vir_est)
    deallocate( E_prim_est)
    deallocate( E_vir_est)
    
    if(inum<n_cont)then
        inum=inum+1
        goto 60
    end if
    
    
    
    
    
    
    
pause   
end program
