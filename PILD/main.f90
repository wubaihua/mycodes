program pild
    use mpi ! 声明使用mpi模块
    implicit real*8(a-h,o-z)
    real*8,allocatable :: p(:),x(:),u(:) ,m4p(:),m4u(:),m4pild(:)
    real*8 mass,est(6)
    real*8,allocatable ::thm_traj(:),u0(:,:)
    integer number(5)
    real*8,allocatable :: corre_fun(:,:,:)
    integer :: mpi_size, mpi_rank, mpi_ierr ! 总进程数、当前进程ID、子程序返回状态
    integer,allocatable :: status(:)


    ! 初始化
    call mpi_init(mpi_ierr)

    ! 获取总进程数和当前进程ID
    call mpi_comm_size(mpi_comm_world, mpi_size, mpi_ierr)
    call mpi_comm_rank(mpi_comm_world, mpi_rank, mpi_ierr)

    
    allocate(status(0:2*mpi_size-1))
    
    
    
    Nbeads=50
    beta=8
    itype=2
    
    Ntraj=40
    x0=0
    mass=1
    
    
    
    
    omega=sqrt(real(Nbeads))/beta
    
    gamma_ad=1E-5
    omega_ad=omega/sqrt(gamma_ad)
    
    igrid=100
    
    ttot_pimd=500
    dt_pimd=0.01
    
   
    
    ttot_pild=14
    dt_pild=0.0001
    
    allocate(thm_traj(Ntraj))
    allocate(u0(Ntraj,Nbeads))
    allocate(p(Nbeads))
    allocate(x(Nbeads))
    allocate(u(Nbeads))
    allocate(m4p(Nbeads))
    allocate(m4u(Nbeads))
    allocate(m4pild(Nbeads))
    allocate(corre_fun(2,int(ttot_pild/dt_pild),Ntraj))
    call pimd_mass(mass,m4u,m4p,Nbeads)
    m4pild=m4u*gamma_ad


    if(mpi_rank==0)then

        do itraj=1,Ntraj
            ! if (mod(itraj, mpi_size) == mpi_rank) then
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

                ! if(mpi_rank==0)then
                !     do irec=1,mpi_size-1
                !         call mpi_recv(itraj_rec,1,mpi_integer,irec,irec+mpi_size,mpi_comm_world,status,mpi_ierr)
                !         call mpi_recv(thm_traj(itraj_rec),1,mpi_double,irec,irec,mpi_comm_world,status,mpi_ierr)
                !         !call mpi_recv(u0(itraj_rec,:),Nbeads,mpi_double,irec,irec+2*mpi_size,mpi_comm_world,status,mpi_ierr)
                !     end do
                ! else
                !     call mpi_send(itraj,1,mpi_integer,0,mpi_rank+mpi_size,mpi_comm_world,mpi_ierr)
                !     call mpi_send(thm_traj(itraj),1,mpi_double,0,mpi_rank,mpi_comm_world,mpi_ierr)
                !     !call mpi_send(u0(itraj,:),Nbeads,mpi_double,0,mpi_rank+2*mpi_size,mpi_comm_world,mpi_ierr)
                ! end if



            ! end if
        end do
        thermmass=sum(thm_traj(:))/Ntraj
        write(*,*) thermmass
        ! stop
    end if

    call mpi_barrier(mpi_comm_world,mpi_ierr)


    if(mpi_rank==0)then
        do irec=1,mpi_size-1
            ! call mpi_send(itraj+irec,1,mpi_integer,irec,irec+mpi_size,mpi_comm_world,mpi_ierr)
            call mpi_send(thermmass,1,mpi_double,irec,irec,mpi_comm_world,mpi_ierr)
            call mpi_send(u0(:,:),Ntraj*Nbeads,mpi_double,irec,irec+2*mpi_size,mpi_comm_world,mpi_ierr)
        end do
    else              
        ! call mpi_recv(itraj_rec,1,mpi_integer,irec,irec+mpi_size,mpi_comm_world,status,mpi_ierr)
        call mpi_recv(thermmass,1,mpi_double,0,mpi_rank,mpi_comm_world,status,mpi_ierr)
        call mpi_recv(u0(:,:),Ntraj*Nbeads,mpi_double,0,mpi_rank+2*mpi_size,mpi_comm_world,status,mpi_ierr)
    end if

    call mpi_barrier(mpi_comm_world,mpi_ierr)
    ! write(*,*) "thermmass=",thermmass
    ! write(*,*) "u0=",u0
    ! stop
    call init_seed(mpi_rank)
    do itraj=1,Ntraj
        if (mod(itraj, mpi_size) == mpi_rank) then
            ! if(mpi_rank==0)then
            !     do irec=1,mpi_size-1
            !         ! call mpi_send(itraj+irec,1,mpi_integer,irec,irec+mpi_size,mpi_comm_world,mpi_ierr)
            !         call mpi_send(thermmass,1,mpi_double,irec,irec,mpi_comm_world,mpi_ierr)
            !         call mpi_send(u0(itraj+irec,:),Nbeads,mpi_double,irec,irec+2*mpi_size,mpi_comm_world,mpi_ierr)
            !     end do
            ! else              
            !     ! call mpi_recv(itraj_rec,1,mpi_integer,irec,irec+mpi_size,mpi_comm_world,status,mpi_ierr)
            !     call mpi_recv(thermmass,1,mpi_double,0,mpi_rank,mpi_comm_world,status,mpi_ierr)
            !     call mpi_recv(u0(itraj,:),Nbeads,mpi_double,0,mpi_rank+2*mpi_size,mpi_comm_world,status,mpi_ierr)
            ! end if
            ! call mpi_barrier(mpi_comm_world,mpi_ierr)
            ! write(*,*) "thermass=",thermass
            ! stop
        
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
            if(mpi_rank==0)then
                do irec=1,mpi_size-1
                    call mpi_recv(itraj_rec,1,mpi_integer,irec,irec+3*mpi_size,mpi_comm_world,status,mpi_ierr)
                    call mpi_recv(corre_fun(:,:,itraj_rec),2*nstep,mpi_double,irec,irec+4*mpi_size,mpi_comm_world,status,mpi_ierr)
                end do
            else
                call mpi_send(itraj,1,mpi_integer,0,mpi_rank+3*mpi_size,mpi_comm_world,mpi_ierr)
                call mpi_send(corre_fun(:,:,itraj),2*nstep,mpi_double,0,mpi_rank+4*mpi_size,mpi_comm_world,mpi_ierr)
            end if
        
        end if
    end do
    
    call mpi_barrier(mpi_comm_world,mpi_ierr)
    
    if(mpi_rank==0)then
        ! write(*,*) corre_fun(1,1,1:10)
        open(23,file="result_traj.dat",status='replace')
        t=0
        do i=1,nstep-1
            t=t+dt_pild
            write(23,*) t,sum(corre_fun(1,i,:))/Ntraj,sum(corre_fun(2,i,:))/Ntraj
        end do
    end if    
    
   
    

    
    
    
 call mpi_finalize(mpi_ierr)   
    
    
! pause   
end program

subroutine init_seed(my_prl)
	integer :: n, ival(8), v(3), i
	integer, allocatable :: seed(:)
	integer, optional :: my_prl
	
	call date_and_time(values=ival)
	v(1) = 101 * ival(8) + 256*ival(7) + mod(ival(8), 103)
    v(2) = ival(6) + 64*ival(5) + ( mod(1993, 2+ival(7) ) + 3 ) * 97 * ival(8) 
    v(3) = ival(3) + 4*ival(2) + 16*ival(1) + ( mod(997, 1+ival(7)) + 5 ) * 101* ival(8)
    
 	call random_seed(size=n)
    allocate(seed(n))
    
    !-- Give the seed an implementation-dependent kick
	call random_seed()
	call random_seed(get=seed)
	
	!-- first Bias
	do i=1, n
    	seed(i) = seed(i) + v(mod(i-1, 3) + 1) + ival(8)
  	enddo
  	  	
  	!-- second Bias
  	if ( n >= size(ival) ) then
        seed( 1 : size(ival) ) = seed( 1 : size(ival) ) + ival(:)
    else
        seed(:) = seed(:) + ival( 1 : n )
    endif
    
    !-- if parallel
    if ( present(my_prl) ) seed(:) = seed(:) + 113 * my_prl
    
  	call random_seed(put=seed)
  	deallocate(seed)
  	!has_initialed = .true.
end subroutine init_seed



