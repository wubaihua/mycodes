! =========================================================
! AIMD program for H2O molecule in ground/excited state
! Using Gaussian 16 for quantum chemistry calcualtion
!
! Complie with def.f90 from wQC
! =========================================================
! Written by Baihua Wu, On August 23, 2020
! Updated by Baihua Wu, On August 23, 2020
! Email: wubaihua@pku.edu.cn
!        wubaihua19@gmail.com
! =========================================================
! ***************************************************************************  
program H2O_AIMD
    use def
    implicit real*8(a-h,o-z)
    type(atomtype) :: atom(3)
    real*8 force(3,3)!1 for atom, 2 for x/y/z
    real*8 NACV(3,3)!1 for atom, 2 for x/y/z
    real*8 ground_eng,excited_eng
    real*8 mass(3),p(3,3)
    character*10 ctraj
    logical alive

    mass=(/15.9994_8,1.00794_8,1.00794_8/)
    atom(:)%name=(/"O ","H ","H "/)
    atom(:)%index=(/8,1,1/)

    atom(:)%x=(/0.0_8, 0.0_8, 0.0_8/)
    atom(:)%y=(/0.0_8, 0.75703900_8, -0.75703900_8/)
    atom(:)%z=(/0.12020700_8, -0.48082800_8, -0.48082800_8/)

    ! atom(:)%x=(/0.0_8, 0.0_8, 0.0_8/)
    ! atom(:)%y=(/0.0_8, 1.24175300_8, -1.24175300_8/)
    ! atom(:)%z=(/-0.00004600_8, 0.00018400_8, 0.00018400_8/)

    atom(:)%x=atom(:)%x*10.0/5.291772108
    atom(:)%y=atom(:)%y*10.0/5.291772108
    atom(:)%z=atom(:)%z*10.0/5.291772108
    atom(:)%charge=(/8,1,1/)

    beta=50
    dt=1
    ttot=20
    nstep=int(ttot/dt)
    Ntraj=1

    

    ! do i=1,3
    !     call box_muller(p(i,1),x2,sqrt(mass(i)/beta),0.0_8)
    !     call box_muller(p(i,2),x2,sqrt(mass(i)/beta),0.0_8)
    !     call box_muller(p(i,3),x2,sqrt(mass(i)/beta),0.0_8)
    ! end do

    do itraj=1,Ntraj
        do i=1,3
            call box_muller(p(i,1),x2,sqrt(mass(i)/beta),0.0_8)
            call box_muller(p(i,2),x2,sqrt(mass(i)/beta),0.0_8)
            call box_muller(p(i,3),x2,sqrt(mass(i)/beta),0.0_8)
        end do
        write(ctraj,"(I10)") itraj
        ! write(*,*) ctraj
        inquire(DIRECTORY="itraj_"//trim(adjustl(ctraj)),exist=alive)
        if(.not.alive)then
            call system("mkdir itraj_"//trim(adjustl(ctraj)))
        end if
        open(90,file="itraj_"//trim(adjustl(ctraj))//"/traj.xyz")
        open(91,file="itraj_"//trim(adjustl(ctraj))//"/energy.dat")
        t_now=0
        do i=1,nstep
            ! call MD_ground_cal(itraj,3,atom,force,ground_eng)
            call MD_excited_cal(3,atom,force,excited_eng,NACV)
            p(:,:)=p(:,:)+force(:,:)*dt/2

            atom(:)%x=atom(:)%x+p(:,1)*dt/mass(:)
            atom(:)%y=atom(:)%y+p(:,2)*dt/mass(:)
            atom(:)%z=atom(:)%z+p(:,3)*dt/mass(:)

            ! call MD_ground_cal(itraj,3,atom,force,ground_eng)
            call MD_excited_cal(3,atom,force,excited_eng,NACV)
            p(:,:)=p(:,:)+force(:,:)*dt/2

            t_now=t_now+dt

            ! call output(90,91,t_now,3,atom,p,mass,ground_eng)
            call output(90,91,t_now,3,atom,p,mass,excited_eng)
            write(*,*) NACV
        end do
        close(90)
        close(91)
    end do

end program


!    ----------------------------------------------------------------------------
!    -----------------------------------------------------------
!    Generate the g16 input file for calculate
!    potential and force of ground state
!    -----------------------------------------------------------
subroutine gen_gau_input_ground(itraj,Natom,atom)
    use def
    implicit real*8(a-h,o-z)
    integer Natom,itraj
    type(atomtype), intent(in) :: atom(Natom)
    character*10 ctraj
    logical alive

    write(ctraj,"(I10)") itraj

    ! inquire(DIRECTORY="itraj_"//trim(adjustl(ctraj)),exist=alive)
    ! if(.not.alive)then
    !     call system("mkdir itraj_"//trim(adjustl(ctraj)))
    ! end if

    open(20,file="itraj_"//trim(adjustl(ctraj))//"/ground.gjf",status="REPLACE")
    write(20,*) "%mem=500MB"
    write(20,*) "%nproc=2"
    write(20,*) "#p b3lyp def2svp force units(au) nosymm"
    write(20,*)
    write(20,*) "generate for AIMD to electronic calculate potential and force."
    write(20,*)
    
    write(20,"(A3)") "0 1"
    do i=1,Natom
        write(20,"(A2,3X,3E18.8)") atom(i)%name, atom(i)%x, atom(i)%y, atom(i)%z
    end do
    write(20,*)

    close(20)




end subroutine

!    ----------------------------------------------------------------------------
!    -----------------------------------------------------------
!    Read the g16 output file to get
!    potential and force of ground state
!    -----------------------------------------------------------
subroutine read_gau_output_ground(itraj,Natom,ground_eng,force)
    implicit real*8(a-h,o-z)
    integer Natom,itraj
    real*8 force(Natom,3),ground_eng
    character*200 c200
    character*20 c20
    character*10 ctraj

    write(ctraj,"(I10)") itraj

    inquire(file="itraj_"//trim(adjustl(ctraj))//"/ground.log",exist=alive)
    if(.not.alive)then
        write(idlog,"(a)") 'Gaussian log file for ground state not found! '
        return
    end if
        
    open(401,file="itraj_"//trim(adjustl(ctraj))//"/ground.log",status='old',position='append')
    backspace(401)
    read(401,"(a)") c200
    if(c200(2:7)=='Normal')then
        write(idlog,"(a)") 'Gaussian ground state end Normally'
    else
        write(idlog,"(a)") 'Gaussian ground state end Error! '
        return
    end if
        
    close(401)

    

    open(31,file="itraj_"//trim(adjustl(ctraj))//"/ground.log")
    do while(.true.)
        read(31,"(a)") c200
        ! write(*,"(a)") c200
        if(index(c200,"SCF Done:")/=0)then
            ! write(*,"(a)") c200
            backspace(31)
            read(31,*) c20,c20,c20,c20,ground_eng
            ! write(*,*)  ground_eng
            exit
        end if
    end do

    do while(.true.)
        read(31,"(a)") c200
        if(index(c200,"Forces (Hartrees/Bohr)")/=0)then
            read(31,"(a)") c200
            read(31,"(a)") c200
            !write(*,*) c200
            do i=1,Natom
                read(31,*) m,n,force(i,:)
                
            end do
            exit
        end if
    end do
    close(31)
end subroutine

!    ----------------------------------------------------------------------------
!    -----------------------------------------------------------
!    ground force calculation in MD process 
!    -----------------------------------------------------------
subroutine MD_ground_cal(itraj,Natom,atom,force,ground_eng)
    use def
    integer Natom,itraj
    type(atomtype), intent(in) :: atom(Natom)
    real*8 force(Natom,3),ground_eng
    character*10 ctraj

    write(ctraj,"(I10)") itraj

    call gen_gau_input_ground(itraj,Natom,atom)
    call system("g16 "//"itraj_"//trim(adjustl(ctraj))//"/ground.gjf")
    call read_gau_output_ground(itraj,Natom,ground_eng,force)

end subroutine

!    ----------------------------------------------------------------------------
!    -----------------------------------------------------------
!    Generate the g16 input file for calculate
!    potential and force of excited state
!    -----------------------------------------------------------
subroutine gen_gau_input_excited(Natom,atom)
    use def
    implicit real*8(a-h,o-z)
    integer Natom
    type(atomtype), intent(in) :: atom(Natom)

    open(21,file="excited.gjf",status="REPLACE")
    write(21,*) "%mem=500MB"
    write(21,*) "%nproc=2"
    write(21,*) "#p b3lyp def2svp force units(au) TD(NAC) nosymm"
    write(21,*)
    write(21,*) "generate for AIMD to electronic calculate potential and force."
    write(21,*)
    
    write(21,"(A3)") "0 1"
    do i=1,Natom
        write(21,"(A2,3X,3E18.8)") atom(i)%name, atom(i)%x, atom(i)%y, atom(i)%z
    end do
    write(21,*)

    close(21)




end subroutine

!    ----------------------------------------------------------------------------
!    -----------------------------------------------------------
!    Read the g16 output file to get
!    potential and force of excited state
!    -----------------------------------------------------------
subroutine read_gau_output_excited(Natom,excited_eng,force,NACV)
    implicit real*8(a-h,o-z)
    integer Natom
    real*8 force(Natom,3),excited_eng,NACV(Natom,3)
    character*200 c200
    character*31 c31

    inquire(file="excited.log",exist=alive)
    if(.not.alive)then
        write(*,"(a)") 'Gaussian log file for excited state not found! '
        return
    end if
        
    open(401,file='excited.log',status='old',position='append')
    backspace(401)
    read(401,"(a)") c200
    if(c200(2:7)=='Normal')then
        write(*,"(a)") 'Gaussian excited state end Normally'
    else
        write(*,"(a)") 'Gaussian excited state end Error! '
        return
    end if
    close(401)

    open(41,file="excited.log")
    do while(.true.)
        read(41,"(a)") c200
        ! write(*,"(a)") c200
        if(index(c200,"Total Energy, E(TD-HF/TD-DFT)")/=0)then
            ! write(*,"(a)") c200
            ! backspace(41)
            ! read(41,"(a,E18.8)") c31,excited_eng
            ! write(*,*)  c31
            id=index(c200,"=")
            c200(1:id)=" "
            read(c200,*) excited_eng
            exit
        end if
    end do

    do while(.true.)
        read(41,"(a)") c200
        if(index(c200,"Forces (Hartrees/Bohr)")/=0)then
            read(41,"(a)") c200
            read(41,"(a)") c200
            !write(*,*) c200
            do i=1,Natom
                read(41,*) m,n,force(i,:)
                
            end do
            exit
        end if
    end do

    do while(.true.)
        read(41,"(a)") c200
        if(index(c200,"Nonadiabatic Coup.")/=0)then
            read(41,"(a)") c200
            read(41,"(a)") c200
            !write(*,*) c200
            do i=1,Natom
                read(41,*) m,n,NACV(i,:)
                
            end do
            exit
        end if
    end do

    
    close(41)
end subroutine

!    ----------------------------------------------------------------------------
!    -----------------------------------------------------------
!    excited force calculation in MD process 
!    -----------------------------------------------------------
subroutine MD_excited_cal(Natom,atom,force,excited_eng,NACV)
    use def
    integer Natom
    type(atomtype), intent(in) :: atom(Natom)
    real*8 force(Natom,3),excited_eng,NACV(Natom,3)

    call gen_gau_input_excited(Natom,atom)
    call system("g16 excited.gjf")
    call read_gau_output_excited(Natom,excited_eng,force,NACV)

end subroutine

!    ----------------------------------------------------------------------------
!    -----------------------------------------------------------
!    write the trajectory file and energy file
!    -----------------------------------------------------------
subroutine output(idxyz,idE,t,Natom,atom,p,mass,potential)
    use def
    implicit real*8(a-h,o-z)
    integer Natom,idxyz,idE
    type(atomtype), intent(in) :: atom(Natom)
    real*8 potential,p(Natom,3),mass(3),t
    real*8 K

    write(idxyz,*) Natom
    write(idxyz,*) "h2o AIMD"
    do i=1,Natom
        write(idxyz,"(A2,3X,3E18.8)") atom(i)%name, atom(i)%x, atom(i)%y, atom(i)%z
    end do

    K=0
    do i=1,3
        K=K+sum(p(:,i)**2/(2*mass(:)))
    end do

    write(idE,"(4E18.8)") t,K,potential,K+potential




end subroutine


!    ----------------------------------------------------------------------------
!    -----------------------------------------------------------
!    Box-MUller method for Gaussian distribution
!    -----------------------------------------------------------
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