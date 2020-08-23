program H2O_AIMD
    use def
    implicit real*8(a-h,o-z)
    type(atomtype) :: atom(3)
    real*8 force(3,3)!1 for atom, 2 for x/y/z
    real*8 ground_eng
    real*8 mass(3),p(3,3)

    mass=(/15.9994_8,1.00794_8,1.00794_8/)
    atom(:)%name=(/"O ","H ","H "/)
    atom(:)%index=(/8,1,1/)
    atom(:)%x=(/0.0_8, 0.0_8, 0.0_8/)
    atom(:)%y=(/0.0_8, 0.75703900_8, -0.75703900_8/)
    atom(:)%z=(/0.12020700_8, -0.48082800_8, -0.48082800_8/)
    atom(:)%x=atom(:)%x*10.0/5.291772108
    atom(:)%y=atom(:)%y*10.0/5.291772108
    atom(:)%z=atom(:)%z*10.0/5.291772108
    atom(:)%charge=(/8,1,1/)

    beta=10
    dt=0.1
    ttot=10
    nstep=int(ttot/dt)

    open(90,file="traj.xyz")
    open(91,file="energy.dat")

    do i=1,3
        call box_muller(p(i,1),x2,sqrt(mass(i)/beta),0.0_8)
        call box_muller(p(i,2),x2,sqrt(mass(i)/beta),0.0_8)
        call box_muller(p(i,3),x2,sqrt(mass(i)/beta),0.0_8)
    end do

    t_now=0
    do i=1,nstep
        call MD_ground_cal(3,atom,force,ground_eng)
        p(:,:)=p(:,:)+force(:,:)*dt/2

        atom(:)%x=atom(:)%x+p(:,1)*dt/mass(:)
        atom(:)%y=atom(:)%y+p(:,2)*dt/mass(:)
        atom(:)%z=atom(:)%z+p(:,3)*dt/mass(:)

        call MD_ground_cal(3,atom,force,ground_eng)
        p(:,:)=p(:,:)+force(:,:)*dt/2

        t_now=t_now+dt

        call output(90,91,t_now,3,atom,p,mass,ground_eng)

    end do





    ! call gen_gau_input_ground(3,atom)
    ! call system("g16 ground.gjf")
    ! call read_gau_output_ground(3,ground_eng,force)

    ! write(*,*) ground_eng
    ! do i=1,3
    !     write(*,*) force(i,:)
    ! end do




end program



subroutine gen_gau_input_ground(Natom,atom)
    use def
    implicit real*8(a-h,o-z)
    integer Natom
    type(atomtype), intent(in) :: atom(Natom)

    open(20,file="ground.gjf",status="REPLACE")
    write(20,*) "%mem=500MB"
    write(20,*) "%nproc=2"
    write(20,*) "#p b3lyp def2svp force units(au)"
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

subroutine read_gau_output_ground(Natom,ground_eng,force)
    implicit real*8(a-h,o-z)
    integer Natom
    real*8 force(Natom,3),ground_eng
    character*200 c200
    character*20 c20

    open(31,file="ground.log")
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

subroutine MD_ground_cal(Natom,atom,force,ground_eng)
    use def
    integer Natom
    type(atomtype), intent(in) :: atom(Natom)
    real*8 force(Natom,3),ground_eng

    call gen_gau_input_ground(Natom,atom)
    call system("g16 ground.gjf")
    call read_gau_output_ground(Natom,ground_eng,force)

end subroutine

subroutine output(idxyz,idE,t,Natom,atom,p,mass,ground_eng)
    use def
    implicit real*8(a-h,o-z)
    integer Natom,idxyz,idE
    type(atomtype), intent(in) :: atom(Natom)
    real*8 ground_eng,p(Natom,3),mass(3),t
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

    write(idE,"(4E18.8)") t,K,ground_eng,K+ground_eng




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