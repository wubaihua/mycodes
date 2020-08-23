program H2O_AIMD
    use def
    implicit real*8(a-h,o-z)
    type(atomtype) :: atom(3)
    real*8 force(3,3)!1 for atom, 2 for x/y/z
    real*8 ground_eng



    atom(:)%name=(/"O ","H ","H "/)
    atom(:)%index=(/8,1,1/)
    atom(:)%x=(/0.0_8, 0.0_8, 0.0_8/)
    atom(:)%y=(/0.0_8, 0.75703900_8, -0.75703900_8/)
    atom(:)%z=(/0.12020700_8, -0.48082800_8, -0.48082800_8/)
    atom(:)%charge=(/8,1,1/)



    call gen_gau_input_ground(3,atom)
    call system("g16 ground.gjf")
    call read_gau_output_ground(3,ground_eng,force)

    write(*,*) ground_eng
    do i=1,3
        write(*,*) force(i,:)
    end do




end program



subroutine gen_gau_input_ground(Natom,atom)
    use def
    implicit real*8(a-h,o-z)
    integer Natom
    type(atomtype), intent(in) :: atom(Natom)

    open(20,file="ground.gjf",status="REPLACE")
    write(20,*) "%mem=500MB"
    write(20,*) "%nproc=2"
    write(20,*) "#p b3lyp def2svp force"
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