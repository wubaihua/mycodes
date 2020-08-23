function Vpot(itype,x)
    implicit real*8(a-h,o-z)
    real*8 x,vpot
    integer itype
    
    select case(itype)
    case(1)
        Vpot=0.5*x**2
    case(2)
        Vpot=0.25*x**4
    end select
    
end function


function dVpot(itype,x)
    implicit real*8(a-h,o-z)
    real*8 x,vpot
    integer itype
    
    select case(itype)
    case(1)
        dVpot=x
    case(2)
        dVpot=x**3
    end select
    
end function


