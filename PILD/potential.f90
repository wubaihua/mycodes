function Vpot(itype,x)
    implicit real*8(a-h,o-z)
    real*8 x,vpot
    integer itype
    
    select case(itype)
    case(1)
        Vpot=0.5*x**2
    case(2)
        Vpot=0.25*x**4
    case(3)
        Vpot=x**2-0.1*x**3+0.1*x**4
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
    case(3)
        dVpot=2*x-0.3*x**2+0.4*x**3
    end select
    
end function

function hessian(itype,x)
    implicit real*8(a-h,o-z)
    real*8 x,hessian
    integer itype
    
    select case(itype)
    case(1)
        hessian=1
    case(2)
        hessian=3*x**2
    case(3)
        hessian=2.0-0.6*x+1.2*x**2
    end select
    
end function


