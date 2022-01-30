!====================================   Bisection_method  ==============================================


!=================================   Definition Of Modules  ==========================================
! Form of Function   :  F(x)=  a1 * x^(n-1)    +  a2 * x^(n-2) +...+ an

! Coefficient_Matrix (1,n)) : Vector         : :/input/: [ a1  a2  ...  an ]

! n                         : Integer Number : :/input/:  Equation's Order ? Polynomial
! a                         : Real    Number : :/input/:  End Point Value  Of The Interval [a,b] Which The Root Of The Function Is Within This Interval.
! b                         : Real    Number : :/input/:  End Ponit Value  Of THe Interval [a,b]
! i                         : counter
! step                      : Integer Number : :/input/:  The Number of Loop's Iterations
! r                         : Real    Number : :/input/:  Accuracy 
! result_x                  : Real    Number : :/output/: Root Of The Function
! Fvalue_1,Fvalue_2         : Variables used in Subroutines 
! j                         : counter


!=========================================   Modules   =================================================
Module Para



end module Para

!==============================================================

module arr

Real*8,Allocatable,dimension(:) :: Coefficient_Matrix


end module arr

!==============================================================

module var

Real*8    :: step,a,b,result_x,r,Fvalue_1,Fvalue_2 ,j
integer :: n,i


end module var
!=========================================================================================================

!============================================ Main Program ===============================================



Program Bisection_method
!Using bisection method to find F(xi)=0 :  a  < xi < b

use arr 
use var

Implicit None

open(101,file='input.txt')
open(102,file='result.txt)
read(101,*),a   
read(101,*),b
read(101,*),step      !The Step to Move From a to b
read(101,*),r         
read(101,*),n         !F(x)=  a1 * x^(n-1)    +  a2 * x^(n-2) +...+ an  :  n=?
allocate(Coefficient_Matrix(n))       ! a1,a2,a3,...,an
read(101,*),Coefficient_Matrix(1:n)
close(101)

!==================== Writing the Results ============================
do j = a,b,step
    call f(Fvalue_1,j)
    call f(Fvalue_2,j+step)
    if ( Fvalue_1*Fvalue_2 <= 0) then
    call function_calculate(result_x,j,j+step)
    write(102,*),' answer = " x" '
    write(102,*),result_x
    end if
end do


close(101)
close(102)


deallocate(Coefficient_Matrix)
pause
end program Bisection_method

!==================================================  End Of Program  ===================================================


!=============================================== Subroutines And Functions =============================================

!=============Subroutine f ============ calculating F(x) value =====================
!=== Subroutine Input : xf
!=== Subroutine Output: Fvalue 

Subroutine f(Fvalue,xf)!========================= calculating F(x) value =============================================

use arr
use var

Implicit None 

Real*8,INTENT(IN) :: xf
Real*8,INTENT(out) :: Fvalue

Fvalue = 0.0
do i=1,n
    Fvalue = Fvalue + Coefficient_Matrix(i)*xf**(n-i)
End do 


END Subroutine f

!==============================================================


!=============Subroutine function_calculate ============ 

!=============== using bisection_method to find f(x)=0 :   j = a_functioncalculate   < x <  j+step = b_functioncalculate ===========================
!=== Subroutine Input : a_functioncalculate,b_functioncalculate
!=== Subroutine Output: _functioncalculate_result 

Subroutine function_calculate(x_functioncalculate_result,a_functioncalculate,b_functioncalculate)
 

use arr
use var

implicit none 
Real*8,INTENT(inout) :: x_functioncalculate_result  
real*8,intent(inout) :: a_functioncalculate,b_functioncalculate
real*8 :: ans,counter_fc,x_value,a_value,b_value

ans=0
counter_fc=1
do while   ( ans == 0)

    x_functioncalculate_result = 0.5 *(a_functioncalculate + b_functioncalculate)

    call f(x_value,x_functioncalculate_result)
    call f(a_value,a_functioncalculate)
    call f(b_value,b_functioncalculate)

    if ( x_value==0 .OR. ((b_functioncalculate-a_functioncalculate)/(2 ** counter_fc)) < 0.5*10**(-r)) then
        ans=1
    else if ( x_value *a_value<=0 )then 
        b_functioncalculate = x_functioncalculate_result
    else if ( x_value *b_value<=0 )then
        a_functioncalculate = x_functioncalculate_result 
    end if 
    counter_fc = counter_fc + 1

end do  

End Subroutine function_calculate
!=============================================== End Of Subroutines And Functions =============================================
