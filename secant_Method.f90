!================================================================ Secant_Method =======================================================================


!=========================================== Definition Of Modules ============================================================================================


!   Form of Function   :  F(x)=  a1 * x^(n-1)    +  a2 * x^(n-2) +...+ an : Polynomial
!   Goal : Finding The Root Of The Function : F(x) = 0

!   n      : :/input/:   : Equation Order= n-1

!   Coefficient_Matrix : :/input/: : Vector/ Number of Arrays=n : ( a1  a2  a3  ....  an )

!   M      : :/input/:   :  Max Number Of loop's Iteration 

!   X0,X1  : :/input/:   : Initial value : This Method starts From two Initial Values.

!   Xn     : :/output/:  : Root Of Function : Result

!   e      : :/input/:   : Error Coresponding TO Accuracy Of The Result
   
!   i,j    : counters

!   zn     :/output/:    (Variable)

!============================================= Modules =====================================================================================================
Module Para

Integer, Parameter :: dp = 8

end module Para

!==============================================================

module arr
use para

Real(Kind = dp),Allocatable,dimension(:) :: Coefficient_Matrix


end module arr

!==============================================================

module var
use para

Real(Kind = dp)   :: e,M,xn,x0,yn,F_result_xn,F_result_yn ,x1,zn
integer :: n,i,j


end module var


!================================================================ Main Program =======================================================================
Program Secant_Method
!F(x)=  a1 * x^(n-1)    +  a2 * x^(n-2) +...+ a0

use para
use arr 
use var

Implicit None


!====== Reading input Data by The Order Below ================

Real(Kind = dp)  :: ans
open(101,file='input.txt')
read(101,*),x0
read(101,*),x1
read(101,*),e    
read(101,*),M   
read(101,*),n         
allocate(Coefficient_Matrix(n))     
read(101,*),Coefficient_Matrix(1:n)
close(101)

open(102,file='result.txt')

xn = x0
yn = x1
ans = 0._dp
j=0



do while ( ans==0 )

call f(F_result_xn ,xn)
call f(F_result_yn ,yn)



zn = yn - (((F_result_yn)*( yn - xn )) / (F_result_yn - F_result_xn ))

if(j>M) then
   write(102,*)," ERROR, DIVERGENt " 
   write(102,*),zn
   ans=1
   else   if ( abs( zn - yn )<= e .OR. F_result_yn == 0 ) then
    write(102,*), " Answer is ; x =  "
    write(102,*),zn

    ans=1
    else
    j=j+1
    xn = yn
    yn = zn
end if

end do
   
 close(102)
 
    deallocate(Coefficient_Matrix)
    
    
pause
end program Secant_Method
!=================================================================== End Of Main Program ======================================================================

 
!================================================================= Subroutines And Functions ===================================================================

!========== Subroutine f ========= calculating F(x) value ===================
! Subroutine inpput : xf
! Subroutine output : Fvalue

Subroutine f(Fvalue,xf) !calculating F(x) value 

use para
use arr
use var


Implicit None 

Real*8,INTENT(IN) :: xf
Real*8,INTENT(out) :: Fvalue

Fvalue = 0.0_dp
do i=1,n
    Fvalue = Fvalue + Coefficient_Matrix(i)*xf**(n-i)
End do 


END Subroutine f
!==================================================================== End Of Subroutines And Functions ===========================================================