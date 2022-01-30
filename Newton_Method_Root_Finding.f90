!=========================================== Newton_Method_Root_Finding =======================================================================================

!=========================================== Definition Of Modules ============================================================================================

!   Form of Function   :  F(x)=  a1 * x^(n-1)    +  a2 * x^(n-2) +...+ an : Polynomial
!   Goal : Finding The Root Of The Function : F(x) = 0
!   n    : :/input/:   : Equation Order= n-1

!   Coefficient_Matrix : :/input/: : Vector/ Number of Arrays=n : ( a1  a2  a3  ....  an )

!   M    : :/input/:   :  Max Number Of loop's Iteration 

!   X0   : :/input/:   : Initial value : This Method starts From An Initial Value.

!   Xn   : :/output/:  : Root Of Function : Result

!   e    : :/input/:   : Error Coresponding TO Accuracy Of The Result

!   Order: :/input/:   : In Newton Method The Root Of The Function Can be Approximated , The order is a Coefficient , used in calculations, Which can Improve 
!                        The Result's Accuracy. For unlnown Orders Use Order=1.

!   i,j  : counters

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

Real(Kind = dp)   :: e,M,xn,x0,yn,F_result,dif_result ,order
integer :: n,i,j


end module var
!==========================================================================================================================================================



!===================================================================== Main Program =======================================================================

Program Newton_Method


use para
use arr 
use var

Implicit None

Real(Kind = dp)  :: ans
open(101,file='input.txt')
read(101,*),x0
read(101,*),order     
read(101,*),e         
read(101,*),M         
read(101,*),n         
allocate(Coefficient_Matrix(n))  

!=============== Reading Coefficient_Matrix ======================     
do i=1,n
read(101,*),Coefficient_Matrix(i)
end do
close(101)

open(102,file='result.txt')


xn = x0
ans = 0._dp
j=0

do while ( ans==0 )

call f(F_result,xn)
call Differential_f(dif_result,xn)



yn = xn - (order*(F_result / dif_result))

if(j>M) then
   write(102,*),"  ERROR,DIVERGENCE  "  
   write(102,*),"  Initial Value " ,xn
   ans=1
   else   if ( abs( yn - xn )<= e .OR. F_result==0 ) then
  write(102,*)," answer =  "
  write(102,*),xn
    ans=1
    else
    j=j+1
    xn = yn
end if

if(ans==1)then
exit
end if

end do  


close(102)    
    deallocate(Coefficient_Matrix)
pause
end program Newton_Method
 
!=================================================================== End Of Main Program ======================================================================
 
!================================================================= Subroutines And Function ===================================================================

!========== Subroutine f ========= calculating F(x) value ===================
! Subroutine inpput : xf
! Subroutine output : Fvalue


Subroutine f(Fvalue,xf) 


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
!====================================================================

!========== Subroutine Differential_f ====== calculating Differential  Of f(x) ===========
! Subroutine inpput : xf
! Subroutine output : dif_f



Subroutine Differential_f(dif_f,xf) 


use para
use arr
use var 

implicit None

Real*8,INTENT(IN)  ::  xf
Real*8,INTENT(OUT)  ::  dif_f

dif_f = 0.0_dp
do i=1,n-1
    dif_f = dif_f +  Coefficient_Matrix(i)*(n-i)*xf**(n-i-1)
End do


End subroutine Differential_f

!============================================================= End Of Subroutines And Function =================================================================
