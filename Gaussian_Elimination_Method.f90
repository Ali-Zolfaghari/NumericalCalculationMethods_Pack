!=========================================== Gaussian_Elimination_Method  ==========================================

!============================================   Definition Of Modules   ==========================================

!   n      :/input/:   : Number of Equations
!   error     : error :)

!   euqations : a11 X1 + a12 X2 + ... + a1n Xn = b1
!               a21 X1 + a22 X2 + ... + a2n Xn = b2
!               .................................
!               an1 X1 + an2 X2 + ... + ann Xn = bn
! =======================    


!    Augmented_Matrix (n,n+1)  : :/input/:
!                                 -                    -
!                                 | a11 a12 ... a1n b1 |
!                                 | .................. |
!							      | an1 an2 ... ann b2 |
!                                 -                    -
! ======     ======     ======     ======      ======      ====== 


!      Vector Result_Matrix   :/output/:   : The final values calculated  

!      Vector A_Matrix        : An Auxiliary Vector for displacement of two lines in Matrix
!       
!      m                      : Coefficient

!      i,j,counter            : counters 

!      w,p,warn               : warning parameters for switching cases 

!===================================================================================================

!========================================= Modules =================================================
Module Para

Integer, Parameter :: dp = 8

end module Para

!==============================================================


module arr
use para


Real(Kind = dp),Allocatable,dimension(:,:) :: Augmented_Matrix
Real(Kind = dp),Allocatable,dimension(:) :: A_Matrix
Real(Kind = dp),Allocatable,dimension(:) :: Result_Matrix

end module arr

!==============================================================

module var
use para

Real(Kind = dp)   :: sum,m
integer :: n,i,j,counter,warn,p,w


end module var


!=============================================  Main Program =====================================

Program Gaussian_Elimination_Method


use para
use arr
use var

Implicit None


open(101,file='input.txt')

Read(101,*),n !n: number of equations

Allocate(Augmented_Matrix(n,n+1))
Allocate(A_Matrix(n+1))
Allocate(Result_Matrix(n))


!=====================Read Augmented_Matrix ============================
do i=1,n
Read(101,*),Augmented_Matrix(i,1:n+1)
enddo



Warn = 0 

Do i=1, n-1 ,1
Do p=i , n ,1
if ( Augmented_Matrix(p,i) /= 0 ) then
exit
else if ( p>= n) then
warn=1
end if
end do


!================== Displacement of two Lines  ============================
if ( p /= i ) then
Do  counter= 1, n+1, 1
A_Matrix(counter) = Augmented_Matrix(p,counter) 
Augmented_Matrix(p,counter)=Augmented_Matrix(i,counter)
Augmented_Matrix(i,counter)=A_Matrix(counter)
End Do
End if  

w=0

do j=i+1,n
  if (Augmented_Matrix(j,i) /= 0 ) then
  w=1 
  end if
end do

if (w==1) then
Do j= i+1 ,n, 1
m = Augmented_Matrix(j,i) / Augmented_Matrix(i,i)
Do counter = 1, n+1,1
Augmented_Matrix(j,counter) = (-m) * Augmented_Matrix(i,counter) + Augmented_Matrix(j,counter)
End Do
End Do
end if


end do
open(102,file='result.txt')

if ( Augmented_Matrix(n,n) == 0) then 
write(102,*)," =====================    THE MATRIX IS NOT REVERSIBLE AND THERE IS NO RESULT   ==================="
pause
stop
End if

!================================== calculating result matrix =====================================
 Result_Matrix=0     

Do i = n, 1, -1
sum =0
Do j = n , i , -1
sum =   Augmented_Matrix(i,j) *  Result_Matrix(j) + sum
End Do

Result_Matrix(i) = ( Augmented_Matrix(i,n+1)  - sum ) / Augmented_Matrix(i,i) 
End Do


do i=1,n
write(102,*), Result_Matrix (i)
enddo 






if ( warn == 1 ) then
write(102,*)," =============     THE MATRIX IS NOT REVERSIBLE AND THERE IS NO RESULT  ====================="
End if
close(102)
close(101)
deallocate (Augmented_Matrix,A_Matrix,Result_Matrix)

End Program Gaussian_Elimination_Method        

!=========================================  End of Program ===============================================