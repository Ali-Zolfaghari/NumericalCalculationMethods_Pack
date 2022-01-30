!======================================== Interpolation_Newton_Method ======================================

!=========================================  Definition Of Variables   ====================================

! n    :/input/:      : number of input data of Function .

! x    :/input/:      : A Real Number Which It's value will be interpolated by the program .

! Matrix              : A Matrix which it's arrays will be used to calculate interpolated values
!                       First coloumn : :/input/: : Xi : input function values
!                       2nd   coloumn : :/input/: : Yi : output function values

! i,j,counter ,p      : counter

! h                   : A constant Value , Diffrence of input function values 

! n_func              : input of The Factorielle Function

! n_fctrl_func        : Result of The Factorielle Function : n! :)

! f_value , f_result  :/output/:  : Inerpolated value corresponding to X 

! w                   : For Switch Cases
!



!===========================================     Modules  ===================================================
Module Para

Integer, Parameter :: dp = 8

end module Para

!==============================================================

module arr
use para

Real(Kind = dp),Allocatable,dimension(:,:) :: Matrix


end module arr

!==============================================================

module var
use para

Real(Kind = dp)   :: sum,mltp,q,h,x,f_result,f_value, n_function, n_fctrl_function, n_fctrl
integer :: n,i,j,m,counter,k,p,w


end module var
!===========================================================================================================

!======================================= Main Program ======================================================
Program Interpolation_Newton_Method

use var
use para 
use arr 
 
 
implicit none


open(101,file='input.txt')
open(102,file='Result.txt')

write(102,*)," n: number of input data of Function "
 read(101,*),n
write(102,*),n

Allocate(Matrix(n,n+1))

Matrix=0


!================== Reading function input ============================================================
write(102,*)," X : function input "
read(101,*),x
write(102,*),x


write(102,*)," Xi : input data "
do i=1,n
read(101,*),Matrix(i,1)
end do 
do i=1,n
write(102,*),Matrix(i,1)
end do 

write(102,*)," Yi : F( Xi ) "
do i=1,n
read(101,*),Matrix(i,2)
end do
do i=1,n
write(102,*),Matrix(i,2)
end do 


h= Matrix(2,1)- Matrix(1,1) ! h is constant 

m=n

do j=3,m
do i=1,n-1

Matrix(i,j) = Matrix(i+1,j-1)- Matrix(i,j-1)

end do
n=n-1 
end do

f_value=0

!=======================   forward divided Method =================================
n=m
do i = 1,m,1

if ( x >= Matrix(i,1) .and. x <= Matrix(i+1,1) ) then

!========================  Calculating q =======================================
q =( x - Matrix(i,1))/ h
exit 
end if
end do

mltp=1.0_dp
F_Result =0.0_dp

do j=2,n+1

do counter=0,j-2,1
if (counter==0) then 
mltp=1.0_dp
else
mltp= mltp * (q - real(counter) +1.0_dp)
end if
end do

call n_Factorielle( real(j-2) , n_fctrl )

F_Result = F_Result + Matrix(i,j)* mltp / (n_fctrl)

end do 

write(102,*)," INTERPOLATED VALUE IS (Forward Divided) : "
Write(102,*),F_Result

!=========================  backward divided =========================================
do i = 1,m,1

if ( x >= Matrix(i,1) .and. x <= Matrix(i+1,1) ) then
q =( Matrix(i,1) -x )/ h
exit 
end if
end do

mltp=1.0_dp
f_value =0.0_dp
p=1
w=0

do j=2,n+1

if(w==1)then 
p=p+1
end if

do counter=0,j-2,1
if (counter==0) then 
mltp=1.0_dp
else
mltp= mltp * (q + real(counter) -1.0_dp)
end if
end do

call n_Factorielle( real(j-2) , n_fctrl )

if (Matrix (i+1,j)/= 0) then
f_value = f_value + Matrix(i+1,j)* mltp / (n_fctrl)
else 
f_value = f_value + Matrix(i+1-p,j)* mltp / (n_fctrl)
w=1
end if

end do


write(102,*)," INTERPOLATED VALUE IS (backward Divided) : "
Write(102,*),f_value

deallocate( matriX)
close(101)
close(102)

End Program Interpolation_Newton_Method

!======================================================== End of Program ====================================================================

!=================================================== Subroutines And Functions ==============================================================

!=== Subroutine n_Factorielle
!=== Subroutine input :  n_func 
!=== Subroutine output:  n_fctrl_func


Subroutine n_Factorielle ( n_func , n_fctrl_func)


 

use para
use arr
use var

Implicit None 

Real(Kind = dp) ,INTENT(INout) :: n_func
Real(Kind = dp) ,INTENT(inout) :: n_fctrl_func

n_fctrl_func= 1

do k=1, (n_func+1),1
 if ( k-1 == 0) then
n_fctrl_func =1
else 
n_fctrl_func = n_fctrl_func * (k-1)
end if
end do

END Subroutine n_Factorielle


!================================================ End of Subroutines And Functions ==============================================================