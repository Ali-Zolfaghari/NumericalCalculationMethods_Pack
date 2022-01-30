!=============================================== Newton Divided Difference Method _ Interpolating Method ===========================================


!====================================================  Definition Of Modules      ===============================================

! X             :/input/:  : The value Which It's Corresponding Value Is Interpolated
! n             :/input/:  : Number of Input data points
! X_Matrix :    :/input/:  : Input data : X(i)
 
! Result_Matrix :  :/input/: A n*n Matrix Which It's First column is Yi_ The Corresponding Value to input data ( Xi )_
!                                            
! i,j           : counter
! f_result      :/output/:  : The Interpolated Value Corresponding to X .

!======================================================== Modules ================================================================
Module Para

Integer, Parameter :: dp = 8

end module Para

!==============================================================

module arr
use para

Real(Kind = dp),Allocatable,dimension(:,:) :: Result_Matrix
Real(Kind = dp),Allocatable,dimension(:) :: X_Matrix



end module arr

!==============================================================

module var
use para

Real(Kind = dp)   ::mltp,x,f_result
integer :: n,i,j


end module var

!======================================================= Main Program =================================================================



Program Newton_Divided_Difference

use para 
use var
use arr

implicit none

open(101,file='input.txt')
open(102,file='result.txt')


!============= Reading function input===================
read(101,*),x
write(102,*)," X : function input"
write(102,*),X

!=========== Reading number of data points =============
read(101,*),n 
allocate(Result_Matrix(n,n))
allocate(X_Matrix(n))


!=============== Reading Input Data ====================
write(102,*)," Xi = input data"
do i=1,n,1
read(101,*),X_Matrix(i)
end do
do i=1,n,1
write(102,*),X_Matrix(i)
end do

Result_Matrix=0
!============= Reading The First Column Of Result_Matrix / Y(i) =============================
write(102,*)," Yi = output data"
do i=1,n,1
read(101,*),Result_Matrix(i,1)
end do
do i=1,n,1
write(102,*),Result_Matrix(i,1)
end do


!===================== Interpolating ===============================
do i=2,n

do j=2,i

Result_Matrix(i,j)= (Result_Matrix(i,j-1) - Result_Matrix(i-1,j-1) )/(X_matrix(i)- X_matrix(i-j+1))

end do

end do

F_Result=0
mltp=1

do i=1,n

if(i==1) then
mltp=1
else 
mltp=mltp* ( x - X_Matrix(i-1))
end if 

F_Result= Result_Matrix(i,i)*mltp +F_Result

end do 

write(102,*)," Interpolated value calculated by Newton Divided Difference Method is :"
write(102,*),F_Result

deallocate(X_Matrix,Result_Matrix)

close(101)
close(102)

End Program Newton_Divided_Difference
!========================================================== End Of Main Program =============================================