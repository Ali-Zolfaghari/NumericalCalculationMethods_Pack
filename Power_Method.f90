!================================================= Power Method / Finding eigenvector and eigenvalue ============================================

!============================================   Definition Of Modules   ==========================================

! A_Matrix : :/input/: : A n*n Matrix Which It's Corresponding eigenvectors and eigenvalues will be calculated.
! n        : :/input/: : The Order Of A_Matrix
! X        : :/input/: : A Vector , An Approximation of The Initial eigenvector 
! number   : :/input/: : The Maximum Number Of Computation's Iterations
! i,j      : counter
! norm_x   : Maximum Array Of Matrix X
! norm_y   : Maximum Array Of Matrix y
! Yp,Xp    : variables
!===================================================== Modules========================================================

Module Para


Integer, Parameter :: dp = 8

end module Para
!========================================================

module arr
use para

Real(Kind = dp),Allocatable,dimension(:,:) :: A_Matrix
Real(Kind = dp),Allocatable,dimension(:) :: X
Real(Kind = dp),Allocatable,dimension(:) :: Y
Real(Kind = dp),Allocatable,dimension(:) :: E



end module arr

!=======================================================
module var
use para

Real(Kind = dp)   :: m , norm_E , error , Xp , Yp , norm_x , sum , norm_y
integer :: n,i,j,k,number


end module var
!==========================================================================================================================

!======================================== Main Program=====================================================================
Program Power_Method


use para
use arr
use var

Implicit None


open(101,file='input.txt')
open(102,file='result.txt')

!Read(101,*),error 

Read(101,*),number 
Read(101,*),n

Allocate(A_Matrix(n,n))
Allocate(X(n)) ! initial  approximation of eigenvector 
Allocate(Y(n))
Allocate(E(n))


do i=1,n
Read(101,*),A_Matrix(i,1:n)
end do
do i=1,n
write(102,*),A_Matrix(i,1:n)
end do

do i=1,n
Read(101,*),X(i)
end do
do i=1,n
write(102,*),x(i)
end do


!=================formation of norm_X==========================

norm_x = maxval(X)

do i=1,n,1
if ( ABS(X(i))== norm_x )then
    Xp=norm_x
    exit
end if
end do
 
do i=1,n,1
X(i)= X(i)/Xp
end do

do k=1,number,1
! ============= Y=A*X ================================================

do i=1,n,1
 sum=0
 do j=1,n,1
 sum=sum+A_Matrix(i,j)*X(j)
 end do
 Y(i)=sum
end do

!===================== formation of norm_Y =========================

norm_Y = maxval(abs(Y))
 

 
do i=1,n,1
if ( ABS(Y(i))== norm_Y )then
    Yp=y(i)
    exit
end if
end do
 
E=0
 
if (Yp==0) then
write(102,*)," eigenvector Corresponding to eigenvalue of '0' "
do i=1,n
write(102,*), X(i)
end do
pause
stop
end if

do i=1,n,1
E(i)= X(i)- Y(i)/Yp
end do

norm_E = maxval( abs(E))

do i=1,n,1
X(i)= Y(i)/Yp
end do



!if ( norm_E < error ) then
end do 

write(102,*)," eigen value is :  "

write(102,*), Yp ! eigen value 
write(102,*), " eigen vaector is : " 

 
do i=1,n
    write(102,*),x(i)
enddo

 

 
close(102)
close(101)
 
 deallocate (A_Matrix,X,Y,E)

 
 
 
 
end program Power_Method

!======================================= End Of Program =====================================================


















