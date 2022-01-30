
!================================================= Gausse Seidel Method ============================================

!============================================   Definition Of Modules   ==========================================

!   n     :/input/:    : Number of Equations
!   error :/input/:    : error :)

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


!   Another form of euqations :  X1 = (1/a11)  * (-(a12 X2 + ... + a1n Xn) + b1) = P12 X2 + ... + P1n Xn + K1
!                                X2 = (1/a22)  * (-(a21 X1 + ... + a2n Xn) + b2) = P21 X1 + ... + P2n Xn + K2
!                                              ........................................
!                                Xn = (1/ann)  * (-(an1 X1 + ... + ann Xn) + bn) = Pn1 X1 + Pn2 X2 + ... + Kn
! ========================


!   Matrix T  (n , n)         :   -                   -
!                                 | 0   P12  ...  P1n |
!                                 | P21  0   ...  P2n |
!                                 | ................. |
!                                 | Pn1 Pn2  ...   0  |
!                                 -                   -

!      Vector X   :/input/:   : The Initial values of Xi which are determined (given) by the user :
!                                    X = ( X1  X2  ...  Xn )t
!       
!      Vector Norm_T          : Summation of Each line of T Matrix.

!      m                      : The Maximum array of Norm_T Vector.

!      i,j                    : counters 
 
!       

!========================================================== Modules =================================================

Module Para

Integer, Parameter :: dp = 8

end module Para

!=============================================================


module arr
use para


Real(Kind = dp),Allocatable,dimension(:,:) :: Augmented_Matrix
Real(Kind = dp),Allocatable,dimension(:,:) :: T,X

Real(Kind = dp),Allocatable,dimension(:) :: Norm_T




end module arr

!==============================================================

module var
use para

Real(Kind = dp)   :: sum1 , m , sum2,error
integer :: n,i,j


end module var

!================================================ Main Program ======================================================


program gausse_seidel_methode

use para
use arr
use var

Implicit None


open(101,file='input.txt')

Read(101,*),n !n: number of equations
Read(101,*),error

Allocate(Augmented_Matrix(n,n+1))
Allocate( T(n,n))


Allocate( X(n,2))
Allocate( Norm_T(n))
 
 
 !======================= Reading Augmented_Matrix ==================
do i=1,n
Read(101,*),Augmented_Matrix(i,1:n+1)
end do
  
  
 ! ===================== Reading Initial Vector X ====================
do i=1,n
read(101,*),x(i,1)
end do


! ====================== Formation of T matrix =======================


Norm_T=0


do i=1,n,1
   do j=1,n,1
	   T(i,j) = - Augmented_Matrix(i,j) / Augmented_Matrix(i,i)
	end do
	t(i,i)=0
end do


! ======================= Formation of Norm_T Vector ===================

do i=1,n,1
    do j = 1,n,1
	Norm_T (i) = T (i,j) + Norm_T(i)
	end do
end do

! ======================= The maximum array of Vector norm T============
m = maxval(Norm_T)
open(102,file='result.txt')
if ( m >= 1 ) then
   write(102,*), "divergent"
   pause 
   stop 
end if



!============================ Calculating Results =======================
error=1
do while(error>0.000001) !  You can Improve The Result's Accuracy by changing This number :)

do i=1,n,1
   sum1 = 0
   do j = 1 , i-1 , 1
        sum1 = sum1 + Augmented_Matrix (i,j) * X(j,2)
	end do 
   sum2 = 0
   do j = i+1 , n,1
        sum2 = sum2 + Augmented_Matrix(i,j) * X(j,1)
   end do

X(i,2) = (Augmented_Matrix( i, n+1 ) -sum1-sum2) / Augmented_Matrix(i ,i) 
end do

error=0
do i=1,n,1
    error = error+(x(i,1)-x(i,2))**2
end do
error=sqrt(error)/n


x(:,1) = X(:,2) 

end do

write(102,*), "============ The Result Vector : ==============="
do i=1,n
    write(102,*),x(i,2)
enddo

close(101)
close(102)

deallocate (Augmented_Matrix,T,X,Norm_T)


End Program gausse_seidel_methode

!==============================================    End of The Program   ====================================================== 






















