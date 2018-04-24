subroutine one_element1(a1,b1,a2,b2,val,S,H,EE,n_bas)
integer i,j,k,l,n_bas
integer a1(6,4),a2(4,4),b1(4,4),b2(6,4)
integer aa1(6,3),aa2(4,3),bb1(4,3),bb2(6,3)
real*8 v,val,val1,val2,nm,norm_aa1,norm_aa2,norm_bb1,norm_bb2
real*8  S(n_bas,n_bas),H(n_bas,n_bas),EE(n_bas,n_bas,n_bas,n_bas)
real*8 over
!--------------------------------------------------------------------------------
! This subroutines compute 
!  part <PHI(i)|H|PHI(j)>  where 
! PHI(i)=phi_11(i)*alpha_1 + phi_12(i)*alpha_2
! PHI(i)=phi_11(j)*alpha_1 + phi_12(j)*alpha_2
! and <alpha_1|alpha_2> = 0 i.e orthogonal
!--------------------------------------------------------------------------------

aa1=a1(:,2:4)   ! storing the components of phi_11(i) 
aa2=a2(:,2:4)   ! storing the components of phi_12(j)


 norm_aa1=0.d0;
 norm_aa2=0.d0;

! Computing norm square of phi_11(i)

 do i=1,6
  do k=1,6
  call norm1(aa1(i,:),aa1(k,:),nm,S,H,EE,n_bas)
  norm_aa1=norm_aa1+a1(i,1)*a1(k,1)*nm/4.d0
  end do
 end do

! Computing norm square of phi_11(j)

 do i=1,4
  do k=1,4
   call norm1(aa2(i,:),aa2(k,:),nm,S,H,EE,n_bas)
   norm_aa2=norm_aa2+a2(i,1)*a2(k,1)*nm*9.d0/4.d0
  end do
 end do

! Computing the <phi_11(i)|H|phi_11(j)>

val1=0
 do i=1,6
   do j=1,4
       call ham(aa1(i,:),aa2(j,:),v,S,H,EE,n_bas)
       val1=val1+a1(i,1)*a2(j,1)*v*3.d0/4.d0
   end do
 end do



bb1=b1(:,2:4);    ! storing the components of phi_12(i)
bb2=b2(:,2:4);    ! storing the components of phi_12(j)


 norm_bb1=0.d0;
 norm_bb2=0.d0;

! Computing norm square of phi_12(i)

 do i=1,4
  do k=1,4
   call norm1(bb1(i,:),bb1(k,:),nm,S,H,EE,n_bas)
   norm_bb1=norm_bb1+b1(i,1)*b1(k,1)*nm*3.d0/4.d0
  end do
 end do

! Computing norm square of phi_12(j)

 do i=1,6
  do k=1,6
  call norm1(bb2(i,:),bb2(k,:),nm,S,H,EE,n_bas)
  norm_bb2=norm_bb2+b2(i,1)*b2(k,1)*nm*3.d0/4.d0
  end do
 end do

! Computing <phi_12(i)|H|phi_12(j)>

val2=0;
  do i=1,4
   do j=1,6
       call ham(bb1(i,:),bb2(j,:),v,S,H,EE,n_bas)
       val2=val2+b1(i,1)*b2(j,1)*v*3.d0/4.d0
   end do
 end do

  over = 0.d0
! Computing overlap between psi_i and psi_j
    do i=1,6
      do k=1,4
         call norm1(aa1(i,:),aa2(k,:),nm,S,H,EE,n_bas)
         over=over+a2(i,1)*a2(k,1)*nm
      end do
    end do
   write(6,*) over

   over = 0.d0

     do i=1,4
       do k=1,6
          call norm1(bb1(i,:),bb2(k,:),nm,S,H,EE,n_bas)
          over=over+b1(i,1)*b1(k,1)*nm
       end do
     end do
     write(6,*) over

! This is more crucial step 
!
!1. This assuming the alpha1 and alpha2 are already orthonormal and normalizing it with norms
!   of individual part
!
   val= (val1+val2)/(sqrt((norm_aa1+norm_bb1))*sqrt((norm_aa2+norm_bb2)))
!   write (6,*) (sqrt((norm_aa1+norm_bb1))*sqrt((norm_aa2+norm_bb2))),norm_aa1/4.0,norm_bb1*3.0/4.0,norm_aa2/4.0,norm_bb2*3.0/4.0
!2. This assuming the alpha1 and alpha2 are already orthonormal and using the direct expression
!   the book by Ruben Paunz
!
! val=(val1+val2)/(6.d0)

!3. This assuming the alpha1 and alpha2 are not normalized and get constructed from conjugate young
!   tabulae
!   alpha1= sqrt(3)[beta*alpha*alpha - alpha*alpha*beta]
!   alpha2= [2*alpha*alpha*beta - beta*lpha*alpha - alpha*beta*alpha ]
!   <alpha1|alpha1> = 6 ; <alpha2|alpha2>=6
!
!  val= (val1*6.d0+val2*3.d0)/(sqrt((norm_aa1*6.d0+norm_bb1*3.d0))*sqrt((norm_aa2*6.d0+norm_bb2*3.d0)))
end subroutine one_element1
