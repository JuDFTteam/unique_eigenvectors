program main
   use  m_unify_zmat
   implicit none
   integer, parameter :: N= 200, lwork = 3*N
   real(kind=8)       :: H(N,N), vecs1(N,N), vecs2(N,N), eig1(N), eig2(N), work(lwork)
   integer :: info, i
   
   ! Hermitian matrix with a few degenerate values
   H = 0.0
   do i = 1,N
      H(i,N-i+1) = floor(sqrt(1.0*i))
   enddo 

   ! hermitianize
   H = 0.5*H +  0.5*transpose(H)

   vecs1 = H 
   vecs2 = H 

   call dsyev("V", "U", N, vecs1, N, eig1, work, lwork, info)
   if(info /= 0) write (*,*) "diag U failed"

   call dsyev("V", "L", N, vecs2, N, eig2, work, lwork, info)
   if(info /= 0) write (*,*) "diag L failed"

   write (*,*) "Before:"
   write (*,*) "Eigenvalue diff norm:", norm2(eig1  - eig2)
   write (*,*) "Eigenvector diff norm:", norm2(vecs1 - vecs2) 
   write (*,*) "Eigenvector maxdiff:", maxval(abs(vecs1 - vecs2))

   call unify_zmat(eig1, vecs1)
   call unify_zmat(eig2, vecs2)

   write (*,*) new_line("a") // " After:"
   !write (*,*) "Eigenvalues:", eig1
   write (*,*) "Eigenvector diff norm:", norm2(vecs1 - vecs2) 
   write (*,*) "Eigenvector maxdiff:", maxval(abs(vecs1 - vecs2))

   call investigate(eig1, vecs1, vecs2)

contains
   subroutine investigate(eigval, eigvec1, eigvec2)
      implicit none 
      real(kind=8), intent(in) :: eigval(:), eigvec1(:,:), eigvec2(:,:)
      integer :: beg_group, end_group, n_g
      integer, allocatable :: groups(:)
      integer :: deg

      groups = make_groups(eigval)
   

      write (*,*) "degeneracy    eigenvalue             norm diff"
      do n_g = 1, size(groups)
         beg_group = sum(groups(1:n_g-1)) + 1
         end_group = sum(groups(1:n_g))

         deg = end_group - beg_group + 1
         write (*,*) deg, eigval(beg_group), norm2(eigvec1(:,beg_group:end_group) - eigvec2(:,beg_group:end_group))
      enddo 
   end subroutine investigate
end program main