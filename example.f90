program main
   use  m_unify_zmat
   implicit none
   integer, parameter :: N=70, lwork = 3*N
   real(kind=8), allocatable :: H(:,:), vecs1(:,:), vecs2(:,:), eig1(:), eig2(:), work(:)
   integer :: info, i, ierr

   allocate(H(N,N), vecs1(N,N), vecs2(N,N), eig1(N), eig2(N), work(lwork))

   ! Hermitian matrix with a few degenerate values
   H = 0.0
   do i = 1,N
      H(i,i) = floor(sqrt(1.0*(i)))
   enddo 
   
   ! hermitianize
   H = 0.5*H +  0.5*transpose(H)
   call shuffle_matrix(H)

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

   call unify_zmat(eig1, vecs1, ierr)
   if(ierr /= 0) stop
   call unify_zmat(eig2, vecs2, ierr)
   if(ierr /= 0) stop

   write (*,*) new_line("a") // " After:"
   !write (*,*) "Eigenvalues:", eig1
   write (*,*) "Eigenvector diff norm:", norm2(vecs1 - vecs2) 
   write (*,*) "Eigenvector maxdiff:", maxval(abs(vecs1 - vecs2))

   call investigate(eig1, vecs1, vecs2)
   !call check_orthogonality(vecs1)
   !call check_eigenvectors(H, vecs1, eig1)
contains
   subroutine shuffle_matrix(H)
      implicit none 
      real(kind=8), intent(inout) :: H(:,:)

      real(kind=8), allocatable  :: P(:,:), tmp(:,:)
      real(kind=8), parameter    :: one = 1.0, zero = 0.0
      real(kind=8)               :: theta
      integer :: i,j, dim 

      dim = size(H,1)

      allocate(P, tmp, mold=H)

      do i = 1,floor(dim/2.0)
         j = N-i+1
         theta = 0.01*i
         call givens(i,j, theta, P)

         ! tmp = H * P
         call dgemm("N", "N", dim, dim, dim, one, H, dim, P, dim, zero, tmp, dim)

         ! H = P^-1 * tmp
         call dgemm("T", "N", dim, dim, dim, one, P, dim, tmp, dim, zero, H, dim)
      enddo
   end subroutine shuffle_matrix

   subroutine givens(i,j,theta, mat)
      implicit none 
      integer, intent(in)         :: i,j
      real(kind=8), intent(in)    :: theta 
      real(kind=8), intent(inout) :: mat(:,:)
      integer :: k

      mat = 0.0

      do k = 1,N 
         mat(k,k) = 1 
      enddo 

      mat(i,i) = cos(theta)
      mat(j,j) = cos(theta)
      mat(j,i) = sin(theta)
      mat(i,j) = - sin(theta)
   
   end subroutine givens

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

   subroutine check_orthogonality(mat)
      implicit none
      real(kind=8), intent(in) :: mat(:,:)
      integer :: i,j 

      do j = 1, size(mat,2)
         do i = 1, j-1 
            if(abs(dot_product(mat(:,i), mat(:,j))) > 1e-8 ) then
               write (*,*) i, j, "not orthogonal", dot_product(mat(:,i), mat(:,j))
            endif
         enddo 
      enddo 

      do i = 1,size(mat,2)
         if(abs(dot_product(mat(:,i), mat(:,i)) - 1) > 1e-8 ) then 
            write (*,*) i, "not normalized"
         endif
      enddo
   end subroutine check_orthogonality

   subroutine check_eigenvectors(A, u, eigval)
      implicit none 
      real(kind=8), intent(in)  :: A(:,:), u(:,:), eigval(:)
      real(kind=8), allocatable :: diff(:)

      real(kind=8), parameter :: one=1.0, zero=0.0 

      do i = 1,size(eigval)
         diff = matmul(A,u(:,i)) - eigval(i)* u(:,i)
         if(norm2(diff) > 1e-10) then
            write (*,*) i, "diff"!, diff 
         endif
      enddo
   end subroutine check_eigenvectors
end program main
