module m_unify_zmat
contains
   subroutine unify_zmat(eigval, z)
      use omp_lib
      implicit none
      real(kind=8), intent(in) :: eigval(:)
      real(kind=8)          :: z(:,:)


      integer :: beg_group, end_group, n_g
      integer, allocatable :: groups(:)

      groups = make_groups(eigval)
   
      do n_g = 1, size(groups)
         beg_group = sum(groups(1:n_g-1)) + 1
         end_group = sum(groups(1:n_g))

         if(end_group <= size(z,2))then 
            call unify_group_operator(beg_group, end_group, z)
         endif
      enddo
   end subroutine unify_zmat

   subroutine rq_decomp(mat, r, q, perm_mat)
      implicit none 
      real(kind=8), intent(inout) :: mat(:,:), r(:,:), q(:,:), perm_mat(:,:) 
      integer      :: dim, lwork, i,j, info, jpvt(size(mat,1))
      real(kind=8) :: tau(size(mat,1)), work(size(mat,1)*3), tmp(size(mat,1), size(mat,2))

      dim   = size(mat,1)
      lwork =size(mat,1)*3
      jpvt  = 0
      if(any(shape(mat) /= shape(r)) .or. any(shape(r) /= shape(mat))) then
         write (*,*) "need square matricies"
      endif

      tmp = transpose(mat)

      call dgeqpf(dim, dim, tmp, dim, jpvt, tau, work, info)
      if(info /= 0) write (*,*) "qr failed"

      r = tmp
      do i = 1,dim
         do j = i+1,3 
            r(j,i) = 0.0
         enddo
      enddo
 
      q = tmp 
      do i =1,3 
         do j=i,3 
            q(i,j) = 0.0
         enddo
      enddo

      call dorgqr(dim, dim, dim, Q, dim, tau, work, lwork, info)
      if(info /= 0) write (*,*) "q build failed"

      write (*,*) "jpvt", jpvt


      perm_mat = permutation_matrix(jpvt)
      Q = transpose(Q)
      R = transpose(R)

   end subroutine rq_decomp 

   function permutation_matrix(jpvt) result(p)
      implicit none 
      integer, intent(in) :: jpvt(:)
      real(kind=8)        :: p(size(jpvt),size(jpvt))
      integer :: i

      p = 0.0 

      do i = 1,size(jpvt)
         p(jpvt(i), i) = 1
      enddo
   end function permutation_matrix

   subroutine unify_group_operator(beg_group, end_group, z)
      implicit none
      integer, intent(in)         :: beg_group, end_group
      real(kind=8), intent(inout) :: z(:,:) 
      real(kind=8), allocatable   :: mtx(:,:), eigval(:), work(:), rq_mat(:,:), tau(:), R(:,:), Q(:,:), tmp(:,:)
      integer :: dim, lwork, info, i, j
      real(kind=8) :: work_size(1)
      real(kind=8), parameter :: zero = 0.0, one = 1.0

      dim = end_group - beg_group + 1 
      allocate(rq_mat(dim, dim), tau(dim), R(dim,dim), Q(dim,dim))
      call make_rq_mat(z(:,beg_group:end_group), rq_mat)
      
      rq_mat = transpose(rq_mat)

      call dgeqrfp(dim, dim, rq_mat, dim, tau, work_size, -1, info)
      if(info /= 0) write (*,*) "problem A"
      lwork = int(work_size(1))
      allocate(work(lwork))
      call dgeqrfp(dim, dim, rq_mat, dim, tau, work, lwork, info)
      if(info /= 0) write (*,*) "problem B"

      deallocate(work)

      q = rq_mat
      r = rq_mat

      do i = 1,dim
         do j = i+1,dim
            r(j,i) = 0.0
         enddo
      enddo
 
      do i =1,3 
         do j=i,dim
            q(i,j) = 0.0
         enddo
      enddo


      call dorgqr(dim, dim, dim, Q, dim, tau, work_size, -1, info)
      if(info /= 0) write (*,*) "problem C"
      lwork = int(work_size(1))
      allocate(work(lwork))
      call dorgqr(dim, dim, dim, Q, dim, tau, work, lwork, info)
      if(info /= 0) write (*,*) "problem D"

      tmp = z(:,beg_group:end_group)
      ! z(:,beg_group:end_group) = matmul(z(:,beg_group:end_group), q)
      call dgemm("N", "N", size(z,1), dim, dim, one, tmp, size(z,1), q, dim, zero, z(:,beg_group:end_group), size(z,1))

   end subroutine unify_group_operator

   subroutine unify_group_lsq(beg_group, end_group, z)
      implicit none
      integer, intent(in)    :: beg_group, end_group
      real(kind=8), intent(inout) :: z(:,:) 
      real(kind=8), allocatable   :: targ(:,:), lhs(:,:), new_basis(:,:)
      real(kind=8)                   :: solution((end_group - beg_group) + 1, (end_group - beg_group) + 1)
      real(kind=8), parameter :: zero = 0.0, one = 1.0
      integer :: nrhs, basis_sz

      nrhs = (end_group - beg_group) + 1

      allocate(targ(size(z,1), nrhs), &
               lhs(size(z,1), nrhs), &
               new_basis(size(z,1), nrhs))
      call set_sin_targ(targ)
      lhs = z(:,beg_group:end_group)

      call leastsq(lhs, targ, solution)
      lhs = z(:,beg_group:end_group)

      basis_sz = size(lhs,1)
      !new_basis = matmul(lhs, solution)
      call dgemm("N", "N", basis_sz, nrhs, nrhs, one, lhs, size(lhs,1), solution, size(solution,1), &
            zero, new_basis, size(new_basis,1))
      !call mod_gram_schmidt(new_basis)

      z(:,beg_group:end_group) = new_basis
   end subroutine unify_group_lsq

   subroutine leastsq(lhs, targ, solution)
      implicit none 
      real(kind=8), intent(inout) :: lhs(:,:), targ(:,:), solution(:,:)
      integer :: m, n, nrhs, lwork, info
      real(kind=8) :: rwork_req(1)
      real(kind=8), allocatable :: rwork(:)


      m = size(lhs, 1)
      n = size(lhs, 2)
      nrhs = size(targ, 2)

      call dgels("N", m, n, nrhs, lhs, m, targ, size(targ,1), rwork_req, -1, info)
      lwork = int(rwork_req(1))

      allocate(rwork(lwork))
      call dgels("N", m, n, nrhs, lhs, m, targ, size(targ,1), rwork, lwork, info)
      if(info /= 0) write (*,*) "leastsq failed"

      solution = targ(:n,:)
   end subroutine leastsq

   subroutine set_sin_targ(targ)
      implicit none 
      real(kind=8), intent(inout) :: targ(:,:)
      integer :: i 
      real, allocatable :: x(:) 

      x = linspace(0.0,4.0, size(targ,1))
      do i = 1, size(targ,2)  
         targ(:,i) = sin(i*x)
      enddo
   end subroutine set_sin_targ

   function proj_r(u,v) result(p)
      implicit none 
      real(kind=8), intent(in) :: u(:), v(:)
      real(kind=8) :: p(size(u)) 

      p = dot_product(u,v) / dot_product(u,u) * u 
   end function

   function proj_c(u,v) result(p)
      implicit none 
      complex, intent(in) :: u(:), v(:)
      complex :: p(size(u)) 

      p = dot_product(u,v) / dot_product(u,u) * u 
   end function

   subroutine mod_gram_schmidt(M)
      implicit none 
      real(kind=8), intent(inout) :: M(:,:)
      
      integer :: k, i 

      do k = 1,size(M,2) 
         do i = 1,k-1 
            M(:,k) = M(:,k) - proj_r(M(:,i), M(:,k)) 
         enddo 

         M(:,k) = M(:,k) / norm2(M(:,k))
      enddo 
   end subroutine mod_gram_schmidt

   function linspace(beg, fin, n_steps) result(res)
      implicit none 
      real, intent(in)    :: beg, fin 
      integer, intent(in) :: n_steps
      real                :: res(n_steps), step

      integer :: i

      step = (fin - beg)  / (n_steps - 1.0)

      do i = 0,n_steps-1 
         res(i+1) = beg + step * i 
      enddo 
   end function linspace

   function make_groups(eigval) result(groups)
      implicit none
      real(8), intent(in)  :: eigval(:)
      integer, allocatable :: groups(:)
      integer              :: color(size(eigval))

      integer :: i, beg, g_cnt, n_groups

      g_cnt = 1
      beg  = 1 

      do while (beg <= size(eigval)) 
         i = beg
         do while (abs(eigval(i) - eigval(beg)) < 1e-8  ) 
            color(i) = g_cnt
            i = i + 1 
            if(i > size(eigval)) exit
         enddo 
         g_cnt = g_cnt + 1 
         beg = i 
      enddo


      n_groups = color(size(eigval))
      allocate(groups(n_groups))

      do i = 1,n_groups 
         groups(i) = count(i == color)
      enddo 
   end function make_groups

   subroutine make_rq_mat(eigvecs, rq_mat)
      implicit none
      real(kind=8), intent(in)                 :: eigvecs(:,:)
      real(kind=8), intent(inout), allocatable :: rq_mat(:,:)
      integer :: n, i, j
      real(kind=8) :: cutoff, lindep
      real(kind=8), allocatable :: tmp(:,:)

      n = size(eigvecs,2)
      cutoff = 0.9 * sqrt(1.0/size(eigvecs,1))
      if(allocated(rq_mat)) deallocate(rq_mat)
      allocate(rq_mat(n,n), tmp(n,n))

      rq_mat = 0.0
      lindep = 0.0  
      
      do i = 1,n 
         if(allocated(tmp)) deallocate(tmp)
         allocate(tmp(i,n))
         tmp = 0.0
         lindep = 0.0

         j = 0
         do while(lindep < cutoff)
            j = j + 1
            if(j > size(eigvecs,1)) then
               write (*,*) "problem with make_rq_mat"
            endif
            
            rq_mat(i,:) = eigvecs(j,:)
            tmp = rq_mat(1:i,:)
            lindep = linear_independency(tmp)
         enddo
      enddo
   end subroutine make_rq_mat

   function linear_independency(mat)
      implicit none
      real(kind=8), intent(inout) :: mat(:,:)
      real(kind=8) :: linear_independency
      
      integer :: info, lwork, iwork(8*minval(shape(mat))), ldmat, m, n
      integer, parameter :: ldu = 1, ldvt = 1

      real(kind=8) :: s(minval(shape(mat))), u(1,1), vt(1,1), work_size(1)
      real(kind=8), allocatable :: work(:)

      ldmat = size(mat, 1)
      m = size(mat,1)
      n = size(mat,2)

      !call dgesdd(jobz,m,n, a,   lda, s, u, ldu, vt, ldvt, work,lwork, iwork, info)
      call dgesdd("N", m, n, mat, ldmat, s, u, ldu, vt, ldvt, work_size, -1, iwork, info)
      
      lwork = int(work_size(1))
      allocate(work(lwork))

      call dgesdd("N", m, n, mat, ldmat, s, u, ldu, vt, ldvt, work, lwork, iwork, info)

      linear_independency = s(size(s))
   end function linear_independency


   subroutine print_mtx(mtx)
      implicit none 
      real(kind=8), intent(in) :: mtx(:,:)
      integer ::i

      do i = 1, size(mtx,1) 
         write (*,*) mtx(i,:)
      enddo
      write (*,*) "####"
   end subroutine print_mtx 
end module m_unify_zmat