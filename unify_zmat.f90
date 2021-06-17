module m_unify_zmat
contains
   subroutine unify_zmat(eigval, z)
      implicit none
      real(kind=8), intent(in) :: eigval(:)
      real(kind=8)          :: z(:,:)

      integer :: beg_group, end_group, n_g
      integer, allocatable :: groups(:)

      groups = make_groups(eigval)
   
      do n_g = 1, size(groups)
         beg_group = sum(groups(1:n_g-1)) + 1
         end_group = sum(groups(1:n_g))

         if(end_group <= size(z,2)) call unify_group(beg_group, end_group, z)
      enddo
   end subroutine unify_zmat

   subroutine unify_group(beg_group, end_group, z)
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
      call mod_gram_schmidt(new_basis)

      z(:,beg_group:end_group) = new_basis
   end subroutine unify_group

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
end module m_unify_zmat