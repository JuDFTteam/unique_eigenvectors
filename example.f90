program main
   use  m_unify_zmat
   implicit none
   integer, parameter :: N= 1000, lwork = 3*N
   real(kind=8)       :: H(N,N), vecs1(N,N), vecs2(N,N), eig1(N), eig2(N), work(lwork)
   integer :: info, i
   
   ! Hermitian matrix with a few degenerate values
   H = 0.0
   do i = 1,(N/2)+1
      H(i,N-i+1) = ceiling(i/2.0)
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

   call unify_zmat(eig1, vecs1)
   call unify_zmat(eig2, vecs2)

   write (*,*) new_line("a") // "After:"
   !write (*,*) "Eigenvalues:", eig1
   write (*,*) "Eigenvector diff norm:", norm2(vecs1 - vecs2) 
   write (*,*) "Eigenvector maxdiff:", maxval(abs(vecs1 - vecs2))
end program main