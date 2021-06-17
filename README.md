# unique_eigenvectors
A set of eigenvectors is transformed into a unique representation.

# Example
The example sets up an hermitian matrix, that is only non-zero on the anti-diagonal and has degenerate eigenvalues. Then it's diagonalized using two different algorithms. 
1) using the upper triangle
2) using the lower triangle

Then the eigenvector matricies are compared before and after applying the unification library:
```
bash-4.4$ ifort  -mkl unify_zmat.f90  example.f90 && ./a.out 
 Before:
 Eigenvalue diff norm:  0.000000000000000E+000
 Eigenvector diff norm:   46.3896540189728     
 Eigenvector maxdiff:   1.41421356237309     
 
 After:
 Eigenvector diff norm:  1.387677879323193E-011
 Eigenvector maxdiff:  9.701184300325849E-012
```