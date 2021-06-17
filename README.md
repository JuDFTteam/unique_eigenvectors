# unique_eigenvectors
A set of eigenvectors is transformed into a unique representation.

# Example
The example sets up an hermitian matrix, that is only non-zero on the anti-diagonal and has degenerate eigenvalues. Then it's diagonalized using two different algorithms. 
1) using the upper triangle
2) using the lower triangle

Then the eigenvector matricies are compared before and after applying the unification library:
```
 Before:
 Eigenvalue diff norm:  0.000000000000000E+000
 Eigenvector diff norm:   20.7364413533277     
 Eigenvector maxdiff:   1.41421356237309     
 
 After:
 Eigenvector diff norm:  4.098685789482567E-008
 Eigenvector maxdiff:  9.685164559236625E-009
 degeneracy    eigenvalue             norm diff
           1  -10.0000000000000       0.000000000000000E+000
          41  -9.50000000000000       3.531725722560567E-008
          29  -9.00000000000000       5.802643958927017E-009
          13  -8.50000000000000       7.439288592879344E-010
           7  -8.00000000000000       1.411707138376848E-009
           8  -7.50000000000000       2.151595122918410E-009
           1  -7.00000000000000       0.000000000000000E+000
           1   7.00000000000000       0.000000000000000E+000
           8   7.50000000000000       2.137558628234475E-009
           7   8.00000000000000       8.705797509887812E-010
          13   8.50000000000000       5.164072213877280E-010
          29   9.00000000000000       3.572111105473027E-009
          41   9.50000000000000       1.932391669180075E-008
           1   10.0000000000000       0.000000000000000E+000

```