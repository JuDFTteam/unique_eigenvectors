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
 Eigenvector diff norm:   17.7763888346312     
 Eigenvector maxdiff:   1.41421356237309     
 
 After:
 Eigenvector diff norm:  3.943535434920788E-008
 Eigenvector maxdiff:  7.662817068376793E-009
 degeneracy    eigenvalue             norm diff
           7  -8.50000000000000       7.623826405426115E-009
          38  -8.00000000000000       2.135776810274473E-008
          14  -7.50000000000000       1.137447650433665E-009
           9  -7.00000000000000       7.222064936932765E-010
           6  -6.50000000000000       2.689916609346051E-010
           1  -6.00000000000000       0.000000000000000E+000
           1   6.00000000000000       0.000000000000000E+000
           6   6.50000000000000       2.027778771164724E-011
           9   7.00000000000000       3.557640927824789E-009
          14   7.50000000000000       1.080178747779877E-009
          38   8.00000000000000       3.200023243658962E-008
           7   8.50000000000000       1.069381635950538E-009
           1   10.0000000000000       0.000000000000000E+000
```

At the end of the example the eigenvectors for each degenerate group are compared. As you can see there is still a high numerical stability for degeneracies smaller around 5, and still performs well upto degeneracies of several 100.