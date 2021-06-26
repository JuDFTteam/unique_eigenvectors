# unique_eigenvectors
A set of eigenvectors is transformed into a unique representation.

# Example
The example sets up an hermitian matrix, that is only non-zero on the anti-diagonal and has degenerate eigenvalues. Then it's diagonalized using two different algorithms. 
1) using the upper triangle
2) using the lower triangle

Then the eigenvector matricies are compared before and after applying the unification routine:
```
 Before:
 Eigenvalue diff norm:   1.7410336091849973E-014
 Eigenvector diff norm:   11.995241560667875     
 Eigenvector maxdiff:   1.9825515426502025     
 
 After:
 Eigenvector diff norm:   4.8937389103031451E-014
 Eigenvector maxdiff:   8.4542058639869201E-015
 degeneracy    eigenvalue             norm diff
           3  0.99999999999999967        3.9639475134322329E-016
           5   1.9999999999999996        1.9358557145744058E-015
           7   2.9999999999999969        5.2559504928307576E-015
           9   3.9999999999999956        8.3543082912107305E-015
          11   4.9999999999999973        1.4658328533276856E-014
          13   5.9999999999999920        4.0361414768320680E-014
          15   6.9999999999999929        1.7600242695908689E-014
           7   7.9999999999999991        1.1826451835715324E-014

```

This algorithm is more stable than previous version, but it depends on the 