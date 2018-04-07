# Racket-Linear-Algebra
A linear algebra module for racket (written before I realized Racket implements one in its math library).
I implemented this while planning the review session for Math 136 and used it as a quick way to check my work.
Vectors should be lists and matrices should be lists of lists.
Minimal error checking is implemented since this was mostly a quick hack.
Don't try to find the determinant of a non-square matrix or do invalid multiplication.

Provided functions:
get-elem: A r c -> Num
  returns the element of A in the rth row and cth column

get-row: A r -> Vec
  returns the rth row of A

get-col: A c -> Vec
  returns the cth column of A

vec-add: v1 v2 -> Vec
  returns the sum of v1 and v2

vec-scmult: c v -> Vec
  returns v multiplied by c
  
dot: v1 v2 -> Num
  returns the dot product of v1 and v2

vec-length: v -> Num
  returns the length of v

cross: v1 v2 -> v
  returns the cross product of v1 and v2
  
proj: v1 v2 -> v
  returns the projection of v1 onto v2

perp: v1 v2 ->
  returns the perpendicular of v1 to v2

matrix-add: A B -> Matrix
  returns the sum of A and B

matrix-scmult c A -> Matrix
  returns A multiplied by c

vec-matrix-mult: A x -> Vec
  returns A multiplied by x
  requires A has the same # of columns as x has entries

matrix-mult: A B -> Matrix
  returns A multiplied by B
  requires A has the same # of columns as B has rows
  
tranpose: A -> Matrix
  returns the transpose of A
  note: A is a list of rows - this will transform A to a list of columns

cofactor: A i j -> Num
  returns the cofactor of A at row i and column j
  requires A is square

determinant: A -> Num
  returns the determinant of A
  requires A is square

identity: dim -> Matrix
  returns the identity matrix with dimension dim

inverse: A -> Matrix
  computes the inverse of A by cofactors

apply-ERO: A ERO -> Matrix
  applies an ERO to A
  an ERO can be:
    '(add index1 c index2): add c times the row at index2 to index1
    '(mult c index): multiply the row at index by c
    '(swap index1 index2): swap the rows at index1 and index

apply-EROs: A (listof ERO) -> Matrix
  applies a list of EROs (in order) to A

elem-matrix: ERO dim -> Matrix
  produces an elementary matrix with dimension dim corresponding to the ERO

solve: A b -> Vec
  solves the system of equations Ax = b using Cramer's rule
  requires A is invertible

characteristic: A -> Matrix
  returns the characteristic polynomial of A
  requires A is square
