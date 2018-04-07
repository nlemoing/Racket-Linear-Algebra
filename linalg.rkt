#lang racket

(provide determinant
         cofactor
         vec-add
         vec-scmult
         vec-matrix-mult
         matrix-add
         matrix-scmult
         matrix-mult
         dot
         vec-length
         transpose
         get-row
         get-col
         get-elem
         inverse
         proj
         perp
         cross
         apply-ERO
         apply-EROs
         elem-matrix
         identity
         solve
         characteristic)

(define (determinant A)
  (local [(define (cofactor-expansion A row col)
            (cond [(empty? row) 0]
                  [else (+ (* (first row)
                              (cond [(= (first row) 0) 0]
                                    [else (cofactor A 0 col)]))
                           (cofactor-expansion A (rest row) (add1 col)))]))
          (define (2d-det A)
            (- (* (first (first A)) (second (second A)))
               (* (first (second A)) (second (first A)))))
          (define (1d-det A)
            (first (first A)))]
    (cond [(= (length A) 1) (1d-det A)]
          [(= (length A) 2) (2d-det A)]
          [else (cofactor-expansion A (first A) 0)])))

(define (cofactor A r c)
  (local [(define (remove-ith-element i lst)
            (cond [(empty? lst) empty]
                  [(= i 0) (rest lst)]
                  [else (cons (first lst)
                              (remove-ith-element (sub1 i) (rest lst)))]))
          (define Aij (remove-ith-element r (map (lambda (x) (remove-ith-element c x)) A)))]
    (* (expt -1 (+ r c))
       (determinant Aij))))           

(define (vec-add x y)
  (foldr (lambda (a b result)
           (cons (+ a b) result))
         empty x y))

(define (vec-scmult c x)
  (map (lambda (y) (* c y)) x))

(define (matrix-add A B)
  (foldr (lambda (a b result)
           (cons (vec-add a b) result))
         empty A B))

(define (matrix-scmult c A)
  (map (lambda (y) (vec-scmult c y)) A))

(define (dot x y)
  (foldr (lambda (a b result)
           (+ (* a b) result))
         0 x y))

(define (vec-length-sqr x)
  (dot x x))

(define (vec-length x)
  (sqrt (vec-length-sqr x)))

(define (vec-matrix-mult A x)
  (foldr (lambda (first rest)
           (cons (dot first x) rest))
         empty A))

(define (matrix-mult A B)
  (transpose (map (lambda (col) (vec-matrix-mult A col)) (transpose B))))

(define (transpose A)
  (local [(define (get-ith-col index A)
            (map (lambda (row) (get-ith-elem index row)) A))]
    (build-list (length (first A)) (lambda (x) (get-ith-col x A)))))

(define (get-ith-elem index row)
  (cond [(empty? row) empty]
        [(= index 0) (first row)]
        [else (get-ith-elem (sub1 index) (rest row))]))

(define (get-row A row)
  (get-ith-elem row A))

(define (get-col A col)
  (get-ith-elem col (transpose A)))

(define (get-elem A row col)
  (get-ith-elem col (get-ith-elem row A)))

;computes inverse by cofactors - will break if A is not n by n
(define (inverse A)
  (local [(define det (determinant A))
          (define (cofactor-matrix A)
            (build-list (length A)
                        (lambda (row)
                          (build-list (length A)
                                      (lambda (col)
                                        (cofactor A row col))))))]
    (cond [(= det 0) false]
          [else (matrix-scmult (/ 1 det) (transpose (cofactor-matrix A)))])))

(define (proj u v)
  (vec-scmult (/ (dot u v)
                 (vec-length-sqr v))
              v))

(define (perp u v)
  (vec-add u (vec-scmult -1 (proj u v))))


(define (cross v w)
  (list (- (* (second v) (third w))
           (* (second w) (third v)))
        (- (* (first w) (third v))
           (* (first v) (third w)))
        (- (* (first v) (second w))
           (* (first w) (second v)))))

(define (apply-ERO A op)
  (local [(define (replace-ith-row index row A)
            (cond [(empty? A) empty]
                  [(= index 0) (cons row (rest A))]
                  [else (cons (first A) (replace-ith-row (sub1 index) row (rest A)))]))
          (define (ERO-mult c index A)
            (local [(define row (vec-scmult c (get-ith-elem index A)))]
              (replace-ith-row index row A)))
          (define (ERO-swap index1 index2 A)
            (local [(define row1 (get-ith-elem index1 A))
                    (define row2 (get-ith-elem index2 A))]
              (replace-ith-row index2 row1 (replace-ith-row index1 row2 A))))
          (define (ERO-add index1 c index2 A)
            (local [(define row1 (get-ith-elem index1 A))
                    (define row2 (get-ith-elem index2 A))
                    (define row (vec-add row1 (vec-scmult c row2)))]
              (replace-ith-row index1 row A)))]
    (cond [(symbol=? (first op) 'add) (ERO-add (second op) (third op) (fourth op) A)]
          [(symbol=? (first op) 'mult) (ERO-mult (second op) (third op) A)]
          [(symbol=? (first op) 'swap) (ERO-swap (second op) (third op) A)])))

(define (apply-EROs A ops)
  (cond [(empty? ops) A]
        [else (apply-EROs (apply-ERO A (first ops)) (rest ops))]))

(define (elem-matrix op dim)
  (apply-ERO (identity dim) op))

(define (identity dim)
  (build-list dim (lambda (x)
                    (build-list dim (lambda (y)
                                      (cond [(= x y) 1]
                                            [else 0]))))))

(define (solve A b)
  (local [(define (solve-inverse A b)
            (vec-matrix-mult (inverse A) b))
          (define (solve-cramer A b)
            (local [(define detA (determinant A))
                    (define (replace-ith-col index A col)
                      (cond [(empty? A) empty]
                            [(= index 0) (cons col (rest A))]
                            [else (cons (first A) (replace-ith-col (sub1 index) (rest A) col))]))]
              (vec-scmult (/ 1 detA) (build-list (length A)
                                                 (lambda (index)
                                                   (determinant (replace-ith-col index (transpose A) b)))))))]
    (solve-cramer A b)))

;we are taking inverse by cofactors
;requires n^2 determinants (although they are (n-1) * (n-1)
;cramer's rule requires n determinants of size n * n
;recall that determinants are recursive and take longer

(define (characteristic A)
  (lambda (x) (determinant (matrix-add A (matrix-scmult (* -1 x) (identity (length A)))))))



; finds the characteristic polynomial of A - returns a 1 argument function

;;Q1

(define I3 (identity 3))

(define op-lst `((swap 0 1)
                 (add 1 -2 0)
                 (add 2 2 0)
                 (add 0 -1 1)
                 (add 2 -2 1)
                 (mult ,(/ -1 6) 2)
                 (add 0 4 2)
                 (add 1 -3 2)))

(define row-ops (map (lambda (x) (elem-matrix x 3)) op-lst))

(define A-1 (foldl (lambda (x y) (matrix-mult x y))
                   I3
                   row-ops))

(define A (foldr (lambda (x y) (matrix-mult (inverse x) y))
                 I3
                 row-ops))
(define row-A (apply-EROs A op-lst))

(define b '(-1 1 0))
(define L '((3 4) (-2 -3)))
(define CL (characteristic L))
(define eigen1 '(2 -1))
(define eigen2 '(-1 1))
(define P (transpose (list eigen1 eigen2)))
(define P-1 (inverse P))

(define v1 '(1 0 -2))
(define v2 '(-1 1 1))
(define p1 (proj v1 v2))
(define p2 (perp v1 v2))
(define n (cross v1 v2))
(define c '(1 1 1))