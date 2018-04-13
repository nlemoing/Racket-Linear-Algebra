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
         solve-cramer
         characteristic
         rank
         REF
         RREF
         REFops
         RREFops)

(define (square-matrix? A)
  (or (empty? A)
      (= (length A) (length (first A)))))

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
    (cond [(not (square-matrix? A)) (error "Cannot take determinant of non-square matrix")]
          [(= (length A) 1) (1d-det A)]
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
  (cond [(not (square-matrix? A)) (error "Matrix must be square to compute inverse")]
        [else
         (local [(define det (determinant A))
                 (define (cofactor-matrix A)
                   (build-list (length A)
                               (lambda (row)
                                 (build-list (length A)
                                             (lambda (col)
                                               (cofactor A row col))))))]
           (cond [(= det 0) false]
                 [else (matrix-scmult (/ 1 det) (transpose (cofactor-matrix A)))]))]))

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

(define (replace-ith-row index row A)
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
    (replace-ith-row index1 row A)))

(define (apply-ERO A op)
  (cond [(symbol=? (first op) 'add) (ERO-add (second op) (third op) (fourth op) A)]
        [(symbol=? (first op) 'mult) (ERO-mult (second op) (third op) A)]
        [(symbol=? (first op) 'swap) (ERO-swap (second op) (third op) A)]))

(define (REF A)
  (row-reduce A 'REF))

(define (REFops A)
  (row-reduce A 'REFops))

(define (RREF A)
  (row-reduce A 'RREF))

(define (RREFops A)
  (row-reduce A 'RREFops))

(define (rank A)
  (foldr (lambda (first rest)
           (cond [(andmap (lambda (x) (= x 0)) first) rest]
                 [else (add1 rest)])) 0 (REF A)))

(define (row-reduce A op)
  (local [(define (row-echelon A RREF? ops?)
            (local [(define (reduce A REF EROs row RREF? ops?)
                      (cond [(or (empty? A)
                                 (empty? (first A)))
                             (cond [(and (not RREF?) (not ops?)) REF]
                                   [(not RREF?) (reverse EROs)]
                                   [else (reduced-row-echelon REF EROs ops?)])]
                            [(andmap (lambda (x) (= x 0))
                                     (map first A)) ; first column has all zeroes
                             (reduce (map rest A) REF EROs row RREF? ops?)]
                            [(= (first (first A)) 0)
                             (local [(define fnz (first-non-zero (map first A)))]
                               (reduce (ERO-swap 0 fnz A)
                                       (ERO-swap row (+ row fnz) REF)
                                       (cons (cons 'swap (cons row (cons (+ row fnz) empty)))
                                             EROs)
                                       row RREF? ops?))]
                            [(not (= (first (first A)) 1))
                             (reduce (ERO-mult (/ 1 (first (first A))) 0 A)
                                     (ERO-mult (/ 1 (first (first A))) row REF)
                                     (cons (cons 'mult (cons (/ 1 (first (first A))) (cons row empty)))
                                           EROs)
                                     row RREF? ops?)]
                            [(ormap (lambda (x) (not (= x 0)))
                                    (rest (map first A))) ; nonzero elements in the column
                             (local [(define fnz (+ 1 (first-non-zero (rest (map first A)))))
                                     (define fnzval (* -1 (foldr (lambda (first result)
                                                                   (cond [(not (= first 0)) first]
                                                                         [else result]))
                                                                 -1
                                                                 (rest (map first A)))))]
                               (reduce (ERO-add fnz fnzval 0 A)
                                       (ERO-add (+ fnz row) fnzval row REF)
                                       (cons (cons 'add (cons (+ fnz row) (cons fnzval (cons row empty))))
                                             EROs)
                                       row RREF? ops?))]
                            [else (reduce (rest (map rest A)) REF EROs (add1 row) RREF? ops?)]))]
              (reduce A A empty 0 RREF? ops?)))
          
          (define (reduced-row-echelon REF EROs ops?)
            (local [(define leading-ones-columns (filter (lambda (x) (not (= x -1)))
                                                         (map first-non-zero REF)))
                    (define (reduce REF EROs ops? ones-lst row)
                      (cond [(empty? ones-lst)
                             (cond [(false? ops?) REF]
                                   [else (reverse EROs)])]
                            [(= (length (filter (lambda (x) (not (= x 0))) (get-col REF (first ones-lst)))) 1)
                             (reduce REF EROs ops? (rest ones-lst) (add1 row))]
                            [else (local [(define col (get-col REF (first ones-lst)))
                                          (define fnz (first-non-zero col))
                                          (define fnzval (* -1 (foldr (lambda (first result)
                                                                        (cond [(not (= first 0)) first]
                                                                              [else result]))
                                                                      -1 col)))]
                                    (reduce (ERO-add fnz fnzval row REF)
                                            (cons (list 'add fnz fnzval row) EROs)
                                            ops? ones-lst row))]))]
              (reduce REF EROs ops? leading-ones-columns 0)))
          
          (define (first-non-zero lst)
            (local [(define (first-non-zero-index lst index)
                      (cond [(empty? lst) -1]
                            [(not (= (first lst) 0)) index]
                            [else (first-non-zero-index (rest lst) (add1 index))]))]
              (first-non-zero-index lst 0)))]
    (cond [(symbol=? op 'REF) (row-echelon A false false)]
          [(symbol=? op 'REFops) (row-echelon A false true)]
          [(symbol=? op 'RREF) (row-echelon A true false)]
          [(symbol=? op 'RREFops) (row-echelon A true true)]
          [else (error "row-reduce: invalid command")])))

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

(define (solve-cramer A b)
  (local [(define detA (determinant A))
          (define (replace-ith-col index A col)
            (cond [(empty? A) empty]
                  [(= index 0) (cons col (rest A))]
                  [else (cons (first A) (replace-ith-col (sub1 index) (rest A) col))]))]
    (cond [(= detA 0) (error "Matrix not invertible, cannot use cramer's rule")]
          [else 
           (vec-scmult (/ 1 detA) (build-list (length A)
                                              (lambda (index)
                                                (determinant (replace-ith-col index (transpose A) b)))))])))

(define (characteristic A)
  (lambda (x) (determinant (matrix-add A (matrix-scmult (* -1 x) (identity (length A)))))))