;; Simple raymarcher

;; SICP Streams

(define-syntax cons-stream
  (syntax-rules ()
                ((x a b) (cons a (delay b)))))

(define stream-car car)

(define (stream-cdr stream)
  (force (cdr stream)))

(define the-empty-stream '())
(define stream-null? null?)

(define (stream-append s1 s2)
  (if (stream-null? s1)
    s2
    (cons-stream (stream-car s1)
                 (stream-append (stream-cdr s1) s2))))

(define (stream-ref s n)
  (if (= n 0)
    (stream-car s)
    (stream-ref (stream-cdr s) (- n 1))))

(define (stream-map proc s)
  (if (stream-null? s)
    the-empty-stream
    (cons-stream 
      (proc (stream-car s))
      (stream-map proc (stream-cdr s)))))

(define (stream-for-each proc s)
  (if (stream-null? s)
    'done
    (begin 
      (proc (stream-car s))
      (stream-for-each proc 
                       (stream-cdr s)))))

(define (display-stream s)
  (stream-for-each (lambda (el) (display el) (newline)) s))

(define (stream-enumerate-interval low high)
  (if (< low high)
      (cons-stream low (stream-enumerate-interval (+ low 1) high))
      the-empty-stream))

(define (stream-zip s1 s2)
  (if (or (stream-null? s1) (stream-null? s2))
      the-empty-stream
      (stream-cons (cons (stream-car s1) (stream-car s2))
                   (stream-zip (stream-cdr s1) (stream-cdr s2)))))

;; Linear algebra

(define (vscale v s)
  (map (lambda (vi) (* s vi)) v))

(define (v+ . vs)
  (apply map (cons + vs)))

;; CSG with signed distance functions

(define (((sdf-combine proc) . sdfs) p)
  (apply proc
         (map
           (lambda (sdf) (sdf p))
           sdfs)))

(define sdf-union (sdf-combine min))
(define sdf-intersection (sdf-combine max))


(define ((sdf-complement sdf) p)
  (- (sdf p)))

(define ((sdf-difference . sdfs) p)
  (if (< (length sdfs) 1)
    (error "sdf-difference needs at least two arguments")
    (apply sdf-intersection
           (cons (car sdfs)
                 (map sdf-complement (cdr sdfs))))))

;; SDF transformations

(define ((sdf-translate sdf x) p)
  (sdf (v+ p x)))

(define ((sdf-scale sdf f) p)
  (* (sdf (vscale p (/ f))) f))


;; Geometric primitives

(define (r2 p)
  (apply + (map * p p)))

(define (r p)
  (sqrt (r2 p)))

(define (sphere p)
  (- (r p) 1.0))

(define (xyplane p) (abs (caddr p)))

(define solid-xyplane caddr) ; Solid plane with normal +z

;; Ray marching

(define MAXDIST2 1e4)
(define EPSILON 1e-4)

(define (march ray-dir sdf p)
  (if (> (r2 p) MAXDIST2)
    '()
    (let ((dist (sdf p)))
      (if (< dist EPSILON)
        p
        (march ray-dir sdf (v+ p (vscale ray-dir dist)))))))

;; Texturing

(define unit+x '(+1 0 0))
(define unit-x '(-1 0 0))
(define unit+y '(0 +1 0))
(define unit-y '(0 -1 0))
(define unit+z '(0 0 +1))
(define unit-z '(0 0 -1))

(define ((sdf-normal sdf) p)
  (let* ((eps EPSILON)
        (eps/2 (/ EPSILON 2)))
  (list
   (/ (- (sdf (v+ p (vscale unit+x eps/2)))
         (sdf (v+ p (vscale unit-x eps/2))))
      eps)
   (/ (- (sdf (v+ p (vscale unit+y eps/2)))
         (sdf (v+ p (vscale unit-y eps/2))))
      eps)
   (/ (- (sdf (v+ p (vscale unit+z eps/2)))
         (sdf (v+ p (vscale unit-z eps/2))))
      eps))))

;; PPM output

(define (make-h-pixel-stream y w)
  (stream-map (lambda (x) (cons x y)) (stream-enumerate-interval 0 w)))

(define (make-pixel-stream-from-row y w h)
  (if (< y h)
    (stream-append (make-h-pixel-stream y w) (make-pixel-stream-from-row (+ y 1) w h))
    the-empty-stream))

(define (make-pixel-stream w h)
  (make-pixel-stream-from-row 0 w h))


(define (pixel-ray-dir w h pixel)
  (let* ((aspect-ratio (/ w h))
         (dir-no-norm (list
                       (* (- (/ (car pixel) w) 0.5) aspect-ratio)
                       (- (/ (cdr pixel) h) 0.5)
                       1)))
    (vscale dir-no-norm (/ (r dir-no-norm)))))

(define (pixel->intersection sdf w h pixel)
  (march (pixel-ray-dir w h pixel) sdf '(0 0 0)))


(define (renderPPM sdf w h filename)
  (call-with-output-file filename
    (lambda (port)
      (let ((pixel-stream (make-pixel-stream w h)))

    ; Header
    (display "P3\n" port)
    (display w port) (display " " port) (display h port) (newline port)
    (display "255\n" port) 

    (stream-for-each (lambda (pixel)
                  (let ((isect (pixel->intersection sdf w h pixel)))
                    (if (null? isect)
                      (display "0 0 0" port)
                      (display "0 255 0" port))
                    (if (= (- w 1) (car pixel))
                      (newline port)
                      (display " " port))))
                pixel-stream)))))


;;;; TEST CODE ;;;;

(define s1 (sdf-translate sphere '(0 0 -0.5)))
(define s2 (sdf-translate sphere '(0 0 +0.5)))
(define lens (sdf-intersection s1 s2))
(define holy-lens (sdf-difference lens (sdf-scale sphere 0.5)))
(define scene (sdf-translate sphere '(0 0 -3)))

(renderPPM scene 640 480 "test.ppm")
(exit)
