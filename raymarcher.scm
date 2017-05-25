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

;; Linear algebra

(define (vscale v s)
  (map (lambda (vi) (* s vi)) v))

(define (v+ . vs)
  (apply map (cons + vs)))

;; CSG 

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

;; Transformation

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

(define (make-h-pixel-stream-from x y w)
  (if (< x w)
    (cons-stream (cons x y) (make-h-pixel-stream-from (+ x 1) y w))
    the-empty-stream))

(define (make-h-pixel-stream y w)
  (make-h-pixel-stream-from 0 y w))

(define (make-pixel-stream-from-row y w h)
  (if (< y h)
    (stream-append (make-h-pixel-stream y w) (make-pixel-stream-from-row (+ y 1) w h))
    the-empty-stream))

(define (make-pixel-stream w h)
  (make-pixel-stream-from-row 0 w h))

(define (pixel-ray-dir w h pixel)
  (let ((dir-no-norm (list
                       (- (/ (car pixel) w) 0.5)
                       (- (/ (cdr pixel) h) 0.5)
                       1)))
    (vscale dir-no-norm (/ (r dir-no-norm)))))

(define ((pixel->intersection w h sdf) pixel)
  (march (pixel-ray-dir w h pixel) sdf '(0 0 0)))


;; PPM output

(define (renderPPM sdf w h)
  (let ((pixel-stream (make-pixel-stream w h)))

    ; Header
    (display "P3\n")
    (display w) (display " ") (display h) (newline)
    (display "255\n") 

    (stream-for-each (lambda (pixel)
                  (let ((isect ((pixel->intersection w h sdf) pixel)))
                    (if (null? isect)
                      (display "0 0 0")
                      (display "0 255 0"))
                    (if (= (- w 1) (car pixel))
                      (newline)
                      (display " "))))
                pixel-stream)))


;;;; TEST CODE ;;;;

(define s1 (sdf-translate sphere '(0 0 -0.5)))
(define s2 (sdf-translate sphere '(0 0 +0.5)))
(define lens (sdf-intersection s1 s2))
(define holy-lens (sdf-difference lens (sdf-scale sphere 0.5)))
(define scene (sdf-translate sphere '(0 0 -3)))

(define pixel-stream (make-pixel-stream 32 20))

(define intersection-stream
  (stream-map (pixel->intersection 32 20 scene)
              pixel-stream))

;(display-stream intersection-stream)

(with-output-to-file "test.ppm"
                     (lambda () (renderPPM scene 640 480)))
(exit)
