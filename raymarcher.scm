;; Simple raymarcher

(define EPSILON 1e-4)
(define GLOBAL-MAXRAYLEN 1e2)
(define PI 3.1415926535897)

;; SICP Streams

(define-syntax cons-stream
  (syntax-rules ()
                ((x a b) (cons a (delay b)))))

(define stream-car car)

(define (stream-cdr stream)
  (force (cdr stream)))

(define the-empty-stream '())
(define stream-null? null?)

(define (stream-append . streams)
  (if (null? streams)
      the-empty-stream
      (if (null? (cdr streams))
          (car streams)
          (let ((first (car streams))
                (rest (cdr streams)))
            
            (if (stream-null? first)
                (apply stream-append rest)
                (cons-stream (stream-car first)
                             (apply stream-append (cons (stream-cdr first) rest))))))))

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
  (stream-for-each (lambda (el) (display el) (display " ")) s))

(define (stream-enumerate-interval low high)
  (if (< low high)
      (cons-stream low (stream-enumerate-interval (+ low 1) high))
      the-empty-stream))

(define (stream-zip s1 s2)
  (if (or (stream-null? s1) (stream-null? s2))
      the-empty-stream
      (cons-stream (list (stream-car s1) (stream-car s2))
                   (stream-zip (stream-cdr s1) (stream-cdr s2)))))

;; Linear algebra

(define unit+x '(+1 0 0))
(define unit-x '(-1 0 0))
(define unit+y '(0 +1 0))
(define unit-y '(0 -1 0))
(define unit+z '(0 0 +1))
(define unit-z '(0 0 -1))
(define vx car)
(define vy cadr)
(define vz caddr)
(define (vxy v) (list (vx v) (vy v)))
(define (vyz v) (list (vy v) (vz v)))
(define (vxz v) (list (vx v) (vz v)))

(define (vscale v s)
  (map (lambda (vi) (* s vi)) v))

(define (v+ . vs)
  (apply map (cons + vs)))

(define (v- . vs)
  (apply map (cons - vs)))

(define (vdot v1 v2)
  (apply + (map * v1 v2)))

(define (vcross v1 v2)
  (let ((v1x (vx v1)) (v1y (vy v1)) (v1z (vz v1))
        (v2x (vx v2)) (v2y (vy v2)) (v2z (vz v2)))
  (list
    (- (* v1y v2z) (* v1z v2y))
    (- (* v1z v2x) (* v1x v2z))
    (- (* v1x v2y) (* v1y v2x)))))

(define (vmag2 p)
  (vdot p p))

(define (vmag p)
  (sqrt (vmag2 p)))

(define (vabs v)
  (map abs v))

(define (vnormalize v)
  (vscale v (/ (vmag v))))

(define (mv* m v)
  (map (lambda (mrow) (vdot mrow v)) m))

(define (Rx theta)
  (let ((costheta (cos theta))
        (sintheta (sin theta)))
    `((1 0 0)
      (0 ,costheta ,(- sintheta))
      (0 ,sintheta ,costheta))))

(define (Ry theta)
  (let ((costheta (cos theta))
        (sintheta (sin theta)))
    `((,costheta 0 ,sintheta)
      (0 1 0)
      (,(- sintheta) 0 ,costheta))))

(define (Rz theta)
  (let ((costheta (cos theta))
        (sintheta (sin theta)))
    `((,costheta ,(- sintheta) 0)
      (,sintheta ,costheta 0)
      (0 0 1))))

;; Operations on signed distance functions

(define (((sdf-combine proc) . sdfs) p)
  (apply proc
         (map (lambda (sdf) (sdf p)) sdfs)))

(define sdf-union (sdf-combine min))
(define sdf-intersection (sdf-combine max))

(define ((sdf-complement sdf) p)
  (- (sdf p)))

(define (sdf-difference . sdfs)
  (if (< (length sdfs) 2)
    (error "sdf-difference needs at least two arguments")
    (sdf-intersection (car sdfs)
                      (sdf-complement (apply sdf-union (cdr sdfs))))))

(define ((sdf-translate sdf x) p)
  (sdf (v- p x)))

(define ((sdf-scale sdf f) p)
  (if (pair? f)
    (* (sdf (map * p (map / f))) (apply min f)) 
    (* (sdf (vscale p (/ f))) f)))

(define (((make-sdf-rotator R) sdf theta) p)
    (sdf (mv* (R (- theta)) p)))

(define sdf-rotate-x (make-sdf-rotator Rx))
(define sdf-rotate-y (make-sdf-rotator Ry))
(define sdf-rotate-z (make-sdf-rotator Rz))

(define (mod x y)
  (- x (* (truncate (/ x y)) y)))

(define ((sdf-tessellate-x sdf x1 x2) p)
  (let* ((dx (- x2 x1))
         (modpx (mod (- (vx p) x1) dx)))
    (sdf (list
          (+ (if (>= modpx 0)
                 modpx
                 (+ modpx dx))
             x1)
          (vy p)
          (vz p)))))

(define (sdf-normal sdf p)
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

;; Geometric primitives

(define (sphere p)
  (- (vmag p) 1.0))

(define solid-xzplane cadr) ; Solid plane with normal +y

(define ((torus minor-rad) p)
  (- (vmag (list (- (vmag (vxy p)) 1.0) (vz p))) minor-rad))

(define (cube p)
  (let ((d (map (lambda (el) (- (abs el) 1.0)) p)))
    (+ (min (apply max d) 0.0)
       (vmag (map
              (lambda (el) (max el 0.0))
              d)))))

;; Operations on textured SDFs

(define ((apply-texture sdf t) p)
  (cons t (sdf p)))

(define ((make-argcomp comp) proc lst)
  (let ((vals (map (lambda (el) (cons (proc el) el)) lst)))
    (cdr (foldl (lambda (prev val)
                  (if (or (null? prev) (comp (car val) (car prev)))
                      val
                      prev))
                '() vals))))

(define argmin (make-argcomp <))
(define argmax (make-argcomp >))

(define ((tsdf-union . tsdfs) p)
  (argmin cdr (map (lambda (tsdf) (tsdf p)) tsdfs)))

(define tsdf-translate sdf-translate)

;; Scene description

(define (make-scene tsdf lights camera)
  (list tsdf lights camera))

(define (scene-tsdf scene)
  (car scene))

(define (scene-lights scene)
  (cadr scene))

(define (scene-camera scene)
  (caddr scene))

;; Cameras

(define (((make-camera location up right) w h) pixel)
  (let* ((aspect-ratio (/ w h))
         (dir-no-norm (list
                       (* (- (/ (car pixel) w) 0.5) aspect-ratio)
                       (- 0.5 (/ (cdr pixel) h))
                       1)))
    (vnormalize dir-no-norm)))

;; Lights

(define (make-light pos col)
  (cons pos col))

(define (light-pos light)
  (car light))

(define (light-col light)
  (cdr light))

;; Ray marching

(define (march ray-dir tsdf p maxRayLen)
  (if (< maxRayLen 0)
    '(0 0 0)
    (let* ((tdist (tsdf p))
           (dist (cdr tdist)))
      (if (< dist EPSILON)
          (let ((tfunc (car tdist))
                (pprime (v- p (vscale ray-dir EPSILON))))
            (tfunc tsdf pprime))
          (march ray-dir tsdf (v+ p (vscale ray-dir dist)) (- maxRayLen dist))))))

;; Shading

(define (color-clip color)
  (map (lambda (c)
         (inexact->exact (round (min 255 (max 0 c)))))
         color))

;; PPM output

(define (make-h-pixel-stream y w)
  (stream-map (lambda (x) (cons x y)) (stream-enumerate-interval 0 w)))

(define (make-pixel-stream-from-row y w h)
  (if (< y h)
    (stream-append (make-h-pixel-stream y w) (make-pixel-stream-from-row (+ y 1) w h))
    the-empty-stream))

(define (make-pixel-stream w h)
  (make-pixel-stream-from-row 0 w h))

(define (make-bitmap scene w h)
  (let* ((pixel-stream (make-pixel-stream w h))
         (dir-stream (stream-map ((scene-camera scene) w h) pixel-stream))
         (ray-colorer (lambda (dir) (march dir (scene-tsdf scene) '(0 0 0)
                                                GLOBAL-MAXRAYLEN)))
         (pixel-color-stream (stream-map ray-colorer dir-stream)))
    (list w h (stream-zip pixel-stream pixel-color-stream))))

(define (bitmap-width bitmap)
  (car bitmap))

(define (bitmap-height bitmap)
  (cadr bitmap))

(define (bitmap-cpixel-stream bitmap)
  (caddr bitmap))

(define (cpixel-color cp)
  (cadr cp))

(define (cpixel-r cp)
  (car (cadr cp)))

(define (cpixel-g cp)
  (cadr (cadr cp)))

(define (cpixel-b cp)
  (caddr (cadr cp)))

(define (cpixel-x cp)
  (car (car cp)))
(define (cpixel-y cp)
  (cdr (car cp)))

(define (renderPPM b filename)
  (call-with-output-file filename
    (lambda (port)
      (display "P3\n" port)
      (display (bitmap-width b) port)
      (display " " port)
      (display (bitmap-height b) port)
      (newline port)
      (display "255\n" port) 

      (stream-for-each (lambda (cp)
                         (display (cpixel-r cp) port)
                         (display " " port)
                         (display (cpixel-g cp) port)
                         (display " " port)
                         (display (cpixel-b cp) port)
                         (if (= (- (bitmap-width b) 1) (cpixel-x cp))
                             (newline port)
                             (display " " port)))
                (bitmap-cpixel-stream b)))))

(define (write-byte int . opt)
  (apply write-char (integer->char int) opt))

(define (write-bytes bytes . opt)
  (if (not (null? bytes))
      (begin
        (apply write-byte (car bytes) opt)
        (apply write-bytes (cdr bytes) opt))))

;; Colours

(define colour-white '(255 255 255))
(define colour-black '(0 0 0))
(define colour-red '(255 0 0))
(define colour-green '(0 255 0))
(define colour-blue '(0 0 255))

;;;; TEST CODE ;;;;

(define scene (make-scene

               (tsdf-translate
                (tsdf-union
                 (apply-texture sphere (lambda (p tsdf) colour-red))
                 (apply-texture (sdf-translate cube '(-2 0 0)) (lambda (p tsdf) colour-green))
                 (apply-texture (sdf-translate cube '(+2 0 0)) (lambda (p tsdf) colour-blue)))
                '(0 0 5))

               (list (make-light '(0 0 0) colour-white))

               (make-camera '(0 0 0) unit+y unit+x)))
               

(define bitmap (make-bitmap scene 640 480))
(renderPPM bitmap "test.ppm")
