;;; benchmarks/siam100_cap.lisp -- per-case time cap for the Maxima halves of
;;; the SIAM 100-digit experiments.  Maxima has no TimeConstrained; Maxima on
;;; SBCL (the build every runner here uses) has sb-ext:with-timeout, which
;;; interrupts the computation and lets the file continue with the next case.
;;; Loaded by siam100_common.mac.  Lisp rather than Maxima because `:lisp`
;;; lines are not accepted inside a batchload'ed file.
(defun $capcall (secs f)
  (handler-case (sb-ext:with-timeout secs (mfuncall f))
    (sb-ext:timeout () '$timeout)
    (error () '$failed)))

;;; Decimal digits a bigfloat carries: its own mantissa width (bits), which is
;;; the fpprec in force when it was made -- not the global fpprec at the time
;;; of the call, which a block-local fpprec has already restored.
(defun $bfdigits (x)
  (max 1 (floor (* (third (car x)) 0.30103))))
