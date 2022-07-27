(define-module (qmckl)
  #:use-module (guix packages)
  #:use-module (gnu packages autotools)
  #:use-module (gnu packages pkg-config)
  #:use-module (gnu packages gcc)
  #:use-module (gnu packages maths)  	;; contains blas, lapack, hdf5
  #:use-module (trexio) 		;; custom trexio module, has to be provided via -L command line option
  #:use-module (guix download)
  #:use-module (guix build-system gnu)
  #:use-module (guix licenses)
  ;; the modules below are additional requirements for building the dev version
  #:use-module (guix git-download)
  #:use-module (gnu packages version-control) ;; this is ONLY needed for cloning the htmlize.el, better to avoid it
  #:use-module (gnu packages python)
  #:use-module (gnu packages emacs)
  #:use-module (gnu packages m4))


(define-public qmckl-hpc-0.2.1
  (package
    (name "qmckl-hpc")
    (version "0.2.1")
    (source (origin
              (method url-fetch)
              (uri (string-append "https://github.com/TREX-CoE/qmckl/releases/download/v" version
				  "/qmckl-" version 
                                  ".tar.gz"))
              (sha256
               (base32
		;; the hash below is produced by guix download <url>
		"18100fd4vp41saxiji734mq5lckjplbnmm1nz29da4azhxzbzki9"
                ))))
    (build-system gnu-build-system)
    (arguments 
      '(#:configure-flags 
	'("--enable-hpc" 
	  "--with-openmp")))
    (inputs 
      `(("trexio", trexio) 
	("gfortran", gfortran)
	("openblas", openblas)
	("lapack", lapack)
	))
    (synopsis "QMCkl: Quantum Monte Carlo Kernel Library")
    (description "The main QMC algorithms are exposed in a simple language and provided a standard API 
		 and tests to enable the development of high-performance QMCkl implementations taking 
		 advantage of modern hardware.")
    (home-page "https://trex-coe.github.io/qmckl/index.html")
    (license bsd-3)))


;; build and installation of this still fails due to some bugs in the tools/build_makefile.py script
;; which fails to parse the git tree
;; Guix does not copy .git folder into the source directory during unpack phase !
(define-public qmckl-master
  (let ((commit "b33d6e523da8528a824bac8bc46a62a79d0dd846")
	(revision "1"))
  (package
    (name "qmckl-master")
    (version (git-version "0.2.3" revision commit))
    (source (origin
              (method git-fetch)
              (uri (git-reference
		     (url "https://github.com/TREX-CoE/qmckl")
		     (commit commit)
		     (recursive? #t))) ;; tried recursive to get htmlize - fails !
              (file-name (git-file-name name version))
              (sha256
               (base32
		;; the hash below is produced by `guix hash -rx .`
		"1a46s14vki0mhx3g5xclhqllkfw3gc0y583a7la1s7ixm5f7zvim"
                ))))
    (build-system gnu-build-system)
    (arguments 
      '(#:configure-flags '("--with-openmp")
	;; ignoring make errors while building htmlized docs is important for the build to pass
	#:make-flags '("-i")
	;;))
	#:phases
	;; this is a workaround to activate QMCKL_DEVEL mode
	(modify-phases %standard-phases
	  (add-after 'unpack 'set_devel 
	    (lambda _ 
	      (mkdir-p ".git")))
	)))
    (inputs
      `(("trexio", trexio) 
	("gfortran", gfortran)
	("openblas", openblas)
	("lapack", lapack)
	("python", python-wrapper)
	("emacs", emacs)
	))
      ;; these inputs are required for autogen.sh to work properly
    (native-inputs
      `(("autoconf", autoconf) 
	("automake", automake)
	("libtool", libtool)
	("pkg-config", pkg-config)
	("git", git) ;; this is ONLY needed for cloning the htmlize.el, better to avoid it
	("m4", m4)
	))
    (synopsis "QMCkl: Quantum Monte Carlo Kernel Library")
    (description "The main QMC algorithms are exposed in a simple language and provided a standard API 
		 and tests to enable the development of high-performance QMCkl implementations taking 
		 advantage of modern hardware.")
    (home-page "https://trex-coe.github.io/qmckl/index.html")
    (license bsd-3))))


(define-public qmckl
  ;; Default version of QMCkl - change this to benchmark installation from Git
  qmckl-hpc-0.2.1)


;; return qmckl variable so that `quix package -f qmckl.scm` works
qmckl
