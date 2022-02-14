;; Thanks to Tobias's answer on Emacs Stack Exchange:
;; https://emacs.stackexchange.com/questions/38437/org-mode-batch-export-missing-syntax-highlighting


(package-initialize)
(add-to-list 'package-archives
             '("gnu" . "https://elpa.gnu.org/packages/"))
(add-to-list 'package-archives
             '("melpa-stable" . "https://stable.melpa.org/packages/"))
(add-to-list 'package-archives
             '("melpa" . "https://melpa.org/packages/"))
(setq package-archive-priorities '(("melpa-stable" . 100)
                                   ("melpa" . 50)
                                   ("gnu" . 10)))


(require 'font-lock)
(setq org-confirm-babel-evaluate nil)
(global-font-lock-mode t)
(setq org-src-fontify-natively t)

(org-babel-do-load-languages
 'org-babel-load-languages
 '(
   (emacs-lisp . t)
   (shell . t)
   (python . t)
   (fortran . t)
   (C . t)
   (org . t)
   (makefile . t)
   ))


; The following is required to compute the file names
(setq top_builddir (or (getenv "top_builddir") "."))
(setq srcdir (or (getenv "srcdir") "."))

(setq src (concat top_builddir "/src/"))
(setq tests (concat top_builddir "/tests/"))
(setq name (file-name-nondirectory (substring buffer-file-name 0 -4)))
(setq f  (concat src name "_f.F90"))
(setq fh_func (concat src name "_fh_func.F90"))
(setq fh_type (concat src name "_fh_type.F90"))
(setq c  (concat src name ".c"))
(setq h_func  (concat src name "_func.h"))
(setq h_type  (concat src name "_type.h"))
(setq h_private_type  (concat src name "_private_type.h"))
(setq h_private_func  (concat src name "_private_func.h"))
(setq c_test  (concat tests "test_" name ".c"))
(setq f_test  (concat tests "test_" name "_f.F90"))
(org-babel-lob-ingest (concat srcdir "/tools/lib.org"))

