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

(require 'cl)
(let* ((required-packages
        '(htmlize
          evil
          org-evil
          org-bullets
          ))
       (missing-packages (remove-if #'package-installed-p required-packages)))
  (when missing-packages
    (message "Missing packages: %s" missing-packages)
    (package-refresh-contents)
    (dolist (pkg missing-packages)
      (package-install pkg)
      (message "Package %s has been installed" pkg))))

(setq backup-directory-alist
      `(("." . ,(concat user-emacs-directory "backups"))))
(setq backup-by-copying t)

(require 'org)
(setq org-format-latex-options (plist-put org-format-latex-options :scale 1.6))

(setq org-hide-leading-stars t)
(setq org-alphabetical-lists t)
(setq org-src-fontify-natively t)
(setq org-src-tab-acts-natively t)
(setq org-src-preserve-indentation t)
(setq org-hide-emphasis-markers nil)
(setq org-pretty-entities nil)
(setq org-confirm-babel-evaluate nil) ;; Do not ask for confirmation all the time!!
(setq python-indent-guess-indent-offset-verbose nil) ;; Remove warning : Can’t guess python-indent-offset 

(org-babel-do-load-languages
 'org-babel-load-languages
 '(
   (emacs-lisp . t)
   (shell . t)
   (python . t)
   (C . t)
   (org . t)
   (makefile . t)
   ))

;; Use python3 instead of python2.7 
(setq org-babel-python-command "python3")

(add-hook 'org-babel-after-execute-hook 'org-display-inline-images)
'(indent-tabs-mode nil)

(require 'evil)
(setq evil-want-C-i-jump nil)
(evil-mode 1)
(global-font-lock-mode t)
(global-superword-mode 1)

(setq line-number-mode 1)
(setq column-number-mode 1)

(evil-select-search-module 'evil-search-module 'evil-search)

(global-set-key (kbd "C-+") 'text-scale-increase)
(global-set-key (kbd "C--") 'text-scale-decrease)


(custom-set-variables
 ;; custom-set-variables was added by Custom.
 ;; If you edit it by hand, you could mess it up, so be careful.
 ;; Your init file should contain only one such instance.
 ;; If there is more than one, they won't work right.
 '(ansi-color-faces-vector
   [default default default italic underline success warning error])
 '(custom-enabled-themes (quote (leuven)))
)
