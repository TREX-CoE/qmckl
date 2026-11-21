#!/bin/bash
# Export QMCkl org-mode docs to HTML.
#
# No htmlize: source blocks are exported as plain <pre class="src src-LANG">
# and highlighted later by build_site.py (Pygments, build time only).
# The "Cannot fontify source block" warnings on stderr are expected and benign.
set -euo pipefail

REPO="${1:?usage: export_docs.sh REPO_DIR}"
cd "$REPO"

count=0
while IFS= read -r f; do
  [[ -z "$f" ]] && continue
  base="${f%.org}"
  echo "export $f"
  # same load order as the repo's tools/build_doc.sh, minus htmlize
  # (config_tangle.el must be loaded AFTER the file: it names tangle
  #  targets from buffer-file-name and ingests the noweb labels from
  #  tools/lib.org; export does not tangle, so no source files are written)
  # org >= 29 htmlizes src blocks at export time if htmlize is loadable
  # (config_tangle.el's (package-initialize) pulls it in from
  # ~/.emacs.d/elpa), emitting inline-CSS <span style="..."> markup;
  # org-html-htmlize-output-type nil is the explicit "plain text" switch.
  # We highlight with Pygments later, so the export must stay plain
  # (build_site.py strips stray fontify spans as a defense)
  emacs --no-init-file --no-site-lisp --quick --batch "org/$f" \
    --load="$REPO/tools/config_tangle.el" \
    --eval '(setq org-html-htmlize-output-type nil)' \
    --eval '(org-html-export-to-html)' \
    2> "org/$base.export-log" \
    || { echo "EXPORT FAILED: $f" >&2; cat "org/$base.export-log" >&2; exit 1; }
  [[ -s "org/$base.html" ]] || { echo "missing output org/$base.html" >&2; exit 1; }
  count=$((count+1))
done < org/table_of_contents

echo "exported $count files"

# Tangle the worked-example programs (literate programming): mo_ortho.c and
# ao_grid.f90 are declared with :tangle in qmckl_examples.org; org-babel-tangle
# resolves the <<qmckl_check_error>> noweb reference. build_site.py copies the
# files into the site and links them from the example cards.
echo "tangle org/qmckl_examples.org"
emacs --no-init-file --no-site-lisp --quick --batch "org/qmckl_examples.org" \
  --eval '(org-babel-tangle)' 2> "org/qmckl_examples.tangle-log"
[[ -s "org/mo_ortho.c" && -s "org/ao_grid.f90" ]] || {
  echo "tangle failed: org/mo_ortho.c / org/ao_grid.f90 missing" >&2
  cat "org/qmckl_examples.tangle-log" >&2; exit 1; }
