#!/bin/bash
# Rebuild the whole QMCkl site.
#
#   build.sh REPO_DIR OUT_DIR
#
# REPO_DIR: a checkout of https://github.com/TREX-CoE/qmckl (org/ sources
#            are read; export + tangle write into REPO_DIR/org/)
# OUT_DIR:  the site root receiving docs/, code/, index.html, css/,
#            fonts/, images/ (locally: the site checkout itself)
#
# Python: $QMCKL_SITE_PYTHON, else tools/venv/bin/python (local uv venv),
# else system python3 (CI: python3 -m pip install -r tools/requirements.txt)
set -euo pipefail

TOOLSDIR="$(cd "$(dirname "$0")" && pwd)"
SITE="$(dirname "$TOOLSDIR")"
REPO="${1:?usage: build.sh REPO_DIR OUT_DIR}"
OUT="${2:?usage: build.sh REPO_DIR OUT_DIR}"

if [[ -n "${QMCKL_SITE_PYTHON:-}" ]]; then
  PY="$QMCKL_SITE_PYTHON"
elif [[ -x "$TOOLSDIR/venv/bin/python" ]]; then
  PY="$TOOLSDIR/venv/bin/python"
else
  PY="$(command -v python3 || true)"
fi
[[ -n "$PY" ]] || { echo "no python found (set QMCKL_SITE_PYTHON)" >&2; exit 1; }
[[ -d "$REPO/org" ]] || { echo "repo not found: $REPO" >&2; exit 1; }
mkdir -p "$OUT"

echo "== 1/4  org -> html (emacs) + tangle worked examples"
"$TOOLSDIR/export_docs.sh" "$REPO"

echo "== 2/4  post-process (MathML, Pygments, layout, links)"
"$PY" "$TOOLSDIR/build_site.py" "$REPO" "$OUT"

echo "== 3/4  assets"
mkdir -p "$OUT/css" "$OUT/fonts" "$OUT/images"
copy_if_changed() { [ "$1" -ef "$2" ] || cp -f "$1" "$2"; }
copy_if_changed "$SITE/css/theme.css" "$OUT/css/theme.css"
copy_if_changed "$TOOLSDIR/404.html" "$OUT/404.html"
: > "$OUT/.nojekyll"
for f in "$SITE/fonts/"*.woff2; do
  copy_if_changed "$f" "$OUT/fonts/$(basename "$f")"
done
[[ -f "$REPO/org/chbrclf.png" ]] && copy_if_changed "$REPO/org/chbrclf.png" "$OUT/images/chbrclf.png"
copy_if_changed "$REPO/QMCkl.png" "$OUT/images/qmckl-logo.png"

echo "== 4/4  done"
echo "serve with:  python3 -m http.server --directory $OUT 8000"
