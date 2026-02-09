#!/usr/bin/env bash

# ------------------------------------
# Check required environment variables
# ------------------------------------

if [ -z ${srcdir} ] ; then
   echo "Error: srcdir environment variable is not defined"
   exit 1
fi

if [ -z ${top_builddir} ] ; then
   echo "Error: top_builddir environment variable is not defined"
   exit 1
fi


set -euo pipefail

# -----------------------------
# patchelf availability
# -----------------------------

if command -v patchelf >/dev/null 2>&1; then
  HAVE_PATCHELF=1
else
  HAVE_PATCHELF=0
  echo
  echo "WARNING: patchelf not found."
  echo "         RPATH will NOT be set."
  echo "         The resulting binaries may require LD_LIBRARY_PATH."
  echo
fi

# -----------------------------
# Metadata
# -----------------------------

PKG=qmckl
VERSION="$VERSION_MAJOR.$VERSION_MINOR.$VERSION_PATCH"
ARCH=$(uname -m)
OS=linux


STAGE=$(mktemp -d)
PREFIX=${STAGE}/${PKG}-${VERSION}-${OS}-${ARCH}

LIBDIR=${PREFIX}/lib
INCDIR=${PREFIX}/include
LICENSEDIR=${PREFIX}/share/licenses

mkdir -p "${LIBDIR}" "${INCDIR}" "${LICENSEDIR}"

# -----------------------------
# Sanity checks
# -----------------------------

test -f src/.libs/libqmckl.so
test -f src/.libs/libqmckl.a

command -v ldd >/dev/null

# -----------------------------
# Install QMCKL libraries
# -----------------------------

libtool --mode=install install \
  src/.libs/libqmckl.so "${LIBDIR}/libqmckl.so"

install -m 644 src/.libs/libqmckl.a "${LIBDIR}/libqmckl.a"

# -----------------------------
# Install headers
# -----------------------------

install -m 644 include/qmckl.h "${INCDIR}/"
install -m 644 include/qmckl_f.F90 "${INCDIR}/"

# -----------------------------
# Dependency embedding helper
# -----------------------------

embed_deps () {
  local so=$1

  ldd "$so" | awk '/=>/ {print $3}' | while read -r dep; do
    [ -z "$dep" ] && continue

    case "$dep" in
      /*)
        base=$(basename "$dep")

        case "$base" in
          libc.so.*|libpthread.so.*|libdl.so.*|ld-linux-*)
            continue
            ;;
        esac

        if [ ! -f "${LIBDIR}/${base}" ]; then
          cp "$dep" "${LIBDIR}/"
        fi
        ;;
    esac
  done
}



# -----------------------------
# Embed dependencies
# -----------------------------

embed_deps "${LIBDIR}/libqmckl.so"

# Recurse once (for OpenBLAS → libgomp etc.)
for so in "${LIBDIR}"/*.so; do
  embed_deps "$so"
done

# -----------------------------
# Fix RPATH for ALL embedded shared objects
# -----------------------------

if [ "${HAVE_PATCHELF}" -eq 1 ]; then
  echo
  echo "Setting RPATH on embedded shared libraries"

  find "${LIBDIR}" -type f -name '*.so*' | while read -r so; do
    echo "  RPATH -> \$ORIGIN : $(basename "$so")"
    patchelf --set-rpath '$ORIGIN' "$so"
  done
else
  echo
  echo "WARNING: patchelf not available."
  echo "         RPATH not set on embedded dependencies."
fi


# -----------------------------
# License files
# -----------------------------

cp "${srcdir}/LICENSE" "${LICENSEDIR}/qmckl.txt"

# -----------------------------
# Validation
# -----------------------------

echo
echo "Final dependency tree (libqmckl.so):"
ldd "${LIBDIR}/libqmckl.so" || true

# -----------------------------
# Tarball
# -----------------------------

TARBALL=${PKG}-${VERSION}-${OS}-${ARCH}.tar.gz
tar -C "${STAGE}" -czf "${top_builddir}/${TARBALL}" \
  "$(basename "${PREFIX}")"

rm -rf "${STAGE}"

echo
echo "Binary distribution created:"
echo "  ${TARBALL}"


