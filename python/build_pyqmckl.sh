
set -e
set -x

source export_files.sh 

cp qmckl.h pyqmckl.i numpy.i src/

cd src/

swig -python -py3 -o pyqmckl_wrap.c pyqmckl.i 

gcc -c -fPIC -I/usr/include/python3.8 ${C_FILES} pyqmckl_wrap.c -ltrexio

gfortran -c -fPIC -I/usr/include/python3.8 ${F_FILES}

gcc -shared ${C_O_FILES} pyqmckl_wrap.o -o _pyqmckl.so

