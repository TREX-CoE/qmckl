""" 
This is the test of the Python API of the QMCkl library.
It is the `bench_mos.c` C code adapted from the `bench`
repo and translated into Python with some modifications.
"""

from os.path import join
import time

import qmckl as pq
from data.data import coord


walk_num  = 100
elec_num  = 158

ITERMAX = 10

ctx = pq.context_create()

try:
    pq.trexio_read(ctx, 'fake.h5')
except RuntimeError:
    print('Error handling check: passed')

fname = join('data', 'Alz_small.h5')

pq.trexio_read(ctx, fname)
print('trexio_read: passed')

pq.set_electron_walk_num(ctx, walk_num)

mo_num = pq.get_mo_basis_mo_num(ctx)
assert mo_num == 404

pq.set_electron_coord(ctx, 'T', coord)

size_max = 5*walk_num*elec_num*mo_num

mo_vgl = pq.get_mo_basis_mo_vgl(ctx, size_max)
assert mo_vgl.size == size_max

start = time.clock_gettime_ns(time.CLOCK_REALTIME)

for _ in range(ITERMAX):
    mo_vgl_in = pq.get_mo_basis_mo_vgl_inplace(ctx, size_max)

end = time.clock_gettime_ns(time.CLOCK_REALTIME)

print(f'Time for the calculation of 1 step : {(end-start)*.000001/ITERMAX} ms')

