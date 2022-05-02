import pyqmckl as pq

from data.data import coord
from os.path import join

ctx = pq.qmckl_context_create()

fname = join('data', 'Alz_small.h5')

rc = pq.qmckl_trexio_read(ctx, fname)
print(pq.qmckl_string_of_error(rc))
