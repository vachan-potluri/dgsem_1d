from _common import *
from change_basis_matrix import *

print_testing_statement("ChangeBasisMatrix")
for N in [2,3]:
    print(f"N = {N}")
    cbm = ChangeBasisMatrix(N)
    print("Matrix:")
    print(cbm.matrix)