from _common import *
from q_gauss_lobatto import *

print_testing_statement("QGaussLobatto")
for npts in [2,3,4,5,6]:
    quad = QGaussLobatto(npts)
    print(f"Number of points: {npts}")
    print(f"Quadrature points: {quad.q_points}")
    print(f"Quadrature weights: {quad.q_weights}\n")