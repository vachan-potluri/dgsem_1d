import common
from q_gauss_lobatto import *

common.print_testing_statement("QGaussLobatto")
for np in [2,3,4,5,6]:
    quad = QGaussLobatto(np)
    print(f"Number of points: {np}")
    print(f"Quadrature points: {quad.q_points}")
    print(f"Quadrature weights: {quad.q_weights}\n")