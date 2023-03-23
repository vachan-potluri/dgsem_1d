from common import *
from dof_handler import *

class DoFVector:
    # Class for storing and operating on a solution vector defined on dofs
    def __init__(self, dh, n_vars=1):
        assert type(dh) == DoFHandler, "Expected a DoFHandler instance as input"
        assert n_vars>=1, "Zero/negative number of variables received"
        self.dof_handler = dh
        self.entries = np.zeros(dh.n_dofs, dtype=object)
        for i in range(dh.n_dofs):
            self.entries[i] = np.zeros(n_vars)