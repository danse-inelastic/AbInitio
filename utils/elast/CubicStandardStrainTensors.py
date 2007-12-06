__doc__ = """This module defines standard forms for cubic strain tensors."""

import numpy as np

# C11 strain
C11 = np.array( [[1,  0,  0],
                 [0,  0,  0],
                 [0,  0,  0]])
# C11-C12 strain (1st form)
C11C12_1 = np.array( [[1,  0,  0],
                      [0, -1,  0],
                      [0,  0, 0]])
# C11-C12 strain (2nd form)
C11C12_2 = np.array( [[1,  0,  0],
                      [0,  1,  0],
                      [0,  0, -2]])
# C44 strain (1st form)
C44_1 = np.array( [[0,  0,  0],
                   [0,  0,  1],
                   [0,  1,  0]])
# C44 strain (2nd form)
C44_2 = np.array( [[0,  1,  1],
                   [1,  0,  1],
                   [1,  1,  0]])

