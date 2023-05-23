from __future__ import absolute_import, division, print_function
"""
Load params for a simulation.
Iris Young, idyoung@lbl.gov
"""

import libtbx.phil as phil

sim_view_phil_scope = phil.parse(
    """
    mosaic_domains = 100
        .type = int
        .expert_level=1
        .help = Number of rotation matrices to model in each simulation.
    """, process_includes=True)

