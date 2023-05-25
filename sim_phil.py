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
    oversampling = 1
        .type = int
        .expert_level=1
        .help = Oversampling of the simulated diffraction in case of e.g.
        .help = very large domain size and wide mosaic angles.
    oscillation_n_steps = 10
        .type = int
        .expert_level=1
        .help = Number of steps to sample in a single oscillation.
        .help = (This number of simulations will be summed to produce a
        .help = simulated rotation image.)
    spectrum_shape = SASE (XFEL)
        .type = str
        .expert_level=1
        .help = Enable overriding default spectrum shape to monochromatic for
        .help = speed, critical when not using GPU, recommended when using
        .help = diffuse scattering.
    """, process_includes=True)

