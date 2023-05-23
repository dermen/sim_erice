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
    rotation {
        sweep_deg = None
            .type = float
            .expert_level=2
            .help = Total range of a sweep in rotation mode (deg).
            .help = If provided, this overrides sweep_n_imgs
        sweep_n_imgs = 100
            .type = int
            .expert_level=1
            .help = Number of images in a sweep.
            .help = sweep_deg will override this param when provided.
        oscillation_n_steps = 10
            .type = int
            .expert_level=1
            .help = Number of steps to sample in a single oscillation.
            .help = (This number of simulations will be summed to produce a
            .help = simulated rotation image.)
            .help = phi_step_size will override this param when provided.
        phi_step_size = None
            .type = float
            .expert_level=2
            .help = Size of a step within an oscillation, i.e. the range in
            .help = phi to be covered by a single simulation.
            .help = If provided, this overrides oscillation_n_steps.
        }

    """, process_includes=True)

