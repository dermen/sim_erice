from __future__ import absolute_import, division, print_function
"""
Load params for a simulation.
Iris Young, idyoung@lbl.gov
"""

import libtbx.phil as phil

sim_view_phil_scope = phil.parse(
    """
    numerical_params {
        domain_size {
            min = 100
                .type = float
                .expert_level=2
            max = 100000
                .type = float
                .expert_level=2
            small_step = 100
                .type = float
                .expert_level=2
            big_step = 1000
                .type = float
                .expert_level=2
            default = 1000
                .type = float
                .expert_level=1
            formatter = "%4.2f"
                .type = str
                .expert_level=3
            label = "Domain size (each edge, Å)"
                .type = str
                .expert_level=3
            row = 0
                .type = int
                .expert_level=3
            column = 0
                .type = int
                .expert_level=3
            }
        mos_ang_deg {
            min = 0.01
                .type = float
                .expert_level=2
            max = 5
                .type = float
                .expert_level=2
            small_step = 0.01
                .type = float
                .expert_level=2
            big_step = 0.1
                .type = float
                .expert_level=2
            default = 0.1
                .type = float
                .expert_level=1
            formatter = "%4.2f"
                .type = str
                .expert_level=3
            label = "Mosaic angle (°)"
                .type = str
                .expert_level=3
            row = 0
                .type = int
                .expert_level=3
            column = 1
                .type = int
                .expert_level=3
            }
        ucell_scale_a {
            min = 0.5
                .type = float
                .expert_level=2
            max = 2.0
                .type = float
                .expert_level=2
            small_step = 0.05
                .type = float
                .expert_level=2
            big_step = 0.10
                .type = float
                .expert_level=2
            default = 1.0
                .type = float
                .expert_level=1
            formatter = "%4.2f"
                .type = str
                .expert_level=3
            label = "Unit cell scale (a)"
                .type = str
                .expert_level=3
            row = 1
                .type = int
                .expert_level=3
            column = 0
                .type = int
                .expert_level=3
            }
        ucell_scale_b {
            min = 0.5
                .type = float
                .expert_level=2
            max = 2.0
                .type = float
                .expert_level=2
            small_step = 0.05
                .type = float
                .expert_level=2
            big_step = 0.10
                .type = float
                .expert_level=2
            default = 1.0
                .type = float
                .expert_level=1
            formatter = "%4.2f"
                .type = str
                .expert_level=3
            label = "Unit cell scale (a)"
                .type = str
                .expert_level=3
            row = 1
                .type = int
                .expert_level=3
            column = 1
                .type = int
            }
        ucell_scale_c {
            min = 0.5
                .type = float
                .expert_level=2
            max = 2.0
                .type = float
                .expert_level=2
            small_step = 0.05
                .type = float
                .expert_level=2
            big_step = 0.10
                .type = float
                .expert_level=2
            default = 1.0
                .type = float
                .expert_level=1
            formatter = "%4.2f"
                .type = str
                .expert_level=3
            label = "Unit cell scale (a)"
                .type = str
                .expert_level=3
            row = 1
                .type = int
                .expert_level=3
            column = 2
                .type = int
            }
        rot_x {
            min = -180
                .type = float
                .expert_level=2
            max = 180
                .type = float
                .expert_level=2
            small_step = 0.1
                .type = float
                .expert_level=2
            big_step = 1.0
                .type = float
                .expert_level=2
            default = 0
                .type = float
                .expert_level=1
            formatter = "%6.2f"
                .type = str
                .expert_level=3
            label = "Rotation (x, °)"
                .type = str
                .expert_level=3
            row = 2
                .type = int
                .expert_level=3
            column = 0
                .type = int
                .expert_level=3
            }
        rot_y {
            min = -180
                .type = float
                .expert_level=2
            max = 180
                .type = float
                .expert_level=2
            small_step = 0.1
                .type = float
                .expert_level=2
            big_step = 1.0
                .type = float
                .expert_level=2
            default = 0
                .type = float
                .expert_level=1
            formatter = "%6.2f"
                .type = str
                .expert_level=3
            label = "Rotation (y, °)"
                .type = str
                .expert_level=3
            row = 2
                .type = int
                .expert_level=3
            column = 1
                .type = int
            }
        rot_z {
            min = -180
                .type = float
                .expert_level=2
            max = 180
                .type = float
                .expert_level=2
            small_step = 0.1
                .type = float
                .expert_level=2
            big_step = 1.0
                .type = float
                .expert_level=2
            default = 0
                .type = float
                .expert_level=1
            formatter = "%6.2f"
                .type = str
                .expert_level=3
            label = "Rotation (z, °)"
                .type = str
                .expert_level=3
            row = 2
                .type = int
                .expert_level=3
            column = 2
                .type = int
            }
        diff_gamma {
            min = 1
                .type = float
                .expert_level=2
            max = 300
                .type = float
                .expert_level=2
            small_step = 1
                .type = float
                .expert_level=2
            big_step = 10
                .type = float
                .expert_level=2
            default = 50
                .type = float
                .expert_level=1
            formatter = "%3.0f"
                .type = str
                .expert_level=3
            label = "Diffuse gamma (Å)"
                .type = str
                .expert_level=3
            row = 3
                .type = int
                .expert_level=3
            column = 0
                .type = int
                .expert_level=3
            }
        diff_sigma {
            min = 0.01
                .type = float
                .expert_level=2
            max = 0.7
                .type = float
                .expert_level=2
            small_step = 0.01
                .type = float
                .expert_level=2
            big_step = 0.05
                .type = float
                .expert_level=2
            default = 0.4
                .type = float
                .expert_level=1
            formatter = "%4.2f"
                .type = str
                .expert_level=3
            label = "Diffuse sigma (Å)"
                .type = str
                .expert_level=3
            row = 3
                .type = int
                .expert_level=3
            column = 1
                .type = int
                .expert_level=3
            }
        diff_aniso {
            min = 0.01
                .type = float
                .expert_level=2
            max = 10
                .type = float
                .expert_level=2
            small_step = 0.1
                .type = float
                .expert_level=2
            big_step = 1
                .type = float
                .expert_level=2
            default = 3
                .type = float
                .expert_level=1
            formatter = "%3.1f"
                .type = str
                .expert_level=3
            label = "Diffuse anisotropy"
                .type = str
                .expert_level=3
            row = 3
                .type = int
                .expert_level=3
            column = 2
                .type = int
                .expert_level=3
            }
        delta_phi {
            min = 0.1
                .type = float
                .expert_level=2
            max = 5
                .type = float
                .expert_level=2
            small_step = 0.05
                .type = float
                .expert_level=2
            big_step = 0.5
                .type = float
                .expert_level=2
            default = 0.25
                .type = float
                .expert_level=1
            formatter = "%3.1f"
                .type = str
                .expert_level=3
            label = "Oscillation width (°)"
                .type = str
                .expert_level=3
            row = 4
                .type = int
                .expert_level=3
            column = 0
                .type = int
                .expert_level=3
            }
        image {
            min = 1
                .type = float
                .expert_level=2
            max = 100
                .type = float
                .expert_level=2
            small_step = 1
                .type = float
                .expert_level=2
            big_step = 10
                .type = float
                .expert_level=2
            default = 1
                .type = float
                .expert_level=1
            formatter = "%3.0f"
                .type = str
                .expert_level=3
            label = "Image no."
                .type = str
                .expert_level=3
            row = 4
                .type = int
                .expert_level=3
            column = 1
                .type = int
                .expert_level=3
            }
        brightness {
            min = 0
                .type = float
                .expert_level=2
            max = 2
                .type = float
                .expert_level=2
            small_step = 0.01
                .type = float
                .expert_level=2
            big_step = 0.1
                .type = float
                .expert_level=2
            default = 0.5
                .type = float
                .expert_level=1
            formatter = "%4.2f"
                .type = str
                .expert_level=3
            label = "Brightness"
                .type = str
                .expert_level=3
            row = 4
                .type = int
                .expert_level=3
            column = 2
                .type = int
                .expert_level=3
            }
        energy {
            min = 6500
                .type = float
                .expert_level=2
            max = 12000
                .type = float
                .expert_level=2
            small_step = 100
                .type = float
                .expert_level=2
            big_step = 1000
                .type = float
                .expert_level=2
            default = 9500
                .type = float
                .expert_level=1
            formatter = "%5.0f"
                .type = str
                .expert_level=3
            label = "Beam energy (eV)"
                .type = str
                .expert_level=3
            row = 5
                .type = int
                .expert_level=3
            column = 0
                .type = int
                .expert_level=3
            }
        bandwidth {
            min = 0.1
                .type = float
                .expert_level=2
            max = 2
                .type = float
                .expert_level=2
            small_step = 0.1
                .type = float
                .expert_level=2
            big_step = 1
                .type = float
                .expert_level=2
            default = 0.3
                .type = float
                .expert_level=1
            formatter = "%3.1f"
                .type = str
                .expert_level=3
            label = "Bandwidth (%)"
                .type = str
                .expert_level=3
            row = 5
                .type = int
                .expert_level=3
            column = 1
                .type = int
                .expert_level=3
            }
    }
    categorical_params {
        spectrum_shape {
            default = Monochromatic *SASE_(XFEL) Gaussian
                .type = choice
                .expert_level=1
            options = Monochromatic
                .type = str
                .multiple = True
                .expert_level=3
            options = SASE (XFEL)
                .type = str
                .multiple = True
                .expert_level=3
            options = Gaussian
                .type = str
                .multiple = True
                .expert_level=3
            label = "Spectrum shape"
                .type = str
                .expert_level=3
            row = 6
                .type = int
                .expert_level=3
            column = 0
                .type = int
                .expert_level=3
            widget = *OptionMenu RadioButtons
                .type = choice
                .expert_level=3
            }
        rotation_mode {
            default = *Stills Rotation
                .type = choice
                .expert_level=1
            options = Stills
                .type = str
                .multiple = True
                .expert_level=3
            options = Rotation
                .type = choice
                .multiple = True
                .expert_level=3
            label = "Experiment mode"
                .type = str
                .expert_level=3
            row = 6
                .type = int
                .expert_level=3
            column = 1
                .type = int
                .expert_level=3
            widget = OptionMenu *RadioButtons
                .type = choice
                .expert_level=3
            }
        diffuse_mode {
            default = *Off On
                .type = choice
                .expert_level=1
            options = On
                .type = str
                .multiple = True
                .expert_level=3
            options = Off
                .type = str
                .multiple = True
                .expert_level=3
            label = "Diffuse scattering"
                .type = str
                .expert_level=3
            row = 6
                .type = int
                .expert_level=3
            column = 2
                .type = int
                .expert_level=3
            widget = OptionMenu *RadioButtons
                .type = choice
                .expert_level=3
            }
        pinned_mode {
            default = *Simulation_only Overlay_with_pinned
                .type = choice
                .expert_level=1
            options = Simulation only
                .type = str
                .multiple = True
                .expert_level=3
            options = Overlay with pinned
                .type = str
                .multiple = True
                .expert_level=3
            label = "Display mode"
                .type = str
                .expert_level=3
            row = 7
                .type = int
                .expert_level=3
            column = 0
                .type = int
                .expert_level=3
            widget = OptionMenu *RadioButtons
                .type = choice
                .expert_level=3
            }
        Fhkl {
            default = On *Off
                .type = choice
                .expert_level=1
            options = On
                .type = str
                .multiple = True
                .expert_level=3
            options = Off
                .type = str
                .multiple = True
                .expert_level=3
            label = "Use structure factors"
                .type = str
                .expert_level=3
            row = 7
                .type = int
                .expert_level=3
            column = 1
                .type = int
                .expert_level=3
            widget = OptionMenu *RadioButtons
                .type = choice
                .expert_level=3
            }
        }
    pdb_file = None
        .type = path
        .expert_level=1
    pdb_id = None
        .type = str
        .expert_level=1
    """)

categorical_params_options = {
    "spectrum_shape": ["Monochromatic", "SASE_(XFEL)", "Gaussian"],
    "rotation_mode": ["Stills", "Rotation"],
    "diffuse_mode": ["Off", "On"],
    "pinned_mode": ["Simulation_only", "Overlay_with_pinned"],
    "Fhkl": ["On", "Off"],
}

