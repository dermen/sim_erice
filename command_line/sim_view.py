from __future__ import absolute_import, division, print_function
"""
Plot the PNG files
Derek Mendez, dermen@lbl.gov

Generate and load simulated images on the fly
Iris Young, idyoung@lbl.gov
# LIBTBX_SET_DISPATCHER_NAME simtbx.sim_view
"""
try:  # py2/3 compat 
    import Tkinter as tk
except ImportError:
    import tkinter as tk

import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')
import pylab as plt
import os, math
from scitbx.matrix import col, sqr

from matplotlib.backends.backend_tkagg import \
    FigureCanvasTkAgg, NavigationToolbar2Tk

import libtbx.load_env
from sim_erice.on_the_fly_simdata import run_simdata, get_SIM, randomize_orientation, set_orientation, sweep
from sim_erice.sim_phil import sim_view_phil_scope
from simtbx.nanoBragg.tst_nanoBragg_multipanel import beam, whole_det
from simtbx.diffBragg import hopper_utils
from sim_erice.local_spectra import spectra_simulation
from iotbx.crystal_symmetry_from_any import extract_from as extract_symmetry_from
from iotbx.pdb.fetch import get_pdb
from cctbx.uctbx import unit_cell
from dxtbx_model_ext import Crystal
from dials.util.options import ArgumentParser
from scitbx import matrix
from scitbx.math import gaussian, euler_angles
from dials.array_family import flex
from libtbx import easy_pickle
from random import randint
import time
from simtbx.diffBragg.utils import ENERGY_CONV, get_laue_group_number

usage="""Usage: libtbx.python sim_view.py [--help] [--config] [--file pdb_filename|--fetch pdbid] [phil_params]

To display the full set of available parameters, use the --config (-c) flag.
To display a comprehensive explanation of the simulator, use the --help (-h) flag.
A model file can be passed with the --file flag instead of the default model, or a model may be fetched from the PDB by supplying the PDB ID after the --fetch flag.
"""

help_message="""SimView: lightweight simulator and viewer for diffraction still images.

%s
Run without any arguments, this program simulates diffraction of a small lysozyme crystal. You may ask the viewer to fetch a different model by entering that PDB ID in the associated text box.

Diffraction from your crystal will be simulated with a set of default parameters for mosaic block size, beam energy and bandpass, etc. that you can change by adjusting each of the dials. Parameters to be adjusted can be selected either with the visual controls or with keyboard shortcuts. Numerical visual controls are made "active" by clicking in the text box or pressing their up or down arrows. The active control is labeled in bold, blue text. The keyboard will also route input to the active control by interpreting a keyboard "up arrow" keypress as the visible "up arrow" button, and the same for the "down arrow". For finer control, you may also enter text directly into the text box. The GUI will respond to new text input only after you press enter. In summary, keyboard shortcuts are as follows:

Up arrow:           increase this parameter (same as pressing the small up arrow next to the control)
Down arrow:         decrease this parameter (same as pressing the small down arrow next to the control)
Shift-up arrow:     increase this parameter a lot
Shift-down arrow:   decrease this parameter a lot
Return:             finalize a changed parameter typed directly into the box

Explanation of parameters:

DomainSize: the small blocks within crystals that are internally perfectly aligned and identical are called mosaic domains. The average size of a domain is one contributing factor determining the sharpness of the Bragg peaks. As the size of the mosaic domain increases, the effect of constructive interference of the scattered X-rays gets stronger, and the peak sharpens.

MosAngDeg: the mosaic domains are misoriented relative to one another by an average angle denoted the mosaic angle, displayed here in degrees. Larger differences in the orientations of the domains produce contributions to Bragg peaks over a wider area on the detector, resulting in more diffuse spots.

Brightness: the overall image scale, which governs the minimum intensity value that will be visible. This is a tradeoff with the maximum value that can be visually differentiated, where everything above this maximum is marked as overflow and colored in yellow in the simulation-only view mode.

Display mode: either "Simulation only" to display a single simulation in cyan, or "Overlay with pinned" to also display a second simulation in the red channel of the same image, generating grayscale where they perfectly match. This can be used to emphasize and explore the effects of changing parameters, by "pinning" a simulation (using the "Update pinned image" button) and then adjusting a parameter.

Unit cell a, b, c: a, b and c are the unit cell axis lengths, and we expose fractional scaling parameters for each of these to simulate e.g. a 10%% larger unit cell in the b axis when setting b_scale to 1.1 -- you should be able to see both the fractional scales in the controls and the resulting cell parameters in the static text to their right. Depending on the symmetry of the crystal, these may all vary independently or may be constrained to scale together, in which case not all unit cell scale controls will be enabled. Larger unit cells produce constructive and destructive interference at smaller angles, resulting in more closely spaced Bragg peaks on the detector.

Missetting angles: in serial crystallography, we are often faced with differences between the true orientation of the crystal and what is determined during indexing. A misorientation of the X or Y axis means a misestimation of which spots are in diffracting conditions, and will affect which spots (or how much of the spots) appear on the detector. Since the Z axis is along the beam direction, misorientation in Z results in rotation of the entire diffraction pattern. In the simulation viewer we expose these misorientations as rotations in x, y and z separately from the crystal orientation matrix which can be updated to a new, random orientation with the "Randomize orientation" button.

Experiment mode: a default mode of stills (serial) crystallography is enabled in the viewer, but you may also switch to rotation mode to simulate a standard experiment with a goniometer. The goniometer axis is identified by a vertical yellow line, and the oscillation width and image number in a sweep can be controlled in the GUI. Advanced users may also want to adjust oscillation step size or total sweep angle with phil parameters when launching the GUI.

Energy/Bandwidth: the energy of the X-ray beam denotes the average energy of an individual pulse. For X-ray free electron laser (XFEL) pulses generated by self-amplified spontaneous emission (SASE), each pulse is spiky and stochastic, with a typical bandwidth of between 10 and 30 eV. (A monochromater may be used to narrow this to around 1 eV, not shown.) The average energy of a series of pulses is matched to the energy entered in this control. You can advance to the next SASE spectrum in the series by pressing "New XFEL pulse." Bandwidth is only relevant to the Gaussian spectrum mode.

Spectrum shape: to simulate images with real SASE spectra, toggle this option to "SASE". To display diffraction produced by a smooth Gaussian function that can be adjusted in energy and bandwidth (illustrative but nonphysical), toggle this to "Gaussian". SASE spectra are stochastic and will intentionally not be modified by the energy and bandwidth controls. The "monochromatic" option is recommended when diffuse scattering is enabled, to offset the greater computational cost.

Use structure factors: when we observe diffraction from a single particle, we are seeing the Fourier transform of the particle, with [a great deal of] added noise. For many copies of a particle all in the same orientation, we can expect the signal to be much stronger. If these aligned particles are also periodically spaced, as in the context of a crystal, they act as a diffraction grating, and we see this signal only at positions of constructive interference (the reciprocal lattice points in diffracting conditions, i.e. Bragg peaks). In order to better see the effects of the above parameters, we can ignore the Fourier transform of the asymmetric unit and pretend all of the Bragg peaks have uniform intensity. This is what we are simulating when we set "Use structure factors" on or off.

Diffuse scattering: this physical effect is a result of variation and imperfection within the crystal and surrounding solvent, and it appears as X-ray scattering between and around the Bragg peaks. Modeling this can be toggled on or off. It is recommended to use the monochromatic spectrum setting when calculating diffuse signal.

Diff_gamma and diff_sigma: parameters describing the diffuse scattering signal produced by long-range correlations. Gamma denotes the correlation length in Ångstroms and sigma squared is the amplitude of the correlated vibrations in square Ångstroms.

Diff_aniso: anisotropic quality to the diffuse scattering. This is arbitrarily assigned directionality, and the adjustible parameter controls the degree to which this effect is observed.

Further notes on the image viewer: matplotlib enables zooming and panning of the figure. We have also added pixel-specific information to the annotation text displayed in the lower right corner, such as the fractional Miller index and resolution at a given position on the detector.
""" % usage

class NumericalParam(object):
    def __init__(self, min, max, small_step, big_step, default,
                 formatter='%4.2f', units_string='', label='', position=None, pos_hivis=None, columnspan=None):
        self.min = min
        self.max = max
        self.sstep = small_step
        self.bstep = big_step
        self.default = default
        self.formatter = formatter
        self.decimals = -1*int(math.log10(self.sstep)-1) # for rounding
        self.units = units_string
        self.label = label
        self.position = position
        self.pos_hivis = pos_hivis
        self.columnspan = columnspan
        self.is_enabled = True
        self.is_active = False
    def register_parent_frame(self, parent_frame):
        self.parent_frame = parent_frame
    def generate_dial(self, extra_logic=lambda: None, hivis=False):
        self.frame = tk.Frame(self.parent_frame)
        self.hivis = hivis
        if hivis:
            row, column = self.pos_hivis
        else:
            row, column = self.position
        self.frame.grid(row=row, column=column,
                        columnspan=self.columnspan,
                        sticky='w', padx=6)
        self.f_label = tk.Label(self.frame, text=self.label)
        self.f_label.pack(side=tk.LEFT, padx=4)
        self.f_units = tk.Label(self.frame, text=self.units)
        self.f_units.pack(side=tk.RIGHT)
        if self.sstep >= 1:
            self.variable = tk.IntVar()
        else:
            self.variable = tk.DoubleVar()
        self.variable.set(self.default)
        self.command = self.make_command(extra_logic=extra_logic)
        self.f_ctrl = tk.Spinbox(self.frame,
                                 from_=self.min,
                                 to=self.max,
                                 increment=self.sstep,
                                 command=self.command,
                                 format=self.formatter,
                                 textvariable=self.variable)
        self.f_ctrl.pack(side=tk.RIGHT)
        self.f_ctrl.bind("<FocusIn>", self.activate)
        self.f_ctrl.bind("<FocusOut>", self.deactivate)
    def make_command(self, extra_logic):
        def on_update_dial(tkevent=None):
            self.activate()
            extra_logic()
        return on_update_dial
    def get_value(self):
        if hasattr(self, 'variable'):
            return self.variable.get()
        else:
            return self.default
    def set_value(self, new_value, callbacks=False):
        self.variable.set(new_value)
        if callbacks:
            self.command()
    def set_small_steps(self):
        self.f_ctrl["increment"] = self.sstep
    def set_big_steps(self):
        self.f_ctrl["increment"] = self.bstep
    def reset(self, callbacks=True):
        self.set_value(self.default, callbacks=callbacks)
    def enable(self):
        for part in (self.f_label, self.f_ctrl, self.f_units):
            part.config(state='normal')
        self.is_enabled = True
    def disable(self):
        for part in (self.f_label, self.f_ctrl, self.f_units):
            part.config(state='disabled')
        self.is_enabled = False
    def activate(self, tkevent=None):
        if hasattr(self, 'is_active') and not self.is_active:
            for part in (self.f_label, self.f_units):
                part.config(font='Helvetica %d bold' % (16 if self.hivis else 10), fg='blue')
            self.is_active = True
            self.f_ctrl.focus_set()
            if hasattr(self, 'handler') and not self.handler.current_param is self:
                self.handler.set_active_param_by_object(self)
    def deactivate(self, tkevent=None):
        for part in (self.f_label, self.f_units):
            part.config(font='Helvetica %d' % (16 if self.hivis else 10), fg='black')
        self.is_active = False

class CategoricalParam(NumericalParam):
    def __init__(self, default, options, label, position=None, pos_hivis=None, columnspan=None):
        self.default = default
        self.options = options
        self.label = label
        self.position = position
        self.pos_hivis = pos_hivis
        self.columnspan = columnspan
        self.is_enabled = True
    def generate_menu(self, extra_logic=lambda: None, hivis=False):
        self.frame = tk.Frame(self.parent_frame)
        if hivis:
            row, column = self.pos_hivis
        else:
            row, column = self.position
        self.frame.grid(row=row, column=column,
                        columnspan=self.columnspan,
                        sticky='w', padx=6)
        self.f_label = tk.Label(self.frame, text=self.label)
        self.f_label.pack(side=tk.LEFT, padx=4)
        self.variable = tk.StringVar()
        self.variable.set(self.default)
        self.command = self.make_command(extra_logic=extra_logic)
        self.f_menu = tk.OptionMenu(self.frame,
                                    self.variable,
                                    *self.options,
                                    command=self.command)
        self.f_menu.pack(side=tk.LEFT)
    def enable(self):
        for part in (self.f_label, self.f_menu):
            part.config(state='normal')
        self.is_enabled = True
    def disable(self):
        for part in (self.f_label, self.f_menu):
            part.config(state='disabled')
        self.is_enabled = False

class RadioParam(CategoricalParam):
    def generate_menu(self, extra_logic=lambda: None, hivis=False):
        self.frame = tk.Frame(self.parent_frame)
        if hivis:
            row, column = self.pos_hivis
        else:
            row, column = self.position
        self.frame.grid(row=row, column=column,
                        columnspan=self.columnspan,
                        sticky='w', padx=6)
        self.f_label = tk.Label(self.frame, text=self.label)
        self.f_label.pack(side=tk.LEFT, padx=4)
        self.intvar = tk.IntVar()
        self.intvar.set(self.options.index(self.default))
        self.command = self.make_command(extra_logic=extra_logic)
        self.option_radios = []
        for i, option in enumerate(self.options):
            radio = tk.Radiobutton(self.frame,
                                   indicatoron=0,
                                   text=option,
                                   variable=self.intvar,
                                   value=i,
                                   command=self.command)
            radio.pack(side=tk.LEFT, padx=4)
            self.option_radios.append(radio)
    def get_value(self):
        return self.options[self.intvar.get()]
    def set_value(self, new_value):
        self.intvar.set(self.options.index(new_value))
    def reset(self, callbacks=True):
        self.intvar.set(self.options.index(self.default))
    def enable(self):
        for part in [self.f_label] + self.option_radios:
            part.config(state='normal')
        self.is_enabled = True
    def disable(self):
        for part in [self.f_label] + self.option_radios:
            part.config(state='disabled')
        self.is_enabled = False

class InfoParam(object):
    def __init__(self, info='', position=None, pos_hivis=None, columnspan=None):
        self.info = info
        self.position = position
        self.pos_hivis = pos_hivis
        self.columnspan = columnspan
    def register_parent_frame(self, parent_frame):
        self.parent_frame = parent_frame
    def generate_info(self, hivis=False):
        self.variable = tk.StringVar()
        self.variable.set(self.info)
        self.f_info = tk.Message(self.parent_frame,
                                 textvariable=self.variable,
                                 width=1200)
        if hivis:
            row, column = self.pos_hivis
        else:
            row, column = self.position
        self.f_info.grid(row=row, column=column,
                        columnspan=self.columnspan,
                        sticky='w', padx=6)
    def set_value(self, new_info):
        if hasattr(self, "variable"):
            self.variable.set(new_info)

class ParamsHandler(object):
    def __init__(self, dict_, track_current_param=False):
        self.all_param_names = []
        self.all_params = []
        for key,value in dict_.items():
            self.__setattr__(key, value)
            self.all_param_names.append(key)
            self.all_params.append(value)
        self.current_param_name = None
        self.current_param = None
    def all_register_parent_frame(self, parent_frame):
        for param in self.all_params:
            param.register_parent_frame(parent_frame)
    def all_register_handler(self):
        for param in self.all_params:
            param.handler = self
    def get_param(self, param_name):
        return self.__getattribute__(param_name)
    def get_param_name(self, param):
        idx = self.all_params.index(param)
        return self.all_param_names[idx]
    def reset_all(self):
        for param_name in self.all_param_names:
            self.get_param(param_name).reset(callbacks=False)
    def get_enabled_param_names(self):
        return [name for name in self.all_param_names if self.get_param(name).is_enabled]
    def set_active_param_by_object(self, new_param):
        if self.current_param is not new_param:
            if self.current_param is not None:
                self.current_param.deactivate()
            self.current_param_name = self.get_param_name(new_param)
            self.current_param = new_param
            self.current_param.activate()
    def set_active_param_by_name(self, new_param_name):
        if self.current_param_name != new_param_name:
            if self.current_param is not None:
                self.current_param.deactivate()
            self.current_param_name = new_param_name
            self.current_param = self.get_param(new_param_name)
            self.current_param.activate()
    def activate_next(self, direction=1):
        if self.current_param_name is None:
            try:
                self.current_param_name = self.get_enabled_param_names()[0]
                self.current_param = self.get_param(self.current_param_name)
                self.current_param.activate()
            except IndexError:
                return
        else:
            enabled_names = self.get_enabled_param_names()
            try:
                next_idx = self.all_param_names.index(self.current_param_name) + 1*direction
                while not self.all_param_names[next_idx] in enabled_names:
                    next_idx += 1*direction
            except IndexError:
                next_idx = 0
            self.current_param.deactivate()
            self.current_param_name = self.all_param_names[next_idx]
            self.current_param = self.get_param(self.current_param_name)
            self.current_param.activate()
    def activate_previous(self):
        self.activate_next(direction=-1)
    def set_small_steps(self):
        for param in self.all_params:
            param.set_small_steps()
    def set_big_steps(self):
        for param in self.all_params:
            param.set_big_steps()

class Button(object):
    def __init__(self, parent_frame, command, label, position, pos_hivis=None, hivis=False):
        self.command = command
        self.button = tk.Button(parent_frame, command=command, text=label)
        if hivis:
            row, column = pos_hivis
        else:
            row, column = position
        self.button.grid(row=row, column=column,
                         sticky='n'+'s'+'e'+'w')
    def press(self):
        self.command()
    def enable(self):
        self.button.configure(state='normal')
    def disable(self):
        self.button.configure(state='disabled')

class TextEntry(object):
    def __init__(self, parent_frame, command, validate_command, label, placeholder_text, position, pos_hivis=None, hivis=False, master=None):
        self.command = self.make_command(command)
        self.validate_command = validate_command
        self.label = label
        self.position = position
        self.pos_hivis = pos_hivis
        self.hivis = hivis
        self.master = master
        self.variable = tk.StringVar()
        self.variable.set(placeholder_text)
        self.frame = tk.Frame(parent_frame)
        if hivis:
            row, column = self.pos_hivis
        else:
            row, column = self.position
        self.frame.grid(row=row, column=column,
                        sticky='w', padx=6)
        self.f_label = tk.Label(self.frame, text=self.label)
        self.f_label.pack(side=tk.LEFT, padx=4)
        self.f_entry = tk.Entry(self.frame,
                                textvariable=self.variable,
                                validatecommand=validate_command)
        self.f_entry.pack(side=tk.LEFT, padx=4)
        self.f_entry.bind("<FocusIn>", self.activate)
        self.f_entry.bind("<FocusOut>", self.deactivate)
        self.is_active=False
    def get_value(self):
        return self.f_entry.get()
    def set_value(self, new_value):
        self.f_entry.delete(0, tk.END)
        self.f_entry.insert(new_value)
    def make_command(self, command):
        def update(tkevent=None):
            text = self.variable.get()
            try:
                self.validate_command(text)
                command(text)
                self.f_entry.configure(fg='black')
            except Exception as e:
                self.f_entry.configure(fg='red')
                self.master._update_status("Failed to load model.")
                raise(e)
        return update
    def activate(self, tkevent=None):
        if self.hivis:
            self.f_label.config(font='Helvetica 16 bold', fg='blue')
        else:
            self.f_label.config(font='Helvetica 10 bold', fg='blue')
        self.is_active = True
    def deactivate(self, tkevent=None):
        if self.hivis:
            self.f_label.config(font='Helvetica 16', fg='black')
        else:
            self.f_label.config(font='Helvetica 10', fg='black')
        self.is_active = False

params_num = ParamsHandler({
    'domain_size':  NumericalParam(min=100,     max=100000, small_step=100,     big_step=1000,  default=1000,   formatter='%4.0f',  units_string=' Å',  label='Domain size (each edge)',    position=(0,0),   pos_hivis=(1,3)),
    'mos_ang_deg':  NumericalParam(min=0.01,    max=5,      small_step=0.01,    big_step=0.1,   default=0.1,    formatter='%4.2f',  units_string='°',   label='Mosaic angle',               position=(0,1),   pos_hivis=(2,3)),
    'ucell_scale_a':NumericalParam(min=0.5,     max=2,      small_step=0.05,    big_step=0.1,   default=1,      formatter='%4.2f',  units_string='',    label='Unit cell scale (a)',        position=(1,0),   pos_hivis=(0,0)),
    'ucell_scale_b':NumericalParam(min=0.5,     max=2,      small_step=0.05,    big_step=0.1,   default=1,      formatter='%4.2f',  units_string='',    label='Unit cell scale (b)',        position=(1,1),   pos_hivis=(0,1)),
    'ucell_scale_c':NumericalParam(min=0.5,     max=2,      small_step=0.05,    big_step=0.1,   default=1,      formatter='%4.2f',  units_string='',    label='Unit cell scale (c)',        position=(1,2),   pos_hivis=(0,2)),
    'rot_x':        NumericalParam(min=-180,    max=180,    small_step=0.1,     big_step=1,     default=0,      formatter='%6.2f',  units_string='°',   label='Rotation (x)',               position=(2,0),   pos_hivis=(2,0)),
    'rot_y':        NumericalParam(min=-180,    max=180,    small_step=0.1,     big_step=1,     default=0,      formatter='%6.2f',  units_string='°',   label='Rotation (y)',               position=(2,1),   pos_hivis=(2,1)),
    'rot_z':        NumericalParam(min=-180,    max=180,    small_step=0.1,     big_step=1,     default=0,      formatter='%6.2f',  units_string='°',   label='Rotation (z)',               position=(2,2),   pos_hivis=(2,2)),
    'diff_gamma':   NumericalParam(min=1,       max=300,    small_step=10,      big_step=30,    default=50,     formatter='%3.0f',  units_string=' Å',  label='Diffuse gamma',              position=(3,1),   pos_hivis=(3,1)),
    'diff_sigma':   NumericalParam(min=0.01,    max=0.7,    small_step=0.01,    big_step=0.05,  default=0.4,    formatter='%4.2f',  units_string=' Å',  label='Diffuse sigma',              position=(3,2),   pos_hivis=(3,2)),
    'diff_aniso':   NumericalParam(min=0.01,    max=5,      small_step=0.1,     big_step=1,     default=3,      formatter='%3.1f',  units_string='',    label='Diffuse anisotropy',         position=(3,3),   pos_hivis=(3,3)),
    'delta_phi':    NumericalParam(min=0.1,     max=5,      small_step=0.05,    big_step=0.5,   default=0.25,   formatter='%3.1f',  units_string='°',   label='Oscillation width',          position=(4,1),   pos_hivis=(4,3)),
    'image':        NumericalParam(min=1,       max=100000, small_step=1,       big_step=10,    default=1,      formatter='%3.0f',  units_string='',    label='Image no.',                  position=(4,2),   pos_hivis=(4,2)),
    'brightness':   NumericalParam(min=0,       max=2,      small_step=0.01,    big_step=0.1,   default=0.7,    formatter='%4.2f',  units_string='',    label='Brightness',                 position=(0,2),   pos_hivis=(5,3)),
    'energy':       NumericalParam(min=6500,    max=12000,  small_step=100,     big_step=1000,  default=9500,   formatter='%5.0f',  units_string=' eV', label='Beam energy',                position=(4,3),   pos_hivis=(5,0)),
    'bandwidth':    NumericalParam(min=0.01,    max=2,      small_step=0.1,     big_step=1,     default=0.3,    formatter='%3.1f',  units_string='%',   label='Bandwidth',                  position=(4,4),   pos_hivis=(5,1)),
}, track_current_param=True)
params_cat = ParamsHandler({
    'spectrum_shape':   CategoricalParam(default='SASE (XFEL)',         options=['Monochromatic', 'Gaussian', 'SASE (XFEL)'],                 label='Spectrum',                 position=(4,5),   pos_hivis=(4,0)),
    'rotation_mode':    CategoricalParam(default='Stills',              options=['Stills', 'Rotation'],                                       label='Experiment mode',          position=(4,0),   pos_hivis=(4,1)),
    'diffuse_mode':     RadioParam(default='Off',                       options=['Off', 'On'],                                                label='Diffuse scattering',       position=(3,0),   pos_hivis=(3,0)),
    'pinned_mode':      RadioParam(default='Simulation only',           options=['Simulation only', 'Overlay with pinned'],                   label='Display mode',             position=(0,3),   pos_hivis=(6,0),  columnspan=2),
    'Fhkl':             RadioParam(default='Off',                       options=['On', 'Off'],                                                label='Use structure factors',    position=(5,0),   pos_hivis=(6,3)),
})
params_info = ParamsHandler({
    'status':   InfoParam(info='Initializing...', position=(5,3), pos_hivis=(9,2), columnspan=2),
    'ucell':    InfoParam(info='', position=(1,3), pos_hivis=(1,0), columnspan=2),
    'sg':       InfoParam(info='', position=(1,5), pos_hivis=(1,2), columnspan=1),
    'recs':     InfoParam(info='It is recommended to use a monochromatic beam for diffuse scattering, for speed.', position=(3,4), pos_hivis=(8,0), columnspan=3),
    'umatrix':  InfoParam(info='', position=(2,3), pos_hivis=(7,1)),
    'shortcuts':InfoParam(info='[SHIFT][up/down arrow key] controls numerical parameters', position=(5,5), pos_hivis=(9,0), columnspan=2),
})

class SimView(tk.Frame):

    def __init__(self, master, params_num, params_cat, params_info, params_hyper, pdbfile, *args, **kwargs):
        tk.Frame.__init__(self, *args, **kwargs)

        self.master = master
        self.params_num = params_num
        self.params_cat = params_cat
        self.params_info = params_info
        self.params_hyper = params_hyper
        self.hivis_mode = self.params_hyper.high_visibility

        self.panel = whole_det[0]
        self.s0 = beam.get_unit_s0()
        fast = self.panel.get_fast_axis()
        slow = self.panel.get_slow_axis()
        offset_orig = (2, -24, -1*self.params_hyper.detector_distance)
        self.panel.set_frame(fast, slow, offset_orig)
        self.beam_center = self.panel.get_beam_centre_px(self.s0)

        self.on_update_pdb(pdbfile, init=True)

        self.SASE_sim = spectra_simulation()
        self.on_update_spectrum(init=True, new_pulse=True, skip_gen_image_data=True)
        self.on_update_diffuse_params(skip_gen_image_data=True)
        self.annotate_hkl = True # display hkl for mouse location

        fsize, ssize = whole_det[0].get_image_size()
        img_sh = 1,ssize, fsize
        self.pfs = hopper_utils.full_img_pfs(img_sh)
        self.fast = self.pfs[1::3].as_numpy_array()
        self.slow = self.pfs[2::3].as_numpy_array()
        self.pan = self.pfs[0::3].as_numpy_array()
        self.img_ref = np.zeros((1, ssize, fsize))
        self.img_sim = np.zeros((1, ssize, fsize))
        self.img_overlay = np.zeros((ssize, fsize, 3))
        self.img_single_channel = np.zeros((ssize, fsize, 3))
        self.start_ori = self.SIM.crystal.dxtbx_crystal.get_U()

        self._set_option_menu()
        self._update_ucell_label(update_sg=True)
        self.on_update_orientation(0,0,0)

        self._init_fig()

        self._pack_canvas()
        self._init_display()
        self._set_visual_defaults()
        self._update_dial(self.params_num.current_param_name)
        self.on_toggle_diffuse_mode(skip_gen_image_data=True)
        self.on_toggle_rotation_mode(skip_gen_image_data=True)
        self.on_toggle_spectrum_shape(skip_gen_image_data=True)
        self._enforce_symmetry_on_controls()
        self._update_status("Ready.")

        self.bind()

    def _update_status(self, new_status):
        self.params_info.get_param("status").set_value(
            "STATUS: " + new_status)
        root.update()

    def instantiate_sims(self):
        if self.params_cat.diffuse_mode == 'On':
            mosaic_domains = self.params_hyper.mosaic_domains_diffuse
        else:
            mosaic_domains = self.params_hyper.mosaic_domains_bragg
        self.SIM = get_SIM(whole_det, beam, self.cryst, self.pdbfile, defaultF=0,
                           oversample=self.params_hyper.oversampling, mosaic_domains=mosaic_domains)
        self.default_amp = np.median(self.SIM.crystal.miller_array.data())
        self.SIM_noSF = get_SIM(whole_det, beam, self.cryst, self.pdbfile, defaultF=self.default_amp, SF=False,
                                oversample=self.params_hyper.oversampling, mosaic_domains=mosaic_domains)
        self.fix_diffBragg_instances()

    def fix_diffBragg_instances(self):
        for SIM in (self.SIM, self.SIM_noSF):
            SIM.D.laue_group_num = self.laue_num
            SIM.D.stencil_size = 1

    def on_update_pdb(self, pdbfile, init=False):
        self._update_status("Loading new PDB model...")
        self.pdbfile = pdbfile
        symmetry = extract_symmetry_from(pdbfile)
        self.laue_num = get_laue_group_number(str(symmetry.space_group_info()))
        sg = str(symmetry.space_group_info())
        if sg == 'P 1':
            raise Exception("Triclinic cells not yet supported.")
        ucell = symmetry.unit_cell()
        fmat = matrix.sqr(ucell.fractionalization_matrix()).transpose()
        self.cryst = Crystal(fmat, sg)
        self.ori = self.cryst.get_U()
        self.instantiate_sims()
        self._make_miller_lookup()  # make a dictionary for faster lookup
        self.xtal = self.SIM.crystal.dxtbx_crystal
        self.ucell = self.xtal.get_unit_cell().parameters()
        self.scaled_ucell = self.ucell
        if not init:
            for axis in ('a','b','c'):
                self.params_num.get_param('ucell_scale_%s' % axis).reset(callbacks=False)
        self.on_update_domain_size(skip_gen_image_data=True)
        self.sg = self.xtal.get_space_group()
        self._check_symmetry()
        if not init:
            self._enforce_symmetry_on_controls()
            self._generate_image_data()
            self._update_ucell_label(update_sg=True)
            self._randomize_orientation(reset=True)

    def _make_miller_lookup(self):
        self._update_status("Generating Miller lookup...")
        ma = self.SIM.crystal.miller_array
        self.amplitude_lookup = {h:val for h,val in zip(ma.indices(), ma.data())}

    def _pack_canvas(self):
        """ embed the mpl figure"""
        self._update_status("Preparing figure canvas...")
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.master) 
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP,fill=tk.BOTH,
            expand=tk.YES)

        self.toolbar = NavigationToolbar2Tk(self.canvas, self.master)
        self.toolbar.update()
        self.canvas.get_tk_widget().pack(side=tk.TOP, 
            fill=tk.BOTH, expand=tk.YES)

    def _set_visual_defaults(self):
        """set defaults for font, size, layout, etc."""
        self._update_status("Setting visual defaults...")
        self.default_font = tk.font.nametofont("TkDefaultFont")
        if self.hivis_mode:
            self.default_font.config(family="Helvetica", size=16)
        else:
            self.default_font.config(family="Helvetica", size=10)
        for categorical_option in self.params_cat.all_params:
            if hasattr(categorical_option, 'f_menu'):
                menu = categorical_option.f_menu
                formatter = menu.nametowidget(menu.menuname)
                formatter.config(font=self.default_font)
        for numerical_option in self.params_num.all_params:
            spinbox = numerical_option.f_ctrl
            spinbox.config(font=self.default_font, width=5)
        self.pdb_entry.f_label.config(font=self.default_font)
        self.pdb_entry.f_entry.config(font=self.default_font)
        for info_param in self.params_info.all_params:
            info_param.f_info.config(font=self.default_font)
        self.params_num.activate_next()

    def _get_miller_index_at_mouse(self, x,y,rot_p):
        t = time.time()
        U = sqr(self.SIM.D.Umatrix)
        xx = col((-1, 0, 0))
        yy = col((0, -1, 0))
        zz = col((0, 0, -1))
        RX = xx.axis_and_angle_as_r3_rotation_matrix(rot_p[0], deg=False)
        RY = yy.axis_and_angle_as_r3_rotation_matrix(rot_p[1], deg=False)
        RZ = zz.axis_and_angle_as_r3_rotation_matrix(rot_p[2], deg=False)
        M = RX * RY * RZ
        rotated_U = M*U
        B = sqr(self.SIM.D.Bmatrix)
        A = rotated_U*B
        A_real = A.inverse()

        xmm, ymm = self.panel.pixel_to_millimeter((x,y))
        s = np.array(self.panel.get_lab_coord((xmm, ymm)))
        s /= np.linalg.norm(s)
        s0 = np.array(self.SIM.beam.nanoBragg_constructor_beam.get_unit_s0())
        wavelen = ENERGY_CONV / self.params_num.energy.get_value()
        q = col( (s-s0) / wavelen )
        hkl_f = np.array(A_real*q) 
        hkl_i = np.ceil(hkl_f - 0.5)
        t = time.time()-t

        return hkl_f, hkl_i

    def _get_diffuse_gamma_portion(self,hkl_f,hkl_i,rot_p):
        U = sqr(self.SIM.D.Umatrix)
        xx = col((-1, 0, 0))
        yy = col((0, -1, 0))
        zz = col((0, 0, -1))
        RX = xx.axis_and_angle_as_r3_rotation_matrix(rot_p[0], deg=False)
        RY = yy.axis_and_angle_as_r3_rotation_matrix(rot_p[1], deg=False)
        RZ = zz.axis_and_angle_as_r3_rotation_matrix(rot_p[2], deg=False)
        M = RX * RY * RZ
        rotated_U = M*U
        B = sqr(self.SIM.D.Bmatrix)
        A = rotated_U*B
        A_real = A.inverse()

        _hkl_f = col(hkl_f)
        _hkl_i = col(hkl_i)

        gamma = sqr((self.diffuse_gamma[0],0,0,
                     0,self.diffuse_gamma[1],0,
                     0,0,self.diffuse_gamma[2]))
        delta_Q = A * (_hkl_f - _hkl_i)
        anisoG_q = gamma*delta_Q
        V_dot_V = anisoG_q.dot(anisoG_q)
        gamma_portion_denom = (1.+V_dot_V * 4.*np.pi*np.pi)
        gamma_portion_denom *= gamma_portion_denom
        gamma_portion = 8.*np.pi*gamma.determinant()/gamma_portion_denom
        return gamma_portion
    
    def _annotate(self):
        self.ax.plot(*self.beam_center,'y+', markersize=20)
        def label_mouse_coords(x, y):
            resol = self.panel.get_resolution_at_pixel(self.s0, (x, y))

            H_str = ""
            if self.annotate_hkl:
                rot_p = (self.params_num.rot_x.get_value() * math.pi/180.,
                         self.params_num.rot_y.get_value() * math.pi/180.,
                         self.params_num.rot_z.get_value() * math.pi/180.)
                hkl_f, hkl_i = self._get_miller_index_at_mouse(x,y,rot_p)
                hkl_dist = np.sqrt(np.sum((hkl_f-hkl_i)**2))

                hkl_key = tuple(hkl_i.astype(int))
                if self.params_cat.Fhkl.get_value() == 'On':
                    if hkl_key in self.amplitude_lookup:
                        amp = self.amplitude_lookup[hkl_key]
                    else:
                        amp = 0
                else:
                    amp = self.default_amp

                hkl_f_str = ", ".join(map(lambda x: "%.2f"%x, hkl_f))
                hkl_i_str = ", ".join(map(lambda x: "%d"%x, hkl_i))
                H_str = "hf,kf,lf=({Hf}) |Hdist|={dH:.3f} |F|={F:.6f} Value={Value:.6f} ".format(
                        Hf=hkl_f_str, Hi=hkl_i_str, F=amp, dH=hkl_dist, Value=self.img_sim[0,int(y),int(x)])
                if self.params_cat.diffuse_mode.get_value() == 'On':
                    gamma_portion = self._get_diffuse_gamma_portion(hkl_f,hkl_i,rot_p)
                    H_str += "G(q)={diff_fac:.6f} ".format(diff_fac=gamma_portion)
            return H_str+"Position ({x}, {y}): {resol:.2f} Å".format(resol=resol, x=int(x), y=int(y))
        self.ax.format_coord = label_mouse_coords

    def _init_fig(self):
        self._update_status("Initializing the matplotlib figure...")
        """ initialize the mpl fig"""
        self.fig = plt.figure(1)
        self.ax = plt.Axes(self.fig, [0,0,1,1])
        self.fig.add_axes(self.ax)
        self.ax.set_axis_off()
        self.ax.set_aspect("equal")
        self.fig.set_size_inches([9.22, 3.8]) 
        self.gonio = []

    def _set_trace(self):
        import pdb
        pdb.set_trace()

    def _set_option_menu(self):
        """create an option menu for selecting params"""
        _options_frame = tk.Frame(self.master)
        _options_frame.pack(side=tk.TOP, expand=tk.NO)
        self.params_num.all_register_parent_frame(_options_frame)
        self.params_num.all_register_handler()
        self.params_cat.all_register_parent_frame(_options_frame)
        self.params_cat.all_register_handler()
        self.params_info.all_register_parent_frame(_options_frame)
        self.params_info.all_register_handler()

        def get_extra_dial_logic(dial_name):
            """set additional responses on updating specific numerical values"""
            if dial_name in ['energy', 'bandwidth']:
                return lambda: self.on_update_spectrum(init=True) # reinit SASE for new energy
            elif dial_name in ['diff_gamma', 'diff_sigma', 'diff_aniso']:
                return self.on_update_diffuse_params
            elif dial_name in ['ucell_scale_a', 'ucell_scale_b', 'ucell_scale_c']:
                axis = dial_name[-1]
                return lambda: self.on_update_ucell(axis)
            elif dial_name == 'domain_size':
                return self.on_update_domain_size
            elif dial_name == 'brightness':
                return self.on_update_normalization
            else:
                return self._generate_image_data
            return dial_specific_logic

        for dial_name in self.params_num.all_param_names:
            self.params_num.get_param(dial_name).generate_dial(extra_logic=get_extra_dial_logic(dial_name), hivis=self.hivis_mode)

        def get_extra_menu_logic(menu_name):
            """set responses for updating menu selections"""
            if menu_name == 'rotation_mode':
                return self.on_toggle_rotation_mode
            elif menu_name == 'spectrum_shape':
                return self.on_toggle_spectrum_shape
            elif menu_name == 'diffuse_mode':
                return self.on_toggle_diffuse_mode
            elif menu_name == 'pinned_mode':
                return self._display
            else:
                return self._generate_image_data

        for menu_name in self.params_cat.all_param_names:
            self.params_cat.get_param(menu_name).generate_menu(extra_logic=get_extra_menu_logic(menu_name), hivis=self.hivis_mode)

        for info_name in self.params_info.all_param_names:
            self.params_info.get_param(info_name).generate_info(hivis=self.hivis_mode)

        # Buttons
        self.new_pulse_button=Button(_options_frame, command=self.on_new_pulse, label="New XFEL pulse", position=(4,6), pos_hivis=(5,2), hivis=self.hivis_mode)
        self.randomize_orientation_button=Button(_options_frame, command=self._randomize_orientation, label="Randomize orientation", position=(2,4), pos_hivis=(7,0), hivis=self.hivis_mode)
        self.update_ref_image_button=Button(_options_frame, command=self._update_pinned, label="Update pinned image", position=(0,5), pos_hivis=(6,2), hivis=self.hivis_mode)
        self.reset_all_button=Button(_options_frame, command=self._reset_all, label="Reset all", position=(0,6), pos_hivis=(7,2), hivis=self.hivis_mode)
        #self.pdb_set_trace_button=Button(_options_frame, command=self._set_trace, label="Enter debugger", position=(5,6), pos_hivis=(7,3), hivis=self.hivis_mode)

        # Text Entry
        def validate(text):
            pdbid = text.strip()
            assert len(pdbid) == 4
            assert pdbid.isalnum()
        def fetch(pdbid):
            try:
                pdbfile = get_pdb(pdbid, "pdb", "rcsb", log=None, format="pdb")
                assert os.path.exists(pdbfile)
            except AssertionError:
                print("Could not fetch requested PDB")
                self._update_status("Failed to fetch requested PDB.")
                return
            try:
                self.on_update_pdb(pdbfile)
            except Exception:
                self._update_status("Failed to load requested PDB. Triclinic cells not yet supported.")
                return
        self.pdb_entry = TextEntry(_options_frame, command=fetch, validate_command=validate, label="PDB ID:", placeholder_text="4bs7", position=(5,1), pos_hivis=(0,3), hivis=self.hivis_mode, master=self)

    def on_toggle_rotation_mode(self, new_mode=None, update_selection=False, skip_gen_image_data=False):
        """enforce monochromatic beam, hide/show rotation specific params"""
        self._update_status("Changing experiment type...")
        if new_mode and update_selection:
            self.params_cat.rotation_mode.set_value(new_mode)
        if self.params_cat.rotation_mode.get_value() == 'Rotation':
            self.on_toggle_spectrum_shape(new_shape='Monochromatic', update_selection=True, skip_gen_image_data=True)
            self.params_cat.spectrum_shape.disable()
            self.params_num.delta_phi.enable()
            self.params_num.image.enable()
            self.gonio.append(self.ax.axvline(x=self.beam_center[0], ymin=0.01, ymax=0.1, color='r'))
            self.gonio.append(self.ax.annotate('  Goniometer axis', self.beam_center, xytext=(self.beam_center[0], self.beam_center[1]+10), color='r'))
        else:
            self.params_cat.spectrum_shape.enable()
            self.params_num.delta_phi.disable()
            self.params_num.image.disable()
            for g in self.gonio:
                g.remove()
                del g
            self.gonio = []
        if not skip_gen_image_data:
            self._generate_image_data()

    def on_toggle_spectrum_shape(self, new_shape=None, update_selection=False, skip_gen_image_data=False):
        """set visibility of other controls dependent on spectrum shape"""
        self._update_status("Changing spectrum type...")
        if new_shape and update_selection:
            self.params_cat.spectrum_shape.set_value(new_shape)
        shape = self.params_cat.spectrum_shape.get_value()
        if shape == 'SASE (XFEL)':
            self.new_pulse_button.enable()
            self.params_num.bandwidth.disable()
        elif shape == 'Gaussian':
            self.new_pulse_button.disable()
            self.params_num.bandwidth.enable()
        elif shape == 'Monochromatic':
            self.new_pulse_button.disable()
            self.params_num.bandwidth.disable()
        self.on_update_spectrum(new_pulse=(new_shape is not None), skip_gen_image_data=skip_gen_image_data)

    def on_toggle_diffuse_mode(self, skip_gen_image_data=False):
        """set visibility of diffuse mode params"""
        self._update_status("Toggling diffuse scattering...")
        if self.params_cat.diffuse_mode.get_value() == 'On':
            for param in [self.params_num.diff_gamma, self.params_num.diff_sigma, self.params_num.diff_aniso]:
                param.enable()
            domains = self.params_hyper.mosaic_domains_diffuse
            spread = 1 if domains > 1 else 0
            self.on_update_diffuse_params(skip_gen_image_data=True)
        else:
            for param in [self.params_num.diff_gamma, self.params_num.diff_sigma, self.params_num.diff_aniso]:
                param.disable()
            domains = self.params_hyper.mosaic_domains_bragg
            spread = 1 if domains > 1 else 0
        # update mosaic domains accd. to hyperparameters, distinct between diffuse and Bragg-only mode
        for SIM in (self.SIM, self.SIM_noSF):
            SIM.crystal.n_mos_domains = domains
            SIM.crystal.mos_spread_deg = spread
            SIM.instantiate_diffBragg(oversample=self.params_hyper.oversampling, device_Id=0, default_F=0)
        self.fix_diffBragg_instances()
        if not skip_gen_image_data:
            self._generate_image_data()

    def on_update_spectrum(self, new_pulse=False, init=False, skip_gen_image_data=False):
        self._update_status("Updating spectrum...")
        shape = self.params_cat.spectrum_shape.get_value()
        if shape == "Gaussian":
            bw = 0.01 * self.params_num.bandwidth.get_value() * self.params_num.energy.get_value() # bandwidth in eV
            gfunc = gaussian.term(1, 4 * math.log(2)/(bw**2)) # FWHM of bw, mu == 0
            self.spectrum_eV = [(energy + self.params_num.energy.get_value(), 1e12 * gfunc.at_x(energy)) \
                                for energy in range(-50,51)]
            self.spectrum_Ang = [(12398./energy, flux) for (energy, flux) in self.spectrum_eV]
        elif shape == "SASE (XFEL)":
            if init or not hasattr(self, "SASE_iter"):
                self.SASE_iter = self.SASE_sim.generate_recast_renormalized_images(
                    energy=self.params_num.energy.get_value(), total_flux=1e12)
            if new_pulse or init or not hasattr(self, "pulse_energies_Ang"):
                self.pulse_energies_Ang, self.flux_list, self.avg_wavelength_Ang = next(self.SASE_iter)
            self.spectrum_Ang = list(zip(self.pulse_energies_Ang, self.flux_list))
        elif shape == "Monochromatic":
            self.spectrum_Ang = [(12398./self.params_num.energy.get_value(), 1e12)] # single wavelength for computational speed
        else:
            raise NotImplementedError("Haven't implemented a spectrum of the requested shape {}".format(shape))
        if not skip_gen_image_data:
            self._generate_image_data()

    def on_new_pulse(self, tkevent=None):
        self.on_update_spectrum(new_pulse=True)

    def on_update_diffuse_params(self, skip_gen_image_data=False):
        """update stored intermediate values"""
        gamma = self.params_num.diff_gamma.get_value()
        sigma = self.params_num.diff_sigma.get_value()
        aniso = self.params_num.diff_aniso.get_value()
        self.diffuse_gamma = (gamma * aniso,
                              gamma,
                              gamma * aniso)
        self.diffuse_sigma = (sigma * sigma * aniso,
                              sigma * sigma,
                              sigma * sigma)

        if not skip_gen_image_data:
            self._generate_image_data()

    def _update_pinned(self, _press=None):
        self._generate_image_data(update_ref=True)

    def _check_symmetry(self):
        """dependent on crystal symmetry: test if axes can be independently adjusted"""
        a, b, c, al, be, ga = self.ucell
        test_c = (a, b, c+10, al, be, ga)
        test_b = (a, b+10, c, al, be, ga)
        self.b_can_scale = True if self.sg.is_compatible_unit_cell(unit_cell(test_b)) else False
        self.c_can_scale = True if self.sg.is_compatible_unit_cell(unit_cell(test_c)) else False
        # if a,b,c remaining: all vary independently
        # if a,b remaining: c scales with a
        # if a,c remaining: b scales with a
        # if a remaining: a,b,c all vary together

    def _enforce_symmetry_on_controls(self):
        """disable ucell params that cannot be independently adjusted"""
        if not self.b_can_scale:
            self.params_num.ucell_scale_b.disable()
        else:
            self.params_num.ucell_scale_b.enable()
        if not self.c_can_scale:
            self.params_num.ucell_scale_c.disable()
        else:
            self.params_num.ucell_scale_c.enable()

    def _update_ucell_label(self, update_sg=False):
        """label current model unit cell and space group"""
        ucell_str = "Unit cell: (%6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f)" % self.scaled_ucell
        self.params_info.ucell.set_value(ucell_str)
        if update_sg:
            sg_str = "Space group: %s" % self.sg.info().symbol_and_number()
            self.params_info.sg.set_value(sg_str)

    def on_update_normalization(self, display=True):
        """update normalization"""
        exponent = self.params_num.brightness.get_value()*2-2
        self.percentile = 100 - 10**exponent
        self._normalize_all_image_data(display=display)

    def _normalize_image_data(self, img_data):
        """scale data to [0,1] where variable %ile intensities and above are set to 1."""
        scale = 1./max(1e-50,np.percentile(img_data, self.percentile))
        scaled_data = img_data * scale
        scaled_data[scaled_data > 1] = 1
        return scaled_data

    def _overloads_image_data(self, img_data):
        """Create an image where values above the maximum in the normalization range are 1 and all other values are 0"""
        scale = 1./max(1e-50,np.percentile(img_data, self.percentile))
        scaled_data = img_data * scale
        scaled_data[scaled_data <= 1] = 0
        scaled_data[scaled_data > 1] = 1
        return scaled_data

    def _normalize_all_image_data(self, display=True):
        """update intensity scaling for the stored pixel data"""
        self._update_status("Normalizing image data...")
        t = time.time()
        self.img_overlay[:,:,0] = self._normalize_image_data(self.img_ref[0]) # red channel
        self.img_overlay[:,:,1] = self._normalize_image_data(self.img_sim[0] + self.img_ref[0]) # green channel (grayscale if identical)
        self.img_overlay[:,:,2] = self._normalize_image_data(self.img_sim[0]) # blue channel
        self.img_single_channel[:,:,0] = 0
        self.img_single_channel[:,:,1] = self._normalize_image_data(self.img_sim[0]) # green channel
        self.img_single_channel[:,:,2] = self._normalize_image_data(self.img_sim[0]) # blue channel
        # Flag overloads using yellow pixels
        img_overloads = self._overloads_image_data(self.img_sim[0])
        self.img_single_channel[:,:,0][img_overloads == 1] = 1
        self.img_single_channel[:,:,1][img_overloads == 1] = 1
        self.img_single_channel[:,:,2][img_overloads == 1] = 0
        t = time.time()-t
        print("Time taken to normalize image data: %8.5f seconds"% t)
        self._update_status("Ready.")
        if display:
            self._display()

    def _generate_image_data(self, update_ref=False, display=True):
        """generate image to match requested params"""
        self._update_status("Generating image data...")
        t = time.time()
        SIM = self.SIM if self.params_cat.Fhkl.get_value() == 'On' else self.SIM_noSF
        diffuse_gamma = self.diffuse_gamma if self.params_cat.diffuse_mode.get_value() == 'On' else None
        diffuse_sigma = self.diffuse_sigma if self.params_cat.diffuse_mode.get_value() == 'On' else None
        if self.params_cat.rotation_mode.get_value() == 'Rotation':
            delta_phi = self.params_num.delta_phi.get_value()
            phi_n_steps = self.params_hyper.oscillation_n_steps
            phi_step = delta_phi/phi_n_steps
            pix = sweep(SIM,
                delta_phi * (self.params_num.image.get_value() - 1), # phi_start
                phi_step, # phi step for simtbx
                delta_phi, # phi range summmed in one image
                self.pfs, self.scaled_ucell,
                tuple(self.ncells),
                (0,
                0,
                self.params_num.rot_z.get_value()*math.pi/180.),
                spectrum=self.spectrum_Ang,
                eta_p=self.params_num.mos_ang_deg.get_value(),
                diffuse_gamma=diffuse_gamma,
                diffuse_sigma=diffuse_sigma,
                oversample=self.params_hyper.oversampling)
        else:
            pix = run_simdata(SIM, self.pfs, self.scaled_ucell,
                tuple(self.ncells),
                (self.params_num.rot_x.get_value()*math.pi/180.,
                self.params_num.rot_y.get_value()*math.pi/180.,
                self.params_num.rot_z.get_value()*math.pi/180.),
                spectrum=self.spectrum_Ang,
                eta_p=self.params_num.mos_ang_deg.get_value(),
                diffuse_gamma=diffuse_gamma,
                diffuse_sigma=diffuse_sigma)
        t = time.time()-t
        self.img_sim[self.pan, self.slow, self.fast] = pix
        if update_ref:
            self.img_ref[self.pan, self.slow, self.fast] = pix
        print("Time taken to update image data: %8.5f seconds"% t)
        if display:
          self._normalize_all_image_data()

    def on_update_domain_size(self, skip_gen_image_data=False):
        """given a target domain size (one side length), determine number of cells to use along each axis"""
        side_length = self.params_num.domain_size.get_value()
        a, b, c, _, _, _ = self.scaled_ucell
        self.ncells = []
        for side in (a, b, c):
            self.ncells.append(max(int(round(side_length/side)), 1))
        if not skip_gen_image_data:
            self._generate_image_data()

    def on_update_ucell(self, axis):
        """scale one or more lengths depending on symmetry"""
        scale = self.params_num.get_param(f'ucell_scale_{axis}'.format()).get_value()
        if axis == "a":
            a = scale * self.ucell[0]
            b = self.scaled_ucell[1] if self.b_can_scale else scale * self.ucell[1]
            c = self.scaled_ucell[2] if self.c_can_scale else scale * self.ucell[2]
        elif axis == "b":
            a = self.scaled_ucell[0]
            b = scale * self.ucell[1]
            c = self.scaled_ucell[2] if self.c_can_scale else scale * self.ucell[2]
        elif axis == "c":
            a = self.scaled_ucell[0]
            b = self.scaled_ucell[1]
            c = scale * self.ucell[2]
        self.scaled_ucell = (a,b,c,*self.ucell[3:6])
        self._update_ucell_label()
        self._generate_image_data()

    def on_update_orientation(self, x, y, z):
        """update U matrix with changes to orientation quaternion"""
        self._update_status("Updating orientation...")
        self.params_info.umatrix.set_value(f"U0 angles: ({x:6.2f}, {y:6.2f}, {z:6.2f})".format())

    def _randomize_orientation(self, _press=None, reset=False):
        ori = randomize_orientation(self.SIM, track_with=self.SIM_noSF, reset=reset)
        self.on_update_orientation(*ori.r3_rotation_matrix_as_x_y_z_angles(deg=True))
        self._generate_image_data(update_ref=True)

    def _init_display(self):
        """initialize the display"""
        self.aximg = self.ax.imshow(self.img_single_channel)
        self._generate_image_data(update_ref=True, display=False)
        self.on_update_normalization(display=False)
        #self.canvas.draw()
        #self._annotate()

    def _display(self):
        """display the current image"""
        if self.params_cat.pinned_mode.get_value() == 'Overlay with pinned':
            imgdata = self.img_overlay
        else:
            imgdata = self.img_single_channel
        self.aximg.set_data(imgdata)
        self.canvas.draw()
        self._annotate()
        self._update_status("Ready.")

    def bind(self):
        """key bindings"""
        # increment or decrement by big steps iff Shift key is held
        self.master.bind_all("<KeyPress-Shift_L>", self._set_big_steps_on_shift)
        self.master.bind_all("<KeyPress-Shift_R>", self._set_big_steps_on_shift)
        self.master.bind_all("<KeyRelease>", self._set_small_steps_on_release_shift)
        #self.master.bind_all("<R>", self._reset)  # reset this dial to default value
        #self.master.bind_all("<Shift-R>", self._reset_all) # reset all controls to default settings
        #self.master.bind_all("<space>", self.on_new_pulse) # repeat the simulation with a new spectrum

        self.master.bind_all("<Return>", self._register_change) # Register an updated value typed into a spinbox
        self.master.bind_all("<KP_Enter>", self._register_change) # Also recognize numpad Enter key

        #self.master.bind_all("<I>", self._toggle_image_mode) # toggle displaying pinned image
        #self.master.bind_all("<F>", self._toggle_Fhkl) # toggle using Fhkl to scale intensities
        #self.master.bind_all("<D>", self._toggle_diffuse_scattering) # toggle diffuses scattering on/off
        #self.master.bind_all("<O>", self._randomize_orientation) # randomize crystal orientation
        #self.master.bind_all("<U>", self._update_pinned) # update pinned image (in red)

    def _register_change(self, tkevent):
        if self.pdb_entry.is_active:
            self.pdb_entry.command()
        else:
            self.params_num.current_param.command()

    def _update_dial(self, new_dial_name):
        self.params_num.set_active_param_by_name(new_dial_name)
        self._display()

    def _set_big_steps_on_shift(self, tkevent):
        #print("registered Shift keypress")
        self.params_num.set_big_steps()

    def _set_small_steps_on_release_shift(self, tkevent):
        if tkevent.keysym in ['Caps_Lock', 'Shift_L', 'Shift_R']:
            #print("Shift key released")
            self.params_num.set_small_steps()

    def _reset(self, tkevent=None):
        self._update_status("Resetting parameter(s)...")
        self.params_num.current_param.reset()

    def _reset_all(self, tkevent=None):
        self.params_cat.reset_all()
        self._randomize_orientation(reset=True)
        self.on_toggle_diffuse_mode(skip_gen_image_data=True)
        self.on_toggle_rotation_mode(skip_gen_image_data=True)
        self.on_toggle_spectrum_shape(skip_gen_image_data=True)
        self.params_num.reset_all()
        self._update_ucell_label()
        self._generate_image_data()

if __name__ == '__main__':
    import sys
    if "-h" in sys.argv or "--help" in sys.argv:
        print(help_message)
        exit()
    if "--config" in sys.argv:
        sys.argv.append("-c") # to recognize --config
    if "-c" in sys.argv:
        # default -c -e5 -a2 for this program
        if not [arg for arg in sys.argv if arg.startswith("-e")]:
            sys.argv.append("-e5")
        if not [arg for arg in sys.argv if arg.startswith("-a")]:
            sys.argv.append("-a2")
    if "--fetch" in sys.argv:
        pos = sys.argv.index("--fetch")
        sys.argv.pop(pos)
        pdbid = sys.argv.pop(pos)
        try:
            pdbfile = get_pdb(pdbid, "pdb", "rcsb", log=sys.stdout, format="pdb")
        except Exception:
            print("Couldn't fetch pdb {}. Try fetching with iotbx.fetch_pdb.".format(pdbid))
            exit()
    elif "--file" in sys.argv:
        pos = sys.argv.index("--file")
        sys.argv.pop(pos)
        pdbfile = sys.argv.pop(pos)
    else:
        pdbfile = libtbx.env.find_in_repositories(
            relative_path="sim_erice/4bs7.pdb",
            test=os.path.isfile)
    if not pdbfile:
        print("Could not load model file. Please supply a valid PDB model with the --file flag or a PDB ID with the --fetch flag.")
        exit()
    parser = ArgumentParser(
        usage=usage,
        phil=sim_view_phil_scope
        )
    params_hyper, options = parser.parse_args(args=sys.argv[1:],
                                              show_diff_phil=True)
    params_cat.spectrum_shape.default = params_hyper.spectrum_shape

    if params_hyper.context == "kokkos":
        os.environ["DIFFBRAGG_USE_KOKKOS"]="1"
    elif "DIFFBRAGG_USE_KOKKOS" in os.environ:
        del os.environ["DIFFBRAGG_USE_KOKKOS"]
    if params_hyper.context == "cuda":
        os.environ["DIFFBRAGG_USE_CUDA"]="1"
    elif "DIFFBRAGG_USE_CUDA" in os.environ:
        del os.environ["DIFFBRAGG_USE_CUDA"]

    from simtbx.diffBragg.device import DeviceWrapper
    with DeviceWrapper(0) as _:
        root = tk.Tk()
        def _close_window():
            root.quit()
            root.destroy()
        root.protocol("WM_DELETE_WINDOW", _close_window)
        
        root.title("SimView")
        
        root.geometry('1920x1140')
        
        frame = SimView(root, params_num, params_cat, params_info, params_hyper, pdbfile)
        
        frame.pack( side=tk.TOP, expand=tk.NO)
        root.mainloop()
#

