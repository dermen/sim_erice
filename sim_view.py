from __future__ import absolute_import, division, print_function
"""
Plot the PNG files
Derek Mendez, dermen@lbl.gov

Generate and load simulated images on the fly
Iris Young, idyoung@lbl.gov
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

from matplotlib.backends.backend_tkagg import \
    FigureCanvasTkAgg, NavigationToolbar2Tk

import libtbx.load_env
from on_the_fly_simdata import run_simdata, get_SIM
from simtbx.nanoBragg.tst_nanoBragg_multipanel import beam, whole_det
from simtbx.diffBragg import hopper_utils
from iotbx.crystal_symmetry_from_any import extract_from as extract_symmetry_from
from dxtbx_model_ext import Crystal
from scitbx import matrix
from dials.array_family import flex
import time

help_message="""SimView: lightweight simulator and viewer for diffraction still images.

Run without any arguments, this program simulates diffraction of a small
lysozyme crystal. You may provide a different pdb file on the command line
instead. The diffraction will be simulated with a set of default parameters
for mosaic block size, beam energy and bandpass, etc. that you can change by
adjusting each of the dials. Parameters to be adjusted can be selected either
with the drop-down menu or with keyboard shortcuts as follows:

Left arrow or p:    previous parameter
Right arrow or n:   next parameter
Up arrow:           increase this parameter
Down arrow:         decrease this parameter
Shift-up arrow:     increase this parameter a lot
Shift-down arrow:   decrease this parameter a lot
Space bar:          simulate a new stochastic XFEL pulse
r:                  reset all parameters to defaults

(Note that matplotlib binds the left arrow key to resetting the field of view,
so both this effect and selection of the previous parameter will take place
if you are zoomed in on an area. The right arrow key can restore the zoom.)

Explanation of parameters:

DomainSize: the small blocks within crystals that are internally perfectly
aligned and identical are called mosaic domains. The average size of a domain
is one contributing factor determining the sharpness of the Bragg peaks. As
the size of the mosaic domain increases, the effect of constructive
interference of the scattered X-rays gets stronger, and the peak sharpens.

MosAngDeg: the mosaic domains are misoriented relative to one another by an
average angle denoted the mosaic angle, displayed here in degrees. Larger
differences in the orientations of the domains produce contributions to
Bragg peaks over a wider area on the detector, resulting in more diffuse
spots.

a, b, c: the parameters of the crystal unit cell. Depending on the symmetry
of the crystal, these may all vary independently or may be constrained to
scale together. Larger unit cells produce constructive and destructive
interference at smaller angles, resulting in more closely spaced Bragg
peaks on the detector.

Missetting angles in degrees: differences between the true orientation of
the crystal and what is determined during indexing. A misorientation of
the X or Y axis means a misestimation of which spots are in diffracting
conditions, and will affect which spots (or how much of the spots) appear
on the detector. Since the Z axis is along the beam direction,
misorientation in Z results in rotation of the entire diffraction pattern.

Energy/Bandwidth: the energy of the X-ray beam denotes the average energy
of an individual pulse. For X-ray free electron laser (XFEL) pulses
generated by self-amplified spontaneous emission (SASE), each pulse is
spiky and stochastic, with a typical bandwidth of between 10 and 30 eV.
A monochromater may be used to narrow this to around 1 eV.

Fhkl: when we observe diffraction from a single particle, we are seeing
the Fourier transform of the particle, with [a great deal of] added
noise. For many copies of a particle all in the same orientation, we can
expect the signal to be much stronger. If these aligned particles are
also periodically spaced, as in the context of a crystal, they act as a
diffraction grating, and we see this signal only at positions of
constructive interference (the reciprocal lattice points in diffracting
conditions, i.e. Bragg peaks). In order to better see the effects of the
above parameters, we can ignore the Fourier transform of the asymmetric
unit and pretend all of the Bragg peaks have uniform intensity. This is
what we are simulating when we display "Fhkl not accounted for."
"""

class SimView(tk.Frame):

    def __init__(self, master, params, pdbfile, *args, **kwargs):
        tk.Frame.__init__(self, *args, **kwargs)

        self.master = master
        self.params = params
        self.dial_names = list(self.params.keys())
        self.current_dial = self.dial_names[0]

        symmetry = extract_symmetry_from(pdbfile)
        sg = str(symmetry.space_group_info())
        ucell = symmetry.unit_cell()
        fmat = matrix.sqr(ucell.fractionalization_matrix())
        cryst = Crystal(fmat, sg)
        self.SIM = get_SIM(whole_det, beam, cryst, pdbfile)
        self.SIM_noSF = get_SIM(whole_det, beam, cryst)
        self.ucell = self.SIM.crystal.dxtbx_crystal.get_unit_cell().parameters()
        self.scaled_ucell = self.ucell
        self._load_params_only()
        self._update_spectrum(new_pulse=True)

        fsize, ssize = whole_det[0].get_image_size()
        img_sh = 1,ssize, fsize
        self.pfs = hopper_utils.full_img_pfs(img_sh)
        self.fast = self.pfs[1::3].as_numpy_array()
        self.slow = self.pfs[2::3].as_numpy_array()
        self.pan = self.pfs[0::3].as_numpy_array()
        self.img = np.zeros((1,ssize, fsize))

        self._set_option_menu()
        self._make_master_label()

        self._init_fig()

        self._pack_canvas()
        self._reset()
        self._display(init=True)

        self.bind()

        # set the option menu first
        self.dial_choice.set(self.current_dial)


    def _pack_canvas(self):
        """ embed the mpl figure"""
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.master) 
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP,fill=tk.BOTH,
            expand=tk.YES)

        self.toolbar = NavigationToolbar2Tk(self.canvas, self.master)
        self.toolbar.update()
        self.canvas.get_tk_widget().pack(side=tk.TOP, 
            fill=tk.BOTH, expand=tk.YES)


    def _init_fig(self):
        """ initialize the mpl fig"""
        self.fig = plt.figure(1)
        self.ax = plt.Axes(self.fig, [0,0,1,1])
        self.fig.add_axes(self.ax)
        self.ax.set_axis_off()
        self.ax.set_aspect("equal")
        self.fig.set_size_inches([9.22, 3.8]) 

    def _set_option_menu(self):
        """create an option menu for selecting params"""
        self.dial_choice = tk.StringVar()
       
        _opt_frame = tk.Frame(self.master)
        _opt_frame.pack(side=tk.TOP, expand=tk.NO)
        
        tk.Label( _opt_frame, text="Adjusting variable: ", 
            font=("Helvetica",18), width=18)\
            .pack(side=tk.LEFT)

        self.param_opt_menu = tk.OptionMenu(_opt_frame, 
            self.dial_choice,
            *self.dial_names,
            command=self._update_dial)
        self.param_opt_menu.pack(side=tk.LEFT)
        self.param_opt_menu.config(width=12)
        self.param_opt_menu.config(font=("Helvetica", 18))

        b=tk.Button(_opt_frame, command=self._reset, text="Reset")
        b.pack(side=tk.LEFT)
        b.config(font=("Helvetica", 18))
        b.config(width=6)

    def _update_dial(self, new_dial):
        if new_dial != self.current_dial:
            self.current_dial = new_dial
            self.dial_choice.set(self.current_dial)
            self._display()

    def _load_params_only(self):
        """load params available only (don't pre-load images). Keys are filenames."""
        self._VALUES = {}
        self._LABELS = {}
        for dial, (_, _, _, _, default) in params.items():
            self._VALUES[dial] = default
            self._LABELS[dial] = self._get_new_label_part(dial, default)

    def _generate_image_data(self):
        """generate image to match requested params"""
        # t = time.time()
        SIM = self.SIM if self._VALUES["Fhkl"] else self.SIM_noSF
        pix = run_simdata(SIM, self.pfs, self.scaled_ucell,
            tuple([(x,x,x) for x in [self._VALUES["MosDom"]]][0]),
            (self._VALUES["RotX"]*math.pi/180.,
            self._VALUES["RotY"]*math.pi/180.,
            self._VALUES["RotZ"]*math.pi/180.),
            spectrum=self.spectrum,
            eta_p=self._VALUES["MosAngDeg"])
        # t = time.time()-t
        self.img[self.pan, self.slow, self.fast] = pix
        self.image = self.img[0]

    def _update_spectrum(self, new_pulse=False):
        if new_pulse:
            self.pulse = [l for l in np.random.normal(0, 1, 10)]
        energy = 12398./self._VALUES["Energy"]
        bw = 0.01*self._VALUES["Bandwidth"]
        scaled_pulse = [l * bw * energy + energy for l in self.pulse]
        self.spectrum = [(energy, 1e12) for energy in scaled_pulse]

    def _update_ucell(self, new_value):
        a,b,c = [x*new_value for x in self.ucell[:3]]
        self.scaled_ucell = (a,b,c,*self.ucell[3:6])

    def _make_master_label(self):
        """label that will be updated for each image"""
        self.master_label = tk.Label(self.master, text=\
"""DomainSize: ____; MosAngleDeg: ____; a,b,c = ____;
Missetting angles in degrees (X,Y,Z) = (____, ____}, ____);
Energy/Bandwidth= ____ / ____; ____""", font="Helvetica 15", width=350)
        self.master_label.pack(side=tk.TOP, expand=tk.NO)

    def _update_label(self):
        self._label = """DomainSize: {mosdom}; MosAngleDeg: {mosang}; a,b,c = {ucell};
Missetting angles in degrees (X,Y,Z) = ({rotx}, {roty}, {rotz});
Energy/Bandwidth= {energy} / {bw}; {Fhkl}""".format(
            mosdom=self._LABELS["MosDom"],
            mosang=self._LABELS["MosAngDeg"],
            ucell=self._LABELS["ucell_scale"],
            energy=self._LABELS["Energy"],
            bw=self._LABELS["Bandwidth"],
            rotx=self._LABELS["RotX"],
            roty=self._LABELS["RotY"],
            rotz=self._LABELS["RotZ"],
            Fhkl=self._LABELS["Fhkl"]
        )

    def _get_new_label_part(self, dial, new_value):
        if dial == "MosDom":
            return "{v}x{v}x{v}".format(v=new_value)
        elif dial == "MosAngDeg":
            return "{:.2f}º".format(new_value)
        elif dial == "ucell_scale":
            a,b,c = self.scaled_ucell[:3]
            return "{a:.2f}, {b:.2f}, {c:.2f}".format(a=a, b=b, c=c)
        elif dial == "Energy":
            return "{:d}".format(new_value)
        elif dial == "Bandwidth":
            return "{:.2f}%".format(new_value)
        elif dial in ["RotX", "RotY", "RotZ"]:
            return "{:+.2f}".format(new_value)
        elif dial == "Fhkl":
            return "Structure factors {}accounted for".format("" if new_value else "not ")

    def _display(self, init=False):
        """display the current image"""

        if init:
            self.aximg = self.ax.imshow(self.image)
        else:
            self.aximg.set_data(self.image)

        self._update_label()
        self.master_label.config(text=self._label)
        self.master_label.config(font=("Courier", 15))

        self.canvas.draw()

    def bind(self):
        """key bindings"""
        self.master.bind_all("<Up>", self._small_step_up)  # increase this parameter
        self.master.bind_all("<Shift-Up>", self._big_step_up) # increase this parameter a lot
        self.master.bind_all("<r>", self._reset)  # reset all params to default settings
        self.master.bind_all("<Down>", self._small_step_down)  # decrease this parameter
        self.master.bind_all("<Shift-Down>", self._big_step_down) # decrease this parameter a lot
        self.master.bind_all("<space>", self._new_pulse) # repeat the simulation with a new spectrum
        
        self.master.bind_all("<Left>", self._prev_dial)
        self.master.bind_all("<p>", self._prev_dial)
        self.master.bind_all("<Right>", self._next_dial)
        self.master.bind_all("<n>", self._next_dial)

    def _next_dial(self, tkevent):
        try:
            new_dial = self.dial_names[self.dial_names.index(self.current_dial) + 1]
            self._update_dial(new_dial)
        except IndexError:
            pass

    def _prev_dial(self, tkevent):
        try:
            new_dial = self.dial_names[self.dial_names.index(self.current_dial) - 1]
            self._update_dial(new_dial)
        except IndexError:
            pass

    def _new_pulse(self, tkevent):
        self._update_spectrum(new_pulse=True)
        self._generate_image_data()
        self._display()

    def _set_new_value(self, dial, new_value):
        self._VALUES[dial] = new_value
        self._LABELS[dial] = sef._get_new_label_part(dial, new_value)
        if dial in ["Energy", "Bandwidth"]:
            self._update_spectrum()
        elif dial == "ucell_scale":
            self._update_ucell(new_value)
        self._generate_image_data()
        self._display()

    def _small_step_up(self, tkevent):
        _, this_max, this_step, _, _ = self.params[self.current_dial]
        this_value = self._VALUES[self.current_dial]
        new_value = this_value + this_step
        if new_value <= this_max:
            self._set_new_value(self.current_dial, new_value)

    def _big_step_up(self, tkevent):
        _, this_max, _, this_step, _ = self.params[self.current_dial]
        this_value = self._VALUES[self.current_dial]
        new_value = this_value + this_step
        if new_value <= this_max:
            self._set_new_value(self.current_dial, new_value)

    def _small_step_down(self, tkevent):
        this_min, _, this_step, _, _ = self.params[self.current_dial]
        this_value = self._VALUES[self.current_dial]
        new_value = this_value - this_step
        if new_value >= this_min:
            self._set_new_value(self.current_dial, new_value)

    def _big_step_down(self, tkevent):
        this_min, _, _, this_step, _ = self.params[self.current_dial]
        this_value = self._VALUES[self.current_dial]
        new_value = this_value - this_step
        if new_value >= this_min:
            self._set_new_value(self.current_dial, new_value)

    def _reset(self, _press=None):
        for dial in self.dial_names:
            default_value = self.params[dial][4]
            self._VALUES[dial] = default_value
            self._LABELS[dial] = self._get_new_label_part(dial, default_value)
        self._generate_image_data()
        self._display(init=True)

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        if "-h" in sys.argv:
            print(help_message)
            exit()
        pdbfile = sys.argv[1]
    else:
        pdbfile = libtbx.env.find_in_repositories(
            relative_path="simtbx/sim_view/data/4bs7.pdb",
            test=os.path.isfile)
        if not pdbfile:
            print("Could not load default model file. Please supply one on the command line.")
            exit()

    # params stored as: [min, max, small_step, big_step, default]
    params = {
        "MosDom":[6, 200, 2, 10, 10],
        "MosAngDeg":[0.0, 1.5, 0.3, 1.0, 0.3001],
        "ucell_scale":[0.9, 1.1, 0.05, 0.1, 1],
        "Energy":[6720, 7180, 10, 30, 6750],
        "Bandwidth":[0, 6, 0.1, 1, 3],
        "RotX": [-180, 180, 0.01, 0.1, 0],
        "RotY": [-180, 180, 0.01, 0.1, 0],
        "RotZ": [-180, 180, 0.01, 0.1, 0],
        "Fhkl":[0, 1, 1, 1, 1]} # binary switch

    root = tk.Tk()
    root.title("SimView")

    root.geometry('940x500')
    #root.resizable(0,0)
    
    frame = SimView(root, params, pdbfile)
    
    frame.pack( side=tk.TOP, expand=tk.YES)
    root.mainloop()
#

