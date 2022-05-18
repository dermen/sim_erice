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

from on_the_fly_simdata import run_simdata, get_SIM
from simtbx.nanoBragg.tst_nanoBragg_multipanel import beam, whole_det
from simtbx.diffBragg import hopper_utils
from iotbx.crystal_symmetry_from_any import extract_from as extract_symmetry_from
from dxtbx_model_ext import Crystal
from scitbx import matrix
from dials.array_family import flex
import time

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
        for dial, (_, _, _, default) in params.items():
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
"""a,b,c = ____; DomainSize: ____; MosAngleDeg: ____;
Missetting angles in degrees (X,Y,Z) = (____, ____}, ____);
Energy/Bandwidth= ____ / ____; ____""", font="Helvetica 15", width=350)
        self.master_label.pack(side=tk.TOP, expand=tk.NO)

    def _update_label(self):
        self._label = """a,b,c = {ucell}; DomainSize: {mosdom}; MosAngleDeg: {mosang};
Missetting angles in degrees (X,Y,Z) = ({rotx}, {roty}, {rotz});
Energy/Bandwidth= {energy} / {bw}; {Fhkl}""".format(
            ucell=self._LABELS["ucell_scale"],
            mosdom=self._LABELS["MosDom"],
            mosang=self._LABELS["MosAngDeg"],
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
        self.master.bind_all("<Up>", self._next)  # next image
        self.master.bind_all("<t>", self._reset)  # default image
        self.master.bind_all("<Down>", self._prev)  # prev image
        self.master.bind_all("<space>", self._new_pulse) # repeat the simulation with a new spectrum
        
        self.master.bind_all("<Left>", self._prev_dial)
        self.master.bind_all("<Right>", self._next_dial)

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

    def _next(self, tkevent):
        _, this_max, this_step, _ = self.params[self.current_dial]
        this_value = self._VALUES[self.current_dial]
        new_value = this_value + this_step
        if new_value <= this_max:
            self._VALUES[self.current_dial] = new_value
            self._LABELS[self.current_dial] = self._get_new_label_part(self.current_dial, new_value)
            if self.current_dial in ["Energy", "Bandwidth"]:
                self._update_spectrum()
            elif self.current_dial == "ucell_scale":
                self._update_ucell(new_value)
            self._generate_image_data()
            self._display()

    def _prev(self, tkevent):
        this_min, _, this_step, _ = self.params[self.current_dial]
        this_value = self._VALUES[self.current_dial]
        new_value = this_value - this_step
        if new_value >= this_min:
            self._VALUES[self.current_dial] = new_value
            self._LABELS[self.current_dial] = self._get_new_label_part(self.current_dial, new_value)
            if self.current_dial in ["Energy", "Bandwidth"]:
                self._update_spectrum()
            elif self.current_dial == "ucell_scale":
                self._update_ucell(new_value)
            self._generate_image_data()
            self._display()

    def _reset(self, _press=None):
        for dial in self.dial_names:
            default_value = self.params[dial][3]
            self._VALUES[dial] = default_value
            self._LABELS[dial] = self._get_new_label_part(dial, default_value)
        self._generate_image_data()
        self._display(init=True)

if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        pdbfile = sys.argv[1]
    else:
        pdbfile = "4bs7.pdb"

    # params stored as: [min, max, step, default]
    params = {
        "MosDom":[6, 30, 2, 10],
        "MosAngDeg":[0.0, 1.5, 0.3, 0.3],
        "ucell_scale":[0.9, 1.1, 0.05, 1],
        "Energy":[6720, 7180, 30, 6750],
        "Bandwidth":[0, 6, 0.5, 3],
        "RotX": [-0.5, 0.5, 0.01, 0],
        "RotY": [-0.5, 0.5, 0.01, 0],
        "RotZ": [-0.5, 0.5, 0.01, 0],
        "Fhkl":[0, 1, 1, 1]} # binary switch

    root = tk.Tk()
    root.title("SimView")

    root.geometry('940x500')
    #root.resizable(0,0)
    
    frame = SimView(root, params, pdbfile)
    
    frame.pack( side=tk.TOP, expand=tk.YES)
    root.mainloop()
#

