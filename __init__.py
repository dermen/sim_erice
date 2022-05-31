import libtbx.load_env
srcdir = libtbx.env.dist_path("sim_erice")
if srcdir is None:
    print("WARNING: please configure sim_erice: libtbx.configure sim_erice")
    exit()
import os.path
modules_dir = os.path.dirname(os.path.abspath(srcdir))
sim_erice = os.path.join(modules_dir,"sim_erice")
