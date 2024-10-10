
import libtbx.load_env
cctbx_src = libtbx.env.dist_path("cctbx")
# Only do this if in a bootstrap env:
if cctbx_src.endswith('modules/cctbx_project/cctbx'):
    srcdir = libtbx.env.dist_path("sim_erice")
    if srcdir is None:
        print("WARNING: please configure sim_erice: libtbx.configure sim_erice")
        exit()
    import os.path
    modules_dir = os.path.dirname(os.path.abspath(srcdir))
    sim_erice = os.path.join(modules_dir,"sim_erice")
