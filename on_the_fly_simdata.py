
from simtbx.diffBragg import utils
from simtbx.nanoBragg import nanoBragg_crystal, nanoBragg_beam, sim_data
from iotbx.crystal_symmetry_from_any import extract_from as extract_symmetry_from
import scitbx
from scitbx.matrix import col, sqr
from math import sin, cos
import numpy as np
from dials.array_family import flex


def get_SIM(dxtbx_det, dxtbx_beam, dxtbx_cryst, Fcalc_pdb=None, defaultF=10, SF=True):
    """

    :param dxtbx_det: detector object
    :param dxtbx_beam: beam object
    :param dxtbx_cryst: crystal object
    :param Fcalc_pdb: pdb file for structure factor gen (if None, flat Fcalc will be used)
    :return: sim_data.SimData instance
    """
    SIM = sim_data.SimData()
    SIM.detector = utils.strip_thickness_from_detector(dxtbx_det)
    SIM.detector = dxtbx_det

    # create nanoBragg crystal
    crystal = nanoBragg_crystal.NBcrystal(init_defaults=False)
    crystal.isotropic_ncells = False
    crystal.dxtbx_crystal = dxtbx_cryst
    crystal.thick_mm = 0.1  # hard code a thickness, will be over-written by the scale

    # mosaic block size
    crystal.Ncells_abc = 10,10,10  # will be update by GUI
    crystal.n_mos_domains = 100
    crystal.mos_spread_deg = 1

    if Fcalc_pdb is not None:
        if SF:
            miller_data = utils.get_complex_fcalc_from_pdb(Fcalc_pdb,
                dmin=1,
                dmax=999,
                wavelength=dxtbx_beam.get_wavelength(),
                k_sol=.5,
                b_sol=50)
            miller_data = miller_data.as_amplitude_array()
        else:
            symmetry = extract_symmetry_from(Fcalc_pdb)
            sg = str(symmetry.space_group_info())
            ucell_p = symmetry.unit_cell()
            miller_data = utils.make_miller_array(sg, ucell_p,defaultF=defaultF)
    else:
        symbol = dxtbx_cryst.get_space_group().info().type().lookup_symbol()
        ucell_p = dxtbx_cryst.get_unit_cell().parameters()
        miller_data = utils.make_miller_array(symbol, ucell_p,defaultF=defaultF)

    crystal.symbol = miller_data.crystal_symmetry().space_group_info().type().lookup_symbol()
    crystal.miller_array = miller_data
    SIM.crystal = crystal

    # create a nanoBragg beam
    beam = nanoBragg_beam.NBbeam()
    beam.size_mm = 0.001
    beam.unit_s0 = dxtbx_beam.get_unit_s0()

    #if spectra_file is not None:
    #    init_spectrum = load_spectra_file(spectra_file, total_flux, spectra_stride, as_spectrum=True)
    #else:
    #    assert total_flux is not None
    init_spectrum = [ (dxtbx_beam.get_wavelength(), 1e12) ]

    beam.spectrum = init_spectrum
    SIM.beam = beam

    # create the diffbragg object, which is the D attribute of SIM
    SIM.panel_id = 0
    SIM.instantiate_diffBragg(oversample=1, device_Id=0, default_F=0)

    return SIM


def get_pfs(all_pid, all_fast, all_slow):
    """
    makes a panel-fast-slow vector for low-level diffBragg method add_diffBragg_spots

    :param all_pid: 1-d array of panel ids (1 per pixel)
    :param all_fast: 1-d array of fast-scan pixel coords
    :param all_slow: 1-d array of slow-scan pixel coords
    :return: pfs vector for run_simdata method
    """
    assert len(all_pid)==len(all_fast)==len(all_slow)
    pan_fast_slow = np.ascontiguousarray((np.vstack([all_pid, all_fast, all_slow]).T).ravel())
    pan_fast_slow = flex.size_t(pan_fast_slow)
    return pan_fast_slow

def randomize_orientation(SIM, seed_rand=32, seed_mersenne=0):
    mersenne_twister = flex.mersenne_twister(seed=seed_mersenne)
    scitbx.random.set_random_seed(seed_rand)
    rand_norm = scitbx.random.normal_distribution(mean=0, sigma=2)
    g = scitbx.random.variate(rand_norm)
    rot = g(1)
    site = scitbx.matrix.col(mersenne_twister.random_double_point_on_sphere())
    ori = site.axis_and_angle_as_r3_rotation_matrix(rot[0],deg=False)
    SIM.crystal.dxtbx_crystal.set_U(ori)
    #SIM.instantiate_diffBragg(oversample=1, device_Id=0, default_F=0)
    SIM.D.Umatrix = ori

def sweep(SIM, phi_start, phistep, osc_deg, *args, **kwargs):
    print("beginning sweep")
    start_ori = SIM.crystal.dxtbx_crystal.get_U()
    SIM.crystal.dxtbx_crystal.rotate_around_origin((0,-1,0), phi_start)
    sum_pix = run_simdata(SIM, *args, **kwargs)
    n_steps = int(osc_deg/phistep)
    for step in range(1, n_steps+1):
        print("step {s} of {n}...".format(s=step, n=n_steps))
        # hard code spindle axis for now
        SIM.crystal.dxtbx_crystal.rotate_around_origin((0,-1,0), phistep)
        SIM.instantiate_diffBragg(oversample=1, device_Id=0, default_F=0)
        pix = run_simdata(SIM, *args, **kwargs)
        sum_pix += pix
    # reset
    print("finished sweep; resetting crystal orientation")
    SIM.crystal.dxtbx_crystal.set_U(start_ori)
    return sum_pix

def run_simdata(SIM, pfs, ucell_p, ncells_p, rot_p, spectrum=None, eta_p=None, G=1,
    diffuse_gamma=None, diffuse_sigma=None):
    """

    :param SIM: sim_data.SimData instance returned by get_SIM method
    :param pfs: pfs vector, flex.size_t , specifies which pixels to simulate
        length is 3xN where N is the number of pixels . See method get_pfs
    :param ucell_p: 6-tuple of unit cell parameters (a,b,c,alpha,beta,gamma) ; Angstrom / degrees
    :param ncells_p: 3-tuple of mosaic block size params, specifies extent of
        mosaic block along unit cell vectors, e.g. 10,5,5, means mosaic block
        is 10 unit cells along unit-cell a direction, and 5 unit cells along b,c
    :param rot_p: 3-tuple of rotational offsets in radians along lab axes.
            The meaning of these numbers is illuminated in the following code
            rot_p[0] is a rotation about (-1,0,0) in radians, and so forth
            THe final perturbation is defined as M=RotX*RotY*RotZ
            >> from scitbx.matrix import sqr col
            >> x = col((-1, 0, 0))
            >> y = col((0, -1, 0))
            >> z = col((0, 0, -1))
            >> RX = x.axis_and_angle_as_r3_rotation_matrix(rot_p[0], deg=False)
            >> RY = y.axis_and_angle_as_r3_rotation_matrix(rot_p[1], deg=False)
            >> RZ = z.axis_and_angle_as_r3_rotation_matrix(rot_p[2], deg=False)
            >> M = RX * RY * RZ
            >> Umat = sqr(dxbtbx_cryst.get_U())
            >> rotated_Umat = M*Umat
            >> dxtbx_cryst.set_U(rotated_Umat)
    :param spectrum: spectrum object list of 2-tuples. each 2-tuple is (wavelength, intensity)
    :param eta_p: float value of rotational mosaicity parameter eta
    :param G: scale factor for bragg peaks (e.g. total crystal volume)
    :param diffuse_gamma: 3-tuple of diffuse scattering param gammma
    :param diffuse_sigma: 3-tuple of diffuse scattering param sigma
    :return: flex array of simulated pixel values (same length as len(lfs) / 3)
    """


    if spectrum is not None:
        SIM.beam.spectrum = spectrum
        SIM.D.xray_beams = SIM.beam.xray_beams

    ucell_man = utils.manager_from_params(ucell_p)
    Bmatrix = ucell_man.B_recipspace   # same as dxtbx crystal .get_B() return value
    SIM.D.Bmatrix = Bmatrix

    if eta_p is not None:
        # NOTE for eta_p we need to also update how we create the SIM instance
        # see hopper_utils for examples (see methods "SimulatorFromExperiment" and "model")
        # FIXME not sure this is right yet
        SIM.update_umats_for_refinement(eta_p)
        SIM.crystal.mos_spread_deg = eta_p

    # Mosaic block
    #    diffuse_params_lookup = {}
    SIM.D.set_ncells_values(ncells_p)

    # diffuse signals
    if diffuse_gamma is not None:
        assert diffuse_sigma is not None
        SIM.D.use_diffuse = True
        SIM.D.diffuse_gamma = diffuse_gamma # gamma has Angstrom units
        SIM.D.diffuse_sigma = diffuse_sigma
    else:
        SIM.D.use_diffuse = False

    SIM.D.raw_pixels_roi *= 0

    ## update parameters:
    # TODO: if not refining Umat, assert these are 0 , and dont set them here
    rotX,rotY,rotZ = rot_p
    SIM.D.set_value(0, rotX)
    SIM.D.set_value(1, rotY)
    SIM.D.set_value(2, rotZ)

    npix = int(len(pfs)/3)
    SIM.D.add_diffBragg_spots(pfs)

    pix = G*SIM.D.raw_pixels_roi[:npix]

    return pix

