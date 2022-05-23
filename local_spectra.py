from __future__ import division
import pickle
import numpy as np
import os
import libtbx.load_env

def get_results():
  data_file = libtbx.env.find_in_repositories(relative_path="sim_erice/spectra.pickle",
            test=os.path.isfile)
  R = pickle.load(open(data_file,"rb"))
  return R

class linear_fit:
  def __init__(self,data):
    self.x = data["expidx"] # the expected index over the 1D spectral distribution (low pass filtered)
    self.y = data["energy"] # the ebeam energy
    #print(len(self.x))
    #print(len(self.y))
    # y = Ap, where A = [[x 1]] and p = [[m], [c]]
    A = np.vstack([self.x, np.ones(len(self.x))]).T
# workaround allows use of non-thread-safe numpy lstsq, even if openMP is enabled elsewhere in the Python program
    import os,omptbx
    workaround_nt = int(os.environ.get("OMP_NUM_THREADS",1))
    omptbx.omp_set_num_threads(1)
    self.m,self.c = np.linalg.lstsq(A,self.y, rcond=-1)[0]
    # says numpy: FutureWarning: `rcond` parameter will change to the default of machine precision times ``max(M, N)`` where M and N are the input matrix dimensions.
    #To use the future default and silence this warning we advise to pass `rcond=None`, to keep using the old, explicitly pass `rcond=-1`.
    omptbx.omp_set_num_threads(workaround_nt)
    # y = mx + c
    # x = (1./m) y - (c/m)

    # Solve the equations without numpy
    # Taken from Bevington, third edition, eqns (6.12,6.13)
    NN = len(self.x)
    from scitbx.array_family import flex
    x_array = flex.double(self.x)
    sumx = flex.sum(x_array)
    sumsq_x = sumx * sumx
    sumx_sq = flex.sum(x_array*x_array)
    sumxy = flex.sum(x_array * flex.double(self.y))
    sumy = flex.sum(flex.double(self.y))
    delta_prime = NN*sumx_sq - sumsq_x
    assert delta_prime != 0.
    mm = (1./delta_prime) * (NN*sumxy - sumx * sumy)
    cc = (1./delta_prime) * (sumx_sq*sumy - sumx * sumxy)
    from libtbx.test_utils import approx_equal
    assert approx_equal(mm,self.m)
    assert approx_equal(cc,self.c)

  def get_residuals(self):
    print(len(self.y))
    calc_idx = (1./self.m)*np.array(self.y) - (self.c/self.m)
    print(len(calc_idx))
    return self.x - calc_idx

class spectra_simulation:
  def __init__(self):
    self.R = get_results()
    self.LF = linear_fit(self.R)
    # get some information to help normalize things
    maxima = []
    bk_subtracted_sum = []
    for image in range(len(self.R["energy"])):
      bk_subtracted_sum.append(np.sum(self.R['spectra'][image]))
      maxima.append(np.max(self.R['spectra'][image]))

    self.max_of_max = max(maxima)
    average_integrated = np.mean(bk_subtracted_sum)
    print("average_integrated",average_integrated)
    self.bk_subtracted_sum = bk_subtracted_sum
    self.average_integrated = average_integrated
    self.NS = len(self.R["spectra"][0]) # number of points in each spectrum
    self.N = len(self.R["spectra"]) # number of events overall
    self.iter_idx = 0

  def generate_recast_renormalized_images(self, energy, total_flux):
    spectrum_fitted_energy = self.LF.m * np.array(range(self.NS)) + self.LF.c
    offset = energy - self.get_average_expected_energy()
    offset_energy = spectrum_fitted_energy + offset
    while 1:
      image = self.iter_idx % self.N
      self.iter_idx += 1

      from scitbx.array_family import flex
      y = flex.double(list(self.R['spectra'][image]))
      ysum = self.bk_subtracted_sum[image]

      expected_energy = self.LF.m * self.R["expidx"][image] + self.LF.c + offset
      print(image,"ebeam = %7.2f eV"%(expected_energy),"%5.1f%% of average pulse intensity"%(100.*
        self.bk_subtracted_sum[image]/self.average_integrated))
      assert expected_energy > 0.
      channel_flux = flex.double(100) # 100 energy channels altogether
      channel_mean_eV = flex.double(range(100)) + energy - 49.5
      eV_to_angstrom = 12398.425
      channel_wavelength = eV_to_angstrom / channel_mean_eV
      for idx in range(len(offset_energy)):
        i_energy = offset_energy[idx]
        channel = int(i_energy - (energy-50))
        if 0 <= channel < 100:
          channel_flux[channel] += self.R['spectra'][image][idx] * total_flux / self.average_integrated
      yield channel_wavelength,channel_flux,eV_to_angstrom / expected_energy

  def get_average_expected_energy(self):
    idx = np.array(self.LF.x)
    fitted_energy = self.LF.m * idx + self.LF.c
    #return np.mean(fitted_energy)
    from scitbx.array_family import flex
    idx_c = flex.double(self.LF.x)
    fitted_energy_c = float(self.LF.m) * idx_c + float(self.LF.c)
    print ("numpy",np.mean(fitted_energy), "flex",flex.mean(fitted_energy_c))
    return flex.mean(fitted_energy_c)

