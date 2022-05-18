from __future__ import division
from matplotlib import pyplot as plt
from LS49.spectra.generate_spectra import spectra_simulation
"""working example:  how the LS49 code can simulate SASA spectra"""

specsim = spectra_simulation()
max_images = 1000
iterator = specsim.generate_recast_renormalized_images(nlimit=max_images,energy=7120.,total_flux=1e12)

for image_no in range(20):
  wavlen, flux, wavelength_A = next(iterator) # list of lambdas, list of fluxes, average wavelength
  # the lambdas are equally spaced in energy, but expressed as wavelength in Angstroms
  # the fluxes are given in units of photons
  plt.bar(wavlen, flux, width=0.0002)
  plt.show()

