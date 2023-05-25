# sim_erice

Bragg simulation and model fitting frontend

## SimView: lightweight simulator and viewer for diffraction still images.

Run without any arguments, this program simulates diffraction of a small lysozyme crystal. You may ask the viewer to fetch a different model by entering that PDB ID in the associated text box. 

Diffraction from your crystal will be simulated with a set of default parameters for mosaic block size, beam energy and bandpass, etc. that you can change by adjusting each of the dials. Parameters to be adjusted can be selected either with the visual controls or with keyboard shortcuts. Numerical visual controls are made "active" by clicking in the text box or pressing their up or down arrows. The active control is labeled in bold, blue text. The keyboard will also route input to the active control by interpreting a keyboard "up arrow" keypress as the visible "up arrow" button, and the same for the "down arrow". For finer control, you may also enter text directly into the text box. The GUI will respond to new text input only after you press enter. In summary, keyboard shortcuts are as follows:

**Up arrow:**           increase this parameter (same as pressing the small up arrow next to the control)

**Down arrow:**         decrease this parameter (same as pressing the small down arrow next to the control)

**Shift-up arrow:**     increase this parameter a lot

**Shift-down arrow:**   decrease this parameter a lot

**Return:**             finalize a changed parameter typed directly into the box

## Explanation of parameters:

**DomainSize:** the small blocks within crystals that are internally perfectly aligned and identical are called mosaic domains. The average size of a domain is one contributing factor determining the sharpness of the Bragg peaks. As the size of the mosaic domain increases, the effect of constructive interference of the scattered X-rays gets stronger, and the peak sharpens.

**MosAngDeg:** the mosaic domains are misoriented relative to one another by an average angle denoted the mosaic angle, displayed here in degrees. Larger differences in the orientations of the domains produce contributions to Bragg peaks over a wider area on the detector, resulting in more diffuse spots.

**Brightness:** the overall image scale, which governs the minimum intensity value that will be visible. This is a tradeoff with the maximum value that can be visually differentiated, where everything above this maximum is marked as overflow and colored in yellow in the simulation-only view mode.

**Display mode:** either "Simulation only" to display a single simulation in cyan, or "Overlay with pinned" to also display a second simulation in the red channel of the same image, generating grayscale where they perfectly match. This can be used to emphasize and explore the effects of changing parameters, by "pinning" a simulation (using the "Update pinned image" button) and then adjusting a parameter. 

**Unit cell a, b, c:** a, b and c are the unit cell axis lengths, and we expose fractional scaling parameters for each of these to simulate e.g. a 10% larger unit cell in the b axis when setting b_scale to 1.1 -- you should be able to see both the fractional scales in the controls and the resulting cell parameters in the static text to their right. Depending on the symmetry of the crystal, these may all vary independently or may be constrained to scale together, in which case not all unit cell scale controls will be enabled. Larger unit cells produce constructive and destructive interference at smaller angles, resulting in more closely spaced Bragg peaks on the detector.

**Missetting angles:** in serial crystallography, we are often faced with differences between the true orientation of the crystal and what is determined during indexing. A misorientation of the X or Y axis means a misestimation of which spots are in diffracting conditions, and will affect which spots (or how much of the spots) appear on the detector. Since the Z axis is along the beam direction, misorientation in Z results in rotation of the entire diffraction pattern. In the simulation viewer we expose these misorientations as rotations in x, y and z separately from the crystal orientation angles x, y and z. The crystal orientation can separately be managed with its own x, y and z controls and the "Randomize orientation" button.

**Experiment mode:** a default mode of stills (serial) crystallography is enabled in the viewer, but you may also switch to rotation mode to simulate a standard experiment with a goniometer. The goniometer axis is identified by a vertical yellow line, and the oscillation width and image number in a sweep can be controlled in the GUI. Advanced users may also want to adjust oscillation step size or total sweep angle with phil parameters when launching the GUI.

**Energy/Bandwidth:** the energy of the X-ray beam denotes the average energy of an individual pulse. For X-ray free electron laser (XFEL) pulses generated by self-amplified spontaneous emission (SASE), each pulse is spiky and stochastic, with a typical bandwidth of between 10 and 30 eV. (A monochromater may be used to narrow this to around 1 eV, not shown.) The average energy of a series of pulses is matched to the energy entered in this control. You can advance to the next SASE spectrum in the series by pressing "New XFEL pulse." Bandwidth is only relevant to the Gaussian spectrum mode.

**Spectrum shape:** to simulate images with real SASE spectra, toggle this option to "SASE". To display diffraction produced by a smooth Gaussian function that can be adjusted in energy and bandwidth (illustrative but nonphysical), toggle this to "Gaussian". SASE spectra are stochastic and will intentionally not be modified by the energy and bandwidth controls. The "monochromatic" option is recommended when diffuse scattering is enabled, to offset the greater computational cost.

**Use structure factors:** when we observe diffraction from a single particle, we are seeing the Fourier transform of the particle, with [a great deal of] added noise. For many copies of a particle all in the same orientation, we can expect the signal to be much stronger. If these aligned particles are also periodically spaced, as in the context of a crystal, they act as a diffraction grating, and we see this signal only at positions of constructive interference (the reciprocal lattice points in diffracting conditions, i.e. Bragg peaks). In order to better see the effects of the above parameters, we can ignore the Fourier transform of the asymmetric unit and pretend all of the Bragg peaks have uniform intensity. This is what we are simulating when we set "Use structure factors" on or off.

**Diffuse scattering:** this physical effect is a result of variation and imperfection within the crystal and surrounding solvent, and it appears as X-ray scattering between and around the Bragg peaks. Modeling this can be toggled on or off. It is recommended to use the monochromatic spectrum setting when calculating diffuse signal.

**Diff_gamma and diff_sigma:** parameters describing the diffuse scattering signal produced by long-range correlations. Gamma denotes the correlation length in Ångstroms and sigma squared is the amplitude of the correlated vibrations in square Ångstroms.

**Diff_aniso:** anisotropic quality to the diffuse scattering. This is arbitrarily assigned directionality, and the adjustible parameter controls the degree to which this effect is observed.

**Further notes on the image viewer:** matplotlib enables zooming and panning of the figure. We have also added pixel-specific information to the annotation text displayed in the lower right corner, such as the fractional Miller index and resolution at a given position on the detector.
