Mike Devereux, 21.9.2020

This folder contains 2 sets of scripts:

#### esp_gdma_grid

- can be used to take a GDMA-format "punch" file and two existing reference Gaussian
  format cube files, one containing the MEP and one containing the electron density
  for the molecule of interest
- for these scripts the GDMA output must have been generated in the same global axis,
  i.e. using the same nuclear coordinates, as the 2 cube files
- the scripts will then generate an MEP cube file across the same grid used in the
  reference cube files using the multipole moments from the GDMA file
- once the multipolar cube file has been generated, a pointwise comparison between
  MEP values using multipoles and MEP values in the reference cube file is performed
- a summary of the RMSE in different regions with boundaries defined by multiples of
  the atomic vdw radii is calculated and presented
- an additional cube file "error.cube" contains the difference between the multipolar
  and reference cube file at each point across the grid, allowing visualization of
  the error in a program capable of rendering cube files acros isodensity surfaces 
  such as gaussview

#### esp_charmm_grid

- the purpose of these scripts is similar to esp_gdma_grid, except they do not 
  require nuclear coordinates in the GDMA file to correspond to coordinates in the
  reference cube files
- this means that, for example, atomic multipole moments calculated in the 
  equilibrium geometry can be used to evaluate the MEP of a different conformation
- this works by converting the GDMA file to CHARMM MTPL format, and using CHARMM's
  local axis system definitions to evaluate the MEP for an arbitrary conformation
- the MEP is evaluated in CHARMM as the interaction between the multipolar molecule
  and a neighboring unit charge. The interaction is evaluated for all grid
  coordinates and written to a file, which is then converted to Gaussian cube format
- as such you must modify the CHARMM files for your system to properly define the
  molecule in CHARMM, using the same atom ordering used in the cube and GDMA files
- to check that you set your system up correctly, you can run one example using 
  multipoles generated with the same nuclear coordinates used to generate the MEP
  cube file. You can then run the esp_gdma_grid scripts and the results should be
  the same for both to within a few tenths of a kcal

#### notes
- as a general reminder when working with these scripts, remember that the gaussian
  cube file typically has coordinates written in Bohr, whereas CHARMM / PDB uses
  Angstrom