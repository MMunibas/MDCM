# MDCM
Machine learning code to fit Minimal Distributed Charge Models (MDCMs)

## Authors:
MDCM scripts and code by Mike Devereux (Michael.Devereux at unibas.ch) except:
* **differential_evolution.f90** by Shun Sakuraba
* **cubefit** original implementation by Oliver Unke, updated and extended by Mike Devereux 
* **punch-to-dcm.f90** original implementation by Shampa Raghunathan, updated by Mike Devereux
* **mtpfit.py** by Oliver Unke
Additional code required for multipole cube generation example:
* **examples/mtpl-cubes/.../tristan-espfit/\*** Tristan Bereau and Christian Kramer


## References:
1. Unke, O. T.; Devereux, M.; Meuwly, M.; Minimal distributed charges:  Multipolar quality at the cost of point charge electrostatics. J. Chem. Phys.2017,147, 161712
1. Devereux, M.; Pezzella, M.; Raghunathan, S.; Meuwly, M.; Polarizable Multipolar Molecular Dynamics Using Distributed Point Charges. Submitted Aug. 2020

## Code Overview:
* Various new and pre-existing utilities are combined. The main fitting routines for MDCM charges are located in bin/src and are compiled to produce (p)cubefit.x.
* The python utility "mtpfit.py" fits high-ranking (*l*=5) atomic multipoles to the molecular electrostatic potential (MEP) across a grid provided in Gaussian cube file format. It also requires a second Gaussian cube file containing the electron density at each corresponding point, in order to define the exclusion zone inside the molecule where points are discarded and the exclusion zone outside a minimum density cut-off where the MEP is typically small.
* The remaining helper scripts provide utilities to assist with file formats etc. and are used by the examples provided

## Compilation:
* Only the MDCM Fortran code requires compilation. Type "make" in the folder that contains the Makefile to create two binaries, *cubefit.x* and *pcubefit.x*. *pcubefit.x* is parallelized using OpenMP, *cubefit.x* is serial code.
* Note that due to the nature of DE, any speedup using the parallelized code can be small and it is usually better to run more fits simultaneously than a single fit with more cores

## Running:
### mtpfit.py
**Fit atomic multipole moments to the molecular electrostatic potential:**
`mtpfit.py -pot $PCUBE -dens $DCUBE -lmax 5 -qtot 0.0`

Options:
* -pot:   MEP cube file
* -dens:  electron density cube file
* -lmax:  maximum rank for atomic multipoles (ditriantapole)
* -qtot:  total molecular charge

### (p)cubefit.x
**Atomic distributed charge fitting:**
`pcubefit.x -greedy $FITTED-MTPL -esp $PCUBE -dens $DCUBE -nacmin $MINCHG -nacmax $MAXCHG -atom $ATOMINDEX -ntry $NTRY -onlymultipoles -v > $OUTFILE`

Options:
* -greedy:  use "greedy" fitting algorithm in differential evolution fitting (recommended), and define $MTPFILE containing fitted atomic multipoles
* -esp: define Gaussian cube file "$PCUBE" containing the molecular electrostatic potential
* -dens: define Gaussian cube file "$DCUBE" containing the molecular electron density
* -nacmin: define lowest number of charges to fit per atom (usually 1)
* -nacmax: define highest number of charges to fit per atom (usually 3-4)
* -atom: define index of atom to be fitted (corresponds to ordering in cube files). Fitting atoms separately allows efficient parallelization of the fitting process.
* -ntry: set number of complete fitting runs. As the fitting code involves making random "mutations" to existing populations of candidate solutions, better results may be obtained by repeating the fitting process a few times and selecting the best result
* -onlymultipoles: state that we want to fit multipole moments only in this step and not atomic charges
* -v: verbose output

**Fragment distributed charge fitting:**
`pcubefit.x -greedy $MTPFILE -esp $PCUBE -dens $DCUBE -ncmin $MINCHG -ncmax $MAXCHG -atom $ATOMLIST -nacmax $MAXCHG -ntry $NTRY -v > $OUTFILE`

Options:
* -greedy:  use "greedy" fitting algorithm in differential evolution fitting (recommended), and define $MTPFILE containing fitted atomic multipoles
* -esp: define Gaussian cube file "$PCUBE" containing the molecular electrostatic potential
* -dens: define Gaussian cube file "$DCUBE" containing the molecular electron density
* -nacmax: the highest number of charges per atom used during atom fitting in the previous step
* -ncmin: define lowest number of charges to fit for this fragment
* -ncmax: define highest number of charges to fit for this fragment
* -atom: define a comma-separated list of atom indices to be fitted without spaces (indices correspond to atom ordering in cube files). Fitting fragments allows efficient fitting of models with too many charges to fit efficiently all at once
* -ntry: set number of complete fitting runs. As the fitting code involves making random "mutations" to existing populations of candidate solutions, better results may be obtained by repeating the fitting process a few times and selecting the best result
* -v: verbose output

**Combining fragment models to build molecular models**
`combine.sh`
Options are set inside the script by defining the shell variables:
* WORKDIR: the root folder of the working directory
* BINDIR: the folder containing the (p)cubefit.x binaries
* FRAGDIR: the folder containing the fitted fragment models
* NFRAG: the number of non-overlapping fragments that were fitted for the molecule
* NFIT: the number of fits that were performed for each fragment
* MINCHGS: the lowest number of charges for the molecular charge model
* MAXCHGS: the maximum number of charges for the molecular charge model

**Molecular MDCM refinement**
`pcubefit.x -xyz $INITIALXYZ $MTPFILE -esp $PCUBE -dens $DCUBE -nacmax $MAXATMCHG -ntry $NTRY -v > $OUTFILE`

Options:
* -xyz: defines the file containing the molecular charge model to be refined (usually the charge model created by combining fitted fragment models), and the file containing the fitted atomic multipoles
* -esp: define Gaussian cube file "$PCUBE" containing the molecular electrostatic potential
* -dens: define Gaussian cube file "$DCUBE" containing the molecular electron density
* -nacmax: the highest number of charges per atom used during atom fitting in the atomic charge fitting step
* -ntry: set number of complete fitting runs. As the fitting code involves making random "mutations" to existing populations of candidate solutions, better results may be obtained by repeating the fitting process a few times and selecting the best result
* -v: verbose output

## Examples
The simplest way to use the code for your own system is to adapt one of the examples provided. They are deliberately quite long-winded and modular to allow easy adaptation to new cases.
### naptha
* This example should act as a template to fit MDCMs for most cases. The basic workflow is:
1. Create a Gaussian-format "cube" file of MEP grid points from a GDMA multipolar output file in "punch" format. This step allows subsequent fitting of an MDCM to GDMA multipoles rather than to a reference ab initio MEP, if ab initio data is to be used then this step can be skipped.
1. Fit new atomic multipoles to the MEP cube file just generated. This step isn't really necessary in the multipolar MEP case but can be easily adapted to fit atomic multipoles to an ab initio MEP instead, so is kept for illustrative purposes.
1. Fit atomic distributed charge models to the ESP generated by the atomic multipoles from the previous step. These charge models will be used as seeds for subsequent DE (Differential Evolution) fitting of fragments
1. DE fitting of MDCMs for fragments or of a whole molecule with the specified number of charges. Scaling with number of charges is bad so a fragment approach allows linear scaling for large molecules at the cost of some accuracy and redundancy of charges in the final model
1. Combination of fitted fragment models with the lowest RMSE to create promising molecular MDCMs
1. Refinement of the assembled molecular models in this case using a simplex algorithm to allow scalability to larger systems. DE refinement would likely yield better results in a reasonable time (several hours) for small systems with up to ca. 18 charges.
* The example is run stepwise using the scripts in the various folders in numerical order
* The final output is an MDCM in the global axis in .xyz file format with fitted charges added in a 5th column after the atomic coordinates in each line
