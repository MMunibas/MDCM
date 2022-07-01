# MDCM-2.0
Updated machine learning code to fit Minimal Distributed Charge Models (MDCMs). Features include:

* Fitting atom-centered multipolar charge models to the molecular ESP, using either constrained Least-Squares fitting (recommended) or Differential Evolution (DE) (not recommended)
* Fitting distributed atomic charge models that can be combined to provide an initial population for subsequent DE fitting of molecules or fragments
* DE fitting of molecular fragments or functional groups to the multipolar ESP, allowing fitting of larger molecules with more charges than would be computationally feasible by fitting all charges in the molecule simultaneously
* Combination of fitted fragment models to create a good initial guess for refinement of a total molecular model
* DE and / or simplex fitting and refinement of molecular charge models to the reference ESP
* Symmetry-constrained fitting to ensure distributed charge models possess the same symmetry as the molecule they describe
* GPU support (still developmental)
* Fitting to multiple conformers simultaneously by defining local axes and two or more reference ESP grids for different conformers. This helps to avoid creating models, in particular for flexible molecules, that are too finely tuned to a single conformer and perform poorly upon conformational change.

## Authors:
MDCM scripts and code by Mike Devereux (Michael.Devereux at unibas.ch) except:
* **differential_evolution.f90** by Shun Sakuraba
* **cubefit** original implementation by Oliver Unke, updated and extended by Mike Devereux
* **symmetry.f90** by Mike Devereux with a contribution from Oliver Unke
* **rmse_gpu.cuf** by Mike Devereux
* **punch-to-dcm.f90** original implementation by Shampa Raghunathan, updated by Mike Devereux
* **mtpfit.py** by Oliver Unke

Additional code required for multipole cube generation example:
* **examples/mtpl-cubes/.../tristan-espfit/\*** Tristan Bereau and Christian Kramer


## References:
1. Unke, O. T.; Devereux, M.; Meuwly, M.; Minimal distributed charges:  Multipolar quality at the cost of point charge electrostatics. J. Chem. Phys.2017,147, 161712
1. Devereux, M.; Pezzella, M.; Raghunathan, S.; Meuwly, M.; Polarizable Multipolar Molecular Dynamics Using Distributed Point Charges. J. Chem. Theory Comput. 2020, 16, 7267

## Code Overview:
* Various new and pre-existing utilities are combined. The main fitting routines for MDCM charges are located in bin/src and are compiled to produce (p)cubefit.x.
* CUDA support was added to produce the binary cudacubefit.x
* The python utility "mtpfit.py" fits high-ranking (*l*=5) atomic multipoles to the molecular electrostatic potential (MEP) across a grid provided in Gaussian cube file format. It also requires a second Gaussian cube file containing the electron density at each corresponding point, in order to define the exclusion zone inside the molecule where points are discarded and the exclusion zone outside a minimum density cut-off where the MEP is typically small.
* The remaining helper scripts provide utilities to assist with file formats etc. and are used by the examples provided

## Compilation:
* Only the MDCM Fortran code requires compilation. Type "make serial" or "make parallel" in the folder that contains the Makefile to create two binaries, *cubefit.x* and *pcubefit.x*. *pcubefit.x* is parallelized using OpenMP, *cubefit.x* is serial code.
* To compile the CUDA code you need to have CUDA-Fortran installed on your computer and in
your PATH. The make command is then "make cuda".
* Note that due to the nature of the fitting procedure it is usually better to run more fits simultaneously than a few fits serially using more cores. For very long fits (many charges and / or ESP fitting points in the grid) more CPU cores or a GPU can give noticeable improvements.

## Running:
### mtpfit.py
**Fit atomic multipole moments to the molecular electrostatic potential:**
`mtpfit.py -pot <cubefile> -dens <cubefile> -lmax <max_rank> -qtot <molecular charge> -fixq <`

Options:
* -pot:   MEP cube file (Gaussian format)
* -dens:  electron density cube file
* -lmax:  maximum rank for atomic multipoles (between 0 and 5, i.e. charge and ditriantapole)
* -qtot:  total molecular charge (a.u., default = 0)
* -fixq:  file containing charges to freeze during fitting (optional)

This code requires a Gaussian-format ESP cube file and an electron density cube file as input to provide reference data for fitting atomic multipole moments. ESP grid points inside the molecule or too far from any atom in the molecule will be excluded. The maximum rank of the fit specifies the highest ranking multipole moments that will be used to describe the ESP. The total charge refers to the total charge of the molecule, as constraints are applied to maintain this value during least-squares fitting.

The fixq option allows the user to specify a file from a previous multipole fit (or any file in similar format), so that the charge terms can be fixed during subsequent multipole fitting. This can be useful if, for example, CHARMM charges are preferred with a purely multipolar correction. The main application, though, is to allow fitting of fragments or functional groups to multiple molecular conformers. As fragments are fitted to the multipolar ESP as reference, multipoles are required for each conformer. This can only work if the total charge of each fragment within the molecule remains fixed for all conformers, so the charges fitted for the first conformer should be applied to the remaining conformers as well.

### (p)cubefit.x
**Atomic distributed charge fitting:**
`pcubefit.x -greedy -mtpfile <fitted_multipole_file> -esp <esp_cube_file> -dens <density_cube_file> -nacmin <min_chgs> -nacmax <max_chgs> -atom <atom_index> -ntry <num_tries> -onlymultipoles -sym -gpu -v > <log_file>`

Options:
* -greedy: use "greedy" fitting algorithm in differential evolution fitting (recommended), 
* -mtpfile: specify file containing atomic multipoles previously fitted with e.g. mtpfit.py
* -esp: specify Gaussian cube file containing the molecular electrostatic potential
* -dens: specify Gaussian cube file containing the molecular electron density
* -nacmin: define lowest number of charges to fit per atom (usually 1)
* -nacmax: define highest number of charges to fit per atom (usually 3-4)
* -atom: define index of atom to be fitted (corresponds to ordering in cube files). Fitting each atom separately allows efficient parallelization of the fitting process.
* -ntry: set number of complete fitting runs. As the fitting code is stochastic, involving random "mutations" to existing populations of candidate solutions, better results can be obtained by repeating the fitting process a few times and selecting the best result
* -onlymultipoles: state that we want to fit atomic charges to multipole moments only in this step and not proceed to full molecular ESP fit
* -sym: state that we want to apply symmetry constraints so that atomic charge models respect the molecular symmetry (optional)
* -gpu: state that we want to run with CUDA support (optional)
* -v: verbose output (optional)

This DE step uses the atomic multipoles from the previous step to construct an "atomic ESP" grid that is used to fit distributed charge models for each atom. We typically fit between 1 and 4 charges per atom, which makes subsequent molecular fitting computationally tractable and still yields molecular charge models of accuracy comparable to a fitted multipole expansion truncated at quadrupole.

The purpose of the atomic charge models is to provide an initial population for subsequent DE fitting that is already focused on promising regions of parameter space. Starting from truly random initial parameter values requires much more exploration of the full-dimensional parameter space when fitting the molecular models, which becomes intractable.

**Fragment distributed charge fitting:**
`pcubefit.x -greedy -mtpfile <fitted_multipole_file> -esp <esp_cube_file> -dens <density_cube_file> [-frames <axis_frames_file> -mtpfile <fitted_multipole_file_2> -esp <esp_cube_file_2> -dens <density_cube_file_2>...] -ncmin <min_chgs> -ncmax <max_chgs> -atom <list> -nacmax <max_chgs> -ntry <num_fits> -converge <cutoff> -gpu -v > $OUTFILE`

Options:
* -greedy: use "greedy" fitting algorithm in differential evolution fitting (recommended)
* -mtpfile: define file containing fitted atomic multipoles
* -esp: define Gaussian cube file containing the molecular electrostatic potential
* -dens: define Gaussian cube file containing the molecular electron density
* -nacmax: the highest number of charges per atom used during atom fitting in the previous (atom-fitting) step
* -ncmin: define lowest number of charges to fit for this fragment
* -ncmax: define highest number of charges to fit for this fragment
* -atom: define a comma-separated list of atom indices to be fitted without spaces (indices correspond to atom ordering in cube files). 
* -ntry: set number of complete fitting runs. As the fitting code involves making random "mutations" to existing populations of candidate solutions, better results may be obtained by repeating the fitting process a few times and selecting the best result
* -frames: if fitting to multiple conformers, reference axis frames must be defined (optional)
* -converge: user-specified convergence criterion (kcal/mol) based on change in RMSE over time (optional)
* -gpu: use a GPU to evaluate the ESP (optional)
* -v: verbose output (optional)

Fitting fragments allows efficient fitting of models with too many charges to fit efficiently all at once. The molecule is divided into fragments by the user, to create fragments that are as large as possible while remaining computationally manageable. The total number of charges to fit for each fragment can be defined as a range, but it is preferable to run a separate fit for each number of charges to allow efficient parallelization.

The atomic multipole file provided is used to generate the fragment ESP from the multipoles of the atoms comprising the fragment. As such, the reference ESP cube file is not really used here to provide reference data.

The initial guess model is created from the combination of atomic charge models from the previous step that are predicted to yield the lowest RMSE when combined for a given number of charges.

Multiple conformers can be used as a reference to make resulting models more robust to conformational change. A set of local reference axes must be provided to allow each candidate model to be transformed to each conformation, using the file defined by the "-frames" flag. The format is e.g.:

ALA3
2    1    4   BO
...
34   33   31  BO

The first line gives the molecule name, subsequent lines contains groups of 3 atom indices that are used to create local axis frames (see https://doi.org/10.1021/ct500511t), the "BO" (or "BI") specfy either the bond-type axis systems described in the reference or a related axis system where the z-axis of the 2nd atom points along the A-B-C bisector, rather than along bond B-A.

After specifying -frames, a separate atomic multipole file (see notes on constraining fragment charge during multipole fitting above), ESP cube file and density cube file is provided for each conformer.

The -converge option can be useful for larger fitting problems. DE fitting of the ESP typically converges quite rapidly to an RMSE within ca. 0.01 kcal/mol of the final solution, but optimization continues as long as the population is still diverse in case a new and better minimum is found. The user can optionally specify looser convergence criteria by specifying a threshold. If the best DE solution improves by less than this threshold over 1000 fitting generations then the process exits early, typically saving a significant amount of time at little cost in accuracy. Note that as long as a diverse fitting population remains this approach potentially misses significantly improved solutions, however.

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
`pcubefit.x -xyz <charges_file> -mtpfile <fitted_multipole_file> -esp <esp_cube_file> -dens <density_cube_file> [-frames <axis_frames_file> -mtpfile <fitted_multipole_file> -esp <esp_cube_file> -dens <density_cube_file> ...] -nacmax <max_chgs> -simplex -ntry $NTRY -sym -gpu -v > <log_file>`

Options:
* -xyz: defines the file containing the molecular charge model to be refined (usually the charge model created by combining fitted fragment models, or from a previous full-molecule fit)
* -mtpfile: the file containing the fitted atomic multipoles
* -esp: define Gaussian cube file containing the molecular electrostatic potential
* -dens: define Gaussian cube file containing the molecular electron density
* -nacmax: the highest number of charges per atom used during atom fitting in the atomic charge fitting step
* -ntry: set number of complete fitting runs. As the fitting code involves making random "mutations" to existing populations of candidate solutions, better results may be obtained by repeating the fitting process a few times and selecting the best result
* -sym: state that we want to apply symmetry constraints so that atomic charge models respect the molecular symmetry (optional)
* -gpu: state that we want to run with CUDA support (optional)
* -frames: if fitting to multiple conformers, reference axis frames must be defined (optional)
* -converge: user-specified convergence criterion (kcal/mol) based on change in RMSE over time (optional)
* -simplex: use simplex refinement only and no DE to make large fitting problems tractable (optional)
* -v: verbose output (optional)

Model refinement is a critical step when fragment fitting, as the individual fragments fitted to multipolar ESPs do not yet yield an optimal solution when combined to describe the molecular ESP. Depending on the size of the total system it may not be feasible to run DE optimization on the entire system, and as we should be reasonably close to minimum after fragment fitting the "-simplex" option is often an effective alternative to create a refined model with many charges.

It is also possible to refine existing models, for example by adding additional conformers with the "-frames" flag and additional -esp, -dens and -mtpfile flags as described for fragment fitting.

The -converge option can also be useful for refinement in improving computational feasibility, as described for fragment fitting.

Finally, this option can be used to fit entire molecules directly to the reference ESP without needing the fragment fitting step for molecules that are sufficiently small / that require sufficiently few charges. In this context the model to be refined is constructed from the atomic charge models that were fitted to the atomic multipolar ESP, as described above.

**Analysis: Creating cube files**
`pcubefit.x -generate [-xyz <charges_file>] [-multipole -mtpfile <fitted_multipole_file>] -esp <esp_cube_file> -dens <density_cube_file> [-frames <axis_frames_file> -esp <esp_cube_file_2> -dens <density_cube_file_2>] -v`

* -xyz: defines the fitted charge or multipoole model that we want to use to generate an ESP cube file
* -multipole: specifies that we want to generate a multipolar ESP from a fitted multipole model instead of using a fitted charge model (optional)
* -esp: defines the Gaussian format cube file containing the reference ESP
* -dens: defines the Gaussian format cube file containing the reference electron density
* -frames: defines the local axis frames if we want to generate a cube file for another conformer (optional)
* -v: verbose output (optional)

This utility is useful for analyzing the quality of models and for visualizing issues with the fitted ESP. By default the molecular geometry and the ESP grid specifications are taken from the supplied cube files, then the requested fitted charge or multipolar model is used to generate the fitted ESP at each point across the reference grid. As such, the fitted grid and the reference grid can subsequently be pointwise compared to evaluate the quality of the model.

There is also the option to supply a local reference axis frames file (see above) with Gaussian cube files representing another conformer using the -frames option with a second -esp and -dens option. In this case the charge (not multipolar!) model that was fitted for the conformer of cube file 1 is transformed using the local reference axes to conformer to, and used to generate the ESP of that conformer using the atomic coordinates and ESP grid specifications in cube file 2. This option can be used to evaluate the preformance of a model fitted for one conformer (or set of conformers) when applied to another conformer, checking how robust the model is to conformational change.

**Analysis: comparing cube files**
`pcubefit.x -analysis -esp <esp_cube_file> -esp2 <esp_cube_file_2> -dens <density_cube_file> -v > <log_file>`

* -esp: Gaussian format cube file containing reference ESP
* -esp2: Gaussian format cube file containing ESP to compare to reference
* -dens: Gaussian format cube file containing electron density
* -v: verbose output (optional)

This mode is used to quantatively compare the ESP in different regions (close, near, far) between two supplied Gaussian-format cube files. This is typically used to compare the ESP from a fitted model to the reference ESP that it was fitted to. In addition to a statistical analysis, a new Gaussian format cube file quantifying the error in ESP across the same grid is created and can be used to visualize spatial distribution of the error in different planes or mapped onto a molecular surface.

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
* The example is run stepwise by adapting and running the scripts in the various folders in numerical order
* The final output is an MDCM in the global axis in .xyz file format with fitted charges added in a 5th column after the atomic coordinates in each line

### dcm
* This example provides basic scripts to convert a GDMA format multipole file to a fixed DCM (not MDCM!) charge arrangement. 
* No fitting is required, the process is an analytical conversion. See https://doi.org/10.1021/ct500511t for details
* In constrast to MDCM, the conversion process is fast (few seconds), but results in redundant charges that will slow down subsequent evaluations of the electrostatic interaction
* The example is run by modifying and executing the dcm-cube-gen.sh script

### mtpl-cubes
* This example uses neither DCM nor MDCM, but is useful for comparing the performance of multipolar models with (M)DCMs in certain scenarios
* The example consists of 2 folders
**esp_gdma_grid** provides a workflow in 1-convert-punch.sh to create a Gaussian format MEP cube file from a GDMA file where both the GDMA file and the reference Gaussian cube files are in the same global axis. The reference Gaussian cube files are used to determine the MEP grid parameters, and the GDMA MEP is compared to the reference cube files after the GDMA multipolar MEP cube file has been generated
**esp_charmm_grid** provides a similar workflow in 1-convert-pun2lpun.sh and 2-run-charmm.sh to generate an MEP cube file from GDMA multipole moments where the GDMA moments are in a different axis system to the reference cube file data. This is useful for example when checking the performance of GDMA multipole moments from the equilibrium geometry of a molecule in describing the MEP around the molecule in an arbitrary conformation from an MD simulation. To achieve this an axis system transformation is necessary, here we use CHARMM's MTPL module for this purpose. Note that the code requires rdkit to generate the CHARMM-format 'lpun' parameter file. Upon completion the GDMA-derived MEP cube file is again compared to the reference cube file to assess the agreement across different regions of the MEP grid.

### symmetry-constraints
* This example shows a symmetry-constrained fit for benzene. Fragment fits are not possible here, but the code automatically restricts fitting to symmetry-unique atomic environments only.
