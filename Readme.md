# OptiPRISMS

Perform inverse analysis to retrieve Crystal Plasticity parameters used for CPFEM simulation through the PRISMS-Plasticity software [[1]](#prisms). This inverse analysis is done using data from in situ tensile tests and Scanning Electron Microscopy Digital Image Correlation (SEM-DIC).

## How it works

This program will optimize CP parameters in order to minimize a given cost function. At each step of the optimization process, it will:

1. generate configuration files for PRISMS-Pasticity;
2. run PRISMS-Pasticity from these files;
3. read results and compare them against experimental data;
4. return a cost function, evidencing *how far* the simulation is, compared to the experiment,
5. remove simulation files.

In order to generate the PRISMS configuration files, this software will parse *template* files in order to retrieve which parameters it shall optimize.

Records for each optimization loop will be stored in a log file, so that the user can track the evolution of optimized parameters and cost functions. In addition, OptiPRISMS will use this file to avoid reruning preexisting simulations, so that the optimization can be stopped then resumed.

## Required materials

- The mesh of the microstructure (in .msh format). One can use [MTEX2Gmsh](https://doriandepriester.github.io/MTEX2Gmsh/) [[2]](#mtex2gmsh) to generate a conforming mesh directly from EBSD data.
- [PRISMS-Plasticity software](https://github.com/prisms-center/plasticity).
- Python 3.6 (or later) with [vtk module](https://pypi.org/project/vtk/).
- Experimental data, consisting in:
    - a macroscopic tensile curve (strain-stress values as a CSV file),
	- SEM-DIC displacement measurements, stored as tabular data.
	
The latter must be formated as column data, in this order:
1. x coordinates of DIC points where DIC measurements are performed,
2. y coordinates of DIC points where DIC measurements are performed,
3. x displacements at 1st step,
3. y displacements at 1st step,
4. correlation coefficients at 1st step,
5. x displacements at 2nd step,
6. y displacements at 2nd step,
7. correlation coefficients at 2nd step

and so on.

## Step-by-step method to run optimization

1. Check out the template files (under the eponymous folder) and give a proper name to every parameter you want to optimize. These names must be precessed by a "$" symbol (e.g. "$a" instead of "a").
2. Be sure that common parameters (e.g. path to mesh file, time increement, etc.) in template files fit with your needs.
3. Edit the [configuration file] to tune optimization-related parameters.
4. With Python3, run the ``OptiPRISMS`` function. E.g.:
````
from OptiPRISMS import runOptim

res = runOptim(config_file='myConfigFile.ini')
````

## Configuration file

This file describes locations of data and parameters for optimization. It divides in sections (some of them are mandatory, other are not).

### Mandatory sections
#### [Initial Guess]

Provide here starting values for the variables you want to optimize. The variables names must be consistent with these in templates for PRISMS parameter files (as detailed [here](step-by-step-method-to-run-optimization)).

#### [Bounds]

- **lower**: list of lower bounds
- **upper**: list of upper bounds

#### [PRISMS]

- **prm file**: path to the main parameter file template
- **latent hardening ratio** : path to latent hardening ratio tempalte
- **path to prisms**= path to PRISMS-Plasticity executable file. This entry is not required if Slurm is used (see [\[Slurm\]] section).

#### [Experimental Data]

- **DIC data**: path to text file with DIC measurements
- **tensile curve**: path to strain-stress values of tensile curve
- **time steps** *(optional)*: increment numbers corresponding to each step of DIC measurements. If not set, they will infered from option ``set Tabular Time Output Table`` in prm template.

#### [Cost Function]

- **weight on tensile curve**: weight to apply to the tensile curve in the overall cost function
- **penalty**: penalty value to raise if the simulation fails

#### [Log File]

- **file path**: path to log file

### Optional sections

#### [Minimize]

Pass here any optional parameter(s) for the ``options`` argument of ``minimize_parallel``. See [the options for ``scipy.optimize.minimize`` with L-BFGS-B method for available arguments](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html#optimize-minimize-lbfgsb).

#### [Minimize parallel]

This section allows to use a parallel implementation of the L-BFGS-B minimizer [[3]](#optim_parallel), through the ``minimize_parallel`` function (from [``optimparallel``](https://github.com/florafauna/optimParallel-python) module).

- **use parallel minimizer**: whether to use ``minimize_parallel``, instead of ``scipy.optimize.minimize`` (default is No)

In addition, any keyword argument normaly passed to the ``parallel`` option of ``minimize_parallel`` can be defined in this section. See [the related documentation](https://github.com/florafauna/optimParallel-python/blob/8bf622be1431ba10fef1d795521a2b1d86307c9d/src/optimparallel.py#L170) for available options.

#### [Slurm]

- **use Slurm**: whether to use the Slurm workload manager. Default is No.
- **batch file**: path to batch file to use for submitting a job running PRISMS-Pasticity. The path to the prm file will be passed as the first argument of this script.

#### [Debug]

- **fake simulations**: If Yes, each simulation is not computed. The related commands are printed instead. Default is No.
- **fake deletions**: Turn off automatic removal of simulation results (step 5 in [How it works] section). Default is No.

## References
<a id="prisms">[1]</a> Yaghoobi et al., (2019). Prisms-plasticity: An open-source crystal plasticity finite element software. Computational Materials Science, 169:109078, https://doi.org/10.1016/j.commatsci.2019.109078

<a id="mtex2gmsh">[2]</a> Depriester et al., (2020). MTEX2Gmsh: a tool for generating 2D meshes from EBSD data. Journal of Open Source Software, 5(52):2094, https://doi.org/10.21105/joss.02094

<a id="optim_parallel">[3]</a> Gerber, F. and Furrer, R. (2019). optimParallel: An R package providing a parallel version of the L-BFGS-B optimization method, The R Journal, 11(1):352–358, http://doi.org/10.32614/RJ-2019-030