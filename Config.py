"""Initial guess and bounds on the parameters to be optimized"""
#  [      q,   h0,        s0,    sinf,    a]
x0=[      2,   150,       20,    150,     3]
bounds=( (1,5),(100,2000),(0,50),(50,250),(1,10) )

"""Maximum number of parallel evaluations"""
max_workers=6
# This number must be lower than the number of CPUs of the 'master' node
# (where run_Opti.py is launched) even if the corresponding simulations are 
# sent to other nodes.

"""Location of log file"""
logfile='Optimization.csv'
# This file will be used to track evolution of the minimizer iterations. In 
# addition, it will be used to avoid redundant calculations, thus allowing to 
# stop and restart the optimization.

"""Path to DIC results"""
DIC_data='ExperimentalResults/DICdata.csv'
# The data must be saved as text file, where colomns are ordered this way:
#   1. x coordinate of nodes (undeformed configuration)
#   2. y coordinate of nodes (undeformed configuration)
#   3. Ux (at 1st step)
#   4. Uy (at 1st step)
#   5. correlation coefficient (at 1st step)
#   6. Ux (at 2nd step)
#   7. Uy (at 2nd step)
#   8. correlation coefficient (at 2nd step)
# and so on.

"""Minimal distance to check if the DIC data are consistent with mesh nodes"""
coincident_distance=1e-5
# This value will be used to detect an error in input data (inconsistent
# locations between DIC measurements and mesh nodes). If the distance between
# one node and the corresponding DIC location is larger than this value, an
# error is raided. Set a negative value to turn the checking routine off.

"""Path to experimental tensile curve"""
tensileCurve=('ExperimentalResults/TensileCurve.csv')

"""Whether the border of the RoI should be considered for computing the
kinematic cost function"""
remove_borders=True

"""Weight applied to the static cost function"""
w_sigma=0.5
# w_sigma=0 means that only the kinematic error will be taken into account
# w_sigma=1 means that only the static error will be taken into account

"""Time steps to be compared with DIC"""
time_steps=(999, 1999, 2999)
# List of time increments where results should be compared against DIC data.
# The values will depend on the time increment used in PRISMS-Plasticity. See
# line 'set Time increments' in Template/prm.prm.

"""Penalty value when simulation fails"""
penalty=5
