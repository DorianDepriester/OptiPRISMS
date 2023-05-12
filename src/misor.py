from triangulate import project_orientation
import numpy as np
import pandas as pd
from orix.quaternion import symmetry
from orix.quaternion.orientation import Orientation
from vtk_utils import create_vtu_from_field

s = symmetry.Oh
fpath = '../example/2.0_150.0_20.0_150.0_3.0/QuadratureOutputs1999.csv'
col_ids=(0,4,5,6,7,8,9)
col_names=('grainID', 'x', 'y', 'z', 'rodx', 'rody', 'rodz')
quadrature =pd.read_csv(fpath, header=None, usecols=col_ids, names=col_names)
z_max=np.max(quadrature['z'])
quadrature_surf = quadrature[quadrature['z']==z_max]

rodrigues = quadrature_surf[['rodx', 'rody', 'rodz']].to_numpy()
mag = np.sqrt(np.sum(rodrigues ** 2, axis=1))
angle = 2*np.arctan(mag)
axis = quadrature_surf[['rodx', 'rody', 'rodz']].to_numpy()
axis = (axis.T / mag).T
o = Orientation.from_axes_angles(-axis, angle, symmetry=s)

xy = quadrature_surf[['x', 'y']].to_numpy()
op = project_orientation(xy, o, np.array([1,1]))

ebsd_file = '../example/ExperimentalResults/EulerAngles_def.txt'
ebsd_data = pd.read_csv(ebsd_file, sep=',')
euler = ebsd_data[['phi1', 'Phi', 'phi2']].to_numpy()
xy_ebsd = ebsd_data[['x', 'y']].to_numpy()
xy_ebsd = xy_ebsd-np.min(xy_ebsd, axis=0)
o_ebsd = Orientation.from_euler(euler, degrees=True, symmetry=s, direction ='crystal2lab')

downsamp = 5
xy_ebsd = xy_ebsd[::downsamp]
o_ebsd = o_ebsd[::downsamp]
xy = xy[::downsamp]
o = o[::downsamp]

op = project_orientation(xy_ebsd, o_ebsd, xy)
#theta = np.ones(shape=xy.shape[0])*np.nan
#for k, opk in enumerate(op):
#    if not any(np.isnan(opk.data[0])):
#        m = Misorientation((o[k].data[0], opk.data[0]), symmetry=(s,s))
#        theta[k] = m.get_distance_matrix(progressbar=False, degrees=True)[1,0]

#not_nan = np.logical_not(np.isnan(theta))
#theta_clean = theta[not_nan]
#xy_clean = xy[not_nan, :]
#create_vtu_from_field('../example/mesh.vtk', xy_clean, theta_clean, 'misorientation1.vtu', 'Misorientation (deg)')

theta = np.ones(shape=xy.shape[0])*np.nan
for k, opk in enumerate(op):
    if not any(np.isnan(opk.data[0])):
        ok = Orientation((o.data[k],op.data[k]),symmetry=o.symmetry)
        theta[k] = ok.get_distance_matrix(degrees=True)[0,1]
not_nan = np.logical_not(np.isnan(theta))
theta_clean = theta[not_nan]
xy_clean = xy[not_nan, :]
create_vtu_from_field('../example/mesh.vtk', xy_clean, theta_clean, 'misorientation5.vtu', 'Misorientation (deg)')