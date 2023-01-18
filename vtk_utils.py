import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import os
import xml.etree.ElementTree as ET
from vtk.numpy_interface import dataset_adapter as dsa
from triangulate import triangular_projection

def readPointsFrom_pvtu(filename):
    reader = vtk.vtkXMLUnstructuredGridReader()
    pts = np.zeros(shape=(0, 2))
    u = np.zeros(shape=(0, 3))
    
    # Retrieve partitions
    if os.path.exists(filename):
        tree = ET.parse(filename)
        root = tree.getroot()
        foldername, filename = os.path.split(filename)
        for piece in root[0].iter('Piece'):
            vtk_file = piece.attrib['Source']
            # Read individual vtk files
            reader.SetFileName('{}/{}'.format(foldername, vtk_file))
            reader.Update()
            data = reader.GetOutput()
            points = data.GetPoints()
            pts_i = vtk_to_numpy(points.GetData())
            u_i = vtk_to_numpy(data.GetPointData().GetArray(0))
            
            # Remove points at z=0 plane and append data
            pts = np.concatenate((pts, pts_i[pts_i[:, 2] == 0., :2]))
            u = np.concatenate((u, u_i[pts_i[:, 2] != 0., :]), axis=0)
            
        # Remove duplicates
        nodes, idpts = np.unique(pts, axis=0, return_index=True)
        return nodes, u[idpts]
    
    else:
        # Return None if pvtu is not found
        return None, None
    
    
def merge_displacement_fields(mesh_vtk, input_pvtu, DIC_data, step):
    # Read mesh file and get node coordinates
    reader = vtk.vtkGenericDataObjectReader()
    reader.SetFileName(mesh_vtk)
    reader.Update()
    mesh = reader.GetOutput()
    meshNew = dsa.WrapDataObject(mesh)
    points_mesh = vtk_to_numpy(meshNew.GetPoints())
    
    # Only keep nodes at the surface of interest
    on_surface = points_mesh[:, 2] == 0.
    point_mesh_surf = points_mesh[on_surface, :2]
    
    # The displacement will be NaN everywhere, except on the surface of interest
    n_pts = len(points_mesh)
    u_SIM_tri = np.ones(shape=(n_pts,3))*np.nan
    nodes, u_SIM = readPointsFrom_pvtu(input_pvtu)
    for i in range(0, len(nodes)):
        # Nodes in mesh are not ordered the same way as in vtu files
        egal = np.argmin((points_mesh[:,0]- nodes[i,0])**2 + (points_mesh[:,1] - nodes[i,1])**2 + points_mesh[:,2]**2)
        u_SIM_tri[egal, :] = u_SIM[i,:]
    
    # Map DIC measurements onto the surface of interest
    data_DIC = np.loadtxt(DIC_data)
    pts_DIC = data_DIC[:,:2]
    u_DIC=data_DIC[:,(3*step-1):(3*step+1)]
    u_DIC_tri = np.zeros(shape=(n_pts,3))
    u_DIC_tri[on_surface, :2], _ = triangular_projection(pts_DIC, u_DIC, point_mesh_surf)
    u_DIC_tri[~on_surface, :] = np.nan
    
    # Map correlation coefficients as well
    C = data_DIC[:, 3*step+1]
    C_tri = np.zeros((n_pts,1))
    C_tri[on_surface,:], _ = triangular_projection(pts_DIC, C, point_mesh_surf)
    C_tri[~on_surface,:] = np.nan
    
    # Error function
    delta_u = np.zeros(shape=(n_pts,3))
    delta_u[:, :2] = u_SIM_tri[:, :2] - u_DIC_tri[:, :2]   
    
    # Save displacement fields in vtu file
    meshNew.PointData.append(u_SIM_tri, "U (FEM)")
    meshNew.PointData.append(u_DIC_tri, "U (DIC)")
    meshNew.PointData.append(delta_u, "Displacement error")
    meshNew.PointData.append(C_tri, "Correlation coefficient")
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName("Displacement_error_{}.vtu".format(step))
    writer.SetInputData(meshNew.VTKObject)
    writer.Write()    
    