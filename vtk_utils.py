import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import os
from xml.etree.ElementTree import parse
from vtk.numpy_interface import dataset_adapter as dsa
from triangulate import triangular_projection


def read_pvtu(filename):
    """
    Read node locations and node displacements from parallel vtu files (pvtu)

    Parameters
    ----------
    filename : str
        Path to pvtu file

    Returns
    -------
    nodes : np.ndarray
        m x 3 array of node coordinates. None if the file is not found.
    u   : np.ndarray
        m x 3 array of node displacements. None if the file is not found.
    """
    reader = vtk.vtkXMLUnstructuredGridReader()
    pts = np.zeros(shape=(0, 2))
    u = np.zeros(shape=(0, 3))
    
    # Retrieve partitions
    if os.path.exists(filename):
        tree = parse(filename)
        root = tree.getroot()
        folder_name, filename = os.path.split(filename)
        for piece in root[0].iter('Piece'):
            vtk_file = piece.attrib['Source']
            # Read individual vtk files
            reader.SetFileName(os.path.join(folder_name, vtk_file))
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
    mesh_new = dsa.WrapDataObject(mesh)
    points_mesh = vtk_to_numpy(mesh_new.GetPoints())
    
    # Only keep nodes at the surface of interest
    on_surface = points_mesh[:, 2] == 0.
    point_mesh_surf = points_mesh[on_surface, :2]
    
    # The displacement will be NaN everywhere, except on the surface of interest
    n_pts = len(points_mesh)
    u_sim_tri = np.ones(shape=(n_pts, 3)) * np.nan
    nodes, u_sim = read_pvtu(input_pvtu)
    for i in range(0, len(nodes)):
        # Nodes in mesh are not ordered the same way as in vtu files
        dist = (points_mesh[:, 0] - nodes[i, 0])**2 + (points_mesh[:, 1] - nodes[i, 1])**2 + points_mesh[:, 2] ** 2
        egal = np.argmin(dist)
        u_sim_tri[egal, :] = u_sim[i, :]
    
    # Map DIC measurements onto the surface of interest
    data_dic = np.loadtxt(DIC_data)
    pts_dic = data_dic[:, :2]
    u_dic = data_dic[:, (3*step-1):(3*step+1)]
    u_dic_tri = np.zeros(shape=(n_pts, 3))
    u_dic_tri[on_surface, :2], _ = triangular_projection(pts_dic, u_dic, point_mesh_surf)
    u_dic_tri[~on_surface, :] = np.nan
    
    # Map correlation coefficients as well
    c = data_dic[:, 3*step+1]
    c_tri = np.zeros((n_pts, 1))
    c_tri[on_surface, :], _ = triangular_projection(pts_dic, c, point_mesh_surf)
    c_tri[~on_surface, :] = np.nan
    
    # Error function
    delta_u = np.zeros(shape=(n_pts, 3))
    delta_u[:, :2] = u_sim_tri[:, :2] - u_dic_tri[:, :2]
    
    # Save displacement fields in vtu file
    mesh_new.PointData.append(u_sim_tri, "U (FEM)")
    mesh_new.PointData.append(u_dic_tri, "U (DIC)")
    mesh_new.PointData.append(delta_u, "Displacement error")
    mesh_new.PointData.append(c_tri, "Correlation coefficient")
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName("Displacement_error_{}.vtu".format(step))
    writer.SetInputData(mesh_new.VTKObject)
    writer.Write()    
    