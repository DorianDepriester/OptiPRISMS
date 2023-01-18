import numpy as np
from scipy.spatial import Delaunay

def triangular_projection(nodes, vect_field, pts):
    """
    Perform triangular projection of displacement fields given on nodes
    onto a set of points, possibly inconsistent with node coordinates.

    Parameters
    ----------
    nodes : numpy.ndarray
        m x 2 table of node coordinates
    vect_field : numpy.ndarray
        m x n array defining the vector field on wants to interpolate
    pts : numpy.ndarray
        p x 2 table of requested points

    Returns
    -------
    u_tri : numpy.ndarray
        p x n array of interpolated field values. The values are np.nan for 
        point not inside the mesh.
    insise_mesh : numpy.ndarray
        Array of bools indicating whether the requested points are in the mesh 

    """
    # Perform Delaunay triangulation and find in which triangle the pts belong to
    tri = Delaunay(nodes)
    triangle_ids = tri.find_simplex(pts)
    triangles = tri.simplices[triangle_ids]
    inside_mesh = triangle_ids != -1
    
    # Only work with points inside the mesh
    n_pts = len(pts)
    pts = pts[inside_mesh]
    ToI = triangles[inside_mesh,:]   # Triangles of Interest
    
    # Compute local coordinates of pts, wrt. their parent triangles
    # Adapted from https://stackoverflow.com/a/57901916/12056867
    A = nodes[ToI[:,0],:]
    B = nodes[ToI[:,1],:]
    C = nodes[ToI[:,2],:]
    v1 = B-A
    v2 = C-A
    v1r=v1.T[:,np.newaxis,:]
    v2r=v2.T[:,np.newaxis,:]
    mat = np.concatenate((v1r,v2r), axis=1)
    inv_mat = np.linalg.inv(mat.T).T
    newp=np.einsum('ijk,jk->ki',inv_mat,(pts-A).T)
    
    # Matrix of basis functions
    phi1 = 1 - np.sum(newp, axis=1, keepdims=True)
    B=np.concatenate((phi1, newp), axis=1)
    
    # The displacement will be nan if the requested point is outside the mesh
    if len(vect_field.shape) == 1:
        # Force field to be at least of column vector
        vect_field = vect_field[:,np.newaxis]
    u_tri = np.ones((n_pts, vect_field.shape[1]))*np.nan
    
    # Project node dispacements with the aid of the matrix of basis function
    u_tri[inside_mesh]=np.einsum('ik,ikj->ij', B, vect_field[ToI,:])

    return u_tri, inside_mesh


