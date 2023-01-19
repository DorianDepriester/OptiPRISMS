import numpy as np
from scipy.spatial import Delaunay


def triangular_projection(nodes, vect_field, pts):
    """
    Given a field known at certain nodes in a 2D plane, estimate the field values at other 2D coordinates by triangular
    projection.

    Parameters
    ----------
    nodes : numpy.ndarray
        m x 2 table of node coordinates
    vect_field : numpy.ndarray
        m x n array defining the n-dimensional vector field at nodes coordinates
    pts : numpy.ndarray
        p x 2 table of coordinates where one wants to evaluate the field

    Returns
    -------
    u_tri : numpy.ndarray
        p x n array of interpolated field values. The values are np.nan for all points not inside the mesh.
    inside_mesh : numpy.ndarray
        Array of bools indicating whether the requested points are in the mesh
    """
    # Perform Delaunay triangulation and find in which triangles the pts belong to
    tri = Delaunay(nodes)
    triangle_ids = tri.find_simplex(pts)
    triangles = tri.simplices[triangle_ids]
    inside_mesh = triangle_ids != -1

    # Only work with points inside the mesh
    pts = pts[inside_mesh]
    toi = triangles[inside_mesh, :]  # Triangles of Interest

    # Compute local coordinates of pts, wrt. their parent triangles
    # Adapted from https://stackoverflow.com/a/57901916/12056867
    origin = nodes[toi[:, 0], :]
    v1 = nodes[toi[:, 1], :] - origin
    v2 = nodes[toi[:, 2], :] - origin
    v1r = v1.T[:, np.newaxis, :]
    v2r = v2.T[:, np.newaxis, :]
    mat = np.concatenate((v1r, v2r), axis=1)
    inv_mat = np.linalg.inv(mat.T).T
    newp = np.einsum('ijk,jk->ki', inv_mat, (pts - origin).T)

    # Matrices of basis functions
    phi1 = 1 - np.sum(newp, axis=1, keepdims=True)
    basis_funcs = np.concatenate((phi1, newp), axis=1)

    if len(vect_field.shape) == 1:
        # Force field to be at least a column vector
        vect_field = vect_field[:, np.newaxis]

    # The displacement will be nan if the requested point is outside the mesh
    u_tri = np.ones((len(pts), vect_field.shape[1])) * np.nan

    # Project node displacements with the aid of the matrix of basis function
    u_tri[inside_mesh] = np.einsum('ik,ikj->ij', basis_funcs, vect_field[toi, :])

    return u_tri, inside_mesh
