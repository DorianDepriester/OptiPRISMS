import numpy as np
from scipy.spatial import Delaunay

def matrix_projection(nodes, pts):
    """
    Compute the projection matrix from Delaunay triangulation.

    Parameters
    ----------
    nodes : numpy.ndarray
        m x 2 table of node coordinates
    pts : numpy.ndarray
        p x 2 table of coordinates where one wants to evaluate the field

    Returns
    -------
    mat : numpy.array
        p x m matrix of projection coefficients. If a requested point is not in a triangle, the corresponding row will
        be NaN.
    inside_mesh : numpy.array
        Array of bools indicating whether the requested points are in the mesh
    """
    # Perform Delaunay triangulation and find in which triangles the pts belong to
    tri = Delaunay(nodes)
    triangle_ids = tri.find_simplex(pts)
    triangles = tri.simplices[triangle_ids]
    inside_mesh = triangle_ids != -1

    # Only work with points inside the mesh
    pts_in = pts[inside_mesh]
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
    newp = np.einsum('ijk,jk->ki', inv_mat, (pts_in - origin).T)

    # Local matrices of basis functions
    phi1 = 1 - np.sum(newp, axis=1, keepdims=True)
    basis_funcs = np.concatenate((phi1, newp), axis=1)

    # Assemble the global matrix
    row_indices = np.repeat(np.arange(toi.shape[0]), toi.shape[1])
    mat_in=np.zeros(shape=(len(toi), len(nodes)))
    mat_in[row_indices,toi.flatten()]=basis_funcs.flatten()
    mat = np.ones(shape=(len(pts), len(nodes)))*np.nan
    mat[inside_mesh, :] = mat_in

    return mat, inside_mesh


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
    mat, inside_mesh = matrix_projection(nodes, pts)
    return mat.dot(vect_field), inside_mesh
