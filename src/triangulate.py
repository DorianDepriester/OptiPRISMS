import numpy as np
from scipy.spatial import Delaunay
from orix.quaternion.orientation import Orientation
from orix.quaternion import Quaternion


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
    # If only one point is requested, treat it as a 2D array anyway
    if len(np.array(pts).shape)==1:
        pts = pts[np.newaxis, :]

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


def project_orientation(ebsd_locations, ebsd_orientations, pts):
    """
    Given of field of orientations, measured at ebsd locations, compute apply triangular projection of requested
    locations to estimate the local orientation at these points.

    Parameters
    ----------
    ebsd_locations : numpy.ndarray
        m x 2 table of coordinate where the orientations are measures (EBSD pixels)
    ebsd_orientations : orix.quaternion.orientation.Orientation
        array of length m of measured orientations
    pts : numpy.array
        p x 2 table of coordinates where one wants to projection the orientations

    Returns
    -------
    orix.quaternion.orientation.Orientation
        array of length p of projected orientations. The orientation is NaN if the requested point is outside the mesh.
    """
    if len(np.array(pts).shape)==1:
        pts = pts[np.newaxis, :]

    mat, inside_mesh = matrix_projection(ebsd_locations, pts)
    mat_red = mat[inside_mesh,:]
    q_mean = np.ones((len(pts),4))*np.nan
    o2 = ebsd_orientations.map_into_symmetry_reduced_zone()
    q = o2.data
    qq = np.einsum('pi,ij,ik->pjk', mat_red, q, q)
    w, v = np.linalg.eig(qq)
    w_max = np.argmax(w, axis=1)
    q_mean_red = v[np.arange(mat_red.shape[0]), :, w_max]
    q_mean[inside_mesh,:] = q_mean_red
    return Orientation(Quaternion(q_mean), ebsd_orientations.symmetry)
