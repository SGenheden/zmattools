"""
Classes to perform geometric calculations
"""
from __future__ import division, print_function, absolute_import

import numpy as np
from parmed.topologyobjects import Bond, Angle
from parmed.geometry import _get_coords_from_atom_or_tuple, _cross

def angle(a1, a2, a3) :
    """
    Calculate the angle between three distances

    Parameters
    ----------
    a1, a2, a3 : parmed.Atom
        the atoms

    Returns
    -------
    float
        the angle in degrees
    """
    return Angle(a1, a2, a3).measure()

def build_xyz(bndatm, angatm, dihedatm, blen, ang, dih) :
    """
    Build Cartesian coordinates from internal coordinates

    Parameters
    ----------
    bndatm : (3x1) ndarray
        the Cartesian coordinates of the atom that makes a bond
    angatm : (3x1) ndarray
        the Cartesian coordinates of the atom that makes an angle
    dihedatm : (3x1) ndarray
        the Cartesian coordinates of the atom that makes a dihedral
    blen : float
        the bond length
    ang : float
        the angle in degrees
    dih : float
        the dihedral in degrees

    Returns
    -------
    (3x1) ndarray
        the Cartesian coordinates of the built atom
    """

    ang = np.radians(ang)
    dih = np.radians(dih)
    blen = blen

    d1 = dihedatm - bndatm
    d2 = angatm   - bndatm

    vy = vnorm(np.cross(d1,d2))
    vx = vnorm(np.cross(d2,vy))
    vz = vnorm(np.cross(vx,vy))

    xbs =  blen * np.sin(ang) * np.cos(dih)
    ybs =  blen * np.sin(ang) * np.sin(dih)
    zbs = -blen * np.cos(ang)

    xyz = np.zeros(3)
    xyz[0] = vx[0]*xbs + vy[0]*ybs + vz[0]*zbs + bndatm[0]
    xyz[1] = vx[1]*xbs + vy[1]*ybs + vz[1]*zbs + bndatm[1]
    xyz[2] = vx[2]*xbs + vy[2]*ybs + vz[2]*zbs + bndatm[2]

    return xyz

def dihedral(a1, a2, a3, a4) :
    """
    Calculates the dihedral angle between four atoms
    Taken from parmed, but modified to get correct sign

    Parameters
    ----------
    a1, a2, a3, a4 : parmed.Atom
        the atoms

    Returns
    -------
    float
        the dihedral in degrees
    """

    p = np.array([_get_coords_from_atom_or_tuple(a1),
              _get_coords_from_atom_or_tuple(a2),
              _get_coords_from_atom_or_tuple(a3),
              _get_coords_from_atom_or_tuple(a4)])
    v1 = p[1] - p[0]
    v2 = p[1] - p[2]
    v3 = p[3] - p[2]
    # Take the cross product between v1-v2 and v2-v3
    v1xv2 = _cross(v1, v2)
    v2xv3 = _cross(v2, v3)
    # Now find the angle between these cross-products
    l1 = np.sqrt(np.dot(v1xv2, v1xv2))
    l2 = np.sqrt(np.dot(v2xv3, v2xv3))
    cosa = np.dot(v1xv2, v2xv3) / (l1 * l2)
    if np.dot(v3, np.cross(v1, v2)) <= 0.0 :
        return np.degrees(np.arccos(cosa))
    else :
        return -np.degrees(np.arccos(cosa))

def distance(a1, a2) :
    """
    Calculate the bond distance between two distances

    Parameters
    ----------
    a1, a2 : parmed.Atom
        the atoms

    Returns
    -------
    float
        the distance in Angstroms
    """
    return Bond(a1, a2).measure()

def majorminor_axis(xyz) :
    """
    Calculates the major and minor axis of a molecule
    
    Parameters
    ----------
    xyz : (nx3) ndarray
        all the Cartesian coordinates

    Returns
    -------
    (1x3) ndarray
        the major axis
    (1x3) ndarray
        the minor axis
    """
    def posneg_vec(vec1, vec2, posvec, negvec, npos, nneg) :

        vlen = vlen3(vec1)*vlen3(vec2)
        vec = (vec1*vec2).sum()
        ang = np.degrees(np.arccos(vec/vlen))
        if ang > 360.0 : ang = ang - 360.0
        if ang < 0.0 : ang = ang + 360.0
        if ang >= -90 and ang < 90 :
            posvec = posvec + vec2
            npos = npos + 1
        else :
            negvec = negvec + vec2
            nneg = nneg + 1
        return posvec, negvec, npos, nneg

    def assign_vec(posvec, negvec, npos, nneg) :

        vec = np.zeros(3)
        if npos > 0 :
            posvec = posvec / float(npos)
        if nneg > 0 :
            negvec = negvec / float(nneg)

        poslen2 = (posvec**2).sum()
        neglen2 = (negvec**2).sum()

        if poslen2 > neglen2 :
            if poslen2 != 0.0 :
                vec = vnorm(posvec)
            else :
                vec[0] = 1.0
        else :
            if neglen2 != 0.0 :
                vec = vnorm(negvec)
            else :
                vec[0] = 1.0
        return vec

    center = xyz.mean(axis=0)
    major = np.zeros(3)
    minor_tmp = np.zeros(3)

    # Major axis
    posvec = np.zeros(3)
    negvec = np.zeros(3)
    npos = 0
    nneg = 0
    for i in range(xyz.shape[0]) :
        vec = xyz[i,:] - center
        posvec,negvec,npos,nneg =  \
            posneg_vec(np.array([1.0,0.0,0.0]),vec, posvec, negvec, npos, nneg)

    major = assign_vec(posvec,negvec,npos,nneg)

    # Minor axis
    posvec = np.zeros(3)
    negvec = np.zeros(3)
    npos = 0
    nneg = 0
    for i in range(xyz.shape[0]) :
        vec = xyz[i,:] - center
        posvec,negvec,npos,nneg = posneg_vec(major,vec,posvec,negvec,npos,nneg,0)

    minor_tmp = assign_vec(posvec,negvec,npos,nneg)
    majlen2 = (major**2).sum()
    minlen2 = (minor_tmp**2).sum()

    if majlen2 < 1E-10 :
        major = np.array([1.0,0.0,0.0])
    elif minlen2 < 1E-10 :
        minor_tmp = np.array([0.0,1.0,0.0])

    if minlen2 > majlen2 :
        temp = np.array(major,copy=True)
        major = np.array(minor_tmp,copy=True)
        minor_tmp = np.array(temp,copy=True)

    if np.abs((major*minor_tmp).sum()) > 0.9 :
        major_swap = np.zeros(3)
        major_swap[0] = major[0]
        major_swap[1] = major[2]
        major_swap[2] = major[1]
        minor = vnorm(np.cross(major,major_swap))
    else :
        minor = vnorm(np.array(minor_tmp,copy=True))
        alpha = (major*minor).sum()
        beta = ((minor**2).sum() - (alpha**2))
        beta = np.sqrt(beta)
        minor = (minor - (alpha*major)) / beta
        minor = vnorm(minor)

    return major,minor

def vnorm(v) :
    """
    Returns a normalised vector
    """
    vlen = np.sqrt((v**2).sum())
    return v / vlen
