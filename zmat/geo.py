"""
Classes to perform geometric calculations
"""
from __future__ import division, print_function, absolute_import

import numpy as np
from parmed.topologyobjects import Atom
from parmed.geometry import _get_coords_from_atom_or_tuple as atom2coords

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
    if dih < 0 :
        dih += 360
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

def internal2cartesian(zvalues, zatoms, make_dummies=False) :
    """
    Transform internal coordinates to Cartesian coordinates

    Updates the zatoms xx, xy and xz properties

    Parameters
    ----------
    zvalues : (nx3) ndarray
        the values of the internal coordinates
    zatoms : list of parmed.Atom
        the atoms with the Cartesian coordinates
    dummies : bool, optional
        make dummy atoms and insert them in the zvalues and zatoms arrays
    """

    if make_dummies :
        xyz0 = np.asarray([atom2coords(atoms[0]) for atoms in zatoms])
        major, minor = majorminor_axis(xyz0)
        center = xyz0.mean(axis=0)

        dm3 = Atom(name="DM3")
        dm3.xx = center[0]+minor[0]
        dm3.xy = center[1]+minor[1]
        dm3.xz = center[2]+minor[2]

        dm2 = Atom(name="DM2")
        dm2.xx = center[0]+major[0]
        dm2.xy = center[1]+major[1]
        dm2.xz = center[2]+major[2]

        dm1 = Atom(name="DM1")
        dm1.xx = center[0]
        dm1.xy = center[1]
        dm1.xz = center[2]

        zatoms[0] = [zatoms[0][0], dm3, dm2, dm1]
        zatoms[1] = [zatoms[1][0], zatoms[1][1], dm3, dm2]
        zatoms[2] = [zatoms[2][0], zatoms[2][1], zatoms[2][2], dm3]

        zvalues[0,0] = distance(zatoms[0][0], zatoms[0][1])
        zvalues[0,1] = angle(zatoms[0][0], zatoms[0][1], zatoms[0][2])
        zvalues[0,2] = dihedral(zatoms[0][0], zatoms[0][1], zatoms[0][2], zatoms[0][3])
        zvalues[1,1] = angle(zatoms[1][0], zatoms[1][1], zatoms[1][2])
        zvalues[1,2] = dihedral(zatoms[1][0], zatoms[1][1], zatoms[1][2], zatoms[1][3])
        zvalues[2,2] = dihedral(zatoms[2][0], zatoms[2][1], zatoms[2][2], zatoms[2][3])

    for i, (zval, atoms) in enumerate(zip(zvalues,zatoms)) :
        xyz = build_xyz(np.asarray(atom2coords(atoms[1])),
                            np.asarray(atom2coords(atoms[2])),
                            np.asarray(atom2coords(atoms[3])),
                            zval[0], zval[1], zval[2])
        atoms[0].xx = xyz[0]
        atoms[0].xy = xyz[1]
        atoms[0].xz = xyz[2]

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

        vlen = np.sqrt((vec1**2).sum())*np.sqrt((vec2**2).sum())
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
        posvec,negvec,npos,nneg = posneg_vec(major,vec,posvec,negvec,npos,nneg)

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
