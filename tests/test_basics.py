
from __future__ import division, print_function, absolute_import

import os
import unittest

import zmat
import numpy as np
import parmed
from parmed import geometry
from parmed.geometry import _get_coords_from_atom_or_tuple as atom2coords

pdb_filename = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..","examples","model.pdb")

class TestZmat(unittest.TestCase):

    def test_generatezmat(self) :
        """ Test the full generation of a zmatrix"""
        saved_zmat = [["C4", "DM3", "DM2", "DM1"],
                      ["C6", "C4", "DM3", "DM2"],
                      ["C5", "C6", "C4", "DM3"],
                      ["H10", "C5", "C4", "C6"],
                      ["C3", "C4", "C6", "C5"],
                      ["C2", "C3", "C4", "C6"],
                      ["H21", "C2", "C3", "C4"],
                      ["H41", "C4", "C6", "C3"],
                      ["C7", "C6", "C4", "C5"],
                      ["H73", "C7", "C6", "C4"],
                      ["H71", "C7", "C6", "H73"],
                      ["H72", "C7", "C6", "H73"],
                      ["H61", "C6", "C4", "C5"],
                      ["H31", "C3", "C4", "C2"],
                      ["H32", "C3", "C4", "C2"],
                      ["H51", "C5", "C4", "C6"],
                      ["H22", "C2", "C3", "H21"],
                      ["H23", "C2", "C3", "H21"]]
        graph = zmat.mol.GraphStructure(verbosity=0)
        graph.initialize(pdb_filename)
        graph.traverse()
        self.assertEqual(graph.zmat, saved_zmat)

    def test_readzmat(self) :
        """ The the functionality to read a z-matrix from disc """
        graph = zmat.mol.GraphStructure(verbosity=0)
        graph.initialize(pdb_filename)
        graph.traverse()
        with open("temp.zmat", "w") as f :
            for row in graph.zmat:
                f.write(" ".join(row)+"\n")

        graph2 = zmat.mol.GraphStructure(verbosity=0)
        with self.assertRaises(Exception) :
            graph2.read_zmat("temp.zmat")
        graph2.initialize(pdb_filename)
        graph2.read_zmat("temp.zmat")
        self.assertEqual(graph2.zmat, graph.zmat)

        os.remove("temp.zmat")

    def test_buildxyz(self) :
        """ Test the geo.build_xyz routine """
        struct = parmed.load_file(pdb_filename)
        a1 = struct.view["@C5"].atoms[0]
        a2 = struct.view["@C4"].atoms[0]
        a3 = struct.view["@C3"].atoms[0]
        a4 = struct.view["@C2"].atoms[0]
        new_xyz = zmat.geo.build_xyz(np.asarray(atom2coords(a3)),
            np.asarray(atom2coords(a2)), np.asarray(atom2coords(a1)),
            np.sqrt(geometry.distance2(a4, a3)),
            geometry.angle(a4, a3, a2), geometry.dihedral(a4, a3, a2, a1))
        self.assertAlmostEqual(new_xyz[0], a4.xx, 3)
        self.assertAlmostEqual(new_xyz[1], a4.xy, 3)
        self.assertAlmostEqual(new_xyz[2], a4.xz, 3)

    def test_internal2cartesian(self) :
        """ Test the geo.internal2cartesian routine """
        graph = zmat.mol.GraphStructure(verbosity=0)

        with self.assertRaises(Exception) :
            _ = graph.zvalues()

        graph.initialize(pdb_filename)
        graph.traverse()
        zval = graph.zvalues()
        xyz0 = np.asarray(graph.structure.coordinates)

        zmat.geo.internal2cartesian(zval, graph.zatoms, make_dummies=True)
        self.assertAlmostEqual(xyz0[0,0], graph.structure.coordinates[0,0], 3)
        self.assertAlmostEqual(xyz0[-1,0], graph.structure.coordinates[-1,0], 3)

        zval[0,2] += 180
        zmat.geo.internal2cartesian(zval, graph.zatoms)
        self.assertNotAlmostEqual(xyz0[0,0], graph.structure.coordinates[0,0], 3)
        self.assertNotAlmostEqual(xyz0[-1,0], graph.structure.coordinates[-1,0], 3)

        zval[0,2] += 180
        zmat.geo.internal2cartesian(zval, graph.zatoms)
        self.assertAlmostEqual(xyz0[0,0], graph.structure.coordinates[0,0], 3)
        self.assertAlmostEqual(xyz0[-1,0], graph.structure.coordinates[-1,0], 3)

if __name__ == '__main__':

    unittest.main()
