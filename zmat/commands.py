"""
Commands for the zmat tool
"""
from __future__ import division, print_function, absolute_import

import argparse
import sys

import numpy as np
from parmed.periodic_table import Element
from parmed.topologyobjects import Bond, Angle, Dihedral
from parmed.geometry import _get_coords_from_atom_or_tuple, _cross

import zmat.mol as mol

class ZmatCommand(object) :

    @staticmethod
    def descr() :
        return "This is generic description"

    @classmethod
    def command_name(cls) :
        return cls.__name__.lower()

    def execute(self) :

        description = self.descr()
        prog = sys.argv[0] + " " + sys.argv[1]
        sys.argv.pop(1)
        argparser = argparse.ArgumentParser(description=description,
                            add_help=True, prog=prog)
        self._add_arguments(argparser)
        self._run(argparser.parse_args())

class GenerateZmat(ZmatCommand) :

    @staticmethod
    def descr() :
        return "Generate a z-matric from a structure"

    @classmethod
    def command_name(cls) :
        return "generate"

    def _add_arguments(self, parser) :
        parser.add_argument('file',help="the structure file")
        parser.add_argument('-v','--verbosity',type=int,help="the verbosity level",default=1)
        parser.add_argument('--gaussian',action="store_true",help="generate Gaussian input",default=False)

    def _run(self, args) :

        graph = mol.GraphStructure(verbosity=args.verbosity)
        graph.initialize(args.file)
        graph.traverse()

        if args.gaussian :
            atoms = [z.split()[0] for z in graph.zmat]
            for i, z in enumerate(graph.zmat):
                zatoms = z.strip().split()
                a1 = graph.structure.view["@%s" % zatoms[0]].atoms[0]
                sys.stdout.write("%s "%Element[a1.atomic_number])
                if i > 0:
                    a2 = graph.structure.view["@%s" % zatoms[1]].atoms[0]
                    bond = Bond(a1, a2)
                    sys.stdout.write("%3d %5.3f" % (atoms.index(zatoms[1]) + 1, bond.measure()))
                if i > 1:
                    a3 = graph.structure.view["@%s" % zatoms[2]].atoms[0]
                    angle = Angle(a1, a2, a3)
                    sys.stdout.write("%3d %8.3f" % (atoms.index(zatoms[2]) + 1, angle.measure()))
                if i > 2:
                    a4 = graph.structure.view["@%s" % zatoms[3]].atoms[0]
                    sys.stdout.write("%3d %8.3f" % (atoms.index(zatoms[3]) + 1, self._dihedral(a1,a2,a3,a4)))
                print("")

    def _dihedral(self, a1, a2, a3, a4) :
        """
        Calculates the dihedral angle between four atoms
        Taken from parmed, but modified to get correct sign
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
