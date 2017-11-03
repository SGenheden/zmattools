"""
Commands for the zmat tool
"""
from __future__ import division, print_function, absolute_import

import argparse
import sys

from parmed.periodic_table import Element
from parmed.topologyobjects import Bond, Angle, Dihedral

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
                    sys.stdout.write("%2d %5.3f" % (atoms.index(zatoms[1]) + 1, bond.measure()))
                if i > 1:
                    a3 = graph.structure.view["@%s" % zatoms[2]].atoms[0]
                    angle = Angle(a1, a2, a3)
                    sys.stdout.write("%2d %8.3f" % (atoms.index(zatoms[2]) + 1, angle.measure()))
                if i > 2:
                    a4 = graph.structure.view["@%s" % zatoms[3]].atoms[0]
                    dihedral = Dihedral(a1, a2, a3, a4)
                    sys.stdout.write("%2d %8.3f" % (atoms.index(zatoms[3]) + 1, dihedral.measure()))
                print("")
