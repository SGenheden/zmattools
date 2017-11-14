"""
Commands for the zmat tool
"""
from __future__ import division, print_function, absolute_import

import argparse
import sys

from parmed.periodic_table import Element

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
            for i, (atoms, zval) in enumerate(zip(graph.zatoms, graph.zvalues())) :
                sys.stdout.write("%s "%Element[atoms[0].element])
                if i > 0 :
                    sys.stdout.write("%3d %5.3f" % (graph.znames.index(atoms[1].name) + 1, zval[0]))
                if i > 1:
                    sys.stdout.write("%3d %8.3f" % (graph.znames.index(atoms[2].name) + 1, zval[1]))
                if i > 2:
                    sys.stdout.write("%3d %8.3f" % (graph.znames.index(atoms[3].name) + 1, zval[2]))
                print("")
