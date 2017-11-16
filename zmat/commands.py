"""
Commands for the zmat tool
"""
from __future__ import division, print_function, absolute_import

import argparse
import sys
import os

import numpy as np
from parmed.periodic_table import Element

import zmat.mol as mol
import zmat.geo as geo

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

class ZmatTraj(ZmatCommand) :

    @staticmethod
    def descr() :
        return "Generate a trajectory with all values of angles and dihedrals"

    @classmethod
    def command_name(cls) :
        return "traj"

    def _add_arguments(self, parser) :
        parser.add_argument('file',help="the structure file")
        parser.add_argument('-z','--zmat',help="the z-mat if it exists")
        parser.add_argument('-d', '--dihedrals',type=int,nargs="+",help="defines the increment, range or sequence of dihedral angles",default=[15])

    def _run(self, args) :

        graph = mol.GraphStructure(verbosity=1)
        graph.initialize(args.file)
        if args.zmat is None :
            graph.traverse()
        else :
            graph.read_zmat(args.zmat)
        zval0 = graph.zvalues()
        geo.internal2cartesian(zval0, graph.zatoms, make_dummies=True)

        if len(args.dihedrals) == 1 :
            dihedrals = range(-180,180+args.dihedrals[0],args.dihedrals[0])
        elif len(args.dihedrals) == 3 :
            dihedrals = range(args.dihedrals[0],args.dihedrals[1]+args.dihedrals[2],
                                args.dihedrals[2])
        else :
            dihedrals = args.dihedrals
        print("Will produce trajectories over %d dihedral angles"%len(dihedrals))
        unique_dihedrals = []
        for i, row in enumerate(graph.zmat[3:],3) :
            dih = row[1]+"-"+row[2]
            dih2 = row[2]+"-"+row[1]
            if dih in unique_dihedrals or dih2 in unique_dihedrals :
                print("Skipping dihedral %d: %s it is non-unique and probably an improper"%(
                    i+1, " ".join(row)))
                continue
            else :
                unique_dihedrals.append(dih)
                print("Producing trajectories for dihedral %d: %s"%(i+1, " ".join(row)))
            zval = np.array(zval0, copy=True)
            with open("%s_dihed_traj%d.xyz"%(os.path.splitext(args.file)[0],
                                                i+1), "w") as f :
                for j, dihedral in enumerate(dihedrals,1) :
                    zval[i, 2] = dihedral
                    geo.internal2cartesian(zval, graph.zatoms)
                    f.write("%d\n%d\n"%(len(graph.zmat),j))
                    for atom in graph.structure.atoms :
                        f.write("%s %8.3f %8.3f %8.3f\n"%(Element[atom.element],
                            atom.xx, atom.xy, atom.xz))
