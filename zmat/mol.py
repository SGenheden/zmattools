"""
Classes to manipulate molecules

Parts of this code is taken from ProtoMS ver 3.0, written by
Richard Bradshaw, Ana Cabedo Martinez, Chris Cave-Ayland, Samuel Genheden
 and Gregory Ross
"""
from __future__ import division, print_function, absolute_import

import sys

import networkx as nx
import parmed
from parmed.structure import Structure

class NodeAtom(object):
    """
    Class to encapsulate an atom its connectivity

    Attributes
    ----------
    bonds : list of strings
      names of all atoms this atom is bonded to
    is_on_loop : dictionary
      a flag for each atom in bonds indicating if the bond is part of a loop
    name : string
      the atom name
    traversed : dictionary
      a flag for each atom in bonds indicating if the bond has been traversed
    zmat : list of strings
      atom names of the atoms defined to this atom in the z-matrix
    """

    def __init__(self):
        self.name = ""
        self.bonds = []
        self.traversed = {}
        self.zmat = []
        self.is_on_loop = {}

    def add_bond(self, atomname):
        """
        Add a bond to this atom
        """
        self.bonds.append(atomname)
        self.is_on_loop[atomname] = False
        self.traversed[atomname] = False

    def sort_bonds(self, metric):
        """
        Sort the bonds based on some metric
        """
        self.bonds = sorted(self.bonds, key=lambda x: metric[x])[::-1]

    def next_bond(self, defined, update=True):
        """
        Find an atom, bonded to this atom that
        has not been defined and where the bond has not been traversed
        """
        next = None
        for bond in self.bonds:
            if not defined[bond] and not self.traversed[bond] and self.is_on_loop[bond]:
                if update:
                    self.traversed[bond] = True
                next = bond
                break

        if next == None:
            for bond in self.bonds:
                if not defined[bond] and not self.traversed[bond] and not self.is_on_loop[bond]:
                    if update:
                        self.traversed[bond] = True
                    next = bond
                    break
        return next

    def backward_bond(self, exclude, defined):
        """
        Find an atom, bonded to this atom that has
        already been defined but is not in the list of excluded atoms
        """
        back = None
        for bond in self.bonds:
            if bond not in exclude and defined[bond]:
                back = bond
                break
        return back

    def improper_dihedral(self, exclude, defined):
        """
        Look for the definition of an improper dihedral:
        see if at least two bonds to this atom has been defined,
        excluding the bond in the exclude list
        """
        a1 = None
        a2 = None
        for bond in self.bonds:
            if bond in exclude:
                continue
            if defined[bond]:
                if a1 == None:
                    a1 = bond
                elif a2 == None:
                    a2 = bond
                    return (a1, a2)
        return (a1, a2)

class GraphStructure(object) :
    """
    Class to encapsulate a structure that is also a graph
    """
    def __init__(self, verbosity=0) :
        self._initialized = False
        self.names = []
        self.names_sorted = []
        self.nodes = {}
        self.graph = None
        self.structure = None
        self.zmat = []
        self.verbosity = verbosity

    def initialize(self, input) :
        """
        Initialize the graph from either a filename or from a parmed.Structure object
        """

        if not isinstance(input, Structure) :
            input = parmed.load_file(input, structure=True)
        self.structure = input

        # Loop over the atoms in the Structure and create NodeAtom objects
        # and a NetworkX graph
        self.names = []
        self.nodes = {}
        self.graph = nx.Graph()
        for atom in input :
            self.names.append(atom.name)
            node = NodeAtom()
            node.name = atom.name
            for partner in atom.bond_partners :
                if partner.name not in self.names :
                    continue
                if partner.name not in node.bonds :
                    node.add_bond(partner.name)
                    self.nodes[partner.name].add_bond(atom.name)
                    self.graph.add_edge(atom.name, partner.name)
            self.nodes[atom.name] = node

        # Mark up the cycles
        for cycle in nx.cycle_basis(self.graph) :
            for name1, name2 in zip(cycle[:-1], cycle[1:]):
                self.nodes[name1].is_on_loop[name2] = True
                self.nodes[name2].is_on_loop[name1] = True
            self.nodes[cycle[0]].is_on_loop[cycle[-1]] = True
            self.nodes[cycle[-1]].is_on_loop[cycle[0]] = True

        # Sort the bonds for each atom by the closeness
        cc = nx.closeness_centrality(self.graph)
        self.names_sorted = sorted(cc, key=cc.get)[::-1]
        if self.verbosity > 1:
            print("%5s %9s"%("Atom", "Closeness"))
        for name in self.names_sorted:
            self.nodes[name].sort_bonds(cc)
            if self.verbosity > 1:
                print("%5s %9.3f"%(name, cc[name]))
        if self.verbosity > 1:
            print("")

        self._initialized = True

    def _define_atom(self, current, previous, defined, terminal):
        """
        Define a new atom in the z-matrix
        """

        defined[current] = True
        # Trivial definition for the first 3 atoms
        if len(self.zmat) == 0:
            self.nodes[current].zmat = "DM3 DM2 DM1".split()
        elif len(self.zmat) == 1 or (len(self.zmat) == 2 and not terminal):
            self.nodes[current].zmat = self.nodes[previous].zmat[:2]
            self.nodes[current].zmat.insert(0, previous)
        # special logical if third atom is terminal, should only occur
        # in molecules like methane
        elif len(self.zmat) == 2 and terminal:
            a2 = self.nodes[previous].backward_bond([current, previous], defined)
            self.nodes[current].zmat = ("%s %s %s" % (previous, a2, "DM3")).split()
        # for the other atoms, let's see if we should define
        # a proper or an improper dihedral
        else:
            a2, a3 = self.nodes[previous].improper_dihedral(current, defined)
            # Make an improper dihedral if the previous atom returned
            # two previously defined atoms
            if a2 != None and a3 != None:
                self.nodes[current].zmat = ("%s %s %s" % (previous, a2, a3)).split()
            # otherwise define a proper dihedral by traversing backwards
            # on the molecular graph, looking for two atoms to form an angle
            # and dihedral
            else:
                a2 = self.nodes[previous].backward_bond([current, previous], defined)
                a3 = self.nodes[a2].backward_bond([previous, a2], defined)
                self.nodes[current].zmat = ("%s %s %s" % (previous, a2, a3)).split()

        if self.verbosity > 1:
            sys.stdout.write(" "+current)
        self.zmat.append(current + " " + " ".join(self.nodes[current].zmat))

    def traverse(self) :
        """
        Traverse the molecular graph and define the z-matrix
        """

        if not self._initialized :
            raise Exception("Function initialize() has not been called!")

        defined = {node:False for node in self.nodes}
        self.zmat = []

        if self.verbosity > 1:
            print("Traversal of the molecular graph:")

        # Start with the most central atom
        branch_atom = self.names_sorted[0]
        terminal_flag = False
        self._define_atom(branch_atom, None, defined, terminal_flag)

        while True:
            # Traverse a branch from an atom with at least two bonds
            previous = branch_atom
            next = self.nodes[branch_atom].next_bond(defined)
            self.nodes[next].traversed[branch_atom] = True
            # Loop until we reach the end of the molecular graph
            while next != None:
                self._define_atom(next, previous, defined, terminal_flag)
                previous = next
                next = self.nodes[previous].next_bond(defined)
                if next != None:
                    self.nodes[next].traversed[previous] = True
                terminal_flag = next == None
            if self.verbosity > 1:
                sys.stdout.write(".")

            # Check if we have more branches to traverse
            if self.nodes[branch_atom].next_bond(defined, update=False) == None:
                if self.verbosity > 1:
                    sys.stdout.write(" : ")
                # Tries to find a new atom to branch off from
                found = False
                for name in self.names_sorted:
                    if defined[name]:
                        continue
                    branch_atom = name
                    terminal_flag = False
                    self._define_atom(branch_atom, self.nodes[name].bonds[0],
                                    defined, terminal_flag)
                    # If the found atom has more than one bonds, we can traverse
                    # its branches
                    if len(self.nodes[branch_atom].bonds) > 1:
                        found = True
                        break
                    else:
                        if self.verbosity > 1:
                            sys.stdout.write("; ")
                # If we could not find any more atoms to branch off from we are
                # done!
                if not found:
                    break
        if self.verbosity > 0:
            print("")
            print("\nZ-matrix:")
            for z in self.zmat:
                print(z)
            print("")
