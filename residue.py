"""
Code should go here
"""

from starting import read_cif
import numpy as np
import networkx as nx
import os

def create_residue(code):
    data_path = ""
    for directory, _, files in os.walk("."):
        if f"{code}.cif" in files:
            data_path = os.path.join(directory, f"{code}.cif")
            return Residue.from_cif(data_path)
    raise FileNotFoundError("No such such cif file with for {code}")

#Custom error for invalid bonds
class BondError(Exception):
    pass

#Custom exception for shape error
class ShapeError(Exception):
    pass

class Residue:
    def __init__(self, name, atoms, bonds,
                 one_letter_code, three_letter_code,
                 residue_type):
        self.name = name
        self.one_letter_code = one_letter_code
        self.three_letter_code = three_letter_code
        self.residue_type = residue_type
        self._atom_ids = atoms["atom_ids"]
        self._symbols = atoms["symbols"]
        self._coordinates = np.array(atoms["coordinates"]).reshape(self.n_atoms, 3)
        self._bonds = bonds
        self._index_map = {atom_name : row_index for row_index, atom_name in enumerate(self.atom_ids)}

    @property
    def n_atoms(self):
        return len(self.atom_ids)

    @property
    def atom_ids(self):
        return self._atom_ids

    @property
    def symbols(self):
        return self._symbols

    #Ensures that bonds are to valid atoms and that the weight of the bond is in string format
    @property
    def bonds(self):
        valid = set(self.atom_ids)
        for bond in self._bonds:
            if bond[0] not in valid or bond[1] not in valid or type(bond[2]) != str:
                raise BondError(f"Following bond is invalid: {bond}")
        return self._bonds

    @property
    def molecule_representation(self):
        molecule_representation = nx.Graph()
        nodes = [(self.atom_ids[i], {'label': self.symbols[i], 'pos': (self.coordinates[i][0], self.coordinates[i][1], self.coordinates[i][2])})
                for i in range(len(self.atom_ids))]
        edges = [(bond[0], bond[1], {'weight': bond[2]}) for bond in self.bonds]
        molecule_representation.add_nodes_from(nodes)
        molecule_representation.add_edges_from(edges)
        return molecule_representation

    #Helper function for determining residue type
    @classmethod
    def _residue_type(cls, type_field):
        if "PEPTIDE" in type_field:
            return "amino_acid"
        else:
            return "nucleic_acid"

    #Helper function for getting initial key for CIF dictionary
    @classmethod
    def _get_initial_key(cls, cif_dict):
        return list(cif_dict.keys())[0]

    #Helper function to combine x, y, z coordinate values for each atom
    @classmethod
    def _get_coordinates(cls, x, y, z):
        #Combine all of the coordinates for each atom as 2-D numpy array
        coordinates = np.empty((len(x), 3))
        for i in range(len(x)):
            coordinates[i][0] = x[i]
            coordinates[i][1] = y[i]
            coordinates[i][2] = z[i]
        return coordinates

    #Helper function to combine atoms_1, atoms_2, and bond order for each atom
    @classmethod
    def _get_bonds(cls, atoms_1, atoms_2, weights):
        return [(atoms_1[i], atoms_2[i], weights[i]) for i in range(len(atoms_1))]

    @classmethod
    def from_cif(cls, filepath):
        #Read CIF file with read_cif from starter.py
        cif_dict = read_cif(filepath)
        #Error handling to isolate FileNotFoundError
        if not cif_dict:
            return 0
        #Extract the initial key from the dictionary
        cif_dict = cif_dict[Residue._get_initial_key(cif_dict)]

        name = cif_dict["_chem_comp.name"]

        one_letter_code = cif_dict["_chem_comp.one_letter_code"]

        three_letter_code = cif_dict["_chem_comp.three_letter_code"]

        residue_type = Residue._residue_type(cif_dict["_chem_comp.type"])

        atom_ids = cif_dict["_chem_comp_atom.atom_id"]

        symbols = cif_dict["_chem_comp_atom.type_symbol"]

        coordinates_x = cif_dict["_chem_comp_atom.pdbx_model_cartn_x_ideal"]
        coordinates_y = cif_dict["_chem_comp_atom.pdbx_model_cartn_y_ideal"]
        coordinates_z = cif_dict["_chem_comp_atom.pdbx_model_cartn_z_ideal"]

        coordinates = Residue._get_coordinates(coordinates_x, coordinates_y, coordinates_z)

        atoms = {
            "atom_ids" : atom_ids,
            "symbols" : symbols,
            "coordinates" : coordinates
        }

        atoms_1 = cif_dict["_chem_comp_bond.atom_id_1"]
        atoms_2 = cif_dict["_chem_comp_bond.atom_id_2"]
        weights = cif_dict["_chem_comp_bond.value_order"]

        #Combine the first atom, the second atom, and weight for each bond as list of tuples
        bonds = Residue._get_bonds(atoms_1, atoms_2, weights)

        #Return a Residue instance that incorporates all of the necessary data parsed
        return Residue(
            name=name,
            atoms=atoms,
            bonds=bonds,
            one_letter_code=one_letter_code,
            three_letter_code=three_letter_code,
            residue_type=residue_type
        )

    #Final Project Additions
    @property
    def coordinates(self):
        return self._coordinates

    @coordinates.setter
    def coordinates(self, new_coordinates):
        new_coordinates = np.array(new_coordinates)
        if new_coordinates.shape != self._coordinates.shape:
            raise ShapeError("Shape of ndarray does not match shape of current coordinates")

        self._coordinates = new_coordinates

    def remove_atom(self, name):
        try:
            index = self._atom_ids.index(name)
        except ValueError:
            print(f"{name} not found in residue")
            #Returns if there are no atoms found with the name but continue with peptide building
            return

        self._atom_ids.pop(index)
        self._symbols.pop(index)
        self._coordinates = np.delete(self._coordinates, index, axis=0)
        for bond in self._bonds[:]:
            if name in bond:
                self._bonds.remove(bond)
        self._index_map = {atom_name : row_index for row_index, atom_name in enumerate(self.atom_ids)}

    #Magic Method One to obtain row_index
    def __getitem__(self, atom_name):
        if atom_name in self._index_map:
            return self._index_map[atom_name]
        else:
            raise KeyError(f"{atom_name} not found in residue")

    @property
    def backbone_vector(self):
        vector = self._coordinates[self._index_map["C"]] - self._coordinates[self._index_map["N"]]
        magnitude = np.linalg.norm(vector)
        backbone_vector = vector / magnitude
        return backbone_vector