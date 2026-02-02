from residue import Residue, create_residue
from mappings import three_to_one, one_to_three
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class NoInputError(Exception):
    pass

#Dictionary for CPK coloring of molecules
cpk_colors = {
    "H" : "white",
    "C" : "black",
    "N" : "blue",
    "O" : "red",
    "F" : "green",
    "Cl" : "green",
    "Br" : "dark red",
    "I" : "dark violet",
    "P" : "orange",
    "S" : "yellow",
    "B" : "beige",
    "He" : "cyan",
    "Ne" : "cyan",
    "Ar" : "cyan",
    "Kr" : "cyan",
    "Xe" : "cyan",
    "Rn" : "cyan",
    "Li" : "violet",
    "Na" : "violet",
    "K" : "violet",
    "Rb" : "violet",
    "Cs" : "violet",
    "Fr" : "violet",
    "Be" : "dark green",
    "Mg" : "dark green",
    "Ca" : "dark green",
    "Sr" : "dark green",
    "Ba" : "dark green",
    "Ra" : "dark green",
    "Ti" : "grey",
    "Fe" : "dark orange",
}

class Peptide:
    def __init__(self, one_letter_sequence):
        self._one_letter_sequence = one_letter_sequence
        self._three_letter_sequence = self._get_three_letter_sequence()
        self._unique_residue_index = {}
        self._bonds = []
        self._peptide_chain = self._get_peptide_chain()

    #Magic Method 2 to return number of amino acids in peptide chain
    def __len__(self):
        return len(self._peptide_chain)

    @property
    def one_letter_sequence(self):
        return self._one_letter_sequence
    @property
    def three_letter_sequence(self):
        return self._three_letter_sequence

    @property
    def unique_residue_index(self):
        return self._unique_residue_index

    @property
    def bonds(self):
        return self._bonds

    @property
    def peptide_chain(self):
        return self._peptide_chain

    def _get_translation_vector(self, prev_res, curr_res):
        target_N = prev_res.coordinates[prev_res["C"]] + (prev_res.backbone_vector * 1.33)
        translation_vector = target_N - curr_res.coordinates[curr_res["N"]]
        return translation_vector

    def _get_three_letter_sequence(self):
        three_letter_sequence = [one_to_three[letter][0] + one_to_three[letter][1:].lower() for letter in self._one_letter_sequence]
        return "-".join(three_letter_sequence)

    def _get_peptide_chain(self):
        peptide_chain = []
        unique_number = 1
        # translating_vector = np.zeros(3)
        # bond_length_vector = np.array([1.33, 0, 0])
        for index, amino_acid in enumerate(self._three_letter_sequence.split("-")):
            unique_name = f"{amino_acid.upper()}_{unique_number}"
            res = create_residue(amino_acid.upper())
            #Special case for first amino acid in chain
            if len(peptide_chain) == 0:
                res.remove_atom("OXT")
                res.remove_atom("HXT")

            #Special case for last amino acid in chain
            elif len(peptide_chain) == len(self._one_letter_sequence) - 1:
                res.remove_atom("H2")
                res.remove_atom("H3")
                prev_res = peptide_chain[-1]
                self._bonds.append([(prev_res, "C"), (res, "N"), "SING"])
                res.coordinates += self._get_translation_vector(prev_res, res)
            else:
                prev_res = peptide_chain[-1]
                res.remove_atom("OXT")
                res.remove_atom("HXT")
                res.remove_atom("H2")
                res.remove_atom("H3")
                self._bonds.append([(prev_res, "C"), (res, "N"), "SING"])
                res.coordinates += self._get_translation_vector(prev_res, res)

            peptide_chain.append(res)
            self._unique_residue_index[unique_name] = (res, index)
            unique_number += 1

        return peptide_chain

    def get_residue(self, index = None, unique_name = None):
        if index is not None:
            return self._peptide_chain[index]
        elif unique_name is not None:
            if unique_name not in self._unique_residue_index:
                raise KeyError(f"{unique_name} not found")
            return self._unique_residue_index[unique_name][0]
        else:
            raise NoInputError("Required to either input index or unique name of targeted residue")

    def visualize_peptide(self, save = False, title = None):
        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(projection="3d")

        all_coordinates = [res.coordinates for res in self._peptide_chain]
        combined_coordinates = np.concatenate(all_coordinates)

        colors = []

        for res in self._peptide_chain:
            for symbol in res.symbols:
                colors.append(cpk_colors[symbol])

        ax.scatter(combined_coordinates[:, 0], combined_coordinates[:, 1], combined_coordinates[:, 2], marker="o", edgecolors='black', facecolors=colors, s=150)

        #Bonds between Residue
        for bond in self._bonds:
            res_1 = bond[0][0]
            res_2 = bond[1][0]
            atom_1 = bond[0][1]
            atom_2 = bond[1][1]

            coordinates = np.array([res_1.coordinates[res_1[atom_1]], res_2.coordinates[res_2[atom_2]]])

            if bond[2] == "SING":
                width = 2
            elif bond[2] == "DOUB":
                width = 6

            ax.plot3D(coordinates[:, 0], coordinates[:, 1], coordinates[:, 2], linewidth=width, color="black")

        #Bonds in Residue
        for res in self._peptide_chain:
            for bond in res.bonds:
                coordinates = np.array([res.coordinates[res[bond[0]]], res.coordinates[res[bond[1]]]])

                if bond[2] == "SING":
                    width = 2
                elif bond[2] == "DOUB":
                    width = 4

                ax.plot3D(coordinates[:, 0], coordinates[:, 1], coordinates[:, 2], linewidth=width, color="black")

        ax.set_xlabel('X (Angstroms)')
        ax.set_ylabel('Y (Angstroms)')
        ax.set_zlabel('Z (Angstroms)')
        ax.set_title(f'{self._one_letter_sequence} Peptide Chain')
        ax.view_init(elev=30, azim=45)
        if not save:
            plt.show()
        else:
            plt.savefig(f"{title}.png")
