from peptide import Peptide

#Creating Peptide by entering one-letter-sequence
one_letter_sequence = "GIGAVLKVLTTGLPALISWIKRKRQQ"
pep = Peptide(one_letter_sequence)

#Accessing Individual Residues
#Determine the index of the letter in the sequence representing the residue you want to access
#Accessing first A or Alanine
res = pep.get_residue(index=3)

#Accessing attributes of residue
print(res.name)
print(res.one_letter_code)
print(res.three_letter_code)
print(res.bonds)

#Getting atom information
#The index represents the information for each atom
print(res.atom_ids[0])
print(res.symbols[0])
print(res.coordinates[0])

#Accessing specific residue in peptide chain
print(pep.get_residue(index=3).one_letter_code)
print(pep.get_residue(unique_name="ALA_4"))

#Visualization
#Use the visualize_peptide method
#Two options, save or show

#Show Only
pep.visualize_peptide()

#Save
# pep.visualize_peptide(save=True, title="melittin")