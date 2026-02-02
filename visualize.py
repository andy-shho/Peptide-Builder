import sys
from peptide import Peptide

peptide = Peptide(sys.argv[1])
peptide.visualize_peptide(save=True, title=sys.argv[2])