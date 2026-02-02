"""
Starting code
"""

import pprint
import json
import os

from gemmi import cif

#Custom exception to handle wrong type of file
class FileExtensionError(Exception):
    pass

# don't forget to unzip your data folder before
# running this.

def read_cif(filepath: str):
    """
    Read a CIF file using gemmi.
    """
    # read CIF with gemmi with error handling
    if os.path.splitext(filepath)[1].lower() != ".cif":
        raise FileExtensionError("File is not a .cif file")

    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Error with file: The file {filepath} was not found; File path may be incorrect")

    doc = cif.read_file(filepath)

    # use as_json in gemmi to convert to json
    # convert from json to Python dict using
    # json.loads
    cif_dict = json.loads(doc.as_json())

    # You now have a dictionary to work with!
    # pprint (from pprint) is "pretty print"
    # will make long dict more readable.
    # pprint.pprint(cif_dict)

    # Consider putting a breakpoint here
    # and examining cif_dict live!
    return cif_dict

if __name__ == "__main__":
    read_cif("data/ALA.cif")