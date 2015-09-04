#!/usr/bin/python
# Author cDejong

# Import 
import sys

# Import RDkit
from rdkit import rdBase
from rdkit import Chem
from rdkit.Chem import MACCSkeys, AllChem
from rdkit.Chem.AtomPairs import Pairs, Torsions
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.ChemicalFeatures import BuildFeatureFactory
from rdkit.Chem import rdMolDescriptors

nbits = 1024
inSmiles = sys.argv[1]
smiles = Chem.MolFromSmiles(inSmiles)

fpdict = {}
fpdict['fcfp6'] = AllChem.GetMorganFingerprintAsBitVect(smiles, 3, useFeatures=True, nBits=nbits)
fpdict['ecfp6'] = AllChem.GetMorganFingerprintAsBitVect(smiles, 3, nBits=nbits)

#Convert to hex for space save, to go back to bin use "bin(int(x, 16))[2:]
for k, v in fpdict.items():
	print k, v.ToBitString()