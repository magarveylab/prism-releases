#!/usr/bin/python
#Author cDejong (Only slight adjustments from original so posting copywrite):
#  Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met: 
#
#     * Redistributions of source code must retain the above copyright 
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following 
#       disclaimer in the documentation and/or other materials provided 
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
#       nor the names of its contributors may be used to endorse or promote 
#       products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#Import 
import sys

#Import RDkit
from rdkit import rdBase
from rdkit import Chem
from rdkit.Chem import MACCSkeys, AllChem
#from rdkit.Avalon import pyAvalonTools
from rdkit.Chem.AtomPairs import Pairs, Torsions
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.ChemicalFeatures import BuildFeatureFactory
#from rdkit.Chem.Pharm2DFP import Pharm2DFPCalculator
from rdkit.Chem import rdMolDescriptors

# implemented fingerprints:
# ECFC0 (ecfc0), ECFP0 (ecfp0), MACCS (maccs), 
# atom pairs (ap), atom pairs bit vector (apbv), topological torsions (tt)
# hashed atom pairs (hashap), hashed topological torsions (hashtt) --> with 1024 bits
# ECFP4 (ecfp4), ECFP6 (ecfp6), ECFC4 (ecfc4), ECFC6 (ecfc6) --> with 1024 bits
# FCFP4 (fcfp4), FCFP6 (fcfp6), FCFC4 (fcfc4), FCFC6 (fcfc6) --> with 1024 bits
# Avalon (avalon) --> with 1024 bits
# RDKit with path length = 5 (rdk5), with path length = 6 (rdk6), with path length = 7 (rdk7)
# 2D pharmacophore (pharm) ?????????????

nbits = 1024
inSmiles = sys.argv[1]
smiles = Chem.MolFromSmiles(inSmiles)

fpdict = {}
fpdict['ecfp0'] = AllChem.GetMorganFingerprintAsBitVect(smiles, 0, nBits=nbits)
fpdict['ecfp2'] = AllChem.GetMorganFingerprintAsBitVect(smiles, 1, nBits=nbits)
fpdict['ecfp4'] = AllChem.GetMorganFingerprintAsBitVect(smiles, 2, nBits=nbits)
fpdict['ecfp6'] = AllChem.GetMorganFingerprintAsBitVect(smiles, 3, nBits=nbits)
fpdict['fcfp2'] = AllChem.GetMorganFingerprintAsBitVect(smiles, 1, useFeatures=True, nBits=nbits)
fpdict['fcfp4'] = AllChem.GetMorganFingerprintAsBitVect(smiles, 2, useFeatures=True, nBits=nbits)
fpdict['fcfp6'] = AllChem.GetMorganFingerprintAsBitVect(smiles, 3, useFeatures=True, nBits=nbits)
fpdict['maccs'] = MACCSkeys.GenMACCSKeys(smiles)
fpdict['ap'] = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(smiles, nBits=nbits)
fpdict['tt'] = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(smiles, nBits=nbits)
fpdict['rdk5'] = Chem.RDKFingerprint(smiles, maxPath=5, fpSize=nbits, nBitsPerHash=2)
fpdict['rdk6'] = Chem.RDKFingerprint(smiles, maxPath=6, fpSize=nbits, nBitsPerHash=2)
fpdict['rdk7'] = Chem.RDKFingerprint(smiles, maxPath=7, fpSize=nbits, nBitsPerHash=2)
#fpdict['avalon'] = fpAvalon.GetAvalonFP(smiles, nbits)

#Convert to hex for space save, to go back to bin use "bin(int(x, 16))[2:]
for k, v in fpdict.items():
	#hexVal = v.ToBinary().encode("hex")
	print k, v.ToBitString()

