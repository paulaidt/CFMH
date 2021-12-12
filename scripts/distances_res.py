#!/usr/bin/env python

from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
import pandas as pd
import numpy as np
import sys
import os
from itertools import repeat
from argparse import ArgumentParser, RawTextHelpFormatter

### ARGUMENTS ###

parser = ArgumentParser(
    description="..."
    )
parser.add_argument(
    '-v','--version', 
    action='version', 
    version='%(prog)s 1.0',
    help="Show program's version number and exit"
    )
parser.add_argument(
    '-i','--input', 
    required=True, 
    help="Path to input file", 
    type=str
    )
parser.add_argument(
    '-o','--output', 
    required=True, 
    help="Path to output file", 
    type=str
    )
parser.add_argument(
    '-d','--distance',
    help="Cutoff disatance for interactions.  Default: 6A",
    type=int, 
    default="6"
    )
parser.add_argument(
    '-s','--cutsasa',
    help="Cutoff for SASA.   Default: 10",
    type=int, 
    default="10"
    )

args = parser.parse_args()


###### MAIN CODE #######
sys.stdout.write("Parsing PDB file...\n")
p = PDBParser(QUIET=1)

sys.stdout.write("Reading structure PDB ...\n")
struct = p.get_structure("ref", args.input)

sys.stdout.write("Calculating SASA for residues ...\n")
sr = ShrakeRupley()
sr.compute(struct, level="R")

model = struct[0]
chainA = model['A']
chainB = model['B']

residA =  list(chainA.get_residues())
residB =  list(chainB.get_residues())

threshold = args.cutsasa

residA_name = []
residB_name = []
distance = []
residA_id = []
residB_id = []

sys.stdout.write("Calculating interaction distances for residues in chain A and B ...\n")

for x in range(len(residA)):
    if residA[x].sasa >= threshold:
        for y in range(len(residB)):
            if residB[y].sasa >= threshold:
                one = residA[x]["CA"].get_coord()
                two = residB[y]["CA"].get_coord()
                d = np.linalg.norm(one-two)
                
                if d < args.distance:
                    # afegir a les llistes
                    residA_name.append(residA[x].get_resname()) # nom
                    residB_name.append(residB[y].get_resname())
                    residA_id.append(residA[x].get_id()[1])
                    residB_id.append(residB[y].get_id()[1])
                    distance.append(d)


tuples = list(zip(residA_name, residA_id, residB_name, residB_id, distance))

df = pd.DataFrame(tuples, columns = ['residA_Name', 'residA_id', 'residB_name', 'residB_id', 'distance'])


f = os.path.basename(args.input)
lst = []
lst.extend(repeat(f, df.shape[0]))
df['file'] = lst



sys.stdout.write("Writing output file...\n")

print(df)

df.to_csv(args.output, header=False)

sys.stdout.write("DONE!\n")
    