#!/usr/bin/env python

from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
import pandas as pd
import numpy as np
import sys
import os
from itertools import repeat
from argparse import ArgumentParser, RawTextHelpFormatter
import mdtraj as md

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
    help="Path to result file", 
    type=str
    )
parser.add_argument(
    '-r','--reference',  
    help="Path to reference file", 
    type=str
    )
parser.add_argument(
    '-x','--input_xtc', 
    required=True, 
    help="Path to input file", 
    type=str
    )
args = parser.parse_args()


###### MAIN CODE #######

#path = os.path.dirname(os.path.realpath(__file__))

try:
    df = pd.read_csv(args.input, header=None, index_col=0)
except:
    sys.stdout.write("NO INTERACTIONS FOUND")
    with open(args.output, 'w') as output:
        output.write("NO INTERACTIONS FOUND")
        output.close()
    sys.exit(0)

df.columns = ['A_name', 'A_id', 'B_name', 'B_id', 'distance', 'pose']

df['interaction_id'] = df['A_id'].astype(str)+"_"+df['B_id'].astype(str)

res = df.groupby(['interaction_id']).size().reset_index(name='counts')

poses = df['pose'].unique()
votes = {}
freq_int = res[res['counts']>1]['interaction_id']

# Voting
for p in poses:
    votes[p] = 0
    p_int = list(df[df['pose'] == p]['interaction_id'])
    for ele in freq_int:
        if ele in p_int:
            votes[p] += 1

pred = sorted(votes.items(), key=lambda x:x[1], reverse=True)[0][0]


program_id, _ = pred.split('.')
program, _ = program_id.split('_')


ref = md.load(args.reference)
xtcfiles = os.path.join(args.input_xtc,program_id+'.xtc')
traj = md.load(xtcfiles, top=ref)

# calculate rmsd from the pose to the ref
rmsd = md.rmsd(traj, ref, 0)

sys.stdout.write("Best pose is: ")
sys.stdout.write(program_id)
sys.stdout.write('\n')
sys.stdout.write("With Root Mean Square Distance = ")
sys.stdout.write(str(rmsd[0])+"\n")

    
with open(args.output, 'w') as output:
    output.write("Best pose is: ")
    output.write(program_id)
    output.write('\n')
    output.write("With Root Mean Square Distance = ")
    output.write(str(rmsd[0])+"\n")
    output.close()
