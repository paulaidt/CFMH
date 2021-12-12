# General imports
import argparse as ap
import mdtraj as md 
import os
import shutil
from multiprocessing import Pool
from functools import partial

def parse_args():
    """
    It parses the command-line arguments.
    Parameters
    ----------
    args : list[str]
        List of command-line arguments to parse
    Returns
    -------
    parsed_args : argparse.Namespace
        It contains the command-line arguments that are supplied by the user
    """
    parser = ap.ArgumentParser(description=" Convert XTC trajectories into PDB files.")
    parser.add_argument("input", type=str,
                        help="Input containing the XTC trajectories.")
    parser.add_argument("reference_file", type=str,
                        help="Path to PDB file to get the topology from.")
    parser.add_argument("output", type=str,
                        help="Output containing the PDB trajectories.")  
    parser.add_argument('--remove',
                        dest="remove",
                        action='store_true',
                        help="Remove original XTC file.")
    parser.add_argument("-c","--n_proc", type=int,
                        help='Number of processor.', default = 1)


    parsed_args = parser.parse_args()
    return parsed_args

def generate_pdb_file(file, topology, remove = False):
    """
    It generates a PDB file with the trajectory of a given XTC file and a reference PDB with the topology, if selected
    it removes the original XTC file.

    Parameters
    ----------
    file : str
        Path to the XTC file. 
    topology : str
        Path to a PDB file that contains the topology of the trajectories.
    remove : bool
        True if the original PDB file has to be removed. Default: False
    """
    structure = md.load(file[1], top = topology)

    #structure.save(os.path.join(file[1].replace('.xtc','.pdb')))

    structure.save(args.output)
    if remove: 
        os.remove(os.path.join(file[1]))

def main(args):
    """
    It converts all the XTC files of a folder as PDB files. 

    Parameters
    ----------
    args : argparse.Namespace
        It contains the command-line arguments that are supplied by the user    
    """
    structures = [args.input]
    #[os.path.join(args.path,file) for file in os.listdir(args.path)
                  #if file.endswith('.xtc')]
    print(structures)

    generate_pdb_file_paral = partial(generate_pdb_file, remove = args.remove, topology = args.reference_file)
    #print(generate_pdb_file_paral)
    
    with Pool(args.n_proc) as p:
        list(p.imap(generate_pdb_file_paral, enumerate(structures)))

if __name__ == '__main__':
    args = parse_args()
    main(args)