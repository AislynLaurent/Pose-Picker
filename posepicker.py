## POSE PICKER
# Developed by: Aislyn Laurent
# Last updated: 22 / 10 / 2020

## IMPORT
# Standard
import os
import csv as csv
from datetime import datetime
# MDTraj
try:
    import mdtraj as md
except ImportError:
    sys.exit('MDTraj (http://mdtraj.org/1.9.3/) module not found. Terminating...')
# Numpy
try:
    import numpy as np 
except ImportError:
    sys.exit('Numpy (https://numpy.org/) module not found. Terminating...')
# Matplotlib
try:
    from matplotlib import pyplot as plt
except ImportError:
    sys.exit('Matplotlib (https://matplotlib.org/) module not found. Terminating...')
# SciKitLearn
try:
    import sklearn as skl
except ImportError:
    sys.exit('One or more SciKitLearn (https://scikit-learn.org/stable/) module(s) not found. Terminating...')
# SKL Extra
from sklearn import cluster
from sklearn import decomposition
from sklearn import metrics
# Scipy
try:
    from scipy.cluster.hierarchy import dendrogram
except ImportError:
    sys.exit('One or more Scipy (https://www.scipy.org/) module(s) not found. Terminating...')

## GET TRAJECTORY
    # Input: File paths for atom types, trajectory files, and the first frame pdb
    # Output: Array containing processed trajectory and atom type array
def get_trajectory(atom_types_input, trajectory_input_folder, first_frame_pdb):
    input_data = []

    # Relevant Atoms
    with open(atom_types_input) as f:
        atom_types = f.read().splitlines()

    ## GET FILES
    dcd_files = []
    for file in os.listdir(trajectory_input_folder):
        if file.endswith('.coor.dcd'):
            dcd_files.append(os.path.join(trajectory_input_folder, file))

    ## SETUP
    # Load trajectories
    trajectory = md.load(dcd_files, top=first_frame_pdb)
    topology = trajectory.topology

    input_data.append(atom_types)
    input_data.append(trajectory)
    input_data.append(topology)

    return(input_data)

## GET ATOM PAIRS
    # Input: MDTraj Trajectory & Topology, name of the substrate, and atom type array
    # Output: Array of pair displacements
    # Notes: Cleans the data for ML. Attempts to narrow down as much as possible
def get_pair_displacements(trajectory, topology, substrate_name, atom_types):
    ## IDENTIFY ATOMS
    # Substrate filter
    print('Finding \'interesting\' atoms...')
    print('\t[1]: Not in the substrate')
    print('\t[2]: Not part of the backbone')
    print('\t[3]: Not in water (HOH)')
    print('\t[4]: Are on the list of atom types')
    substrate = [atom.index for atom in topology.atoms if (atom.residue.name == substrate_name)]
    
    # Atoms of interest filter - atom type in list and not part of the backbone
    interesting_atoms = [atom.index for atom in topology.atoms if (atom.name in atom_types and atom.residue.name != substrate_name and not atom.is_backbone and atom.residue.name != 'HOH')]
    print('\t\t\t\t\tSuccess!')
    print('\t# of \'interesting\' atoms:', end =" ")
    print(np.shape(interesting_atoms))

    # Neighbors - selects everyone nearby frame by frame
    print('Finding substrate neighbours...', end =" ")
    substrate_neighbours = md.compute_neighbors(trajectory, 0.5, substrate, interesting_atoms)
    print('\tSuccess!')

    # Hydrogen atom types
    print('Finding hydrogens...', end =" ")
    h_atom_types = ['h1','h2','h3','h4','h5','ha','hc','hn','ho','hp','hs','HW','hx','H','h1','h2','H3','h4','h5','ha','HA','hc','HK','hn','ho','hp','HP','hs','HT','hw','HW','hx','HZ','HH3']
    upper_h_atom_types = [atom_type.upper() for atom_type in h_atom_types]
    print('\t\t\tSuccess!')

    ## CALCULATE
    # Unique neighbors - makes a 1D array of all the unique neighbour atoms from all frames
    print('Eliminating duplicate atoms...', end =" ")
    unique_neighbor_atoms = []
    unique_pairs = []

    for frame in substrate_neighbours:
        for atom in frame:
            if atom not in unique_neighbor_atoms:
                unique_neighbor_atoms.append(atom)

    unique_neighbor_atoms.sort()
    print('\t\tSuccess!')
    print('\t# of unique atoms:', end =" ")
    print(np.shape(unique_neighbor_atoms))

    # Distance between everyone
    print('Combining atoms into pairs...', end =" ")
    pairs = topology.select_pairs(substrate, unique_neighbor_atoms)
    print('\t\tSuccess!')
    print('\t# of atom pairs:', end =" ")
    print(np.shape(pairs))

    print('Eliminating duplicate pairs...', end =" ")
    for pair in pairs:
        atom1 = topology.atom(pair[0])
        atom2 = topology.atom(pair[1])

        # combine atom names - if both names appear in h_bond_types, it's an H-H pair and can be excluded
        atom_name_tuple = [atom1.name.upper(), atom2.name.upper()]
        h_h_pair = all(atom in upper_h_atom_types for atom in atom_name_tuple)

        # filter pairs in the same reside and H-H
        if atom1.residue != atom2.residue and not h_h_pair:
            unique_pairs.append(pair)
    print('\t\tSuccess!')

    # Who's touching?
    print('Selecting pairs with nearby atoms...', end =" ")
    relevant_pairs = []

    for pair in unique_pairs:
        relevant_pair = filter_pair(trajectory, pair, relevant_pairs)

        if relevant_pair is not False:
            relevant_pairs.append(relevant_pair)

    print('\tSuccess!')
    print('\t# of \'relevant\' atom pairs:', end =" ")
    print(np.shape(relevant_pairs))
    
    # Calculate displacement (vector)
    print('Calculating displacement vectors...', end =" ")
    relevant_displacements = md.compute_displacements(trajectory, relevant_pairs)
    print('\tSuccess!')
    print('\tDisplacement array:', end =" ")
    print(np.shape(relevant_displacements))

    return(relevant_displacements)

## FILTER PAIR
    # Input: MDTraj Trajectory, specific pair, array of all pairs
    # Output: False if the pair is a duplicate or atoms are further then 0.3nm apart, the pair if otherwise
    # Notes: Used to avoid a sad mess of loops - return if it fails the check
def filter_pair(trajectory, pair, relevant_pairs):
    ## Filter atoms that are too far apart

    # If the pair already exists in the relevant_pairs array, return nothing
    if relevant_pairs:
        for relevant_pair in relevant_pairs:
            if pair[0] == relevant_pair[0] and pair[1] == relevant_pair[1]:
                return(False)

    # MDTraj requires a 2D array of pairs to find distances
    pair_2d = []
    pair_2d.append(pair)
    # Compute those distances
    distances = md.compute_distances(trajectory, pair_2d)

    # The distances are in a 2D array, so check each distance value accordingly
    for distance in distances:
        for value in distance:
            if value <= 0.3:
                return(pair)

    return(False)

## COMPRESS DISPLACEMENTS
    # Input: Pair displacement vector array
    # Output: Compressed (3D -> 1D) displacement array
    # Notes: Using the SciKitLearn dimensionality reduction "PCA"
def compress_displacements(pair_displacement):
    compressed_vectors = []

    # For each frame, compress 3D array a single number
    print('Compressing 3D vectors...', end =" ")
    for frame in pair_displacement:
        compressed_frame = decomposition.PCA(n_components=1).fit_transform(frame)
        flattened_compressed_frame = compressed_frame.flatten()
        compressed_vectors.append(flattened_compressed_frame)
    print('\t\tSuccess!')
    print('\tCompressed array:', end =" ")
    print(np.shape(compressed_vectors))

    return(compressed_vectors)

## Model
    # Input: Compressed displacements 
    # Output: 
    # Notes: 
def compile_model(compressed_vectors):
    selected_frames = []
    j_group_small = False

    # Train the model
    print('Training the model...', end =" ")
    model = cluster.AgglomerativeClustering().fit(compressed_vectors)
    print('\t\tSuccess!')

    ## Label
    print('Labeling frames...', end =" ")
    labels = model.labels_

    i = 0
    j = 0
    for label in labels:
        if label == 0:
            i = i+1
        else:
            j = j+1

    print('\t\tSuccess!')
    print('\tGroup labeled 0: '+str(i))
    print('\tGroup labeled 1: '+str(j))

    # Select and print the smaller group of frames
    print('Selecting frames...', end =" ")
    if i > j:
        j_group_small = True

    if j_group_small:
        for index, label in enumerate(labels):
            if label == 1:
                selected_frames.append(index)
    else:
        for index, label in enumerate(labels):
            if label == 0:
                selected_frames.append(index)
    print('\t\tSuccess!')

    print('\nSelected frames (numbering starts at 0):')
    print(selected_frames)
    print('\n')

    ## Assess
    print('Assessing the quality of the model...')
    print('\tSilhouette score (1 is best, 0 is random): ', end =" ")
    print(metrics.silhouette_score(compressed_vectors, labels, metric='euclidean'))

    return(model)

### MAIN
## HOUSEKEEPING
print('\n------------------------------------------------\n Pose Picker Script \t v.0.1\n\n\tWritten By:\t Aislyn Laurent\n\tLast Edited:\t 22-10-2020\n-------------------------------------------------\n')

## USER INPUT
# Data Files
print('Please enter the following information:\n')
trajectory_input_folder = input('Path to the folder containing trajectory files (PDB, XTC, TRR, DCD, binpos, NetCDF or MDTraj HDF5 formats):\n')
atom_types_input = input('Path to a file containing atom types of interest:\n')
first_frame_pdb = input('Path to a PDB file for the first frame:\n')
substrate_name = input('Substrate name:\n')

if trajectory_input_folder and atom_types_input and first_frame_pdb != None:
    print('\t\t\t\t\tSuccess!')
else:
    print('Some error has occured - terminating...')
    sys.exit('Input error - critical data missing.')

print('\n------------------------------------------------\n')

## CLEAN DATA
# Trajectory and atom types
input_data = get_trajectory(atom_types_input, trajectory_input_folder, first_frame_pdb)
atom_types = input_data[0]
trajectory = input_data[1]
topology = input_data[2]

# Get distances
pair_displacement = get_pair_displacements(trajectory, topology, substrate_name, atom_types)

# Compress vectors
compressed_vectors = compress_displacements(pair_displacement)

print('\n------------------------------------------------\n')

## DO ML
sorted_labels = compile_model(compressed_vectors)

print('\n\t\tC O M P L E T E')

print('\n------------------------------------------------\n')