import numpy as np
import pandas as pd

import sys
# Add path of module (temporary measure until pip install works)
sys.path.append('/project/biophys/mouse scRNA-seq/scTOP explorations/scTOP/sctop')
import sctop as top

import h5py

fname = '/projectnb/biophys/mariay/scTOP_manuscript_data/Brain Atlas/Analysis_Zeng_Hippocampus_10X/data/10x_v2/mouse/processed/YaoHippo2020/CTX_Hip_counts_10x.h5'
h5f = h5py.File(fname,'r')
genes = np.array(h5f['data/gene'].asstr())
annotations = pd.read_csv('/projectnb/biophys/mariay/scTOP_manuscript_data/Brain Atlas/Analysis_Zeng_Hippocampus_10X/data/10x_v2/mouse/processed/YaoHippo2020/CTX_Hip_anno_10x/CTX_Hip_anno_10x.csv')
annotations = annotations.set_index('sample_name')

h5_names = np.array(h5f['data/samples']).astype(str)

specifier = 'subclass_label'
class_basis = pd.read_csv('/project/biophys/mouse scRNA-seq/scTOP explorations/scTOP/manuscript/manuscript data/brain_atlas.csv')
training_indices = pd.read_csv('/project/biophys/mouse scRNA-seq/scTOP explorations/scTOP/manuscript/manuscript data/brain_atlas_training_indices.csv')

class_basis = class_basis.set_index('Unnamed: 0')
training_indices = training_indices.loc[:,'0'].values

training_indices_converted = []
for index in training_indices:
    training_indices_converted += [np.where(h5_names == annotations.index[int(index)])[0][0]]
    
# Choose ten populations to make heatmaps
populations = ['Astro', 'Oligo', 'CT SUB', 'CA3', 'Sst', 'L2/3 IT CTX', 'L5 PT CTX', 'Pvalb', 'Vip', 'Lamp5']
population_scores = pd.DataFrame(index=populations, columns=populations)

for population in populations:
    # Random sample of indices for population type, excluding indices used for training data
    current_condition = annotations.loc[:, specifier] == population
    index_options = np.setdiff1d(np.where(current_condition)[0], training_indices)
    random_indices = np.random.choice(index_options, size=min(50000, len(index_options)), replace=False)

    converted_indices = []
    for index in random_indices:
        converted_indices += [np.where(h5_names == annotations.index[index])[0][0]]

    converted_indices = np.array(converted_indices)
    converted_indices.sort()
    
    current_sample = pd.DataFrame(h5f['data/counts'][:,  converted_indices], index = genes)
    processed_sample = top.process(current_sample, average=False)
    current_scores = top.score(class_basis, processed_sample)
    population_scores.loc[population] = current_scores.loc[populations].mean(axis=1)

population_scores.to_csv('/project/biophys/mouse scRNA-seq/scTOP explorations/scTOP/manuscript/manuscript data/brain_atlas_bubbleplot.csv')