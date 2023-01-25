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
    
training_sample_names = h5_names[training_indices_converted]

# Find accuracy scores one batch at a time, so we don't have to load 1.2 million samples at once
batch_size = 50000
sample_index = 0

excluded_types = np.setdiff1d(np.unique(annotations.loc[:, specifier]), class_basis.columns.values)
number_excluded_cells = np.sum([len(np.where(annotations.loc[:, specifier] == cell_type)[0]) for cell_type in excluded_types])
total_cells = len(h5_names) - len(training_sample_names) # total cells in test data (total - training)

included_types_accuracies = {'top1': 0,
                             'top3': 0,
                             'unspecified': 0,
                             'total_cells': total_cells - number_excluded_cells
                            }
excluded_types_accuracies = {'top1': 0,
                             'top3': 0,
                             'unspecified': 0,
                             'total_cells': number_excluded_cells
                            }

predicted_labels = []
predicted_labels_specified = []
true_labels = []

while sample_index < len(h5_names):
    sample_slice = slice(sample_index, sample_index + batch_size) 
    current_sample = pd.DataFrame(h5f['data/counts'][:,  sample_slice], index = genes, columns = h5_names[sample_slice])
    samples_to_drop = np.intersect1d(current_sample.columns, training_sample_names)
    current_sample = current_sample.drop(columns=samples_to_drop)
    
    processed_sample = top.process(current_sample, average=False)
    current_scores = top.score(class_basis, processed_sample)
    del processed_sample
    
    for cell_name, cell_scores in current_scores.items():
#         cell_name = current_sample.columns[cell_label]
        correct_type = annotations.loc[cell_name, specifier]
        true_labels += [correct_type]
       
        sorted_scores = cell_scores.sort_values(ascending=False)
        
        predicted_labels += [sorted_scores.index[0]]

        if cell_scores.max() < 0.1:
            predicted_labels_specified += ['Unspecified']
        else:
            predicted_labels_specified += [sorted_scores.index[0]]
        
        if correct_type in excluded_types:
            if cell_scores.max() < 0.1:
                excluded_types_accuracies['unspecified'] += 1
        else:   
            if sorted_scores.index[0] == correct_type:
                included_types_accuracies['top1'] += 1
            if correct_type in sorted_scores.index[:3]:
                included_types_accuracies['top3'] +=1
            if sorted_scores[0] < 0.1:
                included_types_accuracies['unspecified'] +=1
    
    sample_index += batch_size

# Process the final set of samples

sample_slice = slice(sample_index, len(h5_names)) 
current_sample = pd.DataFrame(h5f['data/counts'][:,  sample_slice], index = genes, columns = h5_names[sample_slice])
samples_to_drop = np.intersect1d(current_sample.columns, training_sample_names)
current_sample = current_sample.drop(columns=samples_to_drop)
    
processed_sample = top.process(current_sample, average=False)
current_scores = top.score(class_basis, processed_sample)
del processed_sample
    
# for cell_label, cell_scores in current_scores.iteritems():
#     cell_name = current_sample.columns[cell_label]
#     correct_type = annotations.loc[cell_name, specifier]
    
#     if correct_type in excluded_types:
#         if cell_scores.max() < 0.1:
#                 excluded_types_accuracies['unspecified'] += 1
#     else:   
#         sorted_scores = cell_scores.sort_values(ascending=False)

#         if sorted_scores.index[0] == correct_type:
#             included_types_accuracies['top1'] += 1
#         if correct_type in sorted_scores.index[:3]:
#             included_types_accuracies['top3'] +=1
#         if sorted_scores[0] < 0.1:
#             included_types_accuracies['unspecified'] +=1

for cell_name, cell_scores in current_scores.items():
#     cell_name = current_sample.columns[cell_label]
    correct_type = annotations.loc[cell_name, specifier]
    true_labels += [correct_type]

    sorted_scores = cell_scores.sort_values(ascending=False)
    
    predicted_labels += [sorted_scores.index[0]]

    if cell_scores.max() < 0.1:
        predicted_labels_specified += ['Unspecified']
    else:
        predicted_labels_specified += [sorted_scores.index[0]]

    if correct_type in excluded_types:
        if cell_scores.max() < 0.1:
            excluded_types_accuracies['unspecified'] += 1
    else:   
        if sorted_scores.index[0] == correct_type:
            included_types_accuracies['top1'] += 1
        if correct_type in sorted_scores.index[:3]:
            included_types_accuracies['top3'] +=1
        if sorted_scores[0] < 0.1:
            included_types_accuracies['unspecified'] +=1

accuracies = pd.DataFrame([included_types_accuracies, excluded_types_accuracies], index = ['included_types', 'excluded_types'])
accuracies.to_csv('/project/biophys/mouse scRNA-seq/scTOP explorations/scTOP/manuscript/manuscript data/brain_atlas_accuracies.csv')

true_and_pred = pd.DataFrame([true_labels, predicted_labels, predicted_labels_specified], index=['true', 'predicted', 'predicted specified'])
true_and_pred.to_csv('/project/biophys/mouse scRNA-seq/scTOP explorations/scTOP/manuscript/manuscript data/brain_atlas_predicted_labels.csv')