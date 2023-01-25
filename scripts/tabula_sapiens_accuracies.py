import numpy as np
import pandas as pd
from scanpy import read_h5ad

import sys
# Add path of module (temporary measure until pip install works)
sys.path.append('/project/biophys/mouse scRNA-seq/scTOP explorations/scTOP/sctop')
import sctop as top

tabula_sapiens_file = read_h5ad('/projectnb/biophys/mariay/scTOP_manuscript_data/TS/TabulaSapiens.h5ad')

# Ignore immune types and other types that aren't specific to organs
types_to_keep = pd.read_csv('/project/biophys/mouse scRNA-seq/scTOP explorations/scTOP/tutorial/manuscript data/ts_types.csv', header=None)
types_to_keep = types_to_keep.loc[:,0].str.replace(' \'|\'|,', '', regex=True)

cell_id_types = tabula_sapiens_file.obs.loc[:,'cell_ontology_class'].astype('string') + " (" + tabula_sapiens_file.obs.loc[:,'organ_tissue'].astype('string') + ")"
keeping_condition = tabula_sapiens_file.obs.loc[:,'cell_ontology_class'].str.contains("|".join(types_to_keep))
kept_obs = tabula_sapiens_file.obs[keeping_condition]
all_types, all_counts = np.unique(kept_obs.loc[:,'cell_ontology_class'].astype('string') + " (" + kept_obs.loc[:,'organ_tissue'].astype('string') + ")", return_counts=True)

# Drop types with counts below 200 cells
ts_types = all_types[all_counts>200]

type_expressions_list = []
training_indices = []

for ts_type in ts_types:
    current_type_indices = np.where(cell_id_types==ts_type)[0]
    
    # Choose n_training cells to use for the basis
    n_training = 200
    random_indices = np.random.choice(current_type_indices, n_training, replace=False)
    training_indices.append(random_indices)
    
    current_raw = pd.DataFrame(tabula_sapiens_file.layers['raw_counts'][random_indices].T.toarray(),
             index = tabula_sapiens_file.var['gene_symbol']
            )
    current_processed = top.process(current_raw, average=True)
    type_expressions_list.append(current_processed)

training_indices = np.concatenate(training_indices)
ts_basis = pd.concat(type_expressions_list, axis=1)
ts_basis.columns = ts_types

all_organs = set(kept_obs.loc[:,'organ_tissue'].astype('string'))

# Exclude lymph nodes because we're ignoring immune cells
all_organs.remove('Lymph_Node')

test_indices = range(tabula_sapiens_file.layers['raw_counts'].shape[0])
test_indices = np.delete(test_indices, training_indices)

for index, index_array in enumerate(np.array_split(test_indices, 20)):
    test_raw = pd.DataFrame(tabula_sapiens_file.layers['raw_counts'][index_array].T.toarray(),
                       index = tabula_sapiens_file.var['gene_symbol'], columns = index_array
                      )
    test_processed = top.process(test_raw)
    del test_raw
    test_projections = top.score(ts_basis, test_processed)
    del test_processed
    test_projections.to_csv('/projectnb/biophys/mariay/scTOP_manuscript_data/TS/test_projections{}.csv'.format(index))
