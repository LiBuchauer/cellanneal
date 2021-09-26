import pandas as pd
import sys
from pathlib import Path

from cellanneal.pipelines import cellanneal_pipe

# extract paths etc from input
celltype_data_path = sys.argv[1]
bulk_data_path = sys.argv[2]
disp_min = float(sys.argv[3])
bulk_min = float(sys.argv[4])
bulk_max = float(sys.argv[5])
maxiter = int(sys.argv[6])
output_path = sys.argv[7]

print("\nWelcome to cellanneal!\n")

""" 1) Import bulk and cell type data """
print('1A. Importing bulk data ...')
try:
    bulk_df = pd.read_csv(bulk_data_path, index_col=0, header=0)
    bulk_names = bulk_df.columns.tolist()
    print('{} bulk samples identified: {}\n'.format(len(bulk_names),
                                                    bulk_names))
except:
    print("""Your bulk data file could not be imported.
    Please check the documentation for format requirements
    and look at the example bulk data file.""")

print('1B. Importing celltype reference data ...')
# import single cell based reference
try:
    celltype_df = pd.read_csv(celltype_data_path, index_col=0, header=0)
    celltypes = celltype_df.columns.tolist()
    print('{} cell types identified: {}\n'.format(len(celltypes), celltypes))
except:
    print("""Your celltype data file could not be imported.
    Please check the documentation for format requirements
    and look at the example celltype data file.""")


# run cellanneal!
cellanneal_pipe(
    celltype_data_path=Path(celltype_data_path),
    celltype_df=celltype_df,
    bulk_data_path=Path(bulk_data_path),  # path object!
    bulk_df=bulk_df,
    disp_min=disp_min,
    bulk_min=bulk_min,
    bulk_max=bulk_max,
    maxiter=maxiter,
    output_path=Path(output_path))
