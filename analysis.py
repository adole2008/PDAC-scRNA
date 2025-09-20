import numpy as np
import pandas as pd
import scanpy as sc

meta = pd.read_csv("datasets/metadata.csv.gz")
type(meta)
meta.columns = ["cell_id", "tissue", "patient", "cell-type"]
print([p for p in meta["tissue"].values])