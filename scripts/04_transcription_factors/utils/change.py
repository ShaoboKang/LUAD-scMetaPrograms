# =============================================================================
# IMPORTANT: Before running this script
# =============================================================================
# This script may contain hard-coded paths. Please review and update
# paths according to your local environment before execution.
# =============================================================================

import os, sys
os.getcwd()
os.listdir(os.getcwd()) 

import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv("/media/desk16/iyun6206/LUAD_project/Result_4/scenic/input/pyscenic_input_counts.csv");
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("/media/desk16/iyun6206/LUAD_project/Result_4/scenic/output/sample.loom",x.X.transpose(),row_attrs,col_attrs);
