"""Metadata-only: list Seurat meta.data columns and cell counts of an object (no expression read). python -I pal_meta_cols.py file.rds"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rds_parse as R, numpy as np
o = R.read_rds(sys.argv[1])
print(type(o), [k for k in o.keys()][:20])
md = o['meta.data']
print('meta.data attrs', list(md.keys()))
names = md['names'].data if 'names' in md else None
print('columns', names)
