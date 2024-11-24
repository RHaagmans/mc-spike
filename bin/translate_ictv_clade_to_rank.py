#!/usr/bin/env python 
import sys
import pandas as pd

ictv_table = sys.argv[1]
clade_rank_table = sys.argv[2]

ictv_df = pd.read_excel(ictv_table, engine="openpyxl")

ictv_clades = ictv_df.loc[:,"Realm":"Genus"]

clade_to_rank = dict()

for rank in ictv_clades.columns:
    print(rank)
    clades = ictv_clades[rank].dropna().unique()
    for clade in clades:
        clade_to_rank[clade] = rank

with open(clade_rank_table, 'w') as table_file:
    clade_rank_tsv = '\n'.join(["{}\t{}".format(clade, rank) for clade, rank in clade_to_rank.items()])
    table_file.write(clade_rank_tsv)


