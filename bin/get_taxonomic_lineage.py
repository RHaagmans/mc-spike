#!/usr/bin/env python 

import sys
import pandas as pd
from collections import defaultdict

table_filename = sys.argv[1]
table_fs = "\t"
taxa_col = 11 # Column containing taxonomix lineage
taxa_sep = ";"# Character separating taxa

ictv_clade_ranks_filename = sys.argv[2]
table_output_filename = sys.argv[3]

# Load ICTV file mapping clade to taxonomic rank
ictv_clade_ranks = {}
with open(ictv_clade_ranks_filename, 'r') as ictv_clade_ranks_file:
    for line in ictv_clade_ranks_file:
        clade, rank = line.strip().split('\t')
        ictv_clade_ranks[clade] = rank

# GeNomad adds herpesviridae, but this is absent from the ICTV VMR
ictv_clade_ranks["Viruses"] = "Superkingdom"
ictv_clade_ranks["Herpesviridae"] = "Family"
ictv_ranks = set(ictv_clade_ranks.values())

table_data = defaultdict(dict)

with open(table_filename, 'r') as genomad_file:
    header = [head.strip() for head in next(genomad_file).split(table_fs)]
    taxa_col_name = header[taxa_col]
    for i, line in enumerate(genomad_file):
        elements = line.split(table_fs)
        table_data[i] = {key: item.strip() for key, item in zip(header, elements)}

for contig, contig_data in table_data.items():
    if contig_data[taxa_col_name] == "Unclassified":
        clades = ["Unclassified" for rank in ictv_ranks]
        ranks = list(ictv_ranks)
    else:
        clades = [clade.strip() for clade in contig_data[taxa_col_name].rstrip('.').split(";") if len(clade) > 0]
        ranks = [ictv_clade_ranks[clade] for clade in clades]
    for rank in ranks:
        table_data[contig].update(zip(ranks, clades))

table_df = pd.DataFrame.from_dict(table_data, orient="index")
table_df.to_csv(table_output_filename, sep=table_fs, index=False)
