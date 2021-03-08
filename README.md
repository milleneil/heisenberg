# heisenberg


Dependencies

python 3.6+
numpy


--- Prep SRA data
1. Download series_matrix.gz files
2. Extract methylation values
	extract_sra_probe_vals.py > temp.tsv
3. Create wide file with columns as probes
	invert_sra_table.py temp.tsv > wide.tsv


