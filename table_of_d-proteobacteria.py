from os.path import dirname, abspath
import pandas as pd
import glob

home = dirname(dirname(abspath(__file__)))
db_folder = home + '/datalocal/d-proteobacteria-protein/' #location of my databases

def get_gcf(filepath):
	db = filepath.split('/')[-1]
	if db.count('_') == 1:
		return a.split('_')[0]
	else:
		return '_'.join(db.split('_', 2)[:2])

gcf_ls = []

for item in glob.glob(db_folder + '*protein.faa'):
	gcf = get_gcf(item)
	gcf_ls.append(gcf)

#open the table from where i'll get the names
key_refseq = pd.read_csv('~/datalocal/bacteria_assembly_summary_refseq_minimal.tsv', sep='\t')
key_genbank = pd.read_csv('~/datalocal/bacteria_assembly_summary_genbank_minimal.tsv', sep='\t')
# headers are ['organism_name', 'infraspecific_name', 'assembly_accession', 'gbrs_paired_asm', 'isolate', 'taxid', 'species_taxid', 'asm_name']

df = pd.DataFrame()
print(len(gcf_ls))
for gcf in gcf_ls:
	if 'GCF' in gcf:
		hit = key_refseq.loc[key_refseq['assembly_accession'] == gcf]
		if len(hit) == 0:
			hit = key_genbank.loc[key_genbank['gbrs_paired_asm'] == gcf]
	elif 'GCA' in gcf:
		hit = key_genbank.loc[key_genbank['assembly_accession'] == gcf]
		if len(hit) == 0:
			hit = key_refseq.loc[key_refseq['gbrs_paired_asm'] == gcf]
	hit = hit.loc[:, ['organism_name', 'infraspecific_name', 'taxid', 'species_taxid']]
	hit.insert(0, 'subject_gcf', gcf)
	df = df.append(hit, ignore_index=True)

print(df.head())
print(len(df))

df_gcf = [i for i in df['subject_gcf']]
for ele in gcf_ls:
	if ele not in df_gcf:
		print(ele)

df.to_csv('d-proteobacteria_for_BBH.tsv', sep='\t', index=False)

