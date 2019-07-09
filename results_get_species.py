import pandas as pd
from Bio import Entrez
Entrez.email = 'td1515@ic.ac.uk'
Entrez.api_key = '41267f8592172caaa22ab00ec006c4330208'

key_refseq = pd.read_csv('~/datalocal/bacteria_assembly_summary_refseq_minimal.tsv', sep='\t')
key_genbank = pd.read_csv('~/datalocal/bacteria_assembly_summary_genbank_minimal.tsv', sep='\t')

#open the results thing and each individual sheet 
bbhresult = 'delta-BBH-results.xlsx' # name of results xlsx
xls = pd.ExcelFile(bbhresult) # 
dic = pd.read_excel(xls, header=0, sheet_name=None) # read the file into a dictionary, keys are worksheets and values are dataframes. None means all worksheets are imported 
print('results xslx is opened')
# list is 'PflA-d-BBH-matches', 'PflA-d-BBH-negs', 'PflB-d-BBH-matches', 'PflB-d-BBH-negs'

def get_entrez_id(accession):
	handle = Entrez.esearch(db='assembly', term=accession, retmax='100000')
	record = Entrez.read(handle)
	handle.close()
	ids = record['IdList']
	return ids

#fetch summary of the assembly identified with the id, 
def get_assembly_summary(ids):
	handle = Entrez.esummary(db='assembly', id=ids, report='full')
	record = Entrez.read(handle, validate=False)
	handle.close()
	summary = record['DocumentSummarySet']['DocumentSummary'][0]
	return summary


def find_species_name(xlsx_dict):
	df = pd.DataFrame()
	for key in xlsx_dict.keys(): 
		for index, gcf in xlsx_dict[key]['subject_gcf'].iteritems():
			if 'GCF' in gcf:
				hit = key_refseq.loc[key_refseq['assembly_accession'] == gcf]
				if hit.empty:
					hit = key_genbank.loc[key_genbank['gbrs_paired_asm'] == gcf]
					if hit.empty:
						print(gcf)
			elif 'GCA' in gcf:
				hit = key_genbank.loc[key_genbank['assembly_accession'] == gcf]
				if hit.empty:
					hit = key_refseq.loc[key_refseq['gbrs_paired_asm'] == gcf]
					if hit.empty:
						print(gcf)
			hit = hit.loc[:, ['organism_name', 'infraspecific_name', 'taxid', 'species_taxid']]
			hit.insert(0, 'subject_gcf', gcf)
			hit['sheet'] = key
			df = df.append(hit, ignore_index=True)
	return df

df = find_species_name(dic)
xl = pd.ExcelWriter('gcf_and_species_names.xlsx')
for key in dic.keys():
	abc = df.loc[df['sheet'] == key]
	abc.to_excel(xl, sheet_name=key, index=False)
xl.save()

		

