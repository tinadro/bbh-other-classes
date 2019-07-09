from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO, Entrez, pairwise2
import pandas as pd
from os.path import dirname, abspath
import glob, sys, csv
import os.path
Entrez.email = 'td1515@ic.ac.uk'
Entrez.api_key = '41267f8592172caaa22ab00ec006c4330208'

# need 2 command line arguments: the initial query .fa sequence, and the database(s) against which the query is to be BLASTed 
# query_fa and database
# database of the initial query organism needs to be manually changed (query_organism_proteins)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# INITIAL BLAST PARAMETERS:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~

query_fa = sys.argv[1]
query_name = query_fa.rsplit(".", 1)[0] # gets the filename without the extension
query_organism_proteins = "/project/home/td1515/datalocal/Cj-subsp-jejuni-81-176-protein/ncbi-genomes-2019-05-03/Cj81-176"
home = dirname(dirname(abspath(__file__)))
database = home + sys.argv[2]
out_format = " '6 qaccver saccver bitscore evalue qlen slen length qstart qend sstart send qcovs pident' "

print(database)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DEFINE HEADERS OF RESULTS TABLES:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

out_columns = out_format[4:-2] # gets only the output format specifiers in out_format
head = out_columns.split(" ") # makes the above string into a list 
gcf_header = "subject_gcf" 
head.append(gcf_header) # attaches 'subject_gcf' to the list. head is now the complete list of headers, used for the results tables. 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GET GCF ASSEMBLY NUMBER TO TRACK STRAIN AND SPECIES OF SUBJECT SEQUENCE:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def get_gcf(a): # gets the GCF_xxxx number from the filename and path of the database 
    db = a.split("/")[-1]
    if db.count("_") == 1:
        return a.split("_")[0]
    else:
        return "_".join(db.split("_", 2)[:2])

gcf = get_gcf(database)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SAVE ALL CREATED FILES HERE:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def makedir(path): # makes a directory if it doesn't yet exist
    try:
        os.mkdir(path)
    except FileExistsError:
        print(path + " exists")

newdir = query_name + "-BBHs/" 
fastadir = newdir + "fa_seq_files/" 

makedir(newdir) # all the results and .faa files will be stored here. 
makedir(fastadir) # the fasta files of the initial best hits will be stored here
print("\n")

# blastp output files are stored in same directory as the database files

#~~~~~~~~~~~~~~~~~~~~~~~~
# INITIAL BLAST SEARCH:
#~~~~~~~~~~~~~~~~~~~~~~~~

output = database + "." + query_name + "_blast-initial.tsv"

# Run BLASTP in command line, according to above settings \
blast_cline = NcbiblastpCommandline(query=query_fa, db=database, outfmt=out_format, out=output, soft_masking='true')
stdout, stderr = blast_cline() 


#~~~~~~~~~~~~~~~~~~~~~~~~
# FILTER BLAST RESULTS:
#~~~~~~~~~~~~~~~~~~~~~~~~

def blast_result_cleanup(blast_file): # function takes blastp tabular output file, retains only hit(s) with best bitscore and evalue
    df = pd.read_csv(blast_file, header=None, sep='\t') # import the blastp output as a pandas.DataFrame, blast_file doesn't have a header
    df.columns = out_columns.strip().split(' ') # naming columns according to -outfmt

    df.sort_values(by=['bitscore', 'evalue'], ascending=[False, True], inplace=True) # sort by bitscore and then e-value

    a = df.loc[df['bitscore'] == df.iat[0, 2]] # object 'a' retains hits with top bitscore
    b = a.loc[a['evalue'] == a.iat[0, 3]] # object 'b' retains hits with top bitscore and e-value

    top_hit_id = b['saccver'] # save the seqID of the top hit(s) into top_hit_id
    top_hit_id = top_hit_id.drop_duplicates() #remove duplicate hits 

    return(b, top_hit_id)

table, top_hit_id = blast_result_cleanup(output)
table[gcf_header] = gcf # adds the gcf number to the blast output table, after filtering it for top hits 

print("BLAST1 HIT:")
print(table) # print the blast output table containing only the top hit(s) 
print('\n')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GET .FA FILES OF TOP BLAST HITS: 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#download .fa file for every top hit of initial blast search
for hit in top_hit_id: # top_hit_id stores the acc.ver numbers of the blast top hits
    handle = Entrez.efetch(db="protein", id=hit, rettype="fasta", retmode="text") # fetch the full sequence of that acc.ver
    record = SeqIO.read(handle, "fasta") # read it out
    handle.close()
    hit_fa = fastadir + hit + ".hit.faa"
    SeqIO.write(record, hit_fa, "fasta") # save it a a file, into the appropriate subdiresctory


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BLAST THE TOP HIT AGAINST GENOME OF INITIAL QUERY ORGANISM:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#take the FASTA(s) I've just made, and blastp against the campy genome
for hit in top_hit_id: # for every FASTA file that was just created
    hit_fa = fastadir + hit + ".hit.faa" # blast query sequence is the FASTA file created above
    out_file = database + "." + query_name + "_blast2-" + hit + ".tsv" # this is the output file 
    blast2_cline = NcbiblastpCommandline(query=hit_fa, db=query_organism_proteins, outfmt=out_format, out=out_file)
    stdout, stderr = blast2_cline() 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE/ OPEN ALL THE FILES IN WHICH RESULTS WILL BE STORED:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

results_faa = open(newdir + query_name + "-BBH_matches.faa", "a+") # makes text file for all bidirectional best hit sequences 
initial_query = SeqIO.read(query_fa, "fasta") # parse the query fasta to get its acc.ver number 

def mktsv_input_header(filename):
	file_exists = os.path.isfile(newdir + query_name + filename) # check if results.tsv already exists
	results_tsv = open(newdir + query_name + filename, "a+") # open (or create) results.tsv
	if not file_exists: #if it didn't exist earlier, add headers
		writer = csv.writer(results_tsv, delimiter='\t')
		writer.writerow(head)
	return(results_tsv)

results_tsv = mktsv_input_header("-BBH_matches.tsv") # make the results table for BBHs
negatives_tsv = mktsv_input_header("-nonBBHs.tsv") # make the results table for initial blast top hits that aren't BBHs


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SAVE THE TOP HIT IF IT'S A BBH:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for file in glob.glob(database + "." + query_name + "_blast2-*.tsv"): # take the reciprtocal blast output file
    table2, top_reciprocal_hits = blast_result_cleanup(file) #filter for only top hits
    table2['query_gcf'] = gcf
    print("BLAST2 HIT:")
    print(table2) 
    print('\n')
    hit = file.split("-")[-1][:-4]
    line = table.loc[table['saccver'] == hit]
    for reciprocal_hit in top_reciprocal_hits: #for each top hit
        if reciprocal_hit == initial_query.id: #if the ID of the top hit is the same as the ID of initial query, append its sequence to results.faa and write its initial blast 
            with open(fastadir+hit+".hit.faa", "r") as hit_file: #open the file titled *top_hit_id*.hit.faa
                hit_seq = hit_file.read() #read the top hit file 
                results_faa.write(hit_seq) #append this top hit file to the results file. 
                results_faa.write("\n")
            line.to_csv(results_tsv, header=False, index=False, sep='\t')
        else:
            line.to_csv(negatives_tsv, header=False, index=False, sep='\t') # if they're not BBHs, append that hit to the negatives file
            
results_faa.close()
results_tsv.close()
negatives_tsv.close()

