# https://www.biostars.org/p/66921/
# https://github.com/nextstrain/fauna/blob/master/vdb/parse.py

# for now, epost only accepts GI value, which have been phased out in
# september 2016
# any younger records cannot be accessed in this way
# crap
# https://ncbiinsights.ncbi.nlm.nih.gov/2016/07/15/ncbi-is-phasing-out-sequence-gis-heres-what-you-need-to-know/


# strategy for now
# tsv + fasta, manually batch download, data wrangle
# see rethinkdb scripts


from Bio import Entrez

fp = "accessions.txt"
l = []
with open(fp, "r") as accessions:
    for accession in accessions:
        l.append(accession.strip("\n"))


# define batch size for download
batchSize = 100
db = "nucleotide"  # "nuccore"
# these are the same thing:
# https://www.ncbi.nlm.nih.gov/nucleotide/X86657
# https://www.ncbi.nlm.nih.gov/nuccore/X86657
Entrez.email = "adrian.viehweger@gmail.com"

# post NCBI query
try:
    search_handle = Entrez.epost(db=db, id='AJ458265')
    search_results = Entrez.read(search_handle)
    webenv, query_key = search_results["WebEnv"], search_results["QueryKey"]
except:
    print("ERROR: Couldn't connect with entrez, please run again")


viruses = []
sequences = []
# fetch all results in batch of batchSize entries at once
for start in range(0, len(l), batchSize):
    # fetch entries in batch
    try:
        handle = Entrez.efetch(
            db=db, rettype="gb", retstart=start, retmax=batchSize,
            webenv=webenv, query_key=query_key)
    except IOError:
        print("ERROR: Couldn't connect with entrez, please run again")
    else:
        print("all good")
#return (viruses, sequences)


"""

def parse_accession_file(self, acc, **kwargs):
    '''
    Parse file for list of accession numbers to be uploaded to vdb
    :return: list of documents(dictionaries of attributes) to upload
    '''
    try:
        handle = open(acc, 'r')
    except IOError:
        raise Exception(acc, "not found")
    else:
        accessions = []
        for acc in handle:
            accessions.append(acc.strip())
    return accessions






    def get_GIs(self, accessions, **kwargs):
        '''
        Use entrez esearch to get genbank identifiers from accession numbers
        '''
        retmax = 10**9
        query = " ".join(accessions)
        handle = Entrez.esearch(db=self.gbdb, term=query, retmax=retmax)
        giList = Entrez.read(handle)['IdList']
return giList




    def get_entrez_viruses(self, giList, **kwargs):
        '''
        Use entrez efetch to get genbank entries from genbank identifiers
        '''
        ## define batch size for download
        batchSize = 100

        # post NCBI query
        try:
            search_handle = Entrez.epost(db=self.gbdb, id=",".join(giList))
            search_results = Entrez.read(search_handle)
            webenv, query_key = search_results["WebEnv"], search_results["QueryKey"]
        except:
            print("ERROR: Couldn't connect with entrez, please run again")

        viruses = []
        sequences = []
        #fetch all results in batch of batchSize entries at once
        for start in range(0,len(giList),batchSize):
            #fetch entries in batch
            try:
                handle = Entrez.efetch(db=self.gbdb, rettype="gb", retstart=start, retmax=batchSize, webenv=webenv, query_key=query_key)
            except IOError:
                print("ERROR: Couldn't connect with entrez, please run again")
            else:
                result = self.parse_gb_entries(handle, **kwargs)
                viruses.extend(result[0])
                sequences.extend(result[1])
return (viruses, sequences)






def parse_gb_entries(self, handle, **kwargs):
    '''
    Go through genbank records to get relevant virus information
    '''
    viruses, sequences = [], []
    SeqIO_records = SeqIO.parse(handle, "genbank")
    for record in SeqIO_records:
        v = {}
        s = {}
        s['source'] = 'genbank'
        s['accession'] = re.match(r'^([^.]*)', record.id).group(0).upper()  # get everything before the '.'?
        s['sequence'] = str(record.seq).lower()
        reference = record.annotations["references"][0]
        if reference.title is not None and reference.title != "Direct Submission":
            s['title'] = reference.title
        else:
            print("Couldn't find reference title for " + s['accession'])
            s['title'] = None
        if reference.authors is not None:
            first_author = re.match(r'^([^,]*)', reference.authors).group(0)
            s['authors'] = first_author + " et al"
        else:
            print("Couldn't parse authors for " + s['accession'])
            s['authors'] = None
            first_author = None
        url = "https://www.ncbi.nlm.nih.gov/nuccore/" + s['accession']
        s['url'] = self.get_doi_url(url, s['title'], first_author)

        record_features = record.features
        for feat in record_features:
            if feat.type == 'source':
                qualifiers = feat.qualifiers
                v['collection_date'] = self.convert_gb_date(qualifiers['collection_date'][0])
                v['country'] = re.match(r'^([^:]*)', qualifiers['country'][0]).group(0)
                if 'strain' in qualifiers:
                    v['strain'] = qualifiers['strain'][0]
                    s['strain'] = qualifiers['strain'][0]
                elif 'isolate' in qualifiers:
                    v['strain'] = qualifiers['isolate'][0]
                    s['strain'] = qualifiers['isolate'][0]
                else:
                    print("Couldn't parse strain name for " + s['accession'])
        v = self.add_virus_fields(v, **kwargs)
        s = self.add_sequence_fields(s, **kwargs)
        viruses.append(v)
        sequences.append(s)
        #self.fix_casing(v)
        #self.fix_boolean(v)
        handle.close()
        print("There were " + str(len(viruses)) + " viruses in the parsed file")
        return (viruses, sequences)



"""