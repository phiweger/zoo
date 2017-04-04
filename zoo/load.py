import Bio.Entrez
import os


def read_accessions(fp):
    with open(fp) as acc_lines:
        return [line.strip() for line in acc_lines]


def accessions_to_gb(accessions, db, batchsize, retmax):
    def batch(sequence, size):
        '''Splitting accession list into batches of <batchsize>.'''
        l = len(accessions)
        for start in range(0, l, size):
            yield sequence[start:min(start + size, l)]

    def extract_records(records_handle):
        '''...'''
        buffer = []
        for line in records_handle:
            if line.startswith("LOCUS") and buffer:
                # yield accession number and record
                yield buffer[0].split()[1], "".join(buffer)
                buffer = [line]
            else:
                buffer.append(line)
        yield buffer[0].split()[1], "".join(buffer)

    def process_batch(accessions_batch):
        '''Main function.'''

        # get GI for query accessions
        query = " ".join(accessions_batch)
        query_handle = Bio.Entrez.esearch(db=db, term=query, retmax=retmax)
        gi_list = Bio.Entrez.read(query_handle)['IdList']

        # get GB files
        search_handle = Bio.Entrez.epost(db=db, id=",".join(gi_list))
        search_results = Bio.Entrez.read(search_handle)
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]

        records_handle = Bio.Entrez.efetch(
            db=db, rettype="gb", retmax=batchsize,
            webenv=webenv, query_key=query_key)
        yield from extract_records(records_handle)
        # yield from, stackoverflow, 9708902

    accession_batches = batch(accessions, batchsize)
    for acc_batch in accession_batches:
        yield from process_batch(acc_batch)


GB_EXT = ".gb"


def write_record(dir, accession, record):
    with open(os.path.join(dir, accession + GB_EXT), "w") as output:
        print(record, file=output)
