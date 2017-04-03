#! /usr/bin/env python3

'''
https://www.biostars.org/p/66921/

Usage:

$ python download_ncbi.py -i data/virome/rna_virome_shi2016.txt \
-d nucleotide -o rna_virome_shi2016
'''


import argparse
import Bio.Entrez
from progressbar import ProgressBar, UnknownLength
import os
import sys


RETMAX = 10**9
GB_EXT = ".gb"


def parse_args(arg_lst):
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True,
                        help="A file with accessions to download")
    parser.add_argument("-d", "--database", type=str, required=True,
                        help="NCBI database ID")
    parser.add_argument("-e", "--email", type=str, required=False,
                        default="some_email@somedomain.com",
                        help="An e-mail address")
    parser.add_argument("-b", "--batch", type=int, required=False, default=100,
                        help="The number of accessions to process per request")
    parser.add_argument("-o", "--output_dir", type=str, required=True,
                        help="The directory to write downloaded files to")

    return parser.parse_args(arg_lst)


def read_accessions(fp):
    with open(fp) as acc_lines:
        return [line.strip() for line in acc_lines]


def accessions_to_gb(accessions, db, batchsize, retmax):
    def batch(sequence, size):
        l = len(accessions)
        for start in range(0, l, size):
            yield sequence[start:min(start + size, l)]

    def extract_records(records_handle):
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


def write_record(dir, accession, record):
    with open(os.path.join(dir, accession + GB_EXT), "w") as output:
        print(record, file=output)


def main(argv):
    args = parse_args(argv)
    accessions = read_accessions(os.path.abspath(args.input))
    op_dir = os.path.abspath(args.output_dir)
    if not os.path.exists(op_dir):
        os.makedirs(op_dir)
    dbase = args.database
    Bio.Entrez.email = args.email
    batchsize = args.batch

    print('Batch size:', batchsize)
    bar = ProgressBar(max_value=UnknownLength)
    counter = 0
    for acc, record in accessions_to_gb(accessions, dbase, batchsize, RETMAX):
        counter += 1
        write_record(op_dir, acc, record)
        bar.update(counter)


if __name__ == "__main__":
    main(sys.argv[1:])
