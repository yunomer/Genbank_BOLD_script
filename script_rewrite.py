import sys
import os
import time
import argparse
import urllib.error
import codecs
from Bio import SeqIO
from Bio import Entrez

Entrez.email = "bob@hotmail.com"

recognition_list_default = ["CO1", "COI", "COXI", "COX1", "MT-CO1", "MT-COX1", "cox", "coi", "cox1", "cytochrome oxidase subunit I", "cytochrome c oxidase subunit 1", "cytochrome oxidase subunit 1", "cytochrome c oxidase subunit I", "CoxI", "Cox1", "co1", "coI"]
header_list_default = ["accession", "collected_by", "collection_date", "gene", "lat_lon", "organism", "sequence", "taxonomy", "notes", "country"]
feature_type_default = ["CDS", "gene"]

parser = argparse.ArgumentParser()
# -rl Recognition List -hl Header List -l Error Log file name -o Output file Type
parser.add_argument('accessions')
parser.add_argument("-a", "--accessions", help="File Containing list of Accession IDs")
parser.add_argument("-r", "--recognition", help="Custom Recognition List")
parser.add_argument("-f", "--feature", help="Custom Feature List")
parser.add_argument("-hl", "--header", help="Custom Headers List and data to pull out")

args = parser.parse_args()


# Function to try to retrieve an annotation from a SeqRecord object, return empty string if not found
def fetch_annotation(string, records):
    output = ""
    for record in range(0, len(records.features)):
        try:
            feature_dict = records.features[record]
            output = feature_dict.qualifiers[string]
        except KeyError:
            pass
    return output


def parser(feature, rec_list, rec, feature_type):
    for key, val in feature.qualifiers.items():
        gene = "No appropriate gene found!"
        sequence = "No sequence found!"
        for value in val:
            if value in rec_list and (feature.type in feature_type):
                gene = value
                sequence = feature.location.extract(rec).seq
                return gene, sequence
        return gene, sequence


# Function to try to retrieve a feature from a SeqRecord object, return empty string if not found
def fetch_feature(string, record):
    output = ""
    try:
        output = record.features[string]
    except KeyError:
        pass
    return output


def grouper(iterable, n, fillvalue=None):
    for i in range(0, len(iterable), n):
        yield iterable[i:i + n]


# Output = Output_file_name, rec_list = recognition_list, logs = logs
def execute(input_file_name, output_file, rec_list, header_list, feature_list):
    output_file = open(output_file, "w")
    # tsv file for data dump based on header list
    tsv_file = open("Result_tsv.tsv", "w")
    # Log file to record Errors
    log_file = open("log_tsv.tsv", "w")
    for header in header_list:
        tsv_file.write(header + "\t")
    tsv_file.write("\n")

    counter = 0

    long_delay = 0
    with codecs.open(input_file_name, 'r') as infile:
        accs = grouper(infile.read().split("\n"), 300, "")
        for chunk in accs:
            print("Processed " + str(counter) + " records...")
            while "" in chunk:
                chunk.remove('')
            acc_list = ",".join(chunk)
            try:
                fetched_gb = Entrez.efetch(db="nucleotide", id=acc_list, rettype="gbwithparts", retmode="text")
                unique_list_error_nonetype = []
                unique_list_found = []
                complete_list = []
                for index, rec in enumerate(SeqIO.parse(fetched_gb, "gb")):
                    complete_list.append(rec.id)
                    for feature in rec.features:
                        ret = parser(feature, rec_list, rec, feature_list)
                        if ret is None:
                            if rec.id not in unique_list_error_nonetype:
                                unique_list_error_nonetype.append(rec.id)
                            break
                        gene = ret[0]
                        sequence = ret[1]
                        if gene != "No appropriate gene found!":
                            output_file.write("> " + rec.id + "|" + str(gene) + "\n" + str(sequence) + "\n")
                            tsv_file.write(rec.id)
                            for header in header_list:
                                fetch_annotation_ret = fetch_annotation(header, rec)
                                if header is "sequence":
                                    fetch_annotation_ret = str(sequence)
                                elif header is "taxonomy":
                                    try:
                                        fetch_annotation_ret = rec.annotations[header]
                                    except IndexError:
                                        fetch_annotation_ret = ""
                                elif len(fetch_annotation_ret) > 0:
                                    fetch_annotation_ret = fetch_annotation_ret[0]
                                else:
                                    fetch_annotation_ret = ""

                                if fetch_annotation_ret is not None:
                                    tsv_file.write(str(fetch_annotation_ret) + "\t")
                                else:
                                    tsv_file.write("\t")
                            tsv_file.write("\n")
                            if rec.id not in unique_list_found:
                                unique_list_found.append(rec.id)
                            counter = counter + 1
                            break
                        else:
                            pass
                output_file.close()
                tsv_file.close()

                # Printing out the IDs with Nonetype error.
                for accession_id in unique_list_error_nonetype:
                    log_file.write(
                        accession_id + "\t" + "No Data found!" + "\n")

                for check in complete_list:
                    if check not in unique_list_found:
                        if check not in unique_list_error_nonetype:
                            log_file.write(
                                check + "\t" + "No Sequence found" + "\n")

                log_file.close()
            except urllib.error.HTTPError:
                if urllib.error.HTTPError.code == 429:
                    time.sleep(5)
                    pass
                else:
                    a = False
            long_delay = long_delay + 1
            if long_delay == 20:
                time.sleep(20)
                long_delay = 0


def gather_status():
    """
    Function that loads all user arguments from console and processes/assigns them to variables
    setting up the environment for data execution

    To run the program, only the input file name is required.
    """
    # If no arguments passed, exit
    if not len(sys.argv) > 1:
        print("Error: No Arguments passed")
        exit(1)

    # list of accessions that need be processed
    if len(sys.argv) == 1:
        input_file_name = sys.argv[1]
    else:
        input_file_name = args.accessions

    # the custom recognition list
    recognition_file = args.recognition

    # the custom feature list
    feature_file = args.feature

    # the custom header list
    header_file = args.header

    # list used to store the recognition list
    recognition_list = []

    # list used to store the recognition list
    feature_list = []

    # list used to store the header list
    header_list = []

    # If user enters a recognition list, use that, else use default list
    if recognition_file is not None:
        try:
            lines = open(recognition_file).readlines()
            for line in lines:
                line = line.rstrip('\n')
                line = line.strip()
                recognition_list.append(line)
        except FileNotFoundError:
            print("Error: File was not found: " + recognition_file)
            exit(1)
    else:
        recognition_list = recognition_list_default

    # If user enters a feature list, use that, else use default list
    if recognition_file is not None:
        try:
            lines = open(feature_file).readlines()
            for line in lines:
                line = line.rstrip('\n')
                line = line.strip()
                feature_list.append(line)
        except FileNotFoundError:
            print("Error: File was not found: " + feature_file)
            exit(1)
    else:
        feature_list = feature_type_default

    if header_file is not None:
        try:
            lines = open(header_file).readlines()
            for line in lines:
                line = line.strip('\n')
                line = line.strip()
                header_list.append(line)
        except FileNotFoundError:
            print("Error: File was not found: " + header_file)
            exit(1)
    else:
        header_list = header_list_default

    input_file_name_edit = os.path.basename(input_file_name)

    output_file_name = str("Result_" + input_file_name_edit.split(".")[0] + ".fasta")

    # Now run the main chuck of code
    execute(input_file_name, output_file_name, recognition_list, header_list, feature_list)


if __name__ == '__main__':
    gather_status()
