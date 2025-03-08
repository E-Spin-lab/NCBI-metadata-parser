import argparse
import pandas as pd
from Bio import Entrez

## FUNCTIONS
def Parse_accession(accessionfile):
    """
    This function processes a text file where each line contains NCBI accession numbers and returns 
    a corresponding list.
    Params:
        accessionfile - str, txt file that contains NCBI accession numbers on each line
    Returns:
        ids - lst, list containing the accession numbers
    """
    with open(accessionfile) as f:
        ids = f.read().split('\n')
    f.close()
    return ids
def NCBI_Fetch(email, ids):
    """
    Using the provided email, this function Connects to the NCBI nucore database and searches using the provided accession numbers
    Params:
        email - str, user's email that has been registered with NCBI
        ids - lst, list containing the accession numbers
    Returns:
        response - ListElement, NCBI formatted each line is a dictionary. each line corresponds to an id
    """
    Entrez.email = email
    handle = Entrez.efetch('nuccore', id=ids, retmode='xml')
    response = Entrez.read(handle)
    return response
def Parse_Metadata_tags(metadata):
    """
    Parse the space separated string to a list.
    Params
        metadata - str, space separated string of genome modifiers
    Returns:
        Simple_Meta - lst, genome modifiers that require a simple function
        Complex_Meta - lst, genome modifiers that require a more complex function
        Missing_Meta - lst, genome modifiers that were misspelt or do not exist
    """
    metadata_list = metadata.split(' ')
    Simple_Meta = []
    Complex_Meta = []
    Missing_Meta = []
    list(map(lambda x: Simple_Meta.append(str(x).lower()) if x.lower() in LIST_01 else None, metadata_list))
    list(map(lambda x: Complex_Meta.append(str(x).lower()) if x.lower() in LIST_02 else None, metadata_list))
    list(map(lambda x: Missing_Meta.append(str(x).lower()) if x.lower() not in LIST_01 and x.lower() not in LIST_02 else None, metadata_list))
    return Simple_Meta, Complex_Meta, Missing_Meta
def Extract_Collected_Metadata_Complex(entry, metadata_item):
    """
    Extract the specified metadata from the NCBI ListElement
    Params
        metadata_item - str, genome modifiers to be extracted
        entry - dict, single line from NCBI retrieved list element, each line corresponds to a single accession
    Returns:
        qualifier - str or lst, if there is no match 'n.a' is returned
    """
    ## Retrieve GBSeq_feature-table corresponding to genome modifiers, ignoring genes and cds etc
    sources = [feature for feature in entry['GBSeq_feature-table'.lower()] if feature['GBFeature_key'] == 'source']
    ## Retrieve specific qualifier
    for source in sources:
        qualifiers = [qual for qual in source['GBFeature_quals'] if qual['GBQualifier_name'] == metadata_item]
        if qualifiers:
            for qualifier in qualifiers:
                return qualifier['GBQualifier_value']
        else:
            return "n.a"

## ARGUMENTS
DESCRIPTION_START = f'NCBI GenBank Metadata Parser.\n Given a list of accession numbers and metadata tags, this script will search NCBI and retrieve the specified \
                    metadata and save to a tsv file. If the accession does not have any information for a particular metadata tag, "n.a" will be returned. For a full list of metadata tags, please see the READ_ME file.'
DESCRIPTION_END =   f'Example:  python NCBI_Metadata_Parser.py --accession_file "PATH/TO/accession_FILE.txt" \
                    --email "my_name@email.com" \
                    --metadata "strain isolate county collection_date"'
parser = argparse.ArgumentParser(description=DESCRIPTION_START, epilog=DESCRIPTION_END)
parser.add_argument('--accession_file', type=str, help="Path to the .txt file that includes NCBI accession numbers. Each accession must be on its line.", required=True)
parser.add_argument('--email', type=str, help='your email registered with NCBI')
parser.add_argument('--metadata', type=str, help="Space separated list of metadata tags. To collect all metadata use --metadata 'all'", required=True)
parser.add_argument('--o', type=str, help="Name of the results file", default='Metadata_Results.txt', required=False)
args = parser.parse_args()

accession_FILE = args.accession_file
EMAIL = args.email
METADATA = args.metadata
OUTFILE = args.o

## define possible metadata tags
LIST_01 = ['gbseq_locus',
'gbseq_length',
'gbseq_strandedness',
'gbseq_moltype',
'gbseq_topology',
'gbseq_division',
'gbseq_update-date',
'gbseq_create-date',
'gbseq_definition',
'gbseq_primary-accession',
'gbseq_accession-version',
'gbseq_other-seqids',
'gbseq_source',
'gbseq_organism',
'gbseq_taxonomy',
'gbseq_references'
]
LIST_02 = ['altitude',
'authority',
'bio_material',
'biotype',
'biovar',
'breed',
'cell_line',
'cell_type',
'chemovar',
'clone',
'collected_by',
'collection_date',
'country',
'cultivar',
'culture_collection',
'dev_stage',
'ecotype',
'forma',
'forma_specialis',
'fwd_primer_name',
'fwd_primer_seq',
'genotype',
'haplogroup',
'haplotype',
'host',
'identified_by',
'isolate',
'isolation source',
'lab_host',
'lat_lon',
'note',
'pathovar',
'pop_variant',
'rev_primer_name',
'rev_primer_seq',
'segment',
'serogroup',
'serotype',
'serovar',
'sex',
'specimen_voucher',
'strain',
'sub_species',
'subclone',
'substrain',
'subtype',
'tissue_lib',
'tissue_type',
'type',
'variety'
]

## Parse metadata
if METADATA.lower() != "all":
    Metadata_lists = Parse_Metadata_tags(metadata = METADATA)
else:
    Metadata_lists = [LIST_01, LIST_02, []]

## Parse accession ids
ids = Parse_accession(accessionfile = accession_FILE)

## Connect to NCBI and collect metadata for each id
print('~~~~~~~~~~~~Connecting to NCBI~~~~~~~~~~~~')
response = NCBI_Fetch(email = EMAIL, ids = ids)

## Create an empty dictionary to store metadata
meta_dict = {}
meta_dict['Accession'] = []
for metadata_item in Metadata_lists[0]:
    meta_dict[metadata_item] = []
for metadata_item in Metadata_lists[1]:
    meta_dict[metadata_item] = [] 

for entry in response:
    ## Extract the Accession
    print(f"~~~~~~~~~~~~Extracting {entry['GBSeq_primary-accession']}~~~~~~~~~~~~")
    meta_dict['Accession'] = entry['GBSeq_primary-accession']
    ## Convert Keys to lowercase
    entry_fmt = dict((k.lower(), v) for k, v in entry.items())
    ## Extract data
    for metadata_item in Metadata_lists[1]:
        result = Extract_Collected_Metadata_Complex(entry = entry_fmt, metadata_item = metadata_item)
        meta_dict[metadata_item].append(result)
    for metadata_item in Metadata_lists[0]:
        try:
            meta_dict[metadata_item].append(entry_fmt[metadata_item])
        except:
            meta_dict[metadata_item].append("n.a")
            Metadata_lists[2].append(metadata_item)

## Convert to panda dataframe and save as a tab separated file
DF_metadata = pd.DataFrame.from_dict(meta_dict)
DF_metadata.to_csv(OUTFILE, sep="\t", index=False)

## Print unrecognized metadata tags
for tag in Metadata_lists[2]:
    print(f"unrecognized metadata tag: {tag}\n")
print("~~~~~~~~~~~~Finished~~~~~~~~~~")