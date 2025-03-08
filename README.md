# Collect and Parse Metadata from NCBI

This repository contains NCBI_Metadata_Parser, a Python script designed to gather metadata from NCBI using a text file composed of accession numbers. Results are saved as tab separated file. If the accession does not have a value for a particular metadata tag, "n.a" will be returned.

## Requirements

The python script utilizes the following python libraries: argparse, biopython, and panda

## Usage

Example

    python NCBI_Metadata_Parser.py --accession_file "PATH/TO/accession_FILE.txt" --email "my_name@email.com" --metadata "strain isolate county collection_date" --o Metadata_Results.txt

### Options

--accession_file (Required) Path to the txt file that contains the NCBI accession numbers. Must be enclosed in quotations if the file path contains spaces.

--email (Required) Your email registered with NCBI

--metadata (Required) Space separated metadata tags that is enclosed in quotations. For example "strain isolate county collection_date" For a full list see the Metadata tags section. Alternatively, the string "all" may be used and all metadata will be returned.

--o (optional) Results File. Default string is Metadata_Results.txt

### Metadata tags

gbseq_locus,
gbseq_length,
gbseq_strandedness,
gbseq_moltype,
gbseq_topology,
gbseq_division,
gbseq_update-date,
gbseq_create-date,
gbseq_definition,
gbseq_primary-accession,
gbseq_accession-version,
gbseq_other-seqids,
gbseq_source,
gbseq_organism,
gbseq_taxonomy,
gbseq_references

altitude,
authority,
bio_material,
biotype,
biovar,
breed,
cell_line,
cell_type,
chemovar,
clone,
collected_by,
collection_date,
country,
cultivar,
culture_collection,
dev_stage,
ecotype,
forma,
forma_specialis,
fwd_primer_name,
fwd_primer_seq,
genotype,
haplogroup,
haplotype,
host,
identified_by,
isolate,
isolation source,
lab_host,
lat_lon,
note,
pathovar,
pop_variant,
rev_primer_name,
rev_primer_seq,
segment,
serogroup,
serotype,
serovar,
sex,
specimen_voucher,
strain,
sub_species,
subclone,
substrain,
subtype,
tissue_lib,
tissue_type,
type,
variety
