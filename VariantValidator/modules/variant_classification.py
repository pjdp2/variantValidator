# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 17:13:08 2021
@author: naomi
@author: Ali
"""

"""
Protein variant format:
    NP_000079.2:p.(Gly197Cys)
    NP_000079.2:p.(G197C)
"""

# Import modules

# Define Ensembl Reference Dictionary
# Object that at the moment simply creates a dictionary and has capacity to
# add further entries. At the moment this is a simple repository which we
# can use in later editions to provide/populate further variant information.



import requests #this is needed to talk to the API
import re  # needed to split the string with multiple delimiters
import json  # needed to create json object


class Ensembl_reference_dict:
    # initiator construct that creates a dictionary attribute when an instance of
    # the class is made.
    def __init__(self):
        self.term_definitions = {}

 # Takes values for each SO descriptor and places in a dictionary using
 # the SO term as a key and displays subsequent information as a list value,
 # for that key.
    def add_entry(self, term, description, SO_number, display_term, impact):
        self.term_definitions[term] = [
            description, SO_number, display_term, impact]


# Calling an instance of the Ensembl_reference_dict class and populating the
# dictionary.
Ensembl_reference = Ensembl_reference_dict()
Ensembl_reference.add_entry(
    "frameshift_variant",
    "A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three",
    "SO:0001589",
    "Frameshift variant",
    "High")
Ensembl_reference.add_entry(
    "stop_gained",
    "A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript",
    "SO:00015872",
    "Stop gained",
    "High")
Ensembl_reference.add_entry(
    "stop_lost",
    "A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript",
    "SO:0001578",
    "Stop lost",
    "High")
Ensembl_reference.add_entry(
    "start_lost",
    "A codon variant that changes at least one base of the canonical start codon",
    "SO:0002012",
    "Start lost",
    "High")
Ensembl_reference.add_entry(
    "synonymous_variant",
    "A sequence variant where there is no resulting change to the encoded amino acid",
    "SO:0001819",
    "Synonymous variant",
    "Low")
Ensembl_reference.add_entry(
    "missense_variant",
    "A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved",
    "SO:0001583",
    "Missense variant",
    "Medium")

# couple of test print statements to access all and specific entries.
# print(Ensembl_reference.term_definitions)
# print(Ensembl_reference.term_definitions['stop_gained'][2])

# Define data
# Will want to replace the variant_accession with a VV input in the long term
#variant_accession = "NP_000079.2:p.(M1G)"
# print(variant_accession)

def make_request(genome_build,variant):
    base_url = 'https://rest.variantvalidator.org/VariantValidator/variantvalidator/'
    content_type = '%3ET/all?content-type=application%2Fjson'
    request = base_url + genome_build + '/' + variant + content_type
    return request

#Requests user input for the genome build
genome_build = input('Please select genome build: (a) GRCh37 or (b) GRCh38 ')
#Prevents user from proceeding without picking correctly 
while not (genome_build == 'a' or genome_build == 'b'):
    print('Genome build not supported.')
    genome_build = input('Please select genome build: (a) GRCh37 or (b) GRCh38 ')

#The above then instructs on which build to use and the choice is stored in the genome variable
if genome_build == "a":
    genome = 'GRCh37'
elif genome_build == "b":
    genome = 'GRCh38'

#NM_000088.3:c.589G
variant_id = input('Please input a RefSeq variant descriptor: ')

while not ('NM_' in variant_id): 
    print('Variant type not supported.')
    variant_id = input('Please input a RefSeq variant descriptor: ')

#Makes the request
request = make_request(genome, variant_id)
print(f'Requesting data from : {request}')
response = requests.get(request)

#Creates a python dictionary object for the returned JSON
response_dictionary = response.json()

#Finds the original variant description from the search to use as a key
keys = list(response_dictionary.keys())
variant_description = keys[0]

#Extracts the variant_accession variables into a smaller dictionary
variant_accession_alternatives = response_dictionary[variant_description]['hgvs_predicted_protein_consequence']

# Extracts a single variant_accession protein descriptor, key can be swapped for alternatives:
# {'lrg_slr': 'LRG_1p1:p.(G197C)', 'lrg_tlr': 'LRG_1p1:p.(Gly197Cys)', 'slr': 'NP_000079.2:p.(G197C)', 'tlr': 'NP_000079.2:p.(Gly197Cys)'} 
variant_accession = variant_accession_alternatives['slr']
variant_accession = str(variant_accession)
print(variant_accession)


# Check variant is formatted correctly
protein_HVGS = re.compile(
    "[N][P][_][0-9]+[\.][0-9]+[:][p][\.][\(][a-zA-Z|*]+[0-9]+[a-zA-Z|*]+[\)]")
HGVS_check = (protein_HVGS.match(variant_accession))
# print(HGVS_check)
if str(HGVS_check) == "None":
    print("Your variant is incorrectly formatted.")
    raise SystemExit(0)

# Split string to get amino acid information
# Note this code would also work to get just the nucleotide variant
variant_accession_split = re.split('[()]', variant_accession)
# print(variant_accession_split)

# define the protein variant
protein_variant = variant_accession_split[1]
# print(protein_variant)

# Use re to split the variant into numbers and letters
number_letter = re.compile("([a-zA-Z]+|[*])([0-9]+)([a-zA-Z]+|[*])")
protein_variant_split = number_letter.match(protein_variant).groups()
print(protein_variant_split)

# Use logic to determine variant type
# This works for three letter and one leter codes
# Edit: added protein_SO_term variable to loop
if protein_variant_split[0] == protein_variant_split[2]:
    #print("Variant is synonymous")
    protein_SO_term = "synonymous_variant"
elif (protein_variant_split[0] != "Ter" or protein_variant_split[0] != "*") and (protein_variant_split[2] == "Ter" or protein_variant_split[2] == "*"):
    #print("Variant is stop gain")
    protein_SO_term = "stop_gained"
elif (protein_variant_split[0] == "Ter" or protein_variant_split[0] == "*") and (protein_variant_split[2] != "Ter" or protein_variant_split[2] != "*"):
    #print("Variant is stop loss")
    protein_SO_term = "stop_lost"
elif protein_variant_split[1] == "1" and (protein_variant_split[0] == "Met" or protein_variant_split[0] == "M")\
        and (protein_variant_split[2] != "Met" or protein_variant_split[2] != "M"):
    #print("Variant is start lost")
    protein_SO_term = "start_lost"
elif (protein_variant_split[0] != "Ter" or protein_variant_split[0] != "*") and (protein_variant_split[2] != "Ter" or protein_variant_split[2] != "*") \
        and (protein_variant_split[1] != "1" and (protein_variant_split[0] != "Met" or protein_variant_split[0] != "M"))\
        and protein_variant_split[0] != protein_variant_split[2]:
    #print("Variant is missense")
    protein_SO_term = "missense_variant"
else:
    #print("Variant type not recognised")
    raise SystemExit(0)

# Open an empty dictionary to store the response
SO_terms_dict = {}

# Add accession + protein SO term to dictionary
SO_terms_dict['Accession'] = variant_accession
SO_terms_dict['SO term'] = protein_SO_term

# Add Ensembl Information to dictionary - impact is currently excluded
SO_terms_dict['SO Information'] = Ensembl_reference.term_definitions[protein_SO_term][0:3]

# Convert the dictionary to a json
SO_terms_output = json.dumps(
    SO_terms_dict,
    sort_keys=True,
    indent=4,
    separators=(
        ',',
        ': '))
print(SO_terms_output)