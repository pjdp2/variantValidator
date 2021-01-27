# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 17:13:08 2021
@author: naomi
@author: Ali

Information in the Ensembl_reference_dict is from Ensembl, available at the following URL:
https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html
The information is the same as the output for sequence ontology terms from the 
Ensembl VEP API

This module assigns sequence ontology terms based on the variant nomenclature. 
v1 - Currently will assign stop_loss, stop_gain, start_loss, missense and
synonymous variants based on a correctly HVGS-formatted protein variant description.

This module has been tested with a string input, and with data obtained from the 
Variant Validator API.
"""

# Import modules
import requests #this is needed to talk to the API
import re  # needed to split the string with multiple delimiters
import json  # needed to create json object

# Define Ensembl Reference Dictionary
# Object that at the moment simply creates a dictionary and has capacity to
# add further entries. At the moment this is a simple repository which we
# can use in later editions to provide/populate further variant information.
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
# Further entries should be added here / code could be used from the Ensembl
# VEP API
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

"""
The code below obtains an input protein variant by querying the Variant Validator API
This can also be replaced with a prompt for a string:
    
variant_accession = input("Please input a RefSeq protein variant: ")
variant_accession = str(variant_accession)

For future integration, this code (end of code section indicated by another long 
comment) should be updated to handle the input object
"""

# Define request to Variant Validator API
def make_request(genome_build,variant):
    base_url = 'https://rest.variantvalidator.org/VariantValidator/variantvalidator/'
    content_type = '%3ET/all?content-type=application%2Fjson'
    request = base_url + genome_build + '/' + variant + content_type
    return request

# Requests user input for the genome build
genome_build = input('Please select genome build: (a) GRCh37 or (b) GRCh38 ')
# Prevents user from proceeding without picking either 'a' or 'b' 
while not (genome_build == 'a' or genome_build == 'b'):
    print('Genome build not supported.')
    genome_build = input('Please select genome build: (a) GRCh37 or (b) GRCh38 ')

# Set genome_build variable based on user input
if genome_build == "a":
    genome = 'GRCh37'
elif genome_build == "b":
    genome = 'GRCh38'

# Example input : NM_000088.3:c.589G
variant_id = input('Please input a RefSeq transcript variant descriptor: ')
# Checks transcript begins 'NM'
while not ('NM_' in variant_id): 
    print('Variant type not supported.')
    variant_id = input('Please input a RefSeq variant descriptor: ')

# Makes the request
request = make_request(genome, variant_id)
print(f'Requesting data from : {request}')
response = requests.get(request)

# Creates a python dictionary object for the returned JSON
response_dictionary = response.json()

# Finds the original variant description from the search to use as a key
keys = list(response_dictionary.keys())
variant_description = keys[0]

# Extracts the variant_accession variables into a smaller dictionary
variant_accession_alternatives = response_dictionary[variant_description]['hgvs_predicted_protein_consequence']

# Extracts a single variant_accession protein descriptor, key can be swapped for alternatives:
variant_accession = variant_accession_alternatives['slr']
variant_accession = str(variant_accession)
print(variant_accession)

# Check variant is formatted correctly
# Note: This was used when protein variant is input as a string, it should not be
# needed if the protein variant is obtained directly from Variant Validator
protein_HVGS = re.compile(
    "[N][P][_][0-9]+[\.][0-9]+[:][p][\.][\(][a-zA-Z|*]+[0-9]+[a-zA-Z|*]+[\)]")
HGVS_check = (protein_HVGS.match(variant_accession))
# print(HGVS_check)
if str(HGVS_check) == "None":
    print("Your variant is incorrectly formatted.")
    raise SystemExit(0)
    
"""
End of variant input code
"""

# NOTE: This code currently applies to protein variants. The same logic could be 
# expanded to process nucleotide variants in the next iteration
# Currently only a single term is assigned to each protein variants. Multiple
# terms would need to be assinged to transcript variants

# Split string to get amino acid information
variant_accession_split = re.split('[()]', variant_accession)

# define the protein variant
protein_variant = variant_accession_split[1]

# Use re to split the variant into numbers and letters
number_letter = re.compile("([a-zA-Z]+|[*])([0-9]+)([a-zA-Z]+|[*])")
protein_variant_split = number_letter.match(protein_variant).groups()
print(protein_variant_split)

# Use logic to determine variant type
# This works for three letter and one letter codes
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
    print("Variant type not recognised")
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


# <LICENSE>
# Copyright (C) 2021 VariantValidator Contributors
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# </LICENSE>