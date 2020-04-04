from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein, HasStopCodon
from Bio.Seq import Seq

import re

sequence_converter = lambda x: x.update({'sequence': Seq(x['sequence'], generic_protein)})

def dna_to_protein(sequence):
	return Seq(sequence, generic_dna).translate()

def protein_string_to_sequence(sequence):
	return Seq(sequence, generic_protein)

def has_stop_codon(protein_seq):
	return isinstance(protein_seq.alphabet, HasStopCodon)

def contains_only_dna_nucleotides(sequence, search=re.compile(r"^[CAGTcagt]+$").search):
	return bool(search(sequence))

def get_user_by_id(user_id):
	from django.contrib.auth.models import User
	return User.objects.filter(id=user_id)[0]
