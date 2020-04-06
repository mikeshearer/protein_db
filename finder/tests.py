from django.test import TestCase

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein

from .utils import (
	dna_to_protein,
	protein_string_to_sequence,
	has_stop_codon,
	contains_only_dna_nucleotides,
	get_user_by_id
	)

# Create your tests here.

class UtilsTestCase(TestCase):
	def test_dna_to_protein(self):
		sequence = "CATCATCAT"
		target = "HHH"
		
		self.assertEqual(target, dna_to_protein(sequence).__str__())

	def test_protein_string_to_sequence(self):
		protein_string = "HHH"

		self.assertIsInstance(
			protein_string_to_sequence(protein_string), Seq
		)
	
		self.assertEqual(
			protein_string_to_sequence(protein_string).__str__(), protein_string
		)

	def test_has_stop_codon(self):
		protein_with_stop_codon = Seq("CATUAGCAT", generic_dna).translate()
		protein_without_stop_codon = Seq("CATCATCAT", generic_dna).translate()

		self.assertTrue(has_stop_codon(protein_with_stop_codon))
		self.assertFalse(has_stop_codon(protein_without_stop_codon))


	def test_contains_only_dna_nucleotides(self):
		correct_sequence = "CATCATCAT"
		incorrect_sequence = "MIKE"

		self.assertTrue(contains_only_dna_nucleotides(correct_sequence))
		self.assertFalse(contains_only_dna_nucleotides(incorrect_sequence))

