from django.db import models
from django.conf import settings

from . import constants


class Protein(models.Model):
	name = models.CharField(max_length=254, default=None)
	source = models.CharField(max_length=254, default=None)
	gene = models.CharField(max_length=254, default=None)
	locus_tag = models.CharField(max_length=254, default=None)
	codon_start = models.CharField(max_length=254, default=None)
	product = models.CharField(max_length=254, default=None)
	protein_id = models.CharField(max_length=254, default=None)
	db_xref = models.CharField(max_length=254, default=None)
	sequence = models.CharField(max_length=constants.MAX_PROTEIN_LENGTH, default=None)
	description = models.CharField(max_length=254, default=-1)

	def __str__(self):
		return self.name

class SearchHistory(models.Model):
	user_id = models.ForeignKey(settings.AUTH_USER_MODEL, on_delete=models.CASCADE)
	dna_sequence = models.CharField(max_length=constants.MAX_DNA_LENGTH, null=True)
	transcribed_sequence = models.CharField(max_length=constants.MAX_PROTEIN_LENGTH, null=True, default=None)
	matched_protein = models.ForeignKey(Protein, on_delete=models.PROTECT, null=True, default=None)
	search_timestamp = models.DateTimeField()
	protein_matched_beginning_index = models.IntegerField(default=None)

