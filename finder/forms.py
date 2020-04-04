from django import forms

from .constants import MAX_DNA_LENGTH

class SearchForm(forms.Form):
	dna_sequence = forms.CharField(label='Sequence', max_length=MAX_DNA_LENGTH)
