from datetime import datetime

from background_task import background
from Bio import SeqIO
from Bio.Seq import Seq


from django.core.cache import cache

from .models import Protein, SearchHistory
from .utils import get_user_by_id

@background(schedule=1)
def sequence_matcher(user_id, input_dna):
	dna = Seq(input_dna)
	transcribed_protein = dna.translate()
	user = get_user_by_id(user_id)
	matched_protein = None

	db = cache.get("proteins")

	# protein_match = next((protein for protein in proteins if protein['sequence'].find(transcribed_protein.__str__()) != -1), None)

	for protein in db:
		match_start_index = protein['sequence'].find(transcribed_protein.__str__())
		if match_start_index != -1:
			matched_protein = Protein.objects.filter(id=protein['id'])[0]
			break;

	s = SearchHistory(
		user_id=user,
		dna_sequence=dna.__str__(),
		transcribed_sequence=transcribed_protein.__str__(),
		matched_protein=matched_protein or None,
		search_timestamp=datetime.now(),
		protein_matched_beginning_index=match_start_index)
	s.save()
