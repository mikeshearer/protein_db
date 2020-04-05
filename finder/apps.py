import subprocess
import os

from django.apps import AppConfig
from django.core.cache import cache

from .utils import sequence_converter

class FinderConfig(AppConfig):
    name = 'finder'

    def ready(self):
    	from .models import Protein
    	proteins = Protein.objects.all().values("id", "sequence")
    	_ = list(map(sequence_converter, proteins))

    	cache.set("proteins", proteins, None)
