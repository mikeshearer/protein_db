from django.contrib import admin

from .models import Protein, SearchHistory

admin.site.register(Protein)
admin.site.register(SearchHistory)
