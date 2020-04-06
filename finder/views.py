from django.contrib.auth import login, authenticate, logout
from django.contrib.auth.forms import UserCreationForm, AuthenticationForm
from django.contrib import messages
from django.core import serializers
from django.shortcuts import render, redirect
from django.http import HttpResponse, HttpResponseRedirect, Http404
from django.urls import reverse
from django.views import generic

from .forms import SearchForm
from .models import Protein, SearchHistory
from .utils import dna_to_protein, has_stop_codon, contains_only_dna_nucleotides
from .worker import sequence_matcher

def _get_username(request):
	username = None
	if request.user.is_authenticated:
		username = request.user.username

	return username

def _get_user_id(request):
	user_id = None
	if request.user.is_authenticated:
		user_id = request.user.id

	return user_id

def _get_user_by_id(user_id):
	return User.objects.filter(id=user_id)[0]

""" Classes """

class ProteinDetailView(generic.DetailView):
	model = Protein
	template_name = 'finder/protein_details.html'
	context_object_name = 'Protein'


class ProteinListView(generic.ListView):
	model = Protein
	template_name = 'finder/proteins.html'
	context_object_name = 'Proteins'


class SearchDetailView(generic.DetailView):
	model = SearchHistory
	template_name = 'finder/search_details.html'
	context_object_name = 'Search History'


class SearchListView(generic.ListView):
	model = SearchHistory
	template_name = 'finder/search_history.html'
	context_object_name = 'Search History'

	def get_queryset(self):
		user_id = _get_user_id(self.request)
		# return SearchHistory.objects.filter(user_id_id=user_id).order_by('-search_timestamp')
		return SearchHistory.objects.all()

	def get_context_data(self, **kwargs):
		context = super(SearchListView, self).get_context_data(**kwargs)
		context['data'] = self.object_list.values()
		return context



def index(request, **kwargs):
	user_id = _get_user_id(request)
	data = []

	if user_id:
		data =  SearchHistory.objects.filter(user_id_id=user_id).order_by('-search_timestamp')

	return render(request, 'finder/index.html', {'data': data})


def signup(request, **kwargs):
	if request.method == "POST":
		form = UserCreationForm(request.POST)
		if form.is_valid():
			user = form.save()
			login(request, user, backend='django.contrib.auth.backends.ModelBackend')
			return redirect('finder:index')
	else:
		form = UserCreationForm()
	return render(request, 'finder/signup.html', {'form': form})


def login_request(request, **kwargs):
	if request.method == "POST":
		form = AuthenticationForm(request=request, data=request.POST)
		if form.is_valid():
			username = form.cleaned_data.get('username')
			password = form.cleaned_data.get('password')
			user = authenticate(username=username, password=password)
			if user is not None:
				login(request, user)
				messages.info(request, f"You're now logged in")
				return redirect('finder:index')
			else:
				messages.error(request, "Invalid username or password")
		else:
			messages.error(request, "Invalid username or password")

	form = AuthenticationForm()
	return render(request, template_name = "finder/login.html", context={"form": form})


def logout_request(request):
	logout(request)
	messages.info(request, "Logged out successfully")
	return redirect('finder:index')


def search_protein_db(request, **kwargs):
	if request.method == "POST":
		form = SearchForm(request.POST)
		user_id = _get_user_id(request)
		if form.is_valid() and user_id:
			if not contains_only_dna_nucleotides(form.cleaned_data['dna_sequence']):
				messages.error(request, "Invalid nucleotides entered")

			sequence_matcher(user_id, form.cleaned_data['dna_sequence'])
			return redirect('finder:index')

		else:
			return redirect('finder:index')
	return redirect('finder:index')
			