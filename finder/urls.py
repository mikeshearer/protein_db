from django.urls import path

from . import views

app_name = 'finder'

urlpatterns = [
	path('', views.index, name='index'),
	path('protein/', views.ProteinListView.as_view(), name='proteins'),
	path('protein/<int:pk>', views.ProteinDetailView.as_view(), name='protein_details'),
	path('search/history', views.SearchListView.as_view(), name='search_history'),
	path('search/<int:pk>', views.SearchDetailView.as_view(), name='search_details'),
	path('search/', views.search_protein_db, name='submit_search'),
	path('signup/', views.signup, name='signup'),
	path('login/', views.login_request, name='login'),
	path('logout/', views.logout_request, name='logout'),
]