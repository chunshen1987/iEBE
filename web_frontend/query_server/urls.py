from django.conf.urls import patterns, url

from query_server import views

urlpatterns = patterns(
    '',
    url(r'$', views.home, name='query'),
    url(r'home/$', views.home, name='query'),
    url(r'query/$', views.query, name='query'),
)
