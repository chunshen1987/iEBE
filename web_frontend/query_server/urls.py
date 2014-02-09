from django.conf.urls import patterns, url
import urllib2

from query_server import views

urlpatterns = patterns(
    '',
    url(r'$', views.home, name='query'),
    url(r'home/$', views.home, name='query'),
    url(r'query/$', views.query, name='query'),
    url(r'query(?<params>.*)', views.query, name='query'),
)
