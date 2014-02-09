from django.conf.urls import patterns, url

from query_server import views

urlpatterns = patterns('',
    url(r'^$', views.query, name='query')
)
