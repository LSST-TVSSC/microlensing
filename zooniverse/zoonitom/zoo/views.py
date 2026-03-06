from django.shortcuts import render
from django.views.generic.detail import DetailView
from django.views.generic.list import ListView

from zoo.models import ZooniverseSubject, ZooniverseClassification, ZooniverseTarget


class ZooniverseTargetDetailView(DetailView):
    model = ZooniverseTarget


class ZooniverseTargetListView(ListView):
    model = ZooniverseTarget
    paginate_by = 100
