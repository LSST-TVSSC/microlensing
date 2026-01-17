from django.db import models

from zoo.client import project
from collections import Counter


class ZooniverseSurvey(models.Model):
    name = models.CharField(max_length=50)

    created = models.DateTimeField(auto_now_add=True)
    updated = models.DateTimeField(auto_now=True)

    def __str__(self):
        return self.name


class ZooniverseTarget(models.Model):
    survey = models.ForeignKey(ZooniverseSurvey, on_delete=models.CASCADE)
    identifier = models.CharField(max_length=128)

    created = models.DateTimeField(auto_now_add=True)
    updated = models.DateTimeField(auto_now=True)

    def __str__(self):
        return self.identifier


class ZooniverseSubject(models.Model):
    subject_id = models.BigIntegerField(unique=True)
    target = models.ForeignKey(
        ZooniverseTarget, on_delete=models.CASCADE, null=True, blank=True
    )

    sequence = models.CharField(
        max_length=50, help_text="Sector, data release, etc.", null=True, blank=True
    )
    data_url = models.URLField()
    start_time = models.DateTimeField(
        null=True, blank=True, help_text="Earliest time in the light curve"
    )
    end_time = models.DateTimeField(
        null=True, blank=True, help_text="Latest time in the light curve"
    )

    metadata = models.JSONField()

    created = models.DateTimeField(auto_now_add=True)
    updated = models.DateTimeField(auto_now=True)

    def talk_url(self):
        return f"https://www.zooniverse.org/projects/{project.slug}/talk/subjects/{self.subject_id}"

    def annotation_totals(self):
        annotations = self.zooniverseclassification_set.values_list(
            "annotation", flat=True
        )
        annotations = [a[0]["value"] for a in annotations]
        return dict(Counter(annotations))


class ZooniverseClassification(models.Model):
    classification_id = models.BigIntegerField(unique=True)
    subject = models.ForeignKey(ZooniverseSubject, on_delete=models.CASCADE)

    user_id = models.BigIntegerField(null=True, blank=True)
    timestamp = models.DateTimeField()
    annotation = models.JSONField()

    created = models.DateTimeField(auto_now_add=True)
    updated = models.DateTimeField(auto_now=True)


class ZooniverseTargetReduction(models.Model):
    """
    Reduced classifications for targets.
    """

    target = models.ForeignKey(ZooniverseTarget, on_delete=models.CASCADE)
    classifications = models.ManyToManyField(ZooniverseClassification)

    reduced_annotations = models.JSONField()

    created = models.DateTimeField(auto_now_add=True)
    updated = models.DateTimeField(auto_now=True)
