import json
import logging

from csv import DictReader

from dateutil.parser import parse as date_parse

from zoo.client import project, workflow
from zoo.models import (
    ZooniverseClassification,
    ZooniverseSubject,
    ZooniverseTarget,
    ZooniverseSurvey,
)

logger = logging.getLogger(__name__)


def generate_subject_export():
    return project.generate_export("subjects")


def generate_classification_export():
    return workflow.generate_export("classifications")


def get_subject_export():
    return project.get_export("subjects").csv_dictreader()


def get_classification_export():
    return workflow.get_export("classifications").csv_dictreader()


def import_classifications():
    """
    Downloads the latest workflow classifications export and creates new ZooniverseClassification
    objects based on it.
    """
    existing_classifications = ZooniverseClassification.objects.all().values_list(
        "classification_id", flat=True
    )
    existing_subjects = ZooniverseSubject.objects.all().values_list(
        "subject_id", flat=True
    )
    for c in get_classification_export():
        classification_id = int(c["classification_id"])
        subject_id = int(c["subject_ids"])
        user_id = c["user_id"]
        if len(user_id) == 0:
            user_id = None
        else:
            user_id = int(user_id)

        if classification_id in existing_classifications:
            continue

        if subject_id not in existing_subjects:
            logger.warning(
                f"Skipping classification {classification_id} for unknown subject {subject_id}"
            )
            continue

        subject = ZooniverseSubject.objects.get(subject_id=subject_id)

        annotation = json.loads(c["annotation"])
        timestamp = date_parse(c["created_at"])

        ZooniverseClassification.objects.create(
            classification_id=classification_id,
            subject=subject,
            user_id=user_id,
            timestamp=timestamp,
            annotation=annotation,
        )


def import_subjects(
    target_identifier=None,
    survey=None,
    survey_identifier=None,
    sequence=None,
    sequence_identifier=None,
):
    """
    Downloads the latest subjects export and creates new ZooniverseSubject objects.

    Options:
        - target_identifier: The metadata key name which gives the target/object ID.
          Any subjects which don't have this metadata key will be skipped.
        - survey: If this and survey_identifier are both provided, filters subjects
          to just the ones in the specified survey. If survey_identifier is not provided,
          assumes all subjects are in the specified survey.
        - survey_identifier: The metadata key name which gives the survey name.
        - sequence: If this and sequence_identifier are both provided, filters subjects
          to just the ones in the specified sequence. Has no effect if sequence_identifier
          is not provided.
        - sequence_identifier: the metadata key name which gives the sequence name (i.e.
          the data release number, sector name, or other grouping).
    """
    if survey is not None:
        survey = ZooniverseSurvey.objects.get_or_create(name=survey)[0]

    existing_subjects = ZooniverseSubject.objects.all()
    if survey is not None:
        existing_subjects = existing_subjects.filter(target__survey=survey)
    existing_subjects = existing_subjects.values_list("subject_id", flat=True)

    count = 0
    for s in get_subject_export():
        if count > 100:
            break
        subject_id = int(s["subject_id"])

        if subject_id in existing_subjects:
            continue

        locations = json.loads(s["locations"])
        metadata = json.loads(s["metadata"])

        if survey_identifier is not None:
            survey_name = metadata.get(survey_identifier, None)
            if survey_name is None:
                continue
            if survey is None:
                survey = ZooniverseSurvey.objects.get_or_create(name=survey_name)[0]
            else:
                if survey.name != survey_name:
                    continue

        target = None
        if target_identifier is not None:
            target_name = metadata.get(target_identifier, None)
            if target_name is None:
                continue
            target = ZooniverseTarget.objects.get_or_create(
                survey=survey, identifier=target_name
            )[0]

        sequence_name = None
        if sequence_identifier is not None:
            sequence_name = metadata.get(sequence_identifier, None)
            if sequence_name is None:
                continue
            if sequence is not None and sequence != sequence_name:
                continue

        ZooniverseSubject.objects.create(
            subject_id=subject_id,
            metadata=s["metadata"],
            data_url=locations["0"],
            target=target,
            sequence=sequence_name,
        )
        count += 1
