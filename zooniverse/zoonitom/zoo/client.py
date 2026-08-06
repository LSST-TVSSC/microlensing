from django.conf import settings

from panoptes_client import Panoptes, Project, Workflow

if (
    settings.ZOONIVERSE_CLIENT_ID
    and settings.ZOONIVERSE_CLIENT_SECRET
    and not Panoptes.client().logged_in
):
    Panoptes.connect(
        client_id=settings.ZOONIVERSE_CLIENT_ID,
        client_secret=settings.ZOONIVERSE_CLIENT_SECRET,
    )


project = Project(settings.ZOONIVERSE_PROJECT_ID)
workflow = Workflow(settings.ZOONIVERSE_WORKFLOW_ID)
