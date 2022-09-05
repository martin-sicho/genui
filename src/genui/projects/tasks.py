"""
tasks

Created by: Martin Sicho
On: 05.09.22, 11:54
"""
from celery import shared_task
from genui.models import helpers
from genui.projects.models import Project


@shared_task(name="createDefaultModels", bind=True)
def createDefaultModels(self, project_id, app):
    project = Project.objects.get(pk=project_id)
    try:
        created = helpers.createDefaultModels(project, app)
    except ModuleNotFoundError:
        return {
            "project": project.id,
            "app": app,
            "models": []
        }
    models = []
    if created:
        models = [x.name for x in created]
        print(f'Created default models {", ".join(models)} from {app} for project {project.name} (owned by {project.owner.username})')

    return {
        "project": project.id,
        "app": app,
        "models": models
    }