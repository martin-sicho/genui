"""
serializers

Created by: Martin Sicho
On: 05-12-19, 12:25
"""
from genui import apps

from genui.models import helpers
from .models import Project
from rest_framework import serializers

# Serializers define the API representation.
class ProjectSerializer(serializers.HyperlinkedModelSerializer):
    owner = serializers.HiddenField(
        default=serializers.CurrentUserDefault()
    )

    class Meta:
        model = Project
        fields = ('id', 'name', 'description', 'created', 'updated', 'owner')
        read_only_fields = ('created', 'updated', 'owner')

    def create(self, validated_data):
        ret = super().create(validated_data)

        for app in apps.all_():
            try:
                created = helpers.createDefaultModels(ret, app)
                if created:
                    print(f'Created default models {", ".join([x.name for x in created])} from {app} for project {ret.name} (owned by {ret.owner.username})')
            except ModuleNotFoundError:
                pass

        return ret


class FilterToProjectMixIn:

    def get_queryset(self):
        queryset = super().get_queryset()
        project = self.request.query_params.get('project_id', None)
        if project is not None:
            queryset = queryset.filter(project__pk=int(project))
        return queryset