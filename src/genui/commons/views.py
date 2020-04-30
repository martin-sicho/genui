"""
views

Created by: Martin Sicho
On: 1/1/20, 4:32 PM
"""


class FilterToProjectMixIn:

    def get_queryset(self):
        queryset = super().get_queryset()
        project = self.request.query_params.get('project_id', None)
        if project is not None:
            queryset = queryset.filter(project__pk=int(project))
        return queryset

class FilterToUserMixIn:
    owner_relation = 'owner'

    def get_queryset(self):
        if self.request.user and not self.request.user.is_anonymous:
            return super().get_queryset().filter(**{
                self.owner_relation: self.request.user
            })
        else:
            return super().get_queryset().none()