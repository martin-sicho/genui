class FilterToUserMixIn:
    owner_relation = 'owner'

    def get_queryset(self):
        if self.request.user and not self.request.user.is_anonymous:
            return super().get_queryset().filter(**{
                self.owner_relation: self.request.user
            })
        else:
            return super().get_queryset().none()