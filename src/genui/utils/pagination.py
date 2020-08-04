"""
pagination

Created by: Martin Sicho
On: 5/6/20, 10:55 AM
"""
from urllib import parse as urlparse

from django.conf import settings
from rest_framework import pagination
from rest_framework.utils.urls import replace_query_param, remove_query_param

class GenuiPagination(pagination.PageNumberPagination):
    page_size = 10

    def __init__(self):
        self.location = settings.GENUI_SETTINGS['HOST_URL'] if 'HOST_URL' in settings.GENUI_SETTINGS else None

    def get_url(self):
        uri = self.request.build_absolute_uri(location=self.location)
        if not self.location:
            query = '&'.join([f'{x}={self.request.query_params[x]}' for x in self.request.query_params])
            return uri + f'?{query}'
        path = self.request.get_full_path()
        return urlparse.urljoin(uri, path)

    def get_next_link(self):
        if not self.page.has_next():
            return None
        page_number = self.page.next_page_number()
        return replace_query_param(self.get_url(), self.page_query_param, page_number)

    def get_previous_link(self):
        if not self.page.has_previous():
            return None
        page_number = self.page.previous_page_number()
        if page_number == 1:
            return remove_query_param(self.get_url(), self.page_query_param)
        return replace_query_param(self.get_url(), self.page_query_param, page_number)

