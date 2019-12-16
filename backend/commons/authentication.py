"""
authentication

Created by: Martin Sicho
On: 16-12-19, 13:18
"""

from rest_framework.authentication import SessionAuthentication, BasicAuthentication

class CsrfExemptSessionAuthentication(SessionAuthentication):

    def enforce_csrf(self, request):
        return  # To not perform the csrf check previously happening

PERMISSIVE_CLASSES = (CsrfExemptSessionAuthentication, BasicAuthentication)


