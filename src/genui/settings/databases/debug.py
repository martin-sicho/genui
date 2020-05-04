"""
debug

Created by: Martin Sicho
On: 5/4/20, 10:46 AM
"""
import os

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql',
        'NAME': 'genui' if not 'POSTGRES_DB' in os.environ else os.environ['POSTGRES_DB'],
        'USER': 'genui' if not 'POSTGRES_USER' in os.environ else os.environ['POSTGRES_USER'],
        'PASSWORD': 'genui' if not 'POSTGRES_PASSWORD' in os.environ else os.environ['POSTGRES_PASSWORD'],
        'HOST': 'localhost',
        'PORT': 5432,
    }
}