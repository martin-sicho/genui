"""
prod

Created by: Martin Sicho
On: 5/4/20, 10:46 AM
"""
import os

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql',
        'NAME': os.environ['POSTGRES_DB'],
        'USER': os.environ['POSTGRES_USER'],
        'PASSWORD': os.environ['POSTGRES_PASSWORD'],
        'HOST': os.environ['POSTGRES_HOST'] if 'POSTGRES_HOST' in os.environ else 'db',
        'PORT': 5432,
    }
}
