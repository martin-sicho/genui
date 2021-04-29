"""
Django settings for genui project.

Generated by 'django-admin startproject' using Django 2.2.7.

For more information on this file, see
https://docs.djangoproject.com/en/2.2/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/2.2/ref/settings/

See https://docs.djangoproject.com/en/2.2/howto/deployment/checklist/
"""

from .genuibase import *

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = 'euws5ei%zq!@0yyo6ta4^e3whylufayu)26th6869x=ljr44=d' if not 'GENUI_BACKEND_SECRET' in os.environ else os.environ['GENUI_BACKEND_SECRET']

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

# Application definition

SITE_ID = 1

INSTALLED_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'polymorphic',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django.contrib.sites',
    'corsheaders',
    'rest_framework',
    'rest_framework.authtoken',
    'allauth',
    'allauth.account',
    'allauth.socialaccount',
    'rest_auth',
    'rest_auth.registration',
    'drf_yasg',
    'django_celery_results',
    'djcelery_model',
    'celery_progress',
    'django_rdkit',
    'rest_framework_extensions',
] + GENUI_SETTINGS['APPS']

MIDDLEWARE = [
    'corsheaders.middleware.CorsMiddleware',
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    # 'django.middleware.csrf.CsrfViewMiddleware', FIXME: this should be handled once we approach more serious production level
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
]

ROOT_URLCONF = 'genui.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [os.path.join(BASE_DIR, 'genui', 'templates')],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

WSGI_APPLICATION = 'genui.wsgi.application'


# Database
# https://docs.djangoproject.com/en/2.2/ref/settings/#databases

# inheriting config should define this
DATABASES = None


# Password validation
# https://docs.djangoproject.com/en/2.2/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {
        'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',
    },
    {
        'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',
    },
]

# Accounts
ACCOUNT_EMAIL_VERIFICATION = 'none' # TODO: if we decide to expose this to the public somehow, we should make this 'mandatory'
ACCOUNT_EMAIL_REQUIRED = True
ACCOUNT_CONFIRM_EMAIL_ON_GET = True
# LOGIN_REDIRECT_URL = '/' # FIXME: should point to a reasonable location

# Internationalization
# https://docs.djangoproject.com/en/2.2/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'UTC'

USE_I18N = False

USE_L10N = False

USE_TZ = True

# Media files
MEDIA_URL = f'{GENUI_SETTINGS["HOST_URL"]}/downloads/'
MEDIA_ROOT = os.path.join(GENUI_SETTINGS['FILES_DIR'], 'media/')


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/2.2/howto/static-files/

STATIC_URL = f'{GENUI_SETTINGS["HOST_URL"]}/static/'
STATIC_ROOT = os.path.join(GENUI_SETTINGS['FILES_DIR'], 'static/')

# rest framework
REST_FRAMEWORK = {
    # will be able to login using the normal Django Framework login views / templates
    'DEFAULT_AUTHENTICATION_CLASSES': (
        'genui.accounts.authentication.CsrfExemptSessionAuthentication',
        'rest_framework.authentication.TokenAuthentication',
    ),

    # Use Django's standard `django.contrib.auth` permissions,
    # or allow read-only access for unauthenticated users.
    'DEFAULT_PERMISSION_CLASSES': [
        'rest_framework.permissions.DjangoModelPermissions'
    ]
}

# swagger docs
SWAGGER_SETTINGS = {
    'SECURITY_DEFINITIONS': {
        'api_key': {
            'type': 'apiKey',
            'in': 'header',
            'name': 'Authorization'
        }
    },
    'LOGIN_URL' : GENUI_SETTINGS['RF_LOGIN_URL'],
    'LOGOUT_URL' : GENUI_SETTINGS['RF_LOGOUT_URL'],
}

# CORS and CSRF
ALLOWED_HOSTS = []
CSRF_TRUSTED_ORIGINS = []
if GENUI_SETTINGS['HOST']:
    ALLOWED_HOSTS.append(GENUI_SETTINGS['HOST'])
    CSRF_TRUSTED_ORIGINS.append(GENUI_SETTINGS['HOST'])
if DOCKER:
    ALLOWED_HOSTS.append(f"{os.environ['GENUI_CONTAINER_PREFIX']}backend")
    CSRF_TRUSTED_ORIGINS.append(f"{os.environ['GENUI_CONTAINER_PREFIX']}backend")

# celery settings
CELERY_BROKER_URL = 'redis://localhost:6379'
if DOCKER:
    CELERY_BROKER_URL = f"redis://{os.environ['GENUI_CONTAINER_PREFIX']}redis:6379"
if 'REDIS_HOST' in os.environ:
    REDIS_PASS = f':{os.environ["REDIS_PASSWORD"]}@' if 'REDIS_PASSWORD' in os.environ else ''
    CELERY_BROKER_URL = f'redis://{REDIS_PASS}{os.environ["REDIS_HOST"]}:6379'
# CELERY_RESULT_BACKEND = 'redis://redis:6379'
CELERY_RESULT_BACKEND = 'django-db'
CELERY_CACHE_BACKEND = 'django-cache'