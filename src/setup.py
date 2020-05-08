"""
setup.py

Created by: Martin Sicho
On: 5/7/20, 12:42 PM
"""
import os

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

about = {}
with open(os.path.join(os.path.dirname(__file__), "genui/about.py")) as fp:
    exec(fp.read(), about)

version = f"{about['__version__']}"
# check if development status and development version are either both empty (final release) or both set (prerelease)
assert (about['__devstatus__'] and about['__devversion__']) or ((about['__devstatus__'] == about['__devversion__']) and (about['__devstatus__'] in (None, '')))
if about['__devstatus__'] and about['__devversion__']:
    version += f".{about['__devstatus__']}{about['__devversion__']}"

setuptools.setup(
    name="genui",
    version=version,
    author="Martin Sicho",
    author_email="martin.sicho@vscht.cz",
    description="Collection of Python modules for the GenUI backend application.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/martin-sicho/genui",
    packages=setuptools.find_packages(include=['genui', 'genui.*']),
    classifiers = [
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    install_requires=[
        'celery>=4.4.0',
        'chembl_webresource_client>=0.10.1',
        'Django>=2.2.8',
        'django-celery-model>=0.2.1',
        'django-celery-results>=1.1.2',
        'django-polymorphic>=2.1.2',
        'djangorestframework>=3.11.0',
        'drf-yasg>=1.17.0',
        'django-rest-auth>=0.9.5',
        'django-allauth>=0.41.0',
        'django-cors-headers>=3.2.1',
        'redis>=3.3.11',
        'requests>=2.22.0'
        'scikit-learn>=0.22.1'
        'joblib>=0.14.1'
        'pandas>=0.25.3'
        'tqdm>=4.41.1'
        'molvs>=0.1.1'
        'psycopg2-binary>=2.8.4',
        'opentsne>=0.3.12',
        'gunicorn>=20.0.4',
        'drugex @ git+https://github.com/martin-sicho/DrugEx.git@feature/api#egg=drugex',
        'celery-progress @ git+https://github.com/czue/celery-progress@f56e19410a3be458de37e38e74f3a720cac5252f#egg=celery-progress',
        'django-rdkit @ git+https://github.com/rdkit/django-rdkit.git@2565d4bf831bf43c997d81b3aa1a39842b7f56b1#egg=django-rdkit',
    ]
)