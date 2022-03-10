"""
setup.py

Created by: Martin Sicho
On: 5/7/20, 12:42 PM
"""
import os

import setuptools

with open("../README.md", "r") as fh:
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
    install_requires=[
        'celery',
        'chembl_webresource_client',
        'Django',
        'django-celery-results',
        'django-polymorphic',
        'djangorestframework',
        'drf-yasg',
        'django-rest-auth',
        'django-allauth',
        'django-cors-headers',
        'redis',
        'requests'
        'scikit-learn'
        'joblib'
        'pandas'
        'psycopg2-binary',
        'opentsne',
        'nvgpu',
        'celery-progress',
        'django-celery-model @ git+https://github.com/martin-sicho/django-celery-model.git@1f606293e959960bdd769e5760140b83f609801a#egg=django-celery-model',
        'django-rdkit @ git+https://github.com/rdkit/django-rdkit.git@v0.3.1#egg=django-rdkit',
        'drugex @ git+https://github.com/martin-sicho/DrugEx.git@feature/api#egg=drugex',
        'chembl_structure_pipeline @ git+https://github.com/chembl/ChEMBL_Structure_Pipeline.git@8d599ad389a458c002be9fd8353a91ebaf370743#egg=chembl_structure_pipeline',
    ]
)