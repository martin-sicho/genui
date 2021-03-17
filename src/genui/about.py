"""
Version and related info.

Created by: Martin Sicho
On: 5/7/20, 2:25 PM
"""

__all__ = [
    "__version__",
    "__devstatus__",
    "__devversion__",
    "get_release_info"
]

# https://setuptools.readthedocs.io/en/latest/setuptools.html#specifying-your-project-s-version
__version__ = "0.0.0"
__devstatus__ = "alpha" # prerelease tag
__devversion__ = "1" # prerelease version

# make sure that the version specification is correct
assert (__devstatus__ and __devversion__) or not (__devstatus__ or __devversion__)

def get_release_info():
    """
    Generate an easy to use dictionary with release information.

    Returns
    -------
    dict
        a dictionary with release information
    """

    version = __version__
    if __devstatus__ and __devversion__:
        version += f".{__devstatus__}{__devversion__}"

    return {
        'version' : __version__,
        'version_full' : version,
        'dev_status' : __devstatus__,
        'dev_version' : __devversion__
    }
