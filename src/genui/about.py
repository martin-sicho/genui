__all__ = [
    "__version__",
    "get_release_info"
]

# https://setuptools.readthedocs.io/en/latest/setuptools.html#specifying-your-project-s-version
__version__ = "0.0.0"

import os

if os.path.exists(os.path.join(os.path.dirname(__file__), "_version.py")):
    from ._version import version
    __version__ = version


def get_release_info():
    """
    Generate an easy-to-use dictionary with release information.

    Returns
    -------
    dict
        a dictionary with release information
    """
    version_split = __version__.split('+')
    version_main = ".".join(version_split[0].split('.')[0:3])


    devstatus = version_split[0].split(".")[-1] if len(version_split) > 1 else None
    devversion = version_split[1].split(".")[-1] if len(version_split) > 1 else None
    return {
        'version': version_main,
        'version_full': __version__,
        'dev_status': devstatus,
        'dev_version': devversion
    }
