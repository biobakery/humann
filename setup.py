import os
from setuptools import setup, find_packages

here = os.path.realpath(os.path.dirname(__file__))
humann2_script = os.path.join(here, "humann2.py")

setup(
    name='humann2src',
    version='0.0.1',
    description='',
    packages=['humann2src'],
    zip_safe=False,
    install_requires=[],
    classifiers=[
        "Development Status :: 3 - Alpha"
    ],
    scripts=[humann2_script],
    entry_points= {
        'distutils.commands': [
            'download = humann2src.utilities:DownloadDBsCommand'
        ]
    }
)
