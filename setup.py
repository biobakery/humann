from setuptools import setup, find_packages

setup(
    name='humann2lib',
    version='0.0.1',
    description='',
    packages=['humann2lib'],
    zip_safe=False,
    install_requires=[],
    classifiers=[
        "Development Status :: 3 - Alpha"
    ],
    entry_points= {
        'distutils.commands': [
            'download = humann2lib.src.utilities:DownloadDBsCommand'
        ],
        'console_scripts': [
            'humann2.py = humann2lib.humann2:main',
        ],
    }
)
