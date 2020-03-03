from setuptools import setup, find_packages

setup(
    name                = 'hier-backbone',
    version             = '0.9',
    description         = 'Extracty hierarchical backbones from a biparite network',
    author              = 'WS Jo, YY Ahn, JP Park'
    author_email        = 'wooseong.jo@gmail.com, yongyeol@gmail.com',
    url                 = 'https://github.com/jeakwon/ccpy',
    download_url        = 'https://github.com/jeakwon/ccpy/archive/0.0.tar.gz',
    install_requires    =  ['networkx'],
    packages            = find_packages(exclude = []),
    keywords            = ['backbone', 'complex networks', 'hierarchy'],
    python_requires     = '>=3',
    package_data        = {},
    zip_safe            = False,
    classifiers         = [
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
)
