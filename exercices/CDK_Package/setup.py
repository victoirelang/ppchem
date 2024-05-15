from setuptools import setup, find_packages

setup(
    name="cdkpackage",
    version="0.1.0",
    author="Victoire",
    author_email="victoirelng18@gmail.com",
    description="A small utility to fetch and display SMILES structures as SVG using an external API.",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url="http://github.com/victoirelang/CDK-Packagegit",
    packages=find_packages(),
    install_requires=[
        "requests",
        "IPython",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)