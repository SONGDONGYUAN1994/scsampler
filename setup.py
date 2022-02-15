import setuptools
with open("README.md", "r") as fh:
    long_description = fh.read()
setuptools.setup(
    name="scsampler", # Replace with your own username
    version="1.0.1",
    author="Dongyuan Song",
    author_email="dongyuansong@ucla.edu",
    description="A package for fast diversity-preserving subsampling of large-scale single-cell transcriptomic data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/SONGDONGYUAN1994/scsampler",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
