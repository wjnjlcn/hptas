import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="hptas",
    version="1.0.0",
    author="Jianan Wang",
    author_email="wjnjlcn@hit.edu.cn",
    description="A library for haplotype phasing using transcriptome information and RNA-seq data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/wjnjlcn/hptas",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    entry_points={
        'console_scripts': [
            'hptas = hptas.program:main'],
    },
)
