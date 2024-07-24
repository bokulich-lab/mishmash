# MISHMASH

MIShMASh (MIcrobiome Sequence and Metadata Availability Standards) is a Python 3 package to download and evaluate sequence data and metadata associated with publications. It enables programmatic access to literature on PubMed Central and to data on the Sequence Read Archive (and other INSDC databases). Users can easily assess the openness of publication data and obtain said data for themselves.

## Installation
Use the package manager `pip` to install MISHMASH and its dependencies.

```bash
pip install git+https://github.com/bokulich-lab/mishmash.git
```

## Usage
MIShMASh provides two main commands to scrape the content of PMC publications and fetch related metadata from the SRA.
 
```shell
scrape get_metadata \
  --email  \
  --accession_ids  \
  --output_file
```
where:
- `email` is your email address (required by NCBI).
- `--accession_ids` is space-separated list of accession IDs to retrieve metadata for.
- `--output_file` specifies the output file to write the retrieved metadata to.

```shell
scrape analyze_pdf \
  --pubmed_central_ids  \
  --output_file
```
where:
- `--pubmed_central_ids` is space-separated list of PubMed Central IDs to perform the pdf analysis on.
- `--output_file` specifies the output file to write the analysis to.

## Contributions
### Pull requests
To set up a development environment, use [Poetry](https://python-poetry.org/).
```console
pip install poetry
poetry install
```
Test the code by running
```console
poetry run pytest
```

## License
MIShMASh is released under a BSD-3-Clause license. See LICENSE for more details.
