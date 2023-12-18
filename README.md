# ORD-MiShMASh project

MiShMASh is a Python 3 package to download and evaluate sequence data and metadata associated with publications.
It enables programmatic access to literature on PubMed Central and to data on the Sequence Read Archive (and other INSDC databases).
Users can easily assess the openness of publication data and obtain said data for themselves.

## Installation
Use the package manager pip to install ord-mishmash, including all dependencies.
```bash
pip install git+
```
## Usage
ord-mishmash provides a couple of 2 main commands to scrape a content of PMC publications  and fetching the related metadata from SRA. 

```shell
scrape get_metadata \
              --email  \
              --accession_ids  \
              --output_file \
             
```
```shell
scrape pdf_analysis \
              --pubmed_central_ids  \
              --output_file \
             
```
where:
- ` --email` is the user email address,
- `--accession_ids` stands for list of acccession ids to be  scraped,
- `--pmc_ids` stands for list of pubmed cetral ids of the publications to be scraped,
- `--output_file` stands for the output file name.


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
ord-mishmash is released under a BSD-3-Clause license. See LICENSE for more details.