# MISHMASH

MIShMASh (MIcrobiome Sequence and Metadata Availability Standards) is a Python 3 package to download and evaluate sequence data and metadata associated with publications. It enables programmatic access to literature on PubMed Central and to data on the Sequence Read Archive (and other INSDC databases). Users can easily assess the openness of publication data and obtain said data for themselves.

## Installation
Use the package manager `pip` to install MISHMASH and its dependencies.

```bash
pip install git+https://github.com/bokulich-lab/mishmash.git
```

## Usage
MIShMASh provides two main commands to evaluate sequence and metadata reporting in publications.

### Evaluate sequencing data reporting
To classify the quality of sequencing data reporting, run `assess_sequences`:

```shell
mishmash assess_sequences \
  --pmc_list  \
  --output_file
```
where:
*`--pmc_list` is space-separated list of PubMed Central IDs to search for associated INSDC accession IDs e.g. the PMC ID **PMC6240460** for [the 2017 article by Naymagon et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6240460/)
* `--output_file` specifies the output file to write the analysis to.

### Evaluate metadata reporting
To retrieve metadata associated with a sequence record from an INSDC (e.g. SRA, DDBJ, ENA) database, run `assess_metadata`:

```shell
mishmash assess_metadata \
  --email  \
  --accession_list  \
  --output_file
```
where:
* `--email` is your email address (required by NCBI).
* `--accession_ids` is a space-separated list of accession IDs to retrieve metadata for. These can be BioProject, BioSample, BioExperiment, or likewise accession IDs that are used within INSDC interfaces e.g. the BioProject ID **PRJNA607574** for the [collection of samples uploaded by the Memorial Sloan Kettering Cancer Center.](https://www.ncbi.nlm.nih.gov/bioproject/?term=(PRJNA607574)%20AND%20bioproject_sra[filter]%20NOT%20bioproject_gap[filter])
* `--output_file` specifies the output file to write the retrieved metadata to.

Optional parameters to `assess_metadata` include:
* `--n_jobs`: an integer value for number of threads in parallelization
* `--verbose`: a flag to print intermediate process outputs to standard output; use in debugging


## Outputs
### `assess_sequences`

This module generates a comma-separated file with the following information:
* PMC ID: Input PubMed Central ID for query
* Sequence Accessibility Badge: Bronze, Silver, or Gold (or "Cannot be determined") as an evaluation of the accessibility of sequencing data from the paper
* INSDC Accessions Numbers: Accession numbers corresponding to the sequencing data uploaded to INSDC databases
* INSDC Database: Database associated with the uploaded sequencing data i.e. SRA, ENA, or DDBJ
* Number of Sequence Records: Total number of sequencing records (INSDC Runs) associated with the input article
* Primer Sequences: If an amplicon-based study, sequences of primers used to amplify variable regions for sequencing; output as a comma-separated string
* Sequencing Method: Probability of sequencing method as either amplicon- or shotgun-based; output as a dictionary
* Includes Code: True/False whether a code repository has been found for the paper
* Code URL: Links to code repositories found in paper; output as a list of strings

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
