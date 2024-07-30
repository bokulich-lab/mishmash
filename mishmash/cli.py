"""
The command-line interface for the Python scraper
"""

import argparse
import nltk
import os

from .fetch_metadata import get_metadata
from .scrape_pdf import analyze_pdf


def install_nltk_punkt_dataset():
    try:
        nltk.data.find("tokenizers/punkt")
    except LookupError:
        nltk.download("punkt")


def main():
    install_nltk_punkt_dataset()

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(required=True)

    md_parser = subparsers.add_parser("get_metadata",
                                      help="Retrieves metadata from INSDC "
                                           "database accession IDs.")
    md_parser.set_defaults(func=get_metadata)
    md_parser.add_argument("--email",
                           help="User email address required for database "
                                "access.",
                           type=str,
                           required=True)
    md_parser.add_argument("--accession_list",
                           n_args="+",
                           help="Space-separated list of INSDC accession IDs "
                                "to retrieve metadata for.",
                           required=True)
    md_parser.add_argument("--n_jobs",
                           help="Number of jobs to run in parallel",
                           default=1,
                           required=False)
    md_parser.add_argument("--output_file",
                           help="File name for output.",
                           type=str,
                           default="output.csv",
                           required=False)
    md_parser.add_argument("--verbose",
                           help="Prints process messages to standard output; "
                                "use for debugging.",
                           action="store_true")

    accession_parser = subparsers.add_parser("get_accessions",
                                             help="From published literature, "
                                                  "retrieves accession IDs "
                                                  "for INSDC datcdabases.")
    accession_parser.set_defaults(func=analyze_pdf)
    accession_parser.add_argument("--pmc_list",
                                  nargs="+",
                                  help="Space-separated list of PubMed Central "
                                       "IDs to evaluate for INSDC accessions.",
                                  required=True)
    accession_parser.add_argument("--output_file",
                                  help="File name for output.",
                                  type=str,
                                  default="output.csv")

    args = parser.parse_args()
    output_df = args.func(args)

    if os.path.exists(args.output_file):
        response = input(
            f"The file '{args.output_file}' already exists. "
            f"Do you want to overwrite it? (y/n): "
        )
        if response.lower() != "y":
            print("Operation aborted by the user.")
            return

    output_df.to_csv(args.output_file)
    print("Results saved to {}".format(args.output_file))


if __name__ == "__main__":
    main()
