"""
The command-line interface for the Python scraper
"""

import argparse
from .metadatafetcher import get_metadata
from .scrape_pdf import pdf_analysis
import nltk
import os


def install_nltk_punkt_dataset():
    try:
        nltk.data.find("tokenizers/punkt")
    except LookupError:
        nltk.download("punkt")


def main():
    install_nltk_punkt_dataset()

    parser = argparse.ArgumentParser()
    parser.add_argument("--output_file",
                        help="File name for output.",
                        type=str,
                        default="output.csv")
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
                           nargs="+",
                           help="Space-separated list of INSDC accession IDs "
                                "to retrieve metadata for.",
                           required=True)

    accession_parser = subparsers.add_parser("get_accessions",
                                             help="From published literature, "
                                                  "retrieves accession IDs "
                                                  "for INSDC databases.")
    accession_parser.set_defaults(func=pdf_analysis)
    accession_parser.add_argument("--pmc_list",
                                  nargs="+",
                                  help="Space-separated list of PubMed Central "
                                       "IDs to evaluate for INSDC accessions.",
                                  required=True)

    args = parser.parse_args()
    output_df = args.func(args)

    if os.path.exists(args.output_file):
        response = input(
            f"The file '{args.output}' already exists. "
            f"Do you want to overwrite it? (y/n): "
        )
        if response.lower() != "y":
            print("Operation aborted by the user.")
            return

    output_df.to_csv(args.output_file)
    print("Results saved to {}".format(args.output_file))


if __name__ == "__main__":
    main()
