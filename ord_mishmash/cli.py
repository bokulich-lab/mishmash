"""
The command-line interface for the PyScraper
"""

import argparse
from .scrape_pdf import scrape_pdf



def main():
    parser = argparse.ArgumentParser(
        description="Python package to evaluate biomedical scietific study papers"
    )
    parser.add_argument(
        "dir_path", type=str,
        help="The local adress of the papers to be evaluated."
    )
    
    args = parser.parse_args()
    file_res = scrape_pdf(args.dir_path)
    print("Result successfully obtained!")

if __name__ == "__main__":
    main()