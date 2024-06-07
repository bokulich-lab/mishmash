"""
The command-line interface for the PyScraper
"""

import argparse
from .metadatafetcher import get_metadata
from .scrape_pdf import pdf_analysis
import nltk
import os 

def installNltkPunktDataSet():
    try:
        nltk.data.find('tokenizers/punkt')
    except LookupError:
        nltk.download('punkt')


def main():
    installNltkPunktDataSet()

    parser = argparse.ArgumentParser(
        description="Python package to evaluate biomedical scietific study papers"
    )
    parser.add_argument('command', choices=['get_metadata','pdf_analysis'])
    parser.add_argument("-e","--email",help = 'User email to use for metadata retrieval.')
    parser.add_argument("-i",'--accession_ids', nargs='+',dest = 'ids', 
                        help = 'Space-separated list of accession IDs to retrieve metadata for.' )
    parser.add_argument('-pmc_ids','--pubmed_central_ids', nargs='+', dest = 'pmc_ids', 
                        help = 'Space-separated list of PubMed Central IDs to perform the pdf analysis on')
    parser.add_argument('-o','--output_file', dest  = 'output', help = 'Output file.')


    
    args = parser.parse_args()
    if args.command == 'get_metadata':
        if not all([args.email,args.ids,args.output ]):       
            parser.error('all parameters -e, -i, -o are required')
        metadata_df = get_metadata(args.email, args.ids)   
        print("Result successfully obtained!")
        if os.path.exists(args.output):
            response = input(f"The file '{args.output}' already exists. Do you want to overwrite it? (y/n): ")
            if response.lower() != 'y':
                print("Operation aborted by the user.")
                return     
        metadata_df.to_csv(args.output)
        print('Metadata saved to {}'.format(args.output))
    else:
        if not all([args.pmc_ids,args.output]):
            parser.error('all parameters -pmc_ids, -o are required')
        df_analysis = pdf_analysis(args.pmc_ids)    
        print("Result successfully obtained!")
        if os.path.exists(args.output):
            response = input(f"The file '{args.output}' already exists. Do you want to overwrite it? (y/n): ")
            if response.lower() != 'y':
                print("Operation aborted by the user.")
                return     
        df_analysis.to_csv(args.output)
        print('Result saved to {}'.format(args.output))
if __name__ == "__main__":
    main()