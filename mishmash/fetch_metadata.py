from .entrezpy_clients._pipelines import _get_run_ids
from .entrezpy_clients._efetch import EFetchAnalyzer
from .scrape_pdf import _check_input_file
import entrezpy.efetch.efetcher as ef

test_ids = ["ERROR"]


def get_metadata(args) -> object:
    """
    Fetch the metadata of corresponding IDs.
    
    Args
    ------
    args : pass
    
    Returns
    -------
    df : dataframe of the metadata collection

    """
    email = args.email
    n_jobs = args.n_jobs

    if args.accession_list:
        accession_list = args.accession_list
    elif args.accession_input_file:
        accession_list = _check_input_file(args.accession_input_file)
    else:
        print("Input accession ID file must be provided via either the "
              "--accession_list or --accession_input_file flag! Please check "
              "your command and try again.")
        exit(1)

    assert isinstance(n_jobs, int)

    run_ids = _get_run_ids(email, accession_list, None, "", n_jobs, "ERROR")

    efetcher = ef.Efetcher(
        "efetcher", email, apikey=None,
        apikey_var=None, threads=n_jobs, qid=None
    )

    metadata_response = efetcher.inquire(
        {
            "db": "sra",
            "id": run_ids,
            "rettype": "xml",
            "retmode": "xml",
            "retmax": len(run_ids),
            "reqsize": 150,
        },
        analyzer=EFetchAnalyzer("ERROR"),
    )
    df = metadata_response.result.metadata_to_df()
    return df
