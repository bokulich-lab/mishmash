from entrezpy_clients._pipelines import _get_run_ids
from entrezpy_clients._efetch import EFetchAnalyzer
import entrezpy.efetch.efetcher as ef

log_level = "ERROR"
nr_jobs = 1
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
    accession_list = args.accession_list

    run_ids = _get_run_ids(email, nr_jobs, accession_list, None, "", log_level)
    efetcher = ef.Efetcher(
        "efetcher", email, apikey=None,
        apikey_var=None, threads=nr_jobs, qid=None
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
        analyzer=EFetchAnalyzer(log_level),
    )
    df = metadata_response.result.metadata_to_df()
    return df
