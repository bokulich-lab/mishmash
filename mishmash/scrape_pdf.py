import re
import requests
import sys
import xmltodict
import pandas as pd

from bs4 import BeautifulSoup, Comment
from collections import Counter
from nltk import sent_tokenize, word_tokenize


project_studies_pattern1 = r"(PRJ(E|D|N)[A-Z][0-9]+)"
project_studies_pattern2 = r"((E|D|S)RP[0-9]{6,})"
biosample_studies_pattern1 = r"(SAM(E|D|N)[A-Z]?[0-9])"
biosample_studies_pattern2 = r"((E|D|S)RS[0-9]{6,})"
experiments_pattern = r"((E|D|S)RX[0-9]{6,})"
runs_pattern = r"((E|D|S)RR[0-9]{6,})"
analysis_pattern = r"((E|D|S)RZ[0-9]{6,})"

# primer_method = r"(16(S|s))"
# metagenomic = r"((M|m)etagenomic)"
patterns = [
    project_studies_pattern1,
    project_studies_pattern2,
    biosample_studies_pattern1,
    biosample_studies_pattern2,
    experiments_pattern,
    runs_pattern,
    analysis_pattern,
]


def _contains_blocking_comment(content) -> bool:
    for element in content(string=lambda text: isinstance(text, Comment)):
        if (
            str(element).strip()
            == "The publisher of this article does not allow downloading "
               "of the full text in XML form."
        ):
            return True
    return False


class PMCScraper:
    def __init__(self, pmc_id):
        """
        Class to scrape a pmc_record.
        
        Inputs
        ------
        pmc_id: `int` PMC record ID.

        """
        self.pmc_id = pmc_id
        self.content = None
        self.core_text = None
        self.accession_tuples = None
        self.sra_records_count = None
        self.sra_record_xmls = None

    def get_xml(self) -> object:
        """
        Retrieve the XML record of the paper text.

        Returns
        -------
        self.content: xml
        """
        if self.content is not None:
            return self.content

        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch." \
              "fcgi?db=pmc&id={}".format(self.pmc_id)
        try:
            r = requests.get(url, 3)
            r.raise_for_status()
        except requests.exceptions.Timeout:
            sys.exit(f"Connection to {url} has timed out. Please retry.")
        except requests.exceptions.HTTPError:
            print(
                f"The download URL {url} is likely invalid.\n", flush=True,
            )
            return 0
        except KeyError:
            print("Key error for: " + url, flush=True)
            return 0
        self.content = BeautifulSoup(r.content, features="xml")
        return self.content

    def contains_blocking_comment(self):
        return _contains_blocking_comment(self.get_xml())

    def get_text(self):
        """
        Retrieve the plain text of the record.

        Returns
        -------
        self.core_text: `str` Article text.
        """
        if self.core_text is not None:
            return self.core_text

        content = self.get_xml()
        if _contains_blocking_comment(content):
            raise RuntimeError(
                "The publisher of this article does not allow downloading "
                "the full text in XML form. Please reevaluate manually!"
            )

        core_text_body = content.find("body")

        core_text_back = content.find("back")
        if core_text_back:
            ref_list = core_text_back.find("ref-list")
            for codetag in ref_list.find_all_next():
                codetag.clear()
            self.core_text = core_text_body.text + core_text_back.text
        else:
            self.core_text = core_text_body.text
        return self.core_text

    def get_accession_tuples(self) -> list:
        """
        Retrieve accession numbers and database names from the record.

        Returns
        -------
        self.accession_tuples: `list` List of tuples containing accession
        numbers and their associated database names.
        """
        if self.accession_tuples is not None:
            return self.accession_tuples

        core_text = self.get_text()
        res = []
        for pattern in patterns:
            matches = re.findall(pattern, core_text)
            if matches:
                res += matches
            self.accession_tuples = res
        return self.accession_tuples

    def get_accession_numbers(self) -> set:
        """
        Retrieve accession numbers from the record.

        Returns
        -------
        Set of accession numbers from the paper.
        """
        accession_set = set(t[0] for t in self.get_accession_tuples())
        return accession_set

    def get_database_names(self) -> list:
        """
        Get the database name from the associated accession IDs.

        Returns
        -------
        `str` of characters indicating the database name.
        """
        char_to_value = {"E": "EMBL-EBI European Nucleotide Archive",
                         "D": "NIG DNA Data Bank of Japan",
                         "S": "NCBI Sequence Read Archive",
                         "N": "NCBI Sequence Read Archive"}
        database_shortcuts = set(t[1] for t in self.get_accession_tuples())
        db_set = set(char_to_value[x] for x in database_shortcuts)
        return ", ".join(list(db_set))

    def get_number_of_records_sra(self) -> int:
        """
        Count the total number of INSDC Run records corresponding to all the
        accession numbers found in the paper.

        Returns
        -------
        self.sra_records_count: `int`
        """
        if self.sra_records_count is not None:
            return self.sra_records_count

        if len(self.get_accession_numbers()) < 1:
            return 0

        else:
            if self.sra_record_xmls is None:
                res_xmls = []
                for n in self.get_accession_numbers():
                    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term={}".format(n)
                    res = requests.get(url)
                    res.raise_for_status()
                    res_xmls.append(xmltodict.parse(res.content))
                self.sra_record_xmls = res_xmls

            total_count = 0
            for record in self.sra_record_xmls:
                total_count += int(record["eSearchResult"]["Count"])
            self.sra_records_count = total_count
            return self.sra_records_count

    @staticmethod
    def _categorize_methods(words):
        amplicon_keywords = {"amplicon", "16s", "marker gene",
                             "marker-gene"}
        metagenomic_keywords = {"metagenomic",
                                "metagenomics", "shotgun"}
        pmwlen = len(amplicon_keywords.intersection(words))
        mwlen = len(metagenomic_keywords.intersection(words))

        if pmwlen > 0 and mwlen == 0:
            return "amplicon_method"
        elif pmwlen == 0 and mwlen > 0:
            return "metagenomics_method"
        elif pmwlen > 0 and mwlen > 0:
            return "both_methods"
        return "no_method"

    def _count_methods(self, sentences):
        sents = Counter()
        for sentence in sentences:
            method = self._categorize_methods(sentence)
            sents[method] += 1
        return sents

    def parse_method(self) -> dict:
        """
        Get the size of each group method (prime_method, metagenomics, both,
        unknown).

        Returns 
        -------
        `dict` providing total counts for each group category.

        """
        sentences = [
            [word.lower() for word in word_tokenize(sentence)]
            for sentence in sent_tokenize(self.get_text())
        ]
        method_dict = dict(self._count_methods(sentences))
        return method_dict

    def get_pcr_primers(self) -> list:
        """
        Get PCR primers.

        Returns
        -------
        `list` of PCR primers as found in the article text.

        """
        pcr_pattern = r"((?:A|G|C|T|N|W|V|M|H){9}(?:A|G|C|T|N|W|V|M|H)+)"
        res = []
        if re.findall(pcr_pattern, str(self.get_text())):
            res += re.findall(pcr_pattern, self.get_text())
        return ", ".join(res)


def analyze_pdf(args) -> dict:
    """
    Gives overview of the paper with respect to the predefined metrics.

    Args
    ----
    args

    Returns
    -------
    df: `pd.DataFrame` summarizing PMC objects.

    """
    pmc_ids = args.pmc_list
    requested_objects = [PMCScraper(id) for id in pmc_ids]
    scrape_objects = [
        x for x in filter(lambda el: not el.contains_blocking_comment(),
                          requested_objects)
    ]
    forbidden_objects = [
        x for x in filter(lambda l: l.contains_blocking_comment(),
                          requested_objects)
    ]
    if len(forbidden_objects) > 0:
        print(
            "Papers represented by following PMC IDs were not fetched; "
            "the publisher of this article does not allow "
            + "downloading of the full text in XML form:"
        )
        for el in forbidden_objects:
            print(el.pmc_id)

    df = pd.DataFrame(
        {
            "PMC ID": [],
            "INSDC Accession Numbers": [],
            "INSDC Database": [],
            "Number of Sequence Records": [],
            "Primer Sequences": [],
            "Sequencing Method": [],
        }
    )
    for el in scrape_objects:
        tmp_df = pd.DataFrame(
            {
                "PMC ID": [el.pmc_id],
                "INSDC Accession Numbers":
                    [", ".join(list(el.get_accession_numbers()))],
                "INSDC Database": [el.get_database_names()],
                "Number of Sequence Records": [el.get_number_of_records_sra()],
                "Primer Sequences": [el.get_pcr_primers()],
                "Sequencing Method": [el.parse_method()],
            }
        )
        df = pd.concat([df, tmp_df])

    return df.set_index("PMC ID", drop=True)
