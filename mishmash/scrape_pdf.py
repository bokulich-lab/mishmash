import importlib.resources
import json
import re
import requests
import sys
import xmltodict

import numpy as np
import pandas as pd

from bs4 import BeautifulSoup, Comment
from collections import Counter
from nltk import sent_tokenize, word_tokenize
from pathlib import Path
from urllib.parse import urlparse


project_studies_pattern1 = r"(PRJ(E|D|N)[A-Z][0-9]{4,6})"
project_studies_pattern2 = r"((E|D|S)RP[0-9]{6,})"
biosample_studies_pattern1 = r"(SAM(E|D|N)[0-9]{8,})"
biosample_studies_pattern2 = r"((E|D|S)RS[0-9]{6,})"
experiments_pattern = r"((E|D|S)RX[0-9]{6,})"
runs_pattern = r"((E|D|S)RR[0-9]{6,})"
analysis_pattern = r"((E|D|S)RZ[0-9]{6,})"
old_sra_pattern = r"((S|D)RA[0-9]{6,}(\.[0-9]+)?)"

patterns = [
    project_studies_pattern1,
    project_studies_pattern2,
    biosample_studies_pattern1,
    biosample_studies_pattern2,
    experiments_pattern,
    runs_pattern,
    analysis_pattern,
    old_sra_pattern
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
        self.method_dict = {}

        # Properties of journal and corresponding author
        self.journal_name = None
        self.publisher_name = None
        self.publish_year = None
        self.institution = None

    def get_xml(self) -> object:
        """
        Retrieve the XML record of the paper text.

        Returns
        -------
        self.content: xml
        """
        if self.content:
            return self.content

        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch." \
              "fcgi?db=pmc&id={}".format(self.pmc_id)
        for i in range(5):
            try:
                r = requests.get(url, 3)
                r.raise_for_status()
                break
            except requests.exceptions.Timeout:
                sys.exit(f"Connection to {url} has timed out. Please retry.")
            except requests.exceptions.HTTPError:
                print(
                    f"The download URL {url} is likely invalid.\n", flush=True,
                )
                continue
            except KeyError:
                print("Key error for: " + url, flush=True)
                return 0

        self.content = BeautifulSoup(r.content, features="xml")
        return self.content

    def contains_blocking_comment(self):
        return _contains_blocking_comment(self.get_xml())

    def get_journal_name(self):
        if not self.journal_name:
            self.get_text()

        return self.journal_name

    def get_publisher_name(self):
        if not self.publisher_name:
            self.get_text()

        return self.publisher_name

    def get_publish_year(self):
        if not self.publish_year:
            self.get_text()

        return self.publish_year

    def get_institution(self):
        if not self.institution:
            self.get_text()

        return self.institution

    def get_text(self):
        """
        Retrieve the plain text of the record.

        Returns
        -------
        self.core_text: `str` Article text.
        """
        if self.core_text:
            return self.core_text

        content = self.get_xml()
        if _contains_blocking_comment(content):
            raise RuntimeError(
                "The publisher of this article does not allow downloading "
                "the full text in XML form. Please reevaluate manually!"
            )

        # Get journal properties
        try:
            self.journal_name = content.find("journal-title").text
        except AttributeError:
            self.journal_name = None           
 
        try:
            self.publisher_name = content.find("publisher-name").text
        except AttributeError:
            self.publisher_name = None

        # Get publish date
        pmc_year_tags = ["pmc-release", "epub", "accepted"]
        for tag in pmc_year_tags:
            pub_date_pmc = content.find_all("pub-date",
                                            {"pub-type": re.compile(tag)})
            if len(pub_date_pmc) > 0:
                self.publish_year = pub_date_pmc[0].year.text
                break

        # Get institutional affiliation of first corresponding author
        cor_list = content.find("contrib", {"contrib-type": "author"})

        try:
            aff_tag = cor_list.find("xref", {"ref-type": "aff"})["rid"]
            inst_list = content.find("aff", {"id": aff_tag}).\
                find_all("institution")
            inst_str = "".join([inst.text for inst in inst_list])
            inst_str = inst_str.rstrip(" ,.;:")
            self.institution = inst_str
        except (IndexError, AttributeError, TypeError):
            pass

        core_text_body = content.find("body")
        core_text_back = content.find("back")
        core_text_front = content.find("front")

        self.core_text = core_text_body.text

        if core_text_back:
            ref_list = core_text_back.find("ref-list")
            for codetag in ref_list.find_all_next():
                codetag.clear()
            self.core_text += core_text_back.text
        if core_text_front:
            self.core_text += core_text_front.text

        return self.core_text

    def get_accession_tuples(self) -> list:
        """
        Retrieve accession numbers and database names from the record.

        Returns
        -------
        self.accession_tuples: `list` List of tuples containing accession
        numbers and their associated database names.
        """
        if self.accession_tuples:
            return self.accession_tuples

        core_text = self.get_text()
        res = []
        for pattern in patterns:
            matches = re.findall(pattern, core_text)
            if matches:
                res += matches
            self.accession_tuples = res

        return list(set(self.accession_tuples))

    def check_non_insdc_db(self) -> str:
        # Checks text for keywords that may denote data upload in non-INSDC
        # databases

        core_text = self.get_text()
        tk_text = sent_tokenize(core_text)

        # Potential keywords, non-case sensitive
        db_name_list = ["figshare", "ega", "european phenome-genome archive",
                        "national genomics data center", "gsa",
                        "genome sequence archive", "ngdc",
                        "china national center for bioinformation",
                        "cncb", "mg-rast", "metagenomic rapid annotations "
                                           "using subsystems technology",
                        "metagenomics rast", "cnsa", "cngb sequence archive",
                        "cngbdb", "china national genebank database"]
        db_name_dict = {"figshare": "Figshare",
                        "ega": "European Phenome-Genome Archive",
                        "european phenome-genome archive":
                            "European Phenome-Genome Archive",
                        "gsa": "China National Center for Bioinformation: "
                               "Genome Sequence Archive",
                        "genome sequence archive": "China National Center for "
                                                   "Bioinformation: Genome "
                                                   "Sequence Archive",
                        "ngdc": "China National Center for Bioinformation: "
                                "National Genomics Data Center",
                        "china national center for bioinformation":
                            "China National Center for Bioinformation",
                        "cncb": "China National Center for Bioinformation:",
                        "mg-rast": "MG-RAST",
                        "metagenomic rapid annotations using subsystems "
                            "technology": "MG-RAST",
                        "metagenomics rast": "MG-RAST",
                        "cnsa": "China National GeneBank Database Sequence "
                                "Archive",
                        "cngb sequence archive": "China National GeneBank "
                                                 "Database Sequence Archive",
                        "cngbdb": "China National GeneBank Database Sequence "
                                  "Archive",
                        "china national genebank database":
                            "China National GeneBank Database Sequence Archive"}
        db_name_re = [fr'(\A|\W)({name})(\W|\Z)' for name in db_name_list]
        db_match_list = [re.search(query, sent, re.IGNORECASE).group(2)
                         for query
                         in db_name_re for sent in tk_text
                         if re.search(query, sent, re.IGNORECASE)]

        url_search_list = ["figshare.com", "ega-archive.org",
                           "ngdc.cncb.ac.cn/gsa", "mg-rast.org",
                           "metagenomics.anl.gov", "db.cngb.org/cnsa"]
        url_search_re = [fr'(\A|\W)({name})(\W|\Z)' for name in
                         url_search_list]
        url_match_list = [re.search(query, sent, re.IGNORECASE).group(2)
                          for query in url_search_re for sent in tk_text
                          if re.search(query, sent, re.IGNORECASE)]

        prep_phrase_list = ["found in", "found at", "deposited in",
                            "deposited into", "deposited on",
                            "accessible at",
                            "available in", "available from", "available on"]
        prep_phrase_re = [fr'(\A|\W)({name})(\W|\Z)' for name in
                          prep_phrase_list]
        prep_match_list = [re.search(query, sent, re.IGNORECASE).group(2)
                           for query in prep_phrase_re for sent in tk_text
                           if re.search(query, sent, re.IGNORECASE)]

        id_word_list = ["accession number", "accession number(s)",
                        "accession numbers",
                        "accession ID", "project ID",
                        "project access number",
                        "ID numbers", r'CRA([0-9]{6})', r'CNP([0-9]{6})',
                        r'[0-9]{7}\.[0-9]', r'mgp[0-9]{5}']
        id_word_re = [fr'(\A|\W)({name})(\W|\Z)' for name in id_word_list]
        id_match_list = [re.search(query, sent, re.IGNORECASE).group(2)
                         for query in id_word_re for sent in tk_text
                         if re.search(query, sent, re.IGNORECASE)]

        # Assume dict insertion order as of Python 3.7+
        match_dict = {"db"   : db_match_list,
                      "url"  : url_match_list,
                      "prep" : prep_match_list,
                      "id"   : id_match_list}

        # For debugging purposes
        match_len_df = pd.DataFrame(data=[[len(d) for d
                                           in match_dict.values()]],
                                    columns=list(match_dict.keys()),
                                    index=[self.pmc_id])

        # To get at least Bronze:
        # >= 1 hit for DB, and # hits (URL + prep + ID) >= 1
        if len(db_match_list) > 0 and (len(url_match_list) + len(
                prep_match_list) + len(id_match_list) > 0):
            db_count = Counter(db_match_list)
            db = db_name_dict[max(db_count).lower()]
            return db

        return None

    def get_accession_numbers(self) -> list:
        """
        Retrieve accession numbers from the record.

        Returns
        -------
        Set of accession numbers from the paper.
        """

        retrieved_tuples = self.get_accession_tuples()
        if retrieved_tuples:
            accession_set = list(set(t[0] for t in
                                     self.get_accession_tuples()))
            return accession_set

        return None

    def get_database_names(self) -> str:
        """
        Get the database name from the associated accession IDs.

        Returns
        -------
        `str` of characters indicating the database name.
        """
        char_to_value = {"E": "EMBL-EBI European Nucleotide Archive",
                         "D": "NIG DNA Data Bank of Japan",
                         "S": "NCBI Sequence Read Archive",
                         "N": "NCBI Sequence Read Archive",
                         ".": "Unknown; old code."}

        retrieved_tuples = self.get_accession_tuples()
        if retrieved_tuples:
            database_shortcuts = set(t[1][0] for t in retrieved_tuples)
            db_set = set(char_to_value[x] for x in database_shortcuts if x in
                         char_to_value.keys())
            return ", ".join(list(db_set))

        return None

    def get_number_of_records_sra(self) -> int:
        """
        Count the total number of INSDC Run records corresponding to all the
        accession numbers found in the paper.

        Returns
        -------
        self.sra_records_count: `int`
        """
        if self.sra_records_count:
            return self.sra_records_count

        retrieved_accession_numbers = self.get_accession_numbers()
        if not retrieved_accession_numbers:
            return 0

        if len(retrieved_accession_numbers) < 1:
            return 0

        # Record count has not yet been processed
        res_xmls = []
        for n in retrieved_accession_numbers:
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/" \
                  "esearch.fcgi?db=sra&term={}".format(n)

            for j in range(5):
                try:
                    res = requests.get(url)
                    res.raise_for_status()
                    break
                except requests.exceptions.Timeout:
                    sys.exit(
                        f"Connection to {url} has timed out. "
                        f"Please retry.")
                except requests.exceptions.HTTPError:
                    print(
                        f"The download URL {url} is likely invalid.\n",
                        flush=True,
                    )
                    continue
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
                             "marker-gene", "ITS", "ITS1", "ITS2"}
        shotgun_keywords = {"metagenomic", "metagenomics", "shotgun",
                            "whole genome", "whole-genome", "genom"}
        pmwlen = len(amplicon_keywords.intersection(words))
        mwlen = len(shotgun_keywords.intersection(words))

        if pmwlen > 0 and mwlen == 0:
            return "amplicon_method"
        elif pmwlen == 0 and mwlen > 0:
            return "shotgun_method"
        elif pmwlen > 0 and mwlen > 0:
            return "both_methods"
        return "no_method"

    def get_method_weights(self):
        self._parse_article_text()

        if not self.method_dict:
            return None

        method_list = ["amplicon_method", "shotgun_method", "both_methods",
                       "no_method"]

        method_dict = {
            method: (self.method_dict[method] if method in
                     self.method_dict.keys() else 0)
            for method in method_list
        }

        no_method_count = method_dict["no_method"]
        total_count = sum(method_dict.values())
        if no_method_count == total_count:
            return None

        notnull_method_count = total_count - no_method_count

        weight_dict = {
            "amplicon": round((method_dict["amplicon_method"] +
                               method_dict["both_methods"] / 2) /
                              notnull_method_count, 2),
            "shotgun": round((method_dict["shotgun_method"] +
                              method_dict["both_methods"] / 2) /
                             notnull_method_count, 2)
        }
        return weight_dict

    def _count_methods(self, sentences):
        sents = Counter()
        for sentence in sentences:
            method = self._categorize_methods(sentence)
            sents[method] += 1
        return sents

    def _parse_article_text(self):
        """
        Get the size of each group method.

        Returns 
        -------
        `dict` providing total counts for each group category.

        """
        sentences = [
            [word.lower() for word in word_tokenize(sentence)]
            for sentence in sent_tokenize(self.get_text())
        ]
        method_dict = dict(self._count_methods(sentences))
        self.method_dict = method_dict
        return

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

    def get_code_links(self):
        url_list = re.findall(r"(https?://\S+)", str(self.get_text()))
        url_list = list(set([url.rstrip(",.;:)]") for url in url_list]))
        code_dict = {"url": None,
                     "has_link": "False"}

        if url_list:
            # Remove URLs that match exclusion criteria i.e. common tools
            with importlib.resources.open_text("mishmash", "urls.json") as file:
                exclude_dict = json.load(file)

            exclude_list = [ex_url.lower() for ex_url in exclude_dict[
                "exclude"]]
            not_in_excl = [url for url in url_list if not
                           any([excl in url.lower() for excl in exclude_list])]

            # Check whether any of these URLs match known repositories
            repo_host_list = ["github", "zenodo", "bitbucket", "figshare",
                              "codeocean"]
            token_list = [urlparse(url) for url in not_in_excl]
            is_repo_match = [any([repo in url.netloc for repo in
                                  repo_host_list]) for
                             url in token_list]

            if any(is_repo_match):
                code_dict["url"] = np.array(url_list)[np.where(
                    is_repo_match)[0]].tolist()
                code_dict["has_link"] = "True"
            else:
                code_dict["has_link"] = "Possible: URL found in paper."

        # If no URLs are found while scraping
        else:
            sentences = [
                [word.lower() for word in word_tokenize(sentence)]
                for sentence in sent_tokenize(self.get_text())
            ]
            repo_keywords = {"github", "zenodo", "bitbucket", "figshare",
                             "code ocean", "codeocean" "repository"}
            repo_match = [any(repo_keywords.intersection(words))
                          for words in sentences]
            if any(repo_match):
                code_dict["has_link"] = "Possible: Repository keywords found " \
                                        "in paper."

        return code_dict


def _check_input_file(inp_file):
    if Path(inp_file).is_file():
        id_df = pd.read_csv(inp_file, header=None)
        id_df.dropna(axis=0, how="any",
                     subset=id_df.columns[0],
                     inplace=True,
                     ignore_index=True)

        if id_df.shape[0] > 0:
            return id_df[id_df.columns[0]].tolist()

        print(f"The provided file does not contain input IDs! Please check "
              f"and try again: {inp_file}")
        exit(1)

    print(f"The provided file path is not valid! Please check and try "
          f"again: {inp_file}")


def analyze_pdf(args,
                pmc_ids: list = None):
    """
    Gives overview of the paper with respect to the predefined metrics.

    Args
    ----
    args
    pmc_ids: :list:

    """
    if args:
        if args.pmc_list:
            pmc_ids = args.pmc_list
        elif args.pmc_input_file:
            pmc_ids = _check_input_file(args.pmc_input_file)
        else:
            print("Input PMC IDs must be provided via either the --pmc_list or "
                  "--pmc_input_file flag! Please check your command "
                  "and try again.")
            exit(1)
    if not pmc_ids:
        print("No input PMC IDs have been detected! "
              "Please check your command and try again.")
        exit(1)

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
            "Papers represented by the following PMC IDs were not fetched; "
            "the publisher of this article does not allow "
            + "downloading of the full text in XML form:"
        )
        for el in forbidden_objects:
            print(el.pmc_id)

    df = pd.DataFrame(
        {
            "PMC ID": [],
            "Sequence Accessibility Badge": [],
            "INSDC Accession Numbers": [],
            "INSDC Database": [],
            "Number of Sequence Records": [],
            "Primer Sequences": [],
            "Sequencing Method Probability": [],
            "Includes Code Repository": [],
            "Code URL": []
        }
    )

    if args and args.include_journal_data:
        df = pd.DataFrame(
            {
                "PMC ID": [],
                "Sequence Accessibility Badge": [],
                "INSDC Accession Numbers": [],
                "INSDC Database": [],
                "Number of Sequence Records": [],
                "Primer Sequences": [],
                "Sequencing Method Probability": [],
                "Includes Code Repository": [],
                "Code URL": [],
                "Publication Year": [],
                "Journal Name": [],
                "Publisher Name": [],
                "First Author Affiliation": []
            }
        )

    for el in scrape_objects:
        pmc_id = el.pmc_id

        insdc_id_list = el.get_accession_numbers()
        if insdc_id_list:
            insdc_id_list = ", ".join(insdc_id_list)

        seq_db = el.get_database_names()
        num_seqs = el.get_number_of_records_sra()
        primer_seqs = el.get_pcr_primers()
        method_prob = el.get_method_weights()
        code_dict = el.get_code_links()
        publish_year = el.get_publish_year()
        journal_name = el.get_journal_name()
        publisher_name = el.get_publisher_name()
        institution = el.get_institution()

        # Alternative check for non-INSDC database hit
        if num_seqs == 0:
            other_db = el.check_non_insdc_db()
            seq_db = other_db

        # Evaluate badge qualifications
        output_badge = "None"
        missing_steps = ""

        # Accession numbers found OR non-INSDC database hit
        if (num_seqs > 0) or (other_db):
            output_badge = "Bronze"  # At minimum

            if method_prob:
                if (method_prob["amplicon"] > method_prob["shotgun"]) and not \
                        primer_seqs:
                    missing_steps += "Primer sequences could not be found! " \
                                     "May require manual review."
                else:
                    output_badge = "Silver"

                    if code_dict["url"]:
                        output_badge = "Gold"
                    else:
                        missing_steps += "Link to code repository could not " \
                                         "be found! May require manual " \
                                         "review."

            else:  # No methods to be found
                missing_steps += "Text does not clearly denote whether an " \
                                 "amplicon or shotgun sequencing paper! May " \
                                 "require manual review."

        else:
            missing_steps += "INSDC accession numbers with corresponding Run " \
                             "IDs could not be found! May require manual " \
                             "review."

        if missing_steps:
            output_badge = f"{output_badge}: {missing_steps}"

        tmp_df = pd.DataFrame(
            {
                "PMC ID": [pmc_id],
                "Sequence Accessibility Badge": [output_badge],
                "INSDC Accession Numbers":
                    [insdc_id_list],
                "Sequence Database": [seq_db],
                "Number of Sequence Records": [num_seqs],
                "Primer Sequences": [primer_seqs],
                "Sequencing Method Probability": [method_prob],
                "Includes Code Repository": [code_dict["has_link"]],
                "Code URL": [code_dict["url"]]
            }
        )

        if args and args.include_journal_data:
            tmp_df = pd.DataFrame(
                {
                    "PMC ID": [pmc_id],
                    "Sequence Accessibility Badge": [output_badge],
                    "INSDC Accession Numbers":
                        [insdc_id_list],
                    "Sequence Database": [seq_db],
                    "Number of Sequence Records": [num_seqs],
                    "Primer Sequences": [primer_seqs],
                    "Sequencing Method Probability": [method_prob],
                    "Includes Code Repository": [code_dict["has_link"]],
                    "Code URL": [code_dict["url"]],
                    "Publication Year": [publish_year],
                    "Journal Name": [journal_name],
                    "Publisher Name": [publisher_name],
                    "First Author Affiliation": [institution]
                }
            )

        df = pd.concat([df, tmp_df])

    return df.set_index("PMC ID", drop=True)
