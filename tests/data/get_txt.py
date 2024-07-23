import requests
from bs4 import BeautifulSoup


s = [7977168]

for i in s:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id={}".format(
        i
    )
    r = requests.get(url)
    content = BeautifulSoup(r.content, features="xml")
    core_text_body = content.find("body")

    core_text_back = content.find("back")
    ref_list = core_text_back.find("ref-list")
    for codetag in ref_list.find_all_next():
        codetag.clear()
    core_text = core_text_body.text + core_text_back.text
    with open("test_sample_8.txt", "w") as f:
        f.write(core_text)
