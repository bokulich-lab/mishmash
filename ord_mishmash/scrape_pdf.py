import sys
import PyPDF2 
from pathlib import Path 
import re
import os
import requests
import xmltodict
from nltk import sent_tokenize, word_tokenize
from collections import Counter
from pdfminer.high_level import extract_pages, extract_text
from pdfminer.layout import  LTTextContainer, LTChar, LTRect, LTFigure



project_studies_pattern1 = r'(PRJ(E|D|N)[A-Z][0-9]+)'
project_studies_pattern2 =  r'((E|D|S)RP[0-9]{6,})'
biosample_studies_pattern1 = r'(SAM(E|D|N)[A-Z]?[0-9])'
biosample_studies_pattern2 =  r'((E|D|S)RS[0-9]{6,})'
experiments_pattern =  r'((E|D|S)RX[0-9]{6,})'
runs_pattern =  r'((E|D|S)RR[0-9]{6,})'
analysis_patern =  r'((E|D|S)RZ[0-9]{6,})'
ref = r'(References)'

primer_metod = r'(16(S|s))'
metagenomic = r'((M|m)etagenomic)'


patterns = [ project_studies_pattern1,
            project_studies_pattern2,
            biosample_studies_pattern1,
            biosample_studies_pattern2,
            experiments_pattern,
            runs_pattern,
            analysis_patern]

class PdfScraper:
    def __init__(self,filename):
        self.filename = filename
        self.reader = PyPDF2.PdfReader(filename)
        self.accession_tuples = None
        self.output_string = None
        self.accession_numbers = None
        self.database_names = None
        self.sra_records_count = None
    
    def get_pages(self) -> str:
        output_string = ''
        count = len(self.reader.pages)
        for i in range(count):
            pageObj =  self.reader.pages[i]
            output_string += pageObj.extract_text()
            self.output_string = output_string
        #print(self.output_string)
        return self.output_string
    def text_extraction(self,element):
        line_text = element.get_text()
    def _is_reference(self,paragraph) -> bool:
        return paragraph.get_text().lstrip().split(' ')[0].lower().rstrip() == 'references'

    def get_pages2(self) -> str:
        page_str = []
        referenceFound = False
        for  page in extract_pages(self.filename):
            if not referenceFound:
                for paragraph in page:
                    if isinstance(paragraph,LTTextContainer):
                        if self._is_reference(paragraph):
                            referenceFound = True
                            break
                        else:
                            page_str.append(paragraph.get_text())
                    else:
                        continue
        return page_str
 

    def get_accession_tuples(self) -> list:
        if self.output_string is None:
            self.get_pages()
        res = []
        for pattern in patterns:
            if re.findall(pattern,self.output_string):
                res += re.findall(pattern,self.output_string)
            self.accession_tuples = res
        return self.accession_tuples
    
    def get_accession_numbers(self) ->list:
        if self.accession_tuples is None:
            self.get_accession_tuples()
        if len(self.accession_tuples) >0: 
            self.accession_numbers = [t[0] for t in self.accession_tuples]
        return self.accession_numbers
    
    def get_database_names(self) ->list:
        if self.accession_tuples is None:
            self.get_accession_tuples()
        if len(self.accession_tuples) >0: 
            self.database_names = [t[1] for t in self.accession_tuples]
        return self.database_names
    
    def get_number_of_records_sra(self) -> int:
        if self.accession_tuples is None:
            self.get_accession_tuples()

        if len(self.get_accession_tuples()) <1:
            self.sra_records_count = 0
            return  self.sra_records_count
        else:
            res_xmls = []
            for n in self.accession_numbers:
                url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term={}'.format(n)
                res = requests.get(url)
                res.raise_for_status()
                res_xmls.append(xmltodict.parse(res.content))
            total_count = 0
            for record in res_xmls:
                total_count += int(record['eSearchResult']['Count'])
                self.sra_records_count = total_count
            return  self.sra_records_count

    def categorize_metods(self,words):
        primer_methods = {'primer','pcr','16s','16 s'}
        metagenomic = {'metagenomic','metagenomics'}
        pmwlen = len(primer_methods.intersection(words))
        mwlen = len(metagenomic.intersection(words))

        if pmwlen > 0 and mwlen ==0:
            return 'primer_metod'
        elif pmwlen == 0 and mwlen > 0:
            return 'metagenomics'
        elif pmwlen > 0 and mwlen > 0:
            return 'both'
        else:
            return 'unknown'
    
    def count_metods(self,sentences):
        sents = Counter()
        for sentence in sentences:
            metod = self.categorize_metods(sentence)
            sents[metod] +=1
        return sents



 


    def parse_metod(self) -> str:
        if self.output_string is None:
            self.get_pages()
        sentences = [[word.lower() for word in word_tokenize(sentence)]
                     for sentence in sent_tokenize(self.output_string)]
        sents = self.count_metods(sentences)
        print(sents)
        return sents
        


        
      



def scrape_pdf(basepath):
     with open("result.txt","w") as file:
        for entry in os.listdir(basepath):
            filename = os.path.join(basepath,entry)
            print(filename)
            c = PdfScraper(filename)
            r = c.get_pages2()
            t = c.parse_metod()
            an = c.get_accession_numbers()
            dn = c.get_database_names()
            file.write('{}   accession number : {}   database : {}\n'.format(filename, an,dn))
 


               
    
def scrape_pdf2():
    print('------------------------------')
    c = PdfScraper('/Users/zsebechle/git/ord_mishmash/ord_mishmash/ord_mishmash/data2/paper7.pdf')
    r = c.get_pages2()
    print(r)
   

if __name__ == "__main__":
  scrape_pdf('/Users/zsebechle/git/ord_mishmash/ord_mishmash/ord_mishmash/data2/')
            
 
