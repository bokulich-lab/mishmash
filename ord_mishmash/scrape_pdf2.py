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
from xml.etree import ElementTree as ET
import urllib.request
from bs4 import BeautifulSoup, Comment

#PRJNA605207
#8751733

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

class PMCscraper:
    def __init__(self,pmc_id):
        ''' Class to scrape a pmc_record.
        
        Inputs
        ------
        pmc_id = int
                PMC record id.

        '''    
        self.pmc_id = pmc_id
        self.core_text = None
        self.accession_tuples = None
        self.accession_numbers = None
        self.database_names = None
        self.sra_records_count = 0
        self.metod_count = None
        self.pcr_primers = None
        self.accession_string = False
        self.content = None
    def get_xml(self) -> object:
        ''' Get the xml record of the text of the paper.

        Returns
        -------
        self.content : xml
        '''
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id={}'.format(self.pmc_id)
        try:
            r = requests.get(url, 3)
            r.raise_for_status()
        except requests.exceptions.Timeout:
            sys.exit(f"Connection to {url} has timed out. Please retry.")
        except requests.exceptions.HTTPError:
            print( f"The download URL:  {url}  is likely invalid.\n",
                  flush=True,
                 )
            return 0
        except KeyError:
            print("Key error for: " + url, flush=True)
            return 0         
        self.content = BeautifulSoup(r.content,features = "xml" )
        return self.content
    
    
    def _contains_blocking_comment(self) -> bool:
   
        if self.content is None:
            self.get_xml()
        
        for element in self.content(string=lambda text: isinstance(text, Comment)):
            if str(element) == 'The publisher of this article does not allow downloading of the full text in XML form.':
                return True
        return False

    class XmlException(Exception):
        pass
    
    def get_text(self):
        ''' Get the plain text of the record.

        Returns
        -------
        self.core_text : str
                        Text of the paper.
        '''
        #if self._contains_blocking_comment():
            #raise XmlException('Not allowed to be scrapped, evaluate manually!')
            
            #return 0
        core_text_body = self.content.find("body")

        core_text_back = self.content.find("back")
        ref_list = core_text_back.find('ref-list')
        for codetag in ref_list.find_all_next():
                codetag.clear() 
        self.core_text = core_text_body.text + core_text_back.text
        return self.core_text
      
    def get_accession_tuples(self) -> list:
        '''Get the accession tuples from the record.

        Returns
        -------
        self.accession_tuples : list 
                            List of tuples - accession numbers, database name.
        '''
        if self.core_text is None:
            self.get_text()
        res = []
        for pattern in patterns:
            if re.findall(pattern, str(self.core_text)):
                res+=re.findall(pattern,self.core_text)
            self.accession_tuples = res
        return self.accession_tuples 

    def get_accession_numbers(self) -> set:
        '''Get the accession numbers from the record.

        Returns
        -------
        self.accession_numbers : set
                          Set of accesion numbers of the paper.
        '''
        if self.accession_tuples is None:
            self.get_accession_tuples()
        if len(self.accession_tuples) >0: 
            self.accession_numbers = set(t[0] for t in self.accession_tuples)
            return self.accession_numbers
        return set()
    
    def get_accession_number_string(self) -> bool:
        '''Get the string: accession
        
        Returns
        -------
        self.accession_string : bool
                        True if there is a match, False otherwise.
        '''
        if self.core_text is None:
            self.get_text()
        if re.search(r'\baccession\b',self.core_text):
            return True
        return False


    def get_database_names(self) ->list:
        ''' Get the database name character.

        Returns
        -------
        self.database_names : list
                    List of characters indicating the database name.
       
        '''
        if self.accession_tuples is None:
            self.get_accession_tuples()
        if len(self.accession_tuples) >0: 
            self.database_names = set(t[1] for t in self.accession_tuples)
        return self.database_names
    
    def get_number_of_records_sra(self) -> int:
        ''' Get total count of records corresponding to the AN.

        Returns
        -------
        self.sra_record_count : int
        '''
        if self.accession_numbers is None:
            self.get_accession_numbers()

        if len(self.accession_numbers) <1:
            self.sra_records_count = 0
            return  self.sra_records_count
        else:
            res_xmls = []
            for n in self.accession_numbers:
                print(n)
                url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term={}'.format(n)
                res = requests.get(url)
                res.raise_for_status()
                res_xmls.append(xmltodict.parse(res.content))
            total_count = 0
            for record in res_xmls:
                total_count += int(record['eSearchResult']['Count'])
                self.sra_records_count = total_count
            return  self.sra_records_count     
    
    def _categorize_metods(self,words):
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
    
    def _count_metods(self,sentences):
        sents = Counter()
        for sentence in sentences:
            metod = self._categorize_metods(sentence)
            sents[metod] +=1
        return sents

    def parse_metod(self) -> Counter:
        '''Get the size of each group metod (prime_metod,metagenomics,both,unknowm).

        Returns 
        -------
        self.metod_count : Counter    

        '''
        if self.core_text is None:
            self.get_text()
        sentences = [[word.lower() for word in word_tokenize(sentence)]
                     for sentence in sent_tokenize(self.core_text)]
        sents = self._count_metods(sentences)
        self.metod_count = sents    
        return sents
        
    def get_pcr_primer(self) -> list:
        '''Get pcr primers.

        Returns
        -------
        self.pcr_primers : list
        '''
        if self.core_text is None:
            self.get_text()
        pcr_pattern = r'((?:A|G|C|T|N|W|V|M|H){9}(?:A|G|C|T|N|W|V|M|H)+)'
        res = []
        if re.findall(pcr_pattern, str(self.core_text)):
                res+=re.findall(pcr_pattern,self.core_text)
        self.pcr_primers = res
        return self.pcr_primers 




               
#7694690, 9714783 ,7802287, 7501109, 7497576 , 7233940
#5464538 , 8156656, 5537870 , 7524393
#8880058, 7977168 (has 3id)
#6017827,8551616(National Microbiology Data Center)
#6033410,8231616,8258384
#8472786,5522834,8751733,8573958,8916702,6887738,3876090,8777273,7286428,6827403,6232675,8300211,5490539,5950546,7524302,6759125,7524390   
#7071329,3682059
#6595049,5647777
#7193032,8794796(2 ids),7197999,9714783,6678981,6627124,9431030

#if __name__ == "__main__":
#  c = PMCscraper(9431030)
#  r= c.get_text()
#  l =c.get_accession_numbers()
#  k = c.get_pcr_primer()
#  print(l)
#  print(k)
