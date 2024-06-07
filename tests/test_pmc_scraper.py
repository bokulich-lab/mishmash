import os
import unittest
from parameterized import parameterized
from pathlib import Path 
from ord_mishmash import PMCscraper,pdf_analysis
from bs4 import BeautifulSoup
import pandas as pd
import xmltodict


THIS_DIR = os.path.dirname(os.path.abspath(__file__))

def fpath(fname):
    return os.path.join(THIS_DIR, fname)

class TestPMCScraper(unittest.TestCase):
    text_file_1 = fpath("data/test_sample_1.txt")
    text_file_2 = fpath("data/test_sample_2.txt")
    text_file_3 = fpath("data/test_sample_3.txt")
    text_file_4 = fpath("data/test_sample_4.txt")
    text_file_5 = fpath("data/test_sample_5.txt")
    text_file_6 = fpath("data/test_sample_6.txt")
    text_file_7 = fpath("data/test_sample_7.txt")
    text_file_8 = fpath("data/test_sample_8.txt")
    xml_file_2 = fpath("data/test_sample_2.xml")
    xml_file_3 = fpath("data/test_sample_3.xml")
    xml_file_9 = fpath("data/test_sample_9.xml")
    xml_file_5 = fpath("data/test_sample_5.xml")
    xml_file_6 = fpath("data/test_sample_6.xml")
    
    @parameterized.expand([(xml_file_2,'\nIntroduction.\n'),
                           (xml_file_9,'\nIntroduction\n\nConclusion\n\n\n\n\n')
                           ])
    def test_get_text(self,data,expected_value):
        a = PMCscraper('id')
        with open(data) as f:
            a.content = BeautifulSoup(f.read(),features = "xml" )
        res = a.get_text()
        self.assertEqual(res,expected_value)

    def test_get_text_exception(self):
        xml_file_1 = fpath("data/test_sample_1.xml")
        a = PMCscraper('id')
        with open(xml_file_1) as f:
            a.content = BeautifulSoup(f.read(),features = "xml" )
        self.assertRaises(RuntimeError, a.get_text)

    @parameterized.expand([(text_file_1,[]),
                           (text_file_7,[('PRJNA605207', 'N')]),
                           (text_file_8,[('PRJNA604899','N'),
                                         ('PRJNA604957','N'),('PRJNA605597','N')])])
    def test_accession_tuples(self,data,expected_res):
        a = PMCscraper('id')
        with open(data) as f:
            
            a.core_text = f.read()
        res = a.get_accession_tuples()
        self.assertEqual(res,expected_res)    
  
    
    @parameterized.expand([(text_file_4),
                           (text_file_5),
                           (text_file_6),
                           (text_file_7)
                           ])

    def test_accession_string(self,data):
        a = PMCscraper('id')
        with open(data) as f:
            a.core_text = f.read()
        res = a.get_accession_number_string()
        self.assertTrue(res, msg= None)
    
    @parameterized.expand([(text_file_1,[]),
                           (text_file_3,['GTGCCAGCMGCCGCGGTAA','GGACTACHVGGGTWTCTAAT']),
                           (text_file_5,['GAACTGCCCACCAACTACAA','CCATCGTGGACAGACATGAA',
                                         'CCAGTGAGTCTTGGCTGACA','TTCAGTGAGGCACAGAATGC']),
                           (text_file_6,['ACTCCTACGGGAGGCAGCAGT','GGACTACNVGGGTWTCTAAT'])
                           ])
    def test_get_pcr_primer(self,data,expected_res):
        a = PMCscraper('id')
        with open(data) as f:
            a.core_text = f.read()
        res = a.get_pcr_primer()
        self.assertEqual(res,expected_res)

    @parameterized.expand([(text_file_1,{'unknown': 13, 'primer_metod': 1,
                                                  'metagenomics': 1}),
                           (text_file_2,{'unknown': 11, 'primer_metod': 2}),
                           (text_file_3,{'unknown': 10, 'primer_metod': 3})
                           ])
    def test_parse_method(self,data,expected_res):
        a = PMCscraper('id')
        with open(data) as f:
            a.core_text = f.read()
        res = a.parse_method()
        self.assertEqual(res,expected_res)
    
    @parameterized.expand([([xml_file_5],163),
                           ([xml_file_5,xml_file_6],173)     
                                ])
    def test_sra_count(self,data,expected_res):
         a = PMCscraper('id')
         a.accession_tuples = [('XXXXX','X')]
         xml_files = []
         for el in data:
            with open(el,'r', encoding='utf-8') as f:
                tmp_file = f.read()
                xml_files.append(xmltodict.parse(tmp_file))
         a.sra_record_xmls = xml_files
         res = a.get_number_of_records_sra()
         self.assertEqual(res,expected_res)   
    
if __name__=='__main__':
    unittest.main()