import os
import unittest
from parameterized import parameterized
from pathlib import Path 
from ord_mishmash import PMCscraper


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



    @parameterized.expand([(text_file_1,[]),
                           (text_file_7,[('PRJNA605207', 'N')]),
                           (text_file_8,[('PRJNA604899','N'),('PRJNA604957','N'),('PRJNA605597','N')])])

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
                           (text_file_5,['GAACTGCCCACCAACTACAA','CCATCGTGGACAGACATGAA','CCAGTGAGTCTTGGCTGACA','TTCAGTGAGGCACAGAATGC']),
                           (text_file_6,['ACTCCTACGGGAGGCAGCAGT','GGACTACNVGGGTWTCTAAT'])
                           ])
    def test_get_pcr_primer(self,data,expected_res):
        a = PMCscraper('id')
        with open(data) as f:
            a.core_text = f.read()
        res = a.get_pcr_primer()
        self.assertEqual(res,expected_res)


    

        


if __name__=='__main__':
    unittest.main()