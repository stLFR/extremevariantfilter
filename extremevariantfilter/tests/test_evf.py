from unittest import TestCase
import extremevariantfilter as evf

class TestCheck_Type(TestCase):
    def test_args(self):
        for arg in ["SNP", "INDEL"]:
            try:
                evf.Check_Type(arg)
            except:
                self.fail("Check_Type encountered an unexpected exception")
    def test_fail(self):
        for arg in ["george", "SNPS", "INDALS", 9]:
            self.assertRaises(ValueError, evf.Check_Type, arg)
            

class TestGet_Name(TestCase):
    def test_path(self):
        inpath = "made/up/path.vcf"
        outpath = "path.filter.vcf"
        self.assertEqual(outpath, evf.Get_Name(inpath))
    def test_path_gz(self):
        inpath = "made/up/path.vcf.gz"
        outpath = "path.filter.vcf"
        self.assertEqual(outpath, evf.Get_Name(inpath))
            