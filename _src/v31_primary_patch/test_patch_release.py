"""Adversarial policy and preservation tests; no corpus mutations."""
import unittest
import uuid
import h5py
import numpy as np
import pandas as pd
from anndata.io import write_elem
import patch_release as p

class Tests(unittest.TestCase):
    def test_policy(self):
        d=pd.DataFrame([
            ['RM','A',False,9000],['MA','A',True,1000],
            ['MA','B',True,9000],['RM','B',True,1000],['RM','B',True,2000],
            ['MA','C',True,1000],['Qiita','C',True,2000],['MBM','C',True,2000],
            ['Qiita','',True,1000],['Qiita','',True,1000],['Qiita','',False,1000],
            ['MA','D',False,9999]],columns=['Database','Run','qc_pass','n_deep_reads'])
        expected=np.zeros(len(d),dtype=bool);expected[[1,4,6,8,9]]=True
        self.assertTrue(np.array_equal(p.select(d),expected))
        self.assertTrue(np.array_equal(p.independent_select(d),expected))

    def setUp(self):
        self.root=p.EXP/'audit/patch_test_fixtures'/uuid.uuid4().hex;self.root.mkdir(parents=True)
        for name in ['a','b']:
            with h5py.File(self.root/name,'w') as f:
                f.create_dataset('matrix',data=[1.,np.nan,2.])
                f.create_group('obs').attrs['column-order']=['qc_pass']
                write_elem(f['obs'],'qc_pass',np.array([True,False,True]))
                f.create_group('uns')
        with h5py.File(self.root/'b','r+') as f:
            write_elem(f['obs'],'analysis_primary',np.array([True,False,True]))
            f['obs'].attrs['column-order']=['qc_pass','analysis_primary']
            write_elem(f['uns'],'analysis_primary_policy','{}')

    def check(self):
        with h5py.File(self.root/'a') as a,h5py.File(self.root/'b') as b: p.preserve(a,b,[])

    def test_allowed_addition_nan(self): self.check()
    def test_matrix_change_rejected(self):
        with h5py.File(self.root/'b','r+') as f:f['matrix'][0]=7
        with self.assertRaises(AssertionError):self.check()
    def test_qc_change_rejected(self):
        with h5py.File(self.root/'b','r+') as f:f['obs/qc_pass'][0]=False
        with self.assertRaises(AssertionError):self.check()
    def test_extra_field_rejected(self):
        with h5py.File(self.root/'b','r+') as f:f['obs/other']=3
        with self.assertRaises(AssertionError):self.check()
    def test_order_change_rejected(self):
        with h5py.File(self.root/'b','r+') as f:f['obs'].attrs['column-order']=['analysis_primary','qc_pass']
        with self.assertRaises(AssertionError):self.check()

if __name__=='__main__':unittest.main()
