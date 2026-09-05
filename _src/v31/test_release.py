"""Small adversarial tests for preservation checker; fixtures stay in experiment."""
import importlib.util
from pathlib import Path
import unittest
import uuid

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem, write_elem

spec = importlib.util.spec_from_file_location('release', Path(__file__).with_name('build_release.py'))
release = importlib.util.module_from_spec(spec)
spec.loader.exec_module(release)


class Checks(unittest.TestCase):
    def setUp(self):
        self.root = release.EXP / 'audit' / 'test_fixtures' / uuid.uuid4().hex
        self.root.mkdir(parents=True)
        for name in ['a.h5', 'b.h5']:
            with h5py.File(self.root / name, 'w') as f:
                f.attrs['encoding'] = 'fixture'
                f.create_dataset('matrix', data=np.array([[1., np.nan], [0., 2.]], dtype='float32'))

    def check(self):
        with h5py.File(self.root / 'a.h5') as a, h5py.File(self.root / 'b.h5') as b:
            rows = []
            release.compare_original(a, b, rows)
            self.assertEqual(len(rows), 1)

    def test_equal_nan_preserved(self):
        self.check()

    def test_changed_value_rejected(self):
        with h5py.File(self.root / 'b.h5', 'r+') as f:
            f['matrix'][1, 1] = 3
        with self.assertRaises(AssertionError):
            self.check()

    def test_extra_original_structure_rejected(self):
        with h5py.File(self.root / 'b.h5', 'r+') as f:
            f['extra'] = 1
        with self.assertRaises(AssertionError):
            self.check()

    def test_attr_change_rejected(self):
        with h5py.File(self.root / 'b.h5', 'r+') as f:
            f.attrs['encoding'] = 'changed'
        with self.assertRaises(AssertionError):
            self.check()

    def test_empty_and_json_category_roundtrip(self):
        with h5py.File(self.root / 'cat.h5', 'w') as f:
            values = ['', 'Human', '[]', '["GOLD"]']
            write_elem(f, 'field', pd.Categorical(values))
            self.assertEqual(list(read_elem(f['field'])), values)


if __name__ == '__main__':
    unittest.main()
