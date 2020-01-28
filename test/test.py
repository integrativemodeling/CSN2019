#!/usr/bin/env python

import unittest
import os
import shutil
import sys
import subprocess
import ihm.reader

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..'))

class Tests(unittest.TestCase):
    def run_modeller_script(self, script_name, model_name, resrng):
        """Run a Modeller script and test the output model"""
        os.chdir(os.path.join(TOPDIR, 'data/'))
        # Run script
        p = subprocess.check_call(["python", script_name, "--test"])
        # Make sure PDB was produced with the requested residue range
        with open('%s.B99990001.pdb' % model_name) as fh:
            pdb_lines = [x for x in fh.readlines() if x.startswith('ATOM')]
        rng = (int(pdb_lines[0][22:26]), int(pdb_lines[-1][22:26]))
        self.assertEqual(rng, resrng)
        
        
    def test_model_CSN1(self):
        """Test generation of full model for CSN1 using Modeller"""
        self.run_modeller_script('./CSN1/CSN1_Modeller.py',
                                 'CSN1', (1, 491))

    def test_model_CSN7B(self):
        """Test generation of comparative model for CSN7B"""
        self.run_modeller_script('./CSN7B/CSN7B_Modeller.py',
                                 'CSN7B', (1, 264))

    def test_simple(self):
        """Test model building"""
        os.chdir(os.path.join(TOPDIR, 'modeling_scripts'))
        p = subprocess.check_call(["python", "smodeling_test.py", "--test"])
        # todo: assert outputs, run analysis


if __name__ == '__main__':
    unittest.main()
