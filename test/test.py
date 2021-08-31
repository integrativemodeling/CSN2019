#!/usr/bin/env python

import unittest
import os
import shutil
import sys
import subprocess
import ihm.reader

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..'))


def make_fake_dcds():
    if os.path.exists('../Ensemble_DCD'):
        shutil.rmtree('../Ensemble_DCD')
    os.mkdir('../Ensemble_DCD')
    for sub in ('A', 'B'):
        for name in ('BMSO', 'BMSO_DHSO', 'CSN', 'CSNn', 'DHSO',
                     'DSSO_BMSO', 'DSSO', 'DSSO_DHSO'):
            with open(os.path.join('..', 'Ensemble_DCD',
                                   '%s_%s.dcd' % (sub, name)), 'w') as fh:
                pass


class Tests(unittest.TestCase):
    def run_modeller_script(self, script_name, model_name, resrng, ali_file):
        """Run a Modeller script and test the output model"""
        # Run script
        p = subprocess.check_call(["python", script_name, "012345",
                                   ali_file, "--test"])
        # Make sure PDB was produced with the requested residue range

        with open('%s.B99990000.pdb' % model_name) as fh:
            pdb_lines = [x for x in fh.readlines() if x.startswith('ATOM')]
        rng = (int(pdb_lines[0][22:26]), int(pdb_lines[-1][22:26]))
        self.assertEqual(rng, resrng)

    def test_model_CSN1(self):
        """Test generation of full model for CSN1 using Modeller"""
        os.chdir(os.path.join(TOPDIR, 'data/CSN1'))
        self.run_modeller_script('CSN1_Modeller.py',
                                 'CSN1', (1, 491), 'CSN1_CSN1.ali')

    def test_model_CSN7B(self):
        """Test generation of comparative model for CSN7B"""
        os.chdir(os.path.join(TOPDIR, 'data/CSN7B'))
        self.run_modeller_script('CSN7B_Modeller.py',
                                 'CSN7B', (1, 264), 'CSN7B_CSN7A.ali')

    def test_simple(self):
        """Test model building"""
        os.chdir(os.path.join(TOPDIR, 'modeling_scripts'))
        p = subprocess.check_call(["python", "smodeling_test.py", "--test"])
        # todo: assert outputs, run analysis

    def test_canonical_mmcif(self):
        """Test generation of canonical mmCIF output"""
        os.chdir(os.path.join(TOPDIR, 'results', 'modeling_scripts'))
        make_fake_dcds()
        if os.path.exists('CSN.cif'):
            os.unlink('CSN.cif')
        # Potentially override methods that need network access
        env = os.environ.copy()
        env['PYTHONPATH'] = os.path.join(TOPDIR, 'test', 'mock') \
                            + ':' + env.get('PYTHONPATH', '')
        p = subprocess.check_call(
            ["python", "create_mmcif_canonical.py"], env=env)
        # Check output file
        self._check_canonical_mmcif_file('CSN.cif')

    def test_non_canonical_mmcif(self):
        """Test generation of non-canonical mmCIF output"""
        os.chdir(os.path.join(TOPDIR, 'results', 'modeling_scripts'))
        make_fake_dcds()
        if os.path.exists('CSNn.cif'):
            os.unlink('CSNn.cif')
        # Potentially override methods that need network access
        env = os.environ.copy()
        env['PYTHONPATH'] = os.path.join(TOPDIR, 'test', 'mock') \
                            + ':' + env.get('PYTHONPATH', '')
        p = subprocess.check_call(
            ["python", "create_mmcif_non_canonical.py"], env=env)
        # Check output file
        self._check_non_canonical_mmcif_file('CSNn.cif')

    def _check_canonical_mmcif_file(self, fname):
        with open(fname) as fh:
            s, = ihm.reader.read(fh)
        self.assertEqual(len(s.citations), 4)
        self.assertEqual(s.citations[0].doi, '10.1073/pnas.1915542117')
        self.assertEqual(len(s.software), 3)
        self.assertEqual(len(s.orphan_starting_models), 8)
        # Should be 1 state
        self.assertEqual(len(s.state_groups), 1)
        state1, = s.state_groups[0]
        # Should be 7 models
        self.assertEqual(sum(len(x) for x in state1), 7)
        # Check # of spheres and atoms in each model
        m = state1[0][0]
        self.assertEqual(len(m._spheres), 2519)
        self.assertEqual(len(m._atoms), 0)
        # Check ensembles
        self.assertEqual([e.num_models for e in s.ensembles],
                         [54702, 132407, 98186, 87368, 243067, 312515, 357350])
        # Seven sets of crosslinks
        self.assertEqual([len(xl.experimental_cross_links)
                         for xl in s.restraints], [74, 141, 40, 31, 91, 79])
        self.assertEqual([len(xl.cross_links)
                         for xl in s.restraints], [74, 140, 40, 28, 91, 78])

    def _check_non_canonical_mmcif_file(self, fname):
        with open(fname) as fh:
            s, = ihm.reader.read(fh)
        self.assertEqual(len(s.citations), 4)
        self.assertEqual(s.citations[0].doi, '10.1073/pnas.1915542117')
        self.assertEqual(len(s.software), 3)
        self.assertEqual(len(s.orphan_starting_models), 8)
        # Should be 1 state
        self.assertEqual(len(s.state_groups), 1)
        state1, = s.state_groups[0]
        # Should be 1 model
        self.assertEqual(sum(len(x) for x in state1), 1)
        # Check # of spheres and atoms in each model
        m = state1[0][0]
        self.assertEqual(len(m._spheres), 2548)
        self.assertEqual(len(m._atoms), 0)
        # Check ensembles
        self.assertEqual([e.num_models for e in s.ensembles], [125750])
        # Seven sets of crosslinks
        self.assertEqual([len(xl.experimental_cross_links)
                         for xl in s.restraints], [86, 186, 75, 34, 107, 131])
        self.assertEqual([len(xl.cross_links)
                         for xl in s.restraints], [86, 183, 75, 33, 107, 131])


if __name__ == '__main__':
    unittest.main()
