from modeller import *
from modeller.automodel import *
import sys, os

log.verbose()
env = environ(rand_seed=int(sys.argv[1]))

# directories for input atom files
env.io.atom_files_directory = ['.', '../']

class MyModel(automodel):
    def special_patches(self, aln):
        self.rename_segments(segment_ids=['A'],
                             renumber_residues=1)

        


a = MyModel(env, 
            alnfile = 'CSN1_CSN1.ali',
            knowns = 'CSNA',
            sequence = 'CSN1',
            assess_methods=(assess.DOPE, assess.GA341, assess.normalized_dope)
            )

a.starting_model= 0
a.ending_model  = 10

if '--test' in sys.argv: a.ending_model = 0
a.make()
