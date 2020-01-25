from modeller import *
from modeller.automodel import *
import sys, os

log.verbose()
env = environ(rand_seed=int(sys.argv[1]))

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

class MyModel(automodel):
    def special_patches(self, aln):
        self.rename_segments(segment_ids=['I'],
                             renumber_residues=1)

        


a = MyModel(env, 
            alnfile = 'CSN7B_CSN7A.ali',
            knowns = 'CSN7A',
            sequence = 'CSN7B',
            assess_methods=(assess.DOPE, assess.GA341, assess.normalized_dope)
            )

a.starting_model= 0
a.ending_model  = 500
a.make()
