import IMP
import RMF
import IMP.atom
import IMP.core
import IMP.algebra
import IMP.container
import IMP.rmf
import os
import sys


import DLFCN as dl
sys.setdlopenflags(dl.RTLD_NOW | dl.RTLD_GLOBAL)

import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking

import IMP.pmi.mmcif

#---------------------------
# Define Input Files
#---------------------------
datadirectory = "../../../data/"
topology_file = datadirectory + "topology_free.txt"

#--------------------------
# Set MC Sampling Parameters
#--------------------------
num_frames = 2
if '--test' in sys.argv:
    num_frames = 20
num_mc_steps = 10
#--------------------------
# Create movers
#--------------------------

# rigid body movement params
rb_max_trans = 1.00
rb_max_rot = 0.01
# flexible bead movement
bead_max_trans = 2.00

#--------------------------------
# Build the Model Representation
#--------------------------------

# Initialize model
m = IMP.Model()

# Create list of components from topology file
topology = IMP.pmi.topology.TopologyReader(topology_file,
                                           pdb_dir=datadirectory,
                                           fasta_dir=datadirectory
                                           )
domains = topology.get_components()

bs = IMP.pmi.macros.BuildSystem(m)

#if '--mmcif' in sys.argv:
# Record the modeling protocol to an mmCIF file
po = IMP.pmi.mmcif.ProtocolOutput(open('CSN.cif', 'w'))
po.system.title = ('Integrative structure of the human CSN complex')
print(po.system.title)
bs.system.add_protocol_output(po)
    
    
# Add publication
#po.system.citations.append(ihm.Citation.from_pubmed_id(11111111))

bs.add_state(topology)
representation, dof = bs.execute_macro(max_rb_trans=rb_max_trans,
                                       max_rb_rot=rb_max_rot,
                                       max_bead_trans=bead_max_trans)


# representation.shuffle_configuration(50)

#--------------------------
# Define Degrees of Freedom
#--------------------------

# Add default mover parameters to simulation
outputobjects = []  # reporter objects (for stat files)
sampleobjects = []  # sampling objects

#-----------------------------------
# Define Scoring Function Components
#-----------------------------------

# Here we are defining a number of restraints on our system.
#  For all of them we call add_to_model() so they are incorporated into scoring
# We also add them to the outputobjects list, so they are reported in stat
# files

# Excluded Volume Restraint
#  To speed up this expensive restraint, we operate it at resolution 20

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "ilan0", sf.evaluate(False)

CSN = []
crs = []

moldict = bs.get_molecules()[0]
for molname in moldict:
    print molname
    for mol in moldict[molname]:

        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(
            mol, scale=2.0, label=molname)
        cr.add_to_model()
        outputobjects.append(cr)
        sampleobjects.append(cr)
        crs.append(cr)

        CSN.append(mol)


print(CSN)
sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "ilan0", sf.evaluate(False)
ev1 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=CSN,
                                                              resolution=10)
ev1.add_to_model()
ev1.set_label('CSN')
outputobjects.append(ev1)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "ilan1", sf.evaluate(False)


# Crosslinks - dataset
#  To use this restraint we have to first define the data format
#  Here assuming that it's a CSV file with column names that may need to change
# Other options include the linker length and the slope (for nudging
# components together)
kw = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
kw.set_unique_id_key("id")
kw.set_protein1_key("prot1")
kw.set_protein2_key("prot2")
kw.set_residue1_key("res1")
kw.set_residue2_key("res2")
kw.set_id_score_key(None)

#############
#############
#####DSS#####
#############
#############

# Medium + High Confidence Intermolecular crosslinks
DSS1 = IMP.pmi.io.crosslink.CrossLinkDataBase(kw)
DSS1.create_set_from_file(datadirectory + 'DSS.Inter.csv')

x_dss1 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=representation,
    CrossLinkDataBase=DSS1,
    length=21,
    label="DSS_Inter",
    resolution=1.0,
    slope=0.02)

x_dss1.rs.set_weight(3.0)
x_dss1.add_to_model()
sampleobjects.append(x_dss1)
outputobjects.append(x_dss1)
dof.get_nuisances_from_restraint(x_dss1)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "ilan3", sf.evaluate(False)

# Medium + High Confidence Intramolecular crosslinks
DSS2 = IMP.pmi.io.crosslink.CrossLinkDataBase(kw)
DSS2.create_set_from_file(datadirectory + 'DSS.Intra.csv')

x_dss2 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=representation,
    CrossLinkDataBase=DSS2,
    length=21,
    label="DSS_Intra",
    resolution=1.0,
    slope=0.02)

x_dss2.rs.set_weight(1.00)
x_dss2.add_to_model()
sampleobjects.append(x_dss2)
outputobjects.append(x_dss2)
dof.get_nuisances_from_restraint(x_dss2)


#############
#############
#####BMS#####
#############
#############

# Medium + High Confidence Intermolecular crosslinks
BMS1 = IMP.pmi.io.crosslink.CrossLinkDataBase(kw)
BMS1.create_set_from_file(datadirectory + 'BMS.Inter.csv')

x_bms1 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=representation,
    CrossLinkDataBase=BMS1,
    length=29,
    label="BMS_Inter",
    resolution=1.0,
    slope=0.02)
x_bms1.rs.set_weight(3.0)
x_bms1.add_to_model()
sampleobjects.append(x_bms1)
sampleobjects.append(x_bms1)
outputobjects.append(x_bms1)
dof.get_nuisances_from_restraint(x_bms1)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "ilan3", sf.evaluate(False)

# Medium + High Confidence Intramolecular crosslinks
BMS2 = IMP.pmi.io.crosslink.CrossLinkDataBase(kw)
BMS2.create_set_from_file(datadirectory + 'BMS.Intra.csv')

x_bms2 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=representation,
    CrossLinkDataBase=BMS2,
    length=29,
    label="BMS_Intra",
    resolution=1.0,
    slope=0.02)

x_bms2.rs.set_weight(1.00)
x_bms2.add_to_model()
sampleobjects.append(x_bms2)
outputobjects.append(x_bms2)
dof.get_nuisances_from_restraint(x_bms2)

Sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "ilan3", sf.evaluate(False)

#############
#############
#####DHS#####
#############
#############

# Medium + High Confidence Intermolecular crosslinks
DHS1 = IMP.pmi.io.crosslink.CrossLinkDataBase(kw)
DHS1.create_set_from_file(datadirectory + 'DHS.Inter.csv')

x_dhs1 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=representation,
    CrossLinkDataBase=DHS1,
    length=21,
    label="DHS_Inter",
    resolution=1.0,
    slope=0.02)
x_dhs1.rs.set_weight(3.0)
x_dhs1.add_to_model()
sampleobjects.append(x_dhs1)
sampleobjects.append(x_dhs1)
outputobjects.append(x_dhs1)
dof.get_nuisances_from_restraint(x_dhs1)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "ilan3", sf.evaluate(False)

# Medium + High Confidence Intramolecular crosslinks
DHS2 = IMP.pmi.io.crosslink.CrossLinkDataBase(kw)
DHS2.create_set_from_file(datadirectory + 'DHS.Intra.csv')

x_dhs2 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=representation,
    CrossLinkDataBase=DHS2,
    length=21,
    label="DHS_Intra",
    resolution=1.0,
    slope=0.02)

x_dhs2.rs.set_weight(1.00)
x_dhs2.add_to_model()
sampleobjects.append(x_dhs2)
outputobjects.append(x_dhs2)
dof.get_nuisances_from_restraint(x_dhs2)

Sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "ilan3", sf.evaluate(False)


###Shuffling for random initial conformations 
IMP.pmi.tools.shuffle_configuration(CSN)
dof.optimize_flexible_beads(200)


#--------------------------
# Monte-Carlo Sampling
#--------------------------
# This object defines all components to be sampled as well as the sampling
# protocol
mc1 = IMP.pmi.macros.ReplicaExchange0(m,
                                root_hier=representation,
                                monte_carlo_sample_objects=dof.get_movers(),
                                output_objects=outputobjects,
                                crosslink_restraints=sampleobjects,
                                monte_carlo_temperature=1.0,

                                simulated_annealing=True,
                                simulated_annealing_minimum_temperature=1.0,
                                simulated_annealing_maximum_temperature=1.5,
                                simulated_annealing_minimum_temperature_nframes=200,
                                simulated_annealing_maximum_temperature_nframes=20,

                                replica_exchange_minimum_temperature=1.5,
                                replica_exchange_maximum_temperature=2.5,

                                number_of_best_scoring_models=1,
                                monte_carlo_steps=num_mc_steps,
                                number_of_frames=num_frames,
                                global_output_directory="output")

mc1.execute_macro()

def fix_rmf_file(original_rmf, molnames, tmpd):
    # Put the molecules in original_rmf in the given order, and return the
    # fixed RMF file name
    m = IMP.Model()
    rh = RMF.open_rmf_file_read_only(original_rmf)
    h = IMP.rmf.create_hierarchies(rh, m)
    state, = h[0].get_children()
    children = state.get_children()

    names = {}
    for c in children:
        state.remove_child(0)
        names[c.get_name()] = c

    for mn in molnames:
        state.add_child(names[mn])

    rmf_file = os.path.join(tmpd, 'new.rmf')
    rh = RMF.create_rmf_file(rmf_file)
    IMP.rmf.add_hierarchies(rh, h)
    IMP.rmf.save_frame(rh)
    del rh
    return rmf_file



import ihm.dumper
import ihm.format
import ihm.location
import ihm.representation
import ihm.startmodel
import ihm.dataset
import ihm.protocol
import ihm.analysis
import ihm.model
import ihm.restraint
import ihm.geometry
import ihm.cross_linkers

import tempfile
import shutil

print(po.system.title)

# Correct number of output models to account for multiple runs
protocol = po.system.orphan_protocols[-1]
protocol.steps[-1].num_models_end = 5250000

# Get last protocol in the file
protocol = po.system.orphan_protocols[-1]
