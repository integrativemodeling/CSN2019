#######################################
#                                     #
#       Creating mmCIF file           #
#       For the CSN complexe          #
#          Ilan E Chemmama            #
#       ilan - Salilab - UCSF         #
#                                     #
# Modified from the script written by #
#       Ignacia Echeverria            #
#   iecheverria - Salilab - UCSF      #
#       ignacia@salilab.org           #
#                                     #
#######################################

import IMP
import RMF
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.mmcif

import numpy as np
import sys

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




########################
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


########################
models = {1: {'analysis_dir': '../',
              'n_models':7500000,
              'n_cluster':125750,
              'model_precision':22,
              'dcd_file':'',
              'description': 'BMSO+DHSO+DSSO'}}
sys.argv = ['','--mmcif']

exec(open('scripts_canonical.py').read())
print('print', po.asym_units)


#po.system.title = ('Genetic interaction mapping informs integrative determination of biomolecular assembly structures')
for m, info in models.items():
    #bs.system.add_protocol_output(po)

    #rg = ihm.restraint.RestraintGroup(all_rest)
    #po.system.restraint_groups.append(rg)

    # Correct number of output models to account for multiple runs
    print(po.system)
    protocol = po.system.orphan_protocols[-1]
    protocol.steps[-1].num_models_end = info['n_models']

    # Get last protocol in the file
    protocol = po.system.orphan_protocols[-1]
    # State that we filtered the 200000 frames down to one cluster of
    # N models:
    analysis = ihm.analysis.Analysis()
    protocol.analyses.append(analysis)
    analysis.steps.append(ihm.analysis.ClusterStep(
        feature='RMSD', num_models_begin=info['n_models'],
        num_models_end=info['n_cluster']))


    e = po._add_simple_ensemble(analysis.steps[-1],
                                name="Ensemble generated using "+str(info['description'])+" datasets ",
                                num_models=info['n_cluster'],
                                drmsd=info['model_precision'], num_models_deposited=1,
                                localization_densities={}, ensemble_file=None)
    
    # Add the model from RMF
    import tempfile
    import shutil

    tmpd = tempfile.mkdtemp()
    centroid = fix_rmf_file(info["analysis_dir"]+'/cluster_center_model.rmf3',
                            moldict,
                            tmpd)

    rh = RMF.open_rmf_file_read_only(centroid)
    IMP.rmf.link_hierarchies(rh, [representation])
    IMP.rmf.load_frame(rh, RMF.FrameID(0))
    del rh
    model = po.add_model(e.model_group)

    # Add localization densities
    # Look up the ihm.AsymUnit corresponding to a PMI component name
    for asym in po.asym_units:
        name = asym.split('.')[0]
        fname = str(info["analysis_dir"])+'/LPD_'+name+'.mrc'
        print('fname', fname)
        loc = ihm.location.OutputFileLocation(fname)
        den = ihm.model.LocalizationDensity(file=loc, asym_unit=po.asym_units[asym])
        # Add to ensemble
        e.densities.append(den)

    
po.flush()
    
    
