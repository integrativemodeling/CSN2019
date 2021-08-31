######################################
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
models = {1: {'analysis_dir': '../Cop9_DSS_BMS_DHS_BioRep/prefilter_2/',
              'n_models':5250000,
              'n_cluster':54702,
              'nA':12886,
              'nB':16227,
              'model_precision':16,
              'dcd_file':'CSN',
              'description': 'BMSO+DHSO+DSSO',
              'xls':'DSSO_DHSO_BMSO'},
          2: {'analysis_dir': '../Cop9_BMS_DHS/prefilter/',
              'n_models':5250000,
              'n_cluster':132407,
              'nA':15909,
              'nB':12318,
              'model_precision':22,
              'dcd_file':'BMSO_DHSO',
              'description': 'BMSO+DHSO',
              'xls':'DHSO_BMSO'},
          3: {'analysis_dir': '../Cop9_DSS_DHS/prefilter/',
              'n_models':5250000,
              'n_cluster':98186,
              'nA':13663,
              'nB':15486,
              'model_precision':24,
              'dcd_file':'DSSO_DHSO',
              'description': 'DHSO+DSSO',
              'xls':'DSSO_DHSO'},
          4: {'analysis_dir': '../Cop9_DSS_BMS/prefilter/',
              'n_models':5250000,
              'n_cluster':87368,
              'nA':13902,
              'nB':15481,
              'model_precision':27,
              'dcd_file':'DSSO_BMSO',
              'description': 'BMSO+DSSO',
              'xls':'DSSO_BMSO'},
          5: {'analysis_dir': '../Cop9_DSS/prefilter/',
              'n_models':5250000,
              'n_cluster':243067,
              'nA':11137,
              'nB':9382,
              'model_precision':27,
              'dcd_file':'DSSO',
              'description': 'DSSO',
              'xls':'DSSO'},
          6: {'analysis_dir': '../Cop9_DHS/prefilter/',
              'n_models':5250000,
              'n_cluster':312515,
              'nA':14175,
              'nB':13364,
              'model_precision':29,
              'dcd_file':'DHSO',
              'description': 'DHSO',
              'xls':'DHSO'},
          7: {'analysis_dir': '../Cop9_BMS/prefilter/',
              'n_models':5250000,
              'n_cluster':357350,
              'nA':7793,
              'nB':14049,
              'model_precision':37,
              'dcd_file':'BMSO',
              'description': 'BMSO',
              'xls':'BMSO'}}


sys.argv = ['','--mmcif']

exec(open('script_canonical.py').read())
print('print', po.asym_units)

po.system.citations.append(ihm.Citation.from_pubmed_id(32034103))
#po.system.title = ('Genetic interaction mapping informs integrative determination of biomolecular assembly structures')

#add uniprot links

Uniprot={'CSN1.0':'Q13098',
         'CSN2.0':'P61201',
         'CSN3.0':'Q9UNS2',
         'CSN4.0':'Q9BT78',
         'CSN5.0':'Q92905',
         'CSN6.0':'Q7L5N1',
         'CSN7.0':'Q9H9Q2',
         'CSN8.0':'Q99627'}
         #'CSN9.0':'Q8WXC6'}
for prot, entry in Uniprot.items():
     ref = ihm.reference.UniProtSequence.from_accession(entry)
     po.asym_units[prot].entity.references.append(ref)

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

    loc = ihm.location.WorkflowFileLocation(path='script_canonical.py',
                                            details="Main modeling script for the human CSN complex")

    po.system.locations[0] = loc
    #po.system.locations.append(loc)
    
    e = po._add_simple_ensemble(analysis.steps[-1],
                                name="Cluster 0"+str(info['description']), num_models=info['n_cluster'],
                                drmsd=info['model_precision'], num_models_deposited=1,
                                localization_densities={},
                                ensemble_file=None)
    for sample in ('A', 'B'):
        dcd_location = ihm.location.OutputFileLocation(path='../Ensemble_DCD/%s_%s.dcd' % (sample, info['dcd_file']),
                                                       details="Coordinates of subsample %s structures in the largest cluster" % sample)
        ss = ihm.model.IndependentSubsample(
            name='Cluster 0 subsample %s' % sample,
            num_models=info['n%s'%sample],
            file=dcd_location)
        e.subsamples.append(ss)
                                                       
    
    # Add the model from RMF
    import tempfile
    import shutil

    tmpd = tempfile.mkdtemp()
    centroid = fix_rmf_file(
        '../Localization_Probability_Densities/IntegrativeStructure.CSN/Structure_'+str(info['xls'])+'/cluster_center_model.rmf3',
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
        fname = '../Localization_Probability_Densities/IntegrativeStructure.CSN/Structure_'+str(info['xls'])+'/'+str(info['xls'])+'.'+name+'.LPD.mrc'
        print('fname', fname)
        loc = ihm.location.OutputFileLocation(fname)
        den = ihm.model.LocalizationDensity(file=loc, asym_unit=po.asym_units[asym])
        # Add to ensemble
        e.densities.append(den)

    from ihm import cross_linkers

    for r in po.system.restraints:
        if hasattr(r, 'linker_type') and r.linker_type == 'DSS_Inter':
            r.linker_type = 'DSSO'
        elif hasattr(r, 'linker_type') and r.linker_type == 'DSS_Intra':
            r.linker_type = 'DSSO'
        elif hasattr(r, 'linker_type') and r.linker_type == 'DHS_Inter':
            r.linker_type = 'DHSO'
        elif hasattr(r, 'linker_type') and r.linker_type == 'DHS_Intra':
            r.linker_type = 'DHSO'
        elif hasattr(r, 'linker_type') and r.linker_type == 'BMS_Inter':
            r.linker_type = 'BMSO'
        elif hasattr(r, 'linker_type') and r.linker_type == 'BMS_Intra':
            r.linker_type = 'BMSO'

        elif hasattr(r, 'linker') and r.linker.auth_name == 'DSS_Inter':
            r.linker = cross_linkers.dsso
        elif hasattr(r, 'linker') and r.linker.auth_name == 'DSS_Intra':
             r.linker = cross_linkers.dsso
        elif hasattr(r, 'linker') and r.linker.auth_name == 'DHS_Inter':
             r.linker = cross_linkers.dhso
        elif hasattr(r, 'linker') and r.linker.auth_name == 'DHS_Intra':
             r.linker = cross_linkers.dhso
        elif hasattr(r, 'linker') and r.linker.auth_name == 'BMS_Inter':
             r.linker = cross_linkers.bmso
        elif hasattr(r, 'linker') and r.linker.auth_name == 'BMS_Intra':
             r.linker = cross_linkers.bmso


        # Point to repositories where files are deposited
        scriptrepos = ihm.location.Repository(
            doi="10.5281/zenodo.3827934",
            root=".",
            top_directory="modeling_scripts",
            url="https://zenodo.org/record/3827934/files/modeling_scripts.zip")
        datarepos = ihm.location.Repository(
            doi="10.5281/zenodo.3827934",
            root="../data",
            top_directory="data",
            url="https://zenodo.org/record/3827934/files/data.zip")

        densityrepos = ihm.location.Repository(
            doi="10.5281/zenodo.3827934",
            root='../Localization_Probability_Densities',
            top_directory="Localization_Probability_Densities",
            url="https://zenodo.org/record/3827934/files/Localization_Probability_Densities.zip")

        dcdrepos = ihm.location.Repository(
            doi="10.5281/zenodo.3827934",
            root="../Ensemble_DCD",
            top_directory="Ensemble_DCD",
            url="https://zenodo.org/record/3827934/files/Ensemble_DCD.zip")

        po.system.update_locations_in_repositories([scriptrepos, datarepos, densityrepos, dcdrepos])

po.flush()

