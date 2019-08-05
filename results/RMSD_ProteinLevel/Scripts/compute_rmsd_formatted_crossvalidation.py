import os
import sys

import IMP
import RMF

import IMP.algebra
import IMP.atom
import IMP.rmf
import IMP.pmi

import IMP.pmi.analysis
import IMP.pmi.output


import pandas as pd
import itertools
import random
import math
import numpy as np

import multiprocessing as mp


class RMSDCalculator(object):

    def __init__(self,
                 clustering_directory='../',
                 cluster='0',
                 align_to=[],
                 number_of_models=0):

        self.clustering_dir = clustering_directory
        self.cluster = cluster
        self.number_of_models = int(number_of_models)
        self.align_to = align_to
        self.selection = []

        self.manager = mp.Manager()
        self.lock = self.manager.Lock()

        self.all_coords = self.manager.dict()
        self.ali_coords = self.manager.dict()

        self.nproc = 16
        self.nprod = 4
        self.RMSD_pairs = self.manager.dict()
            
        rmfs_cluster = self._get_rmfs_cluster()
        if 0 < self.number_of_models:
            rmfs_cluster = random.sample(rmfs_cluster, self.number_of_models)

        self.structures = rmfs_cluster
        self.get_rbs_from_centroid()
        self._read_all_coords_RBs_ignacia(rmfs_cluster)

    def compute_RMSD_all_versus_centroid(self):
        rmf3 = self.clustering_dir + '/cluster.' + \
            str(self.cluster) + '/cluster_center_model.rmf3'
        
        particles1, reference1 = self._get_RBs_particle_coordinates(rmf3)
        #ntot = len(self.all_coords.values())
        ntot = len(self.structures)
        print(ntot)
        
        rmfs_cluster = list(zip(range(len(self.structures)), self.structures))
        ND = int(np.ceil(len(rmfs_cluster)/float(self.nproc)))

        rmfs_dict = {}
        for k in range(self.nproc - 1):
            rmfs_dict[k] = rmfs_cluster[(k*ND):(k*ND+ND)]
	rmfs_dict[self.nproc-1] = rmfs_cluster[((self.nproc-1)*ND):(len(rmfs_cluster))]

        if self.nproc > 1:

	    output = mp.Queue()
            processes = [mp.Process(target=self._compute_RMSD_mp,
                                    args=(particles1, reference1, rmfs_dict[x])) for x in range(self.nproc)]
            for p in processes:
                p.start()

            for p in processes:
                p.join()


                
        elif self.nproc == 1:
            self._compute_RMSD_mp(particles1, reference1, rmfs_dict[0])
            
            
        self.all_versus_centroid = 1
        self.write_RMSD_output()

    def _compute_RMSD_mp(self, particles1, reference1, rmfs_dict_sub):

        for i, structure in rmfs_dict_sub:
            particles2, reference2 = self._get_RBs_particle_coordinates(structure)
            self.compute_RMSD(particles1, reference1, particles2, reference2)


    def write_RMSD_output(self):
        file_matrix = 'RMSD_mean_pairs_all_versus_centroid.csv'
        file_dev_matrix = 'RMSD_dev_pairs_all_versus_centroid.csv'
        file_index = 'RBs_indexes_RMSD_calculation_all_versus_centroid.dat'
        file_values = 'RMSD_pairs_all_versus_centroid_'


        M = pd.DataFrame(0, columns = self.rb_components.keys(), index = self.rb_components.keys())
        for pair, values in  self.RMSD_pairs.items():
            #print(pair, values)
            M.loc[pair[0], pair[1]] = np.mean(values)
            M.loc[pair[1], pair[0]] = M.loc[pair[0], pair[1]]
            np.savetxt(file_values+str(pair[0])+'_'+str(pair[1])+'.csv', np.array(values))
        print(M)
              
        M.to_csv(file_matrix)

        out = open(file_index, 'w')
        for k, v in self.rb_components.items():
            out.write(str(k)+'\t'+str(v).strip('[]')+'\n')
        out.close()
    
    def _get_RBs_particle_coordinates(self, rmf3):

        model=IMP.Model()
        h = IMP.pmi.analysis.get_hiers_from_rmf(model,0,rmf3)[0]
        
        particles_dict = {}
        reference_dict = {}

        for rb, sel in self.rb_components.items():
            particles = []
            for comp in self.rb_components[rb]:
                if len(self.selection) == 0 or (comp[0] in self.selection):
                    r1 = int(comp[1].split('-')[0])
                    r2 = int(comp[1].split('-')[1])
                    sel = IMP.atom.Selection(h,
                                             molecule=comp[0],
                                             residue_indexes = range(r1,r2),
                                             resolution=1)
                    particles += sel.get_selected_particles()
               
            coords = [np.array(IMP.core.XYZ(p).get_coordinates()) for p in particles]
            particles_dict[rb] = coords
            
        particles = []
        for ref in self.align_to:
            sel = IMP.atom.Selection(h,
                                     molecule = ref[0],
                                     residue_indexes = range(int(ref[1]), int(ref[2])),
                                     resolution = 1)
            particles += sel.get_selected_particles()
        coords = [np.array(IMP.core.XYZ(p).get_coordinates()) for p in particles]
        reference_dict['ref'] = coords
        
        del model
        del h
        return particles_dict, reference_dict
    
    def _get_rmfs_cluster(self):

        '''
        Read identities of model in clusters and then
        get the names of all rmfs
        TODO: write unique class to read the 
        '''

        rmfs_dic = {}
        for line in open(self.clustering_dir+'/Identities_A.txt', 'r'):
            vals = line.split()
            vals = line.replace('/scratch/', '/scratch/Cop9/').split()
            rmfs_dic[vals[1]] = vals[0]
        for line in open(self.clustering_dir+'/Identities_B.txt', 'r'):
            vals = line.split()
            vals = line.replace('/scratch/', '/scratch/Cop9/').split()
            rmfs_dic[vals[1]] = vals[0]
            
        # Read rmfs in cluster.0
        rmfs = []
        for mod in open(self.clustering_dir+'/cluster.'+str(self.cluster)+'.sample_A.txt','r'):
            vals = mod.split()[0]
            try:
                rmfs.append(rmfs_dic[vals])
            except:
                print('Model missing: ', vals)
        for mod in open(self.clustering_dir+'/cluster.'+str(self.cluster)+'.sample_B.txt','r'):
            vals = mod.split()[0]
            try:
                rmfs.append(rmfs_dic[vals])
            except:
                print('Model missing: ', vals)
                
        return rmfs

    def _read_coords_RBs(self, numbered_rmfs, a):
        for i, file in numbered_rmfs:
            particles = self._get_RBs_particle_coordinates_0(file)
            reference = self._get_RBs_particle_coordinates_1(file)
            self.all_coords[i] = particles
            self.ali_coords[i] = reference

            
    def get_rbs_from_centroid(self):
        rmf3 = self.clustering_dir+'/cluster.'+str(self.cluster)+'/cluster_center_model.rmf3'
        self._read_representation(rmf3)

    def _read_representation(self, rmf3):
        m = IMP.Model()
        rh_ref = RMF.open_rmf_file_read_only(rmf3)
        h_ref = IMP.rmf.create_hierarchies(rh_ref, m)[0]
        IMP.rmf.load_frame(rh_ref, RMF.FrameID(0))

        self.components = []
        self.struct_components = {}
        self.flex_components = {}
        rb_components = {}

        n_struct, n_flex = 0, 0
        for state in h_ref.get_children():
            for component in state.get_children():
                if len(self.selection) == 0 or (component.get_name() in self.selection):
                    n_residues_prot = 0
                    self.components.append(component.get_name())
                    
                    leaves = IMP.atom.get_leaves(component)
                    residues = len(leaves)
                    residue_range = str(1)+'-'+str(residues)
                    # n_res = residues[-1] - residues[0]
                    rb_components[str(component.get_name())] = [(component.get_name() ,residue_range)]
                    
        print("Ilan: ", rb_components)
        # self.struc_coverage = 100.0* n_struct/(n_struct + n_flex)
        self.struc_coverage = 100.0
        self.rb_components = {i:v for i, (k,v) in enumerate(rb_components.items())}
        
        for i,j in itertools.combinations_with_replacement(self.rb_components.keys(),2):
            self.RMSD_pairs[(i,j)] = []
    
        print('Coverage', self.struc_coverage)
        print('Number of rigid-bodies', len(self.rb_components))

        del m

    def compute_RMSD(self, particles1, reference1, particles2, reference2):
        import pyRMSD.RMSDCalculator
        calculator_name = "QCP_SERIAL_CALCULATOR"
        
        if set(particles1.keys()) != set(particles2.keys()):
            raise ValueError('ERROR: rmfs have different RBs')
        
        else:

            for rb1, rb2 in itertools.combinations_with_replacement(particles1.keys(), 2):
                cc = np.array([particles1[rb1] + particles1[rb2],
                               particles2[rb1] + particles2[rb2]])
            
                cr = np.array([reference1['ref'],
                               reference2['ref']])

                calculator = pyRMSD.RMSDCalculator.RMSDCalculator(calculator_name,
                                                                  fittingCoordsets=cr,
                                                                  calculationCoordsets=cc)
                
                rmsd = calculator.pairwise(0, 1, get_superposed_coordinates = False)
                with self.lock:
                    self.RMSD_pairs[(rb1, rb2)] = self.RMSD_pairs[(rb1, rb2)] + [rmsd]


    def _read_all_coords_RBs_ignacia(self, rmfs_cluster):
        rmfs_cluster = list(zip(range(len(rmfs_cluster)), rmfs_cluster))
        ND = int(np.ceil(len(rmfs_cluster)/float(self.nproc)))
        rmfs_dict = {}

        for k in range(self.nproc-1):
            rmfs_dict[k] = rmfs_cluster[(k*ND):(k*ND+ND)]
        rmfs_dict[self.nproc-1] = rmfs_cluster[((self.nproc-1)*ND):(len(rmfs_cluster))]
        output = mp.Queue()
        processes = [mp.Process(target=self._read_coords_RBs,
                                args=(rmfs_dict[x], 0)) for x in range(self.nproc)]
        for p in processes:
            p.start()
        for p in processes:
            p.join()
            
    def _get_RBs_particle_coordinates_0(self, rmf3):
        model=IMP.Model()
        h = IMP.pmi.analysis.get_hiers_from_rmf(model,0,rmf3)[0]

        particles_dict = {}

        for rb, sel in self.rb_components.items():
            particles = []
            for comp in self.rb_components[rb]:
                if len(self.selection) == 0 or (comp[0] in self.selection):
                    r1 = int(comp[1].split('-')[0])
                    r2 = int(comp[1].split('-')[1])
                    sel = IMP.atom.Selection(h,
                                             molecule=comp[0],
                                             residue_indexes = range(r1,r2),
                                             resolution=1)
                    particles += sel.get_selected_particles()

            coords = [np.array(IMP.core.XYZ(p).get_coordinates()) for p in particles]
            particles_dict[rb] = coords
        del model
        return particles_dict

    def _get_RBs_particle_coordinates_1(self, rmf3):
        model=IMP.Model()
        h = IMP.pmi.analysis.get_hiers_from_rmf(model,0,rmf3)[0]

        reference_dict = {}
        particles = []
        for ref in self.align_to:
            sel = IMP.atom.Selection(h,
                                     molecule = ref[0],
                                     residue_indexes = range(int(ref[1]), int(ref[2])),
                                     resolution = 1)
            particles += sel.get_selected_particles()
        coords = [np.array(IMP.core.XYZ(p).get_coordinates()) for p in particles]
        reference_dict['ref'] = coords

        del model
        return reference_dict

    def compute_RMSD_all_versus_centroid_ignacia(self):
        rmf3 = self.clustering_dir + '/cluster.' + \
            str(self.cluster) + '/cluster_center_model.rmf3'

        particles1 = self._get_RBs_particle_coordinates_0(rmf3)
        reference1 = self._get_RBs_particle_coordinates_1(rmf3)

        ntot = len(self.all_coords.values())
        print(ntot)

        number = range(ntot)
        ND = int(np.ceil(len(number)/float(self.nprod)))
        
        output = mp.Queue()
        processes = [mp.Process(target=self._compute_RMSD_2,
                                args=(particles1, reference1, number[ND * x : ND * x + ND])) for x in range(self.nprod)]
        for p in processes:
            p.start()
        for p in processes:
            p.join()

        '''
        for i, coordinates_all in enumerate(self.all_coords.values()):
            particles2 = coordinates_all
            reference2 = self.ali_coords.values()[i]
            self.compute_RMSD(particles1, reference1, particles2, reference2)
            if i % 100 == 0:
                print(str(i) + ' of ' + str(ntot))
        '''
        self.all_versus_centroid = 1
        self.write_RMSD_output()

    def _compute_RMSD_2(self, particles1, reference1, number):

        import pyRMSD.RMSDCalculator
        os.system("taskset -p 0xfffff %d" % os.getpid())
        
        calculator_name = "QCP_OMP_CALCULATOR"
        
        
        for i in number:
            print('Ilan: ', i)
            reference2 = self.ali_coords.values()[i]
            particles2 = self.all_coords.values()[i]
            
            if set(particles1.keys()) != set(particles2.keys()):
                raise ValueError('ERROR: rmfs have different RBs')

            else:

                for rb1, rb2 in itertools.combinations_with_replacement(particles1.keys(), 2):
                    cc = np.array([particles1[rb1] + particles1[rb2],
                                   particles2[rb1] + particles2[rb2]])

                    cr = np.array([reference1['ref'],
                                   reference2['ref']])

                    calculator = pyRMSD.RMSDCalculator.RMSDCalculator(calculator_name,
                                                                      fittingCoordsets=cr,
                                                                      calculationCoordsets=cc)
                    calculator.setNumberOfOpenMPThreads(int(self.nprod))
                    
                    rmsd = calculator.pairwise(0,
                                               1,
                                               get_superposed_coordinates = False)
                    with self.lock:
                        self.RMSD_pairs[(rb1, rb2)] =  self.RMSD_pairs[(rb1, rb2)] + [rmsd]
        
alignment = [('CSN1',431,461),
             ('CSN2',417,443),
             ('CSN3',368,401),
             ('CSN4',365,406),
             ('CSN5',296,333),
             ('CSN6',271,316),
             ('CSN7',163,212),
             ('CSN8',194,209)]

dirs = ['BMS', 'DHS', 'DSS', 'BMS_DHS', 'DSS_DHS', 'DSS_BMS']
for dir in dirs:
    os.mkdir('dRMS_CSNc_vs_%s' % dir)
    os.chdir('dRMS_CSNc_vs_%s' % dir) 

    D = RMSDCalculator('/scratch/Cop9/Cop9_%s/prefilter/' % dir,
                       0,
                       align_to = alignment,
                       number_of_models = 7500)

    D.compute_RMSD_all_versus_centroid()
    os.chdir('..')

