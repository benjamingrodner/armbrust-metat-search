"""
Filename: cl.tree_trim.py
Author: Benjamin Grodner (https://github.com/benjamingrodner)
Date: 2025-02-05
Description:
    Given contigs mapped to gene annotations, taxonomic annotations, and sample read counts,
    trim the taxonomic tree such that each node in the tree fulfils some threshold
    such as number of contigs, number of genes, presence across multiple samples, etc.

    The input is a csv file in tidy format (https://tidyr.tidyverse.org/articles/tidy-data.html).

    Save the tree and the new mapping of contig to taxon.
"""
###################################################################################################
# Imports
###################################################################################################

import fn_metat_files as fnf
from collections import defaultdict
from ete4 import NCBITaxa
from time import time
import numpy as np



###################################################################################################
# Object
###################################################################################################


        

class TreeTrim:
    """
    Class to manges functions for trimming trees based on thresholding at each node.      
    """
    
    def __init__(
            self, 
            fn_tidycsv, 
            col_contig='contig', 
            col_ko='KO',
            col_taxon='taxon',
            col_estcounts='estcounts',
            col_sample='fn_sample_counts'
            ):
        """
        Initializes a TreeTrim instance

        Params:

        
        Returns:
            None
        """
        # Params
        self.fn_tidycsv = fn_tidycsv
        self.col_contig = col_contig
        self.col_ko = col_ko
        self.col_taxon = col_taxon
        self.col_estcounts = col_estcounts
        self.col_sample = col_sample

        # Builtin vars
        self.ncbi = NCBITaxa()  # ete4 access to ncbi database
        self.dict_name_filtfunc = {
            "mean_nkos": self.get_mean_nkos,
            "nkos_in_gt_minsamples": self.get_nkos_in_gt_minsamples,
        }  # Functions used to filter tree

        # Execute functions
        print('Loading dicts...')
        fnf.getmem()
        t0 = time()
        self._get_dicts()  # Load from tidytable
        t1 = time()
        el = (t1 - t0)
        elm = el // 60
        els = el % 60
        print(f"Time to load dicts: {elm} min {els} sec")
        fnf.getmem()

        self._get_set_taxids()  # All taxids present
        self.tree = self.ncbi.get_topology(self.set_taxids)  # Build tree
        self._get_dict_taxon_contigs()  # invert taxon contig dict
        self.nsamples = len(next(iter(self.dict_contig[self.col_estcounts].values())))  # number of samples


    def _get_dicts(self):
        """
        Convert the tidytable csv into nested dictionaries
        """
        # Load the csv file as a dict mapping contigs to other values
        col_key = self.col_contig
        cols_key2value = [self.col_ko, self.col_taxon]
        cols_key2list = [self.col_estcounts, self.col_sample]
        self.dict_contig = fnf.tidytable_to_dict(
            fn = self.fn_tidycsv, 
            col_key=col_key,
            cols_key2value=cols_key2value,
            cols_key2list=cols_key2list, 
            floats=[self.col_estcounts]
            )
        # # Convert taxon and ko dicts to independent dicts
        # self.dict_contig_taxon = {}
        # self.dict_contig_ko = {}
        # self.dict_contig_estcounts = {}
        # self.dict_contig_samples = {}
        # for c, d in self.dict_contig.items():
        #     self.dict_contig_taxon[c] = d[self.col_taxon]
        #     self.dict_contig_ko[c] = d[self.col_ko]
        #     self.dict_contig_estcounts[c] = d[self.col_estcounts]
        #     self.dict_contig_samples[c] = d[self.col_sample]

        

    def _get_set_taxids(self):
        """
        Get a set of taxids from the input dictionary that are present in the ncbi taxonomy database.
        """
        taxids_all = []
        for t in self.dict_contig[self.col_taxon].values():
            t = t
            if int(t):
                try:
                    self.ncbi.get_lineage(t)
                    taxids_all.append(t)
                except:
                    pass
        self.set_taxids = set(taxids_all)


    def _get_taxon_contigs_full_lineage(self):
        """
        Get the full lineage for each contig and then build a dictionary mapping
        taxon to a list of all contigs within that taxonomy.  
        """
        self.dict_taxfull_contigs = defaultdict(list)
        for c, t in self.dict_contig[self.col_taxon].items():
            if t in self.set_taxids:
                lin = self.ncbi.get_lineage(t[0])
                for l in lin:
                    self.dict_taxfull_contigs[str(l)].append(c)


    def _get_dict_taxon_contigs(self):
        """
        Invert the key and value for dict_contig_taxon and build a list of contigs mapped at each
        taxonomic level (as long as it's in the ncbi database).
        """
        self.dict_taxon_contigs = defaultdict(list)
        for c, tid in self.dict_contig[self.col_taxon].items():
            if tid in self.set_taxids:
                self.dict_taxon_contigs[tid].append(c)


    # filter value calculation
    def get_mean_nkos(self, contigs):
        """
        Given a set of contigs within a taxon, count the number of KOs with nonzero read counts
        in each sample and take the mean.

        Params:
            contigs : list[str]
                List of contig ids

        Returns:
            float
                Mean number of KOs across samples
        """
        dict_sam_kos = defaultdict(set)
        for c in contigs:
            ko = self.dict_contig[self.col_ko][c]
            estcounts = self.dict_contig[self.col_estcounts][c]
            for sam, ec in enumerate(estcounts):
                if ec:
                    dict_sam_kos[sam].add(ko)
        nkos = []
        for _, kos in dict_sam_kos.items():
            _nkos = len(kos) if kos else 0
            nkos.append(_nkos)
        return np.mean(nkos) if nkos else 0
    

    def get_nkos_in_gt_minsamples(self, contigs, minsamples):
        """
        Given a set of contigs within a taxon, count the number of genes with nonzero read counts
        in more than a minimum number of samples. This value then used to 

        Params:
            contigs : list[str]
                List of contig ids
            minsamples : int
                Minimum number of samples a gene has to be present in to be counted

        Returns:
            float
                Number of genes in greater than the minimum number of samples
        """

        dict_ko_sam_estcounts = defaultdict(lambda: 
            defaultdict(float)
        )
        for c in contigs:
            ko = self.dict_contig[self.col_ko][c]
            sams = self.dict_contig[self.col_sample][c]
            estcounts = self.dict_contig[self.col_estcounts][c]
            for s, ec in zip(sams, estcounts):
                dict_ko_sam_estcounts[ko][s] += float(ec)
        nkos = 0
        for ko, dict_sam_estcounts in dict_ko_sam_estcounts.items():
            ko_in_nsam = 0
            for _, ec in dict_sam_estcounts.items():
                ko_in_nsam += bool(ec)
            if ko_in_nsam >= minsamples:
                nkos += 1
        return nkos            


    def trim_tree(self, filt_func_name, thresh, **fncargs):
        """
        Trim the tree such that the contigs assigned to each node have some measured value passing
        a threshold. 

        Produces self.treetrim and self.dict_taxtrim_contigs

        Params: 
            filt_func_name : str
                Specify which function to use to measure the contigs at a node. Current options are
                'mean_nkos', 'nkos_in_gt_minsamples'.
            thresh : float
                Threshold value used to decide whether to keep or trim a node
            **fnargs : various
                Keyword arguments to pass to the function specified in filt_func_name
        
        Returns:
            None
        """
        # Remove tree levels below filter value
        dict_taxtrim_contigs = defaultdict(list)
        # Traverse from leaves to root
        for n in self.tree.traverse(strategy='postorder'):
            tid = n.name
            # Get any new contigs from children
            if tid in dict_taxtrim_contigs:
                contigs = dict_taxtrim_contigs[tid]
            # Otherwise get old contigs for node
            else:
                contigs = self.dict_taxon_contigs[tid]
            # skip if root node
            if not n.is_root:
                # if the mean nkos is less than desired
                func = self.dict_name_filtfunc[filt_func_name]
                filtvalue = func(contigs, **fncargs)
                if filtvalue < thresh:
                    # nkos = len(set([self.dict_contig[self.col_ko][c] for c in contigs]))
                    # print(f'{tid} trimmed...with {len(contigs)} contigs and {nkos} kos')
                    # Then add these contigs to the parent taxon id
                    tid_parent = next(n.ancestors()).name
                    if tid_parent in dict_taxtrim_contigs:
                        dict_taxtrim_contigs[tid_parent] += contigs
                    else:
                        contigs_parent = self.dict_taxon_contigs[tid_parent]
                        dict_taxtrim_contigs[tid_parent] += contigs_parent + contigs
                    # Remove existing taxid mapping
                    if tid in dict_taxtrim_contigs:
                        del dict_taxtrim_contigs[tid]
                else:
                    dict_taxtrim_contigs[tid] = contigs
            else:
                dict_taxtrim_contigs[tid] = contigs
        # Store tax -> contig mapping
        self.dict_taxtrim_contigs = dict_taxtrim_contigs
        self.dict_contig_taxtrim = {}
        for taxtrim, contigs in dict_taxtrim_contigs.items():
            for c in contigs:
                self.dict_contig_taxtrim[c] = taxtrim
        # Store new tree
        taxids_trim = list(dict_taxtrim_contigs.keys())
        self.treetrim = self.ncbi.get_topology(taxids_trim)