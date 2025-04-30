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
            "nkos_in_gt_minsamples": self.get_nkos_in_gt_minsamples,
            "nkos_in_gt_minbatches": self.get_nkos_in_gt_minbatches,
        }  # Functions used to filter tree
            # "mean_nkos": self.get_mean_nkos,

        # Execute functions
        print(f'Loading dicts for {fn_tidycsv}...')
        fnf.getmem()
        t0 = time()
        self._get_dicts()  # Load from tidytable
        t1 = time()
        el = (t1 - t0)
        elm = el // 60
        els = el % 60
        print(f"Time to load dicts from {fn_tidycsv}: {elm} min {els} sec")
        fnf.getmem()

        self._get_set_taxids()  # All taxids present
        self.tree = self.ncbi.get_topology(self.set_taxids)  # Build tree
        self._get_dict_taxon_contigs()  # invert taxon contig dict
        # self.nsamples = len(next(iter(self.dict_contig[self.col_sample].values())))  # number of samples

    def _get_dicts(self):
        """
        Convert the tidytable csv into nested dictionaries
        """
        # Load the csv file as a dict mapping contigs to other values
        kwargs = {
            'col_key': self.col_contig,
            'cols_key2value': [self.col_ko, self.col_taxon],
            'cols_key2list': [self.col_estcounts, self.col_sample],
            'floats': [self.col_estcounts]
        }
        self.dict_contig = fnf.tidytable_to_dict(
            fn = self.fn_tidycsv, 
            **kwargs
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
       #  Deprecated  # 

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
        in at least a minimum number of samples. 

        Params:
            contigs : list[str]
                List of contig ids
            minsamples : int
                Minimum number of samples a gene has to be present in to be counted

        Returns:
            float
                Number of genes in greater than the minimum number of samples
            dict[dict[float]]
                mapping gene to a dict mapping samples to estcounts summed across contigs                
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
        return nkos, dict_ko_sam_estcounts

    def get_nkos_in_gt_minbatches(self, contigs, dict_sample_batch, minsamples=3, minbatches=2):
        """
        Given a set of contigs within a taxon, count the number of genes with nonzero read counts
        in at least a minimum number of samples across a minimum number of batches. 

        Params:
            contigs : list[str]
                List of contig ids
            dict_sample_batch : dict[str]
                Dictionary mapping sample name to batch name
            minsamples : int
                Minimum number of samples a gene has to be present in (in a batch) for a batch
                to be counted
            minbatches : int
                Minimum number of batches a gene has to be present in to be counted

        Returns:
            float
                Number of genes in >= minsamples in >= minbatches
            dict[dict[float]]
                mapping gene to a dict mapping samples to estcounts summed across contigs

        """    
        # Build output dict, sum ec counts across contigs for each gene
        dict_ko_sam_estcounts = defaultdict(lambda: 
            defaultdict(float)
        )
        for c in contigs:
            ko = self.dict_contig[self.col_ko][c]
            sams = self.dict_contig[self.col_sample][c]
            estcounts = self.dict_contig[self.col_estcounts][c]
            for s, ec in zip(sams, estcounts):
                batch = dict_sample_batch[s]
                dict_ko_sam_estcounts[ko][s] += float(ec)
        # Count the number of kos
        nkos = 0
        for ko, dict_sam_estcounts in dict_ko_sam_estcounts.items():
            # Get the number of samples this ko is in for each batch
            ko_in_nsam = defaultdict(float)
            for s, ec in dict_sam_estcounts.items():
                b = dict_sample_batch[s]
                ko_in_nsam[b] += bool(ec)
            # Count the number of batches
            nbatch = 0
            for _, nsam in ko_in_nsam.items():
                # Only count batches that meet the criteria
                if nsam >= minsamples:
                    nbatch += 1
            # Only count kos that meet the criteria
            if nbatch >= minbatches:
                nkos += 1
        return nkos, dict_ko_sam_estcounts

    def trim_tree(self, filt_func_name, thresh, **fncargs):
        """
        Trim the tree such that the contigs assigned to each node have some measured value passing
        a threshold. 

        Produces self.treetrim and self.dict_taxtrim_contigs

        Params: 
            filt_func_name : str
                Specify which function to use to measure the contigs at a node. Current options are
                'nkos_in_gt_minsamples' and 'nkos_in_gt_minbatches'.
            thresh : float
                Threshold value used to decide whether to keep or trim a node
            **fncargs : various
                Keyword arguments to pass to the function specified in filt_func_name
        
        Returns:
            None
        """
        # Remove tree levels below filter value
        dict_taxtrim_contigs = defaultdict(list)
        # Traverse from leaves to root
        dict_taxtrim_ko_sam_estcounts = {}
        for n in self.tree.traverse(strategy='postorder'):
            tid = n.name
            # Get any new contigs from children
            if tid in dict_taxtrim_contigs:
                contigs = dict_taxtrim_contigs[tid]
            # Otherwise get old contigs for node
            else:
                contigs = self.dict_taxon_contigs[tid]
            # Get the filter function value
            func = self.dict_name_filtfunc[filt_func_name]
            filtvalue, dict_ko_sam_estcounts = func(contigs, **fncargs)
            # skip if root node
            if not n.is_root:
                # if the mean nkos is less than desired
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
                    dict_taxtrim_ko_sam_estcounts[tid] = dict_ko_sam_estcounts
            else:
                dict_taxtrim_contigs[tid] = contigs
                dict_taxtrim_ko_sam_estcounts[tid] = dict_ko_sam_estcounts
        # Store tensor dict
        self.dict_taxtrim_ko_sam_estcounts = dict_taxtrim_ko_sam_estcounts
        # Store tax -> contig mapping
        self.dict_taxtrim_contigs = dict_taxtrim_contigs
        self.dict_contig_taxtrim = {}
        for taxtrim, contigs in dict_taxtrim_contigs.items():
            for c in contigs:
                self.dict_contig_taxtrim[c] = taxtrim
        # Store new tree
        taxids_trim = list(dict_taxtrim_contigs.keys())
        self.treetrim = self.ncbi.get_topology(taxids_trim)
