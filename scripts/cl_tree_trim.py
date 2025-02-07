"""
Filename: cl.tree_trim.py
Author: Benjamin Grodner (https://github.com/benjamingrodner)
Date: 2025-02-04
Description:
    Given a list of contigs mapped to gene annotations, taxonomic annotations, and sample read counts,
    trim the taxonomic tree such that each node in the tree fulfils some threshold
    such as number of contigs, number of genes, presence across multiple samples, etc.

    Save the tree and the new mapping of contig to taxon.
"""
###################################################################################################
# Imports
###################################################################################################

import functions.fn_metat_files as fnf
from collections import defaultdict
from ete4 import NCBITaxa
import numpy as np

###################################################################################################
# Object
###################################################################################################

class TreeTrim:
    """
    Class to manges functions for trimming trees based on thresholding at each node.      
    """
    
    def __init__(self, dict_contig_taxon, dict_contig_ko, dict_contig_estcounts):
        """
        Initializes a TreeTrim instance

        Params:
            dict_contig_taxon : (dict[str, str])
                Dictionary mapping contig name to taxonomic id anotation
            dict_contig_ko : (dict[str, str])
                Dictionary mapping contig name to gene id annotation
            dict_contig_estcounts : (dict[str, list[float]])
                Dictionary mapping contig name to a list of read counts for each sample.
                The number of samples must be the same for each contig (each list the same length).

        Returns:
            None
        """
        # Params
        self.dict_contig_taxon = dict_contig_taxon
        self.dict_contig_ko = dict_contig_ko
        self.dict_contig_estcounts = dict_contig_estcounts

        # Builtin vars
        self.ncbi = NCBITaxa()  # ete4 access to ncbi database
        self.dict_name_filtfunc = {
            "mean_nkos": self.get_mean_nkos,
            "nkos_in_gt_minsamples": self.get_nkos_in_gt_minsamples,
        }  # Functions used to filter tree

        # Execute functions
        self._get_set_taxids()  # All taxids present
        self.tree = self.ncbi.get_topology(self.set_taxids)  # Build tree
        self._get_dict_taxon_contigs()  # invert taxon contig dict
        self.nsamples = len(next(iter(self.dict_contig_estcounts.values())))  # number of samples


    def _get_set_taxids(self):
        """
        Get a set of taxids from the input dictionary that are present in the ncbi taxonomy database.
        """
        taxids_all = []
        for t in self.dict_contig_taxon.values():
            t = t[0]
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
        for c, t in self.dict_contig_taxon.items():
            if t[0] in self.set_taxids:
                lin = self.ncbi.get_lineage(t[0])
                for l in lin:
                    self.dict_taxfull_contigs[str(l)].append(c)


    def _get_dict_taxon_contigs(self):
        """
        Invert the key and value for dict_contig_taxon and build a list of contigs mapped at each
        taxonomic level (as long as it's in the ncbi database).
        """
        self.dict_taxon_contigs = defaultdict(list)
        for c, t in self.dict_contig_taxon.items():
            tid = t[0]
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
            ko = self.dict_contig_ko[c]
            estcounts = self.dict_contig_estcounts[c]
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

        dict_ko_estcounts = defaultdict(lambda: np.zeros(self.nsamples))
        for c in contigs:
            ko = self.dict_contig_ko[c]
            estcounts = self.dict_contig_estcounts[c]
            estcounts = np.array(estcounts)
            dict_ko_estcounts[ko] += estcounts
        nkos = 0
        for ko, estcounts in dict_ko_estcounts.items():
            ko_in_nsam = estcounts.astype(bool).sum()
            if ko_in_nsam >= minsamples:
                nkos += 1
        return nkos            


    def trim_tree(self, filt_func_name, thresh, **fnargs):
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
                contigs = self.dict_tax_contigs[tid]
            # skip if root node
            if not n.is_root:
                # if the mean nkos is less than desired
                func = self.dict_name_filtfunc[filt_func_name]
                filtvalue = func(contigs, **fnargs)
                if filtvalue < thresh:
                    # Then add these contigs to the parent taxon id
                    tid_parent = next(n.ancestors()).name
                    if tid_parent in dict_taxtrim_contigs:
                        dict_taxtrim_contigs[tid_parent] += contigs
                    else:
                        contigs_parent = self.dict_tax_contigs[tid_parent]
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
        # Store new tree
        taxids_trim = list(dict_taxtrim_contigs.keys())
        self.treetrim = self.ncbi.get_topology(taxids_trim)