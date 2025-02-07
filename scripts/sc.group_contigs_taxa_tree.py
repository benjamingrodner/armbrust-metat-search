
"""
Script Name: sc.group_contigs_taxa_tree.py
Author: Benjamin Grodner (https://github.com/benjamingrodner)
Date: 2025-02-04
Description:
    Given a list of contigs mapped to gene annotations, taxonomic annotations, and sample read counts,
    trim the taxonomic tree such that each node in the tree fulfils some threshold
    such as number of contigs, number of genes, presence across multiple samples, etc.

    Save the tree and the new mapping of contig to taxon.

Usage:
    python sc.group_contigs_taxa_tree.py -m <fn_metat> -t <fn_targets> -k <name_key> -v <name_value> -c <columns_line_number> -ik <idx_key> -ivl <idx_value> -o <output_fn> --verbose

"""

###################################################################################################
# Imports
###################################################################################################
from collections import defaultdict
import argparse

###################################################################################################
# Functions
###################################################################################################

def trim_tree(tree, dict_tax_contigs, thresh, fn_contigs2filtvalue, **fn_args):
    # Remove tree levels below filter value
    dict_taxtrim_contigs = defaultdict(list)
    # Traverse from leaves to root
    for n in tree.traverse(strategy='postorder'):
        tid = n.name
        # Get any new contigs from children
        if tid in dict_taxtrim_contigs:
            contigs = dict_taxtrim_contigs[tid]
        # Otherwise get old contigs for node
        else:
            contigs = dict_tax_contigs[tid]
        # skip if root node
        if not n.is_root:
            # if the mean nkos is less than desired
            filtvalue = fn_contigs2filtvalue(contigs, **fn_args)
            if filtvalue < thresh:
                # Then add these contigs to the parent taxon id
                tid_parent = next(n.ancestors()).name
                if tid_parent in dict_taxtrim_contigs:
                    dict_taxtrim_contigs[tid_parent] += contigs
                else:
                    contigs_parent = dict_tax_contigs[tid_parent]
                    dict_taxtrim_contigs[tid_parent] += contigs_parent + contigs
                # Remove existing taxid mapping
                if tid in dict_taxtrim_contigs:
                    del dict_taxtrim_contigs[tid]
            else:
                dict_taxtrim_contigs[tid] = contigs
        else:
            dict_taxtrim_contigs[tid] = contigs
    return dict_taxtrim_contigs

def parse_arguments():
    """
    Parse command-line arguments.
    
    Returns:
        argparse.Namespace: Parsed arguments as an object.
    """
    parser = argparse.ArgumentParser()
    # # Check if there is a stdin 
    # if sys.stdin.isatty():
    #     # If not add input dir from arguments
    #     parser.add_argument(
    #         '-id', '--input_dir', type=str, required=True, 
    #         help="Path to the input directory, can also use stdin with no '-id' argument."
    #     )
    parser.add_argument(
        '-m', '--fn_metat', type=str, required=True,
        help="Path to metaT file."
    )
    parser.add_argument(
        '-mt', '--fn_metat_tar', nargs='?', const='',
        help="If the metat file is a tarball, what is the name of the file within the tarball to extract."
    )
    parser.add_argument(
        '-t', '--fn_targets', type=str, required=True,
        help="Path to a file with a list of targets to find."
    )

    parser.add_argument(
        '-k', '--name_key', type=str, nargs='?', const='',
        help="Column name in metat file to be used for dictionary keys."
    )
    parser.add_argument(
        '-vl', '--name_value', type=str, nargs='?', const='',
        help="Column name in metat file to be used for dictionary values."
    )
    parser.add_argument(
        '-c', '--columns_line_number', type=int, nargs='?', const=1,
        help="Line in the file to search for column names. Default is line 1 (counting starts at 1 not 0)."
    )
    parser.add_argument(
        '-ik', '--idx_key', type=int, nargs='?', const=0,
        help="Column name in metat file to be used for dictionary keys."
    )
    parser.add_argument(
        '-ivl', '--idx_value', type=int, nargs='?', const=1,
        help="Column name in metat file to be used for dictionary values."
    )

    parser.add_argument(
        '-o', '--output_fn', type=str, default="", 
        help="Path to the output file. Default is 'output/{fn_metat}-{fn_targets}-dict-{name_key}-{name_value}.json'."
    )
    parser.add_argument(
        '--verbose', action='store_true', 
        help="Enable verbose mode for detailed logging."
    )
    
    args = parser.parse_args()

    # # Check if there is stdin
    # if not sys.stdin.isatty():
    #     # Read all lines from stdin
    #     input_dir = sys.stdin.read() 
    # else: 
    #     # If not add input dir from args  
    #     input_dir = args.input_dir

    return(args)

###################################################################################################
# Main
###################################################################################################

def main():
    """
    Main function to execute script logic.
    """

    # Load taxa file
    with open(fn_dict_contig_tax, 'r') as f:
        dict_contig_tax = json.load(f)

    # get dict ko contigs
    with open(fn_dict, 'r') as f:
        dict_ko_contigs = json.load(f)
    
    # full tax list
    taxids_all = []
    for t in dict_contig_tax.values():
        t = t[0]
        if int(t):
            try:
                ncbi.get_lineage(t)
                taxids_all.append(t)
            except:
                pass

    set_taxids = set(taxids_all)    

    # get tree
    tree = ncbi.get_topology(set_taxids)

    # Get dict taxid -> contig list (full lineage)
    dict_tax_contigs = defaultdict(list)
    for c, t in dict_contig_tax.items():
        if t[0] in set_taxids:
            lin = ncbi.get_lineage(t[0])
            for l in lin:
                dict_tax_contigs[l].append(c)

    # dict contig to ko
    dict_contig_ko = {}
    for ko, contigs in dict_ko_contigs.items():
        for c in contigs:
            dict_contig_ko[c] = ko

    # # filenames read counts for each sample
    # dir_counts = '/scratch/bgrodner/iron_ko_contigs/metat_search_results/dicts_iron_KO_contig/dicts_contig_count/tables_norm_count/NPac.G1PA.bf100.id99.aa.best.Kofam.incT30.csv.gz-iron_KOs.txt-tables_norm_count'
    # glob_fn = f'{dir_counts}/*.tsv-table_norm_count.csv'
    # fns = glob.glob(glob_fn)

    # # lineage = [str(leaf)] + [n.name for n in tree[str(leaf)].ancestors()]
    # dict_lin_nkos = {}
    # dict_lin_ko_cnts = defaultdict(lambda: defaultdict(list))
    # # Load dict_tax_contigs, dict_contig_ko

    # # Get number of kos for each lineage
    # for n in tree.traverse():
    #     l = n.name
    #     contigs = dict_tax_contigs[int(l)]
    #     kos = [dict_contig_ko[c] for c in contigs]
    #     # add number of kos to lineage info
    #     set_kos = set(kos)
    #     dict_lin_nkos[l] = len(set_kos)

    # # Get KO counts for each sample
    # for fn in tqdm(fns):
    #     if '_3um' in fn:
    #         # Load dict_contig_count
    #         df = pd.read_csv(fn)
    #         ctg = df.iloc[:,0].values
    #         cnt = df.iloc[:,1].values
    #         dict_ctg_cnt = dict(zip(ctg, cnt))
    #         # Count each lineage
    #         # n_kos = 0
    #         # mean_frac_fill = 0
    #         for n in tree.traverse():
    #             l = n.name
    #             # if (n_kos < thresh_n_kos) or (mean_frac_fill < thresh_mff):
    #             contigs = dict_tax_contigs[int(l)]
    #             # counts and kos
    #             dict_ko_cnt = defaultdict(lambda: 0)
    #             for c in contigs:
    #                 ko = dict_contig_ko[c]
    #                 c_ = re.sub(r'_\d+$','',c)  # manage 6tr filenaming
    #                 dict_ko_cnt[ko] += dict_ctg_cnt[c_]
    #             # Add ko count to lineage info
    #             for ko, cnt in dict_ko_cnt.items():
    #                 dict_lin_ko_cnts[l][ko].append(cnt)'

    # dict contig -> est counts
    dir_counts = '/scratch/bgrodner/iron_ko_contigs/metat_search_results/dicts_iron_KO_contig/dicts_contig_count'
    glob_fn = f'{dir_counts}/*G1PA*3um*.json'
    fns = glob.glob(glob_fn)
    dict_contig_estcounts = defaultdict(list)
    contigs_all = list(dict_contig_tax.keys())
    samples = []
    for fn in fns:
        if '_3um' in fn:
            samples.append(re.search(r'(?<=.tar.gz).+(?=\.tsv-)', fn)[0])
            # Load counts
            with open(fn, 'r') as f:
                dict_ctg_cnt = json.load(f)
            for c in contigs_all:
                c_ = re.sub(r'_\d+$','',c)
                # map contig to counts
                dict_contig_estcounts[c].append(float(dict_ctg_cnt[c_][0]))


    # dict mappint taxid to max mean n samples across kos
    dict_lin_maxmeannsam = {}
    for n in tree.traverse():
        l = n.name
        n_sam = []
        for ko, counts in dict_lin_ko_cnts[l].items():
            ns = np.array(counts).astype(bool).sum()
            n_sam.append(ns)
        max_n_sam = np.max(n_sam)
        mean_n_sam = np.round(np.mean(n_sam), 1)
        dict_lin_maxmeannsam[l] = [max_n_sam, mean_n_sam]

    # How many kos are in each sample across the nodes?
    dict_tax_sam_kos = defaultdict(lambda: defaultdict(set))
    for n in tree_trim.traverse():
        contigs = dict_taxtrim_contigs[n.name]
        for c in contigs:
            ko = dict_contig_ko[c]
            estcounts = dict_contig_estcounts[c]
            for sam, ec in enumerate(estcounts):
                if ec:
                    dict_tax_sam_kos[n.name][sam].add(ko)

    # filter value calculation
    def get_mean_nkos(contigs):
        dict_sam_kos = defaultdict(set)
        for c in contigs:
            ko = dict_contig_ko[c]
            estcounts = dict_contig_estcounts[c]
            for sam, ec in enumerate(estcounts):
                if ec:
                    dict_sam_kos[sam].add(ko)
        nkos = []
        for _, kos in dict_sam_kos.items():
            nkos.append(len(kos))
        return np.mean(nkos)

    filt_mean_nkos = 30




    # SAve tree
    # get dict contig to taxon-trim
    # Save dict    
    return 

if __name__ == "__main__":
    main()
