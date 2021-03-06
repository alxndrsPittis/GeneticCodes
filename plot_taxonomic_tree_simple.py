import sys
import os
from sys import exit
from argparse import ArgumentParser, ArgumentError
from collections import OrderedDict

from ete3 import PhyloTree, TreeStyle, NodeStyle, faces, NCBITaxa

__DESCRIPTION__ = """
Plotting of NCBI taxonomy and species trees
and a phylogenetic profile associated to them.
"""

def remove_single_inter_nodes(tNCBI):
    for d in tNCBI.get_descendants():
        if len(d.children) == 1:
            d.delete(prevent_nondicotomic=False)
    return tNCBI


if __name__ == "__main__":
    parser = ArgumentParser(description=__DESCRIPTION__)

    parser.add_argument("-i", "--infile", dest="infile",
                        type=str,
                        required=True,
                        help="""full path of taxids infile""")

    parser.add_argument("--newick", dest="newick",
                        type=str,
                        help="""
                        The newick of a species tree + code2taxid in infile
                        """)    

    parser.add_argument("--mode", dest="mode",
                        type=str,
                        help="""
                        "r" rectangular or "c" circular mode
                        """)

    parser.add_argument("--inter", dest="inter",
                        action="store_true",
                        help="""
                        If true, intermediate nodes are kept
                        """)

    parser.add_argument("--taxoncolors", dest="taxoncolors",
                        type=str,
                        help="""
                        path of color dictionary
                        """)

    parser.add_argument("--save", dest="save",
                        type=str,
                        help="""
                        If provided, tree rendered to file
                        """)

    parser.add_argument("--no_internal_names", dest="no_internal_names",
                        action="store_true",
                        help="""
                        If true, internal names are not plotted
                        """)

    parser.add_argument("--no_intermediate_nodes", dest="no_intermediate_nodes",
                        action="store_true",
                        help="""
                        If true, internal nodes are removed
                        """)    

    parser.add_argument("--no_names", dest="no_names",
                        action="store_true",
                        help="""
                        If true, leaf names are not plotted
                        """)    

    parser.add_argument("--ultrametric", dest="ultrametric",
                        type=int,
                        help="""
                        If true, tree is transformed to ultrametic                        
                        """)

    parser.add_argument("--swap", dest="swap",
                        type=str,
                        default="all",
                        help="""
                        If string=='', all branches are swaped. Otherwise node taxids shall be defined - for visual inspection
                        """)        

    parser.add_argument("--profile", dest="profile",
                        type=str,
                        help="""
                        Phylogenetic profile, tab delimited, taxid per line, column per COG, if T/F presence/absence, else number
                        """)

    parser.add_argument("--bubbles", dest="bubbles",
                        type=str,
                        help="""
                        Plots bubbles in nodes based on taxid 2 value
                        tab file.
                        """)    

    parser.add_argument("--mcl_clusters", dest="mcl_clusters",
                        type=str,
                        help="""
                        MCL output file, cluster per line
                        """)

    parser.add_argument("--mcl_cluster_colors", dest="mcl_cluster_colors",
                        type=str,
                        help="""
                        color per line
                        """)    

    parser.add_argument("--hmmsearch_profile", dest="hmmsearch_profile",
                        type=str,
                        help="""
                        profile from hmmsearch outfile
                        """)    

    parser.add_argument("--seqid2taxid", dest="seqid2taxid",
                        type=str,
                        help="""
                        seqID to taxid dictionary
                        """)

    parser.add_argument("--colorbar", dest="colorbar",
                        action="store_true",
                        help="""
                        Colorbar for the heatmap with matplotlib
                        """)

    parser.add_argument("--colorbar_save", dest="colorbar_save",
                        type=str,
                        help="""
                        save path of Colorbar for the heatmap with matplotlib
                        """)        
    

    args = parser.parse_args()
    infile = args.infile
    mode = args.mode
    newick = args.newick

    if newick:
        t = PhyloTree(args.newick)      
        species2taxid = dict([ line.split()[0], line.strip().split()[1] ] for line in open(infile))
        taxids = set(species2taxid.values())
    else:
        ncbi = NCBITaxa()
        taxids = set([ line.strip() for line in open(infile) ])


    if args.taxoncolors:
        taxon2color = dict([int(line.split()[0]), line.split()[1]] for line in open(args.taxoncolors))

    tNCBI = ncbi.get_topology(taxids, intermediate_nodes=True)
    tNCBI = tNCBI.search_nodes(name="2759")[0]
    ncbi.annotate_tree(tNCBI, taxid_attr="name")
    tax2node = dict([node.taxid, node] for node in tNCBI.traverse())

    if args.no_intermediate_nodes:
        for node in tNCBI.get_descendants():
            if len(node.children) == 1:
                node.delete()
        if len(tNCBI.children) == 1:
            tNCBI = tNCBI.children[0]
    
    tax2node = {}
    for node in tNCBI.traverse():
        tax2node[node.taxid] = node
        if args.taxoncolors:
            if node.taxid in taxon2color:
                node.add_feature("bgcolor", taxon2color[node.taxid])

    if args.swap:
        if args.swap == "all":
            tNCBI.swap_children()
            for node in tNCBI.get_descendants():
                if not node.is_leaf():
                    node.swap_children()
        else:
            taxids2swap == args.swap.split()
            for taxid in taxids2swap:
                tax2node[taxid].swap_children()
        
    if newick:
        for leaf in t.iter_leaves():
            leaf.add_feature("bgcolor", tax2node[species2taxid[leaf.name]].bgcolor)
            leaf.add_feature("taxid", species2taxid[leaf.name])
            leaf.add_feature("sci_name", tax2node[species2taxid[leaf.name]].sci_name)
            
        for node in t.traverse():
            if not node.is_leaf():                
                leaves = [tax2node[leaf.taxid] for leaf in node.iter_leaves()]
                common = tNCBI.get_common_ancestor(leaves)
                node.add_feature("taxid", common.taxid)
                node.add_feature("sci_name", common.sci_name)
                for n in node.traverse():
                    n.add_feature("bgcolor", common.bgcolor)
                
    else:
        pass
        """
        for d in tNCBI.get_descendants():
            if len(d.children) == 1:
                d.delete()
        """
    if args.ultrametric:
        print "ultrametric tree"
        tNCBI.convert_to_ultrametric(tree_length=args.ultrametric)

    if args.bubbles:
        taxid2value = dict([line.split("\t")[0], float(line.split("\t")[1])] for line in open(args.bubbles))
        max_bubble = max(taxid2value.values())
        
    def layout(node):
        node.img_style['size'] = 0
        node.img_style['vt_line_width'] = 0
        node.img_style['hz_line_width'] = 0
        if "bgcolor" in node.features:
            node.img_style['bgcolor'] = node.bgcolor
        if node.is_leaf():
            if not args.no_names:
                name = faces.TextFace(node.sci_name, fsize=12, fstyle='italic')
                faces.add_face_to_node(name, node, 0, aligned=True)
            fake = faces.TextFace(" ")
            fake.background.color = "white"
            faces.add_face_to_node(fake, node, 1, aligned=True) # fake
        else:
            if not args.no_internal_names and node.get_distance(tNCBI, topology_only=True) < 3:
                name = faces.TextFace(node.sci_name, fsize=12, fstyle='italic')
                faces.add_face_to_node(name, node, 0, position='branch-top')

            
    S = TreeStyle()
    #S.allow_face_overlap = True
    S.show_leaf_name = False
    #S.scale = 200
    #S.draw_aligned_faces_as_table = True
    #S.aligned_table_style = 0
    #S.min_leaf_separation = 1
    if args.mode == 'r':
        S.mode = 'r'
    elif args.mode == 'c':
        S.mode = 'c'

    if args.save:
        tNCBI.render(file_name=args.save, layout=layout, tree_style=S)
    else:
        print "showing"
        tNCBI.show(layout=layout, tree_style=S)


        
