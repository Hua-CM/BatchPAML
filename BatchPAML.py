# -*- coding: utf-8 -*-
# @Time : 2021/7/15 15:10
# @Author : Zhongyi Hua
# @FileName: BatchPAML.py
# @Usage:
# @Note:
# @E-mail: njbxhzy@hotmail.com

import argparse
import sys
from Bio import Phylo
from tempfile import mktemp
from paml import CodemlExecutor, batch_calculate
from gadget import unroot_tree
from interface import ortho_calculate
from os import remove


def getArgs():
    parser = argparse.ArgumentParser(
        description='This package was for calculating positive selection in batches of families using PAML')
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-o', '--out', required=True,
                        help='<file_path> output path')
    parent_parser.add_argument('-c', '--codeml', default=None,
                        help='<bin_path> codeml bin path (if not in ENVIRONMENT PATH)')
    parent_parser.add_argument('-d', '--chi2', default=None,
                        help='<bin_path> chi2 bin path (if not in ENVIRONMENT PATH)')
    parent_parser.add_argument('-@', '--threads', default=2, type=int,
                        help='<int> The number of threads to use (default: 2)')
    parent_parser.add_argument('-i', '--icode', default=1, type=int,
                        help='<int> The codon table number according to NCBI (default: 1)')

    sub_parser = parser.add_subparsers(dest='universe/orthofinder', title='Available')
    sub_parser.required = True

    universe_parser = sub_parser.add_parser('universe', help='use meta file for batchPAML', parents=[parent_parser])
    universe_parser.add_argument('-m', '--meta', required=True,
                                 help='<file_path> The meta file')
    universe_parser.add_argument('-t', '--tree', required=True,
                                 help='<file_path> tree path (tree in newick format)')
    universe_parser.set_defaults(subcmd='universe')

    ortho_parser = sub_parser.add_parser('orthofinder', help='use meta file for batchPAML', parents=[parent_parser])
    ortho_parser.add_argument('-f', '--fasta', required=True,
                              help='<file_path> The cds fasta file (All species in one)')
    ortho_parser.add_argument('-g', '--groups', required=True,
                              help='<file_path> The Orthogroups.tsv')
    ortho_parser.add_argument('-s', '--single', required=True,
                              help='<file_path> The Orthogroups_SingleCopyOrthologues.txt')
    ortho_parser.add_argument('-t', '--tree', required=True,
                              help='<file_path> Species tree with marked foreground branch')
    ortho_parser.add_argument('-m', '--macse', required=True,
                              help='<file_path> MACSE java bin')
    ortho_parser.set_defaults(subcmd='orthofinder')

    args = parser.parse_args()
    return args


def main(args):
    tree = Phylo.read(args.tree, 'newick')
    tmp_tree_path = mktemp()
    # unroot tree
    unroot_tree(tree)
    Phylo.write(tree, tmp_tree_path, 'newick')
    # wrap executor
    if args.codeml is None:
        if sys.platform == 'linux':
            args.codeml = 'codeml'
        elif sys.platform == 'win32':
            args.codeml='codeml.exe'
    if args.chi2 is None:
        if sys.platform == 'linux':
            args.chi2 = 'chi2'
        elif sys.platform == 'win32':
            args.codeml = 'chi2.exe'
    codeml_ins = CodemlExecutor(args.codeml, args.chi2)
    if args.subcmd == 'universe':
        family_list = batch_calculate(args.meta, tmp_tree_path, codeml_ins, args.threads, args.icode)
    if args.subcmd == 'orthofinder':
        family_list = ortho_calculate(args.macse, codeml_ins, args.fasta,
                                      args.single, args.groups, tmp_tree_path, args.threads, args.icode)
    with open(args.out, 'w') as f_out:
        f_out.write('\n'.join(['\t'.join(_) for _ in family_list]))
    remove(tmp_tree_path)


if __name__ == '__main__':
    main(getArgs())
