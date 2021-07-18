# -*- coding: utf-8 -*-
# @Time : 2021/7/16 14:01
# @Author : Zhongyi Hua
# @FileName: gadget.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

import os
import re
from collections import defaultdict


def del_dir(_dir):
    for r, d, f in os.walk(_dir):
        for files in f:
            os.remove(os.path.join(r, files))
        os.removedirs(r)


def read_fasta(_fasta):
    with open(_fasta, 'r') as _f_in:
        _fasta = list(line.rstrip() for line in _f_in)
    d_fasta = defaultdict(str)
    for _line in _fasta:
        if re.search(r'^>', _line):
            scaffold_id = _line.split()[0].replace('>', '')
            continue
        d_fasta[scaffold_id] += _line
    return d_fasta


def convertfasta2pml(_fa_aligned, _paml):
    """
    convert aligned sequences from fasta file to PAML file
    :param _fa_aligned: input aligned fasta file path
    :param _paml: output PAML format file path
    :return: None
    """
    d_fasta = read_fasta(_fa_aligned)
    # remove stop codon at the end
    for _id, _seq in d_fasta.items():
        if re.match('(TAG)|(TAA)|(TGA)', _seq[-3:]):
            d_fasta[_id] = _seq[:-3]
    # for feature remove frameshift sequence
    """
    tmp_fasta = d_fasta.copy()
    for _id, _seq in tmp_fasta.items():
        if '!' in _seq:
            d_fasta.pop(_id)
    """
    # check aligned
    _tmp_list = [len(_seq) for _seq in d_fasta.values()]
    assert max(_tmp_list) == min(_tmp_list), \
        'The sequence length in aligned fasta is different or Some with stop codon, some not?'
    # reformat
    _out_list = [_id+'\n'+_seq+'\n' for _id, _seq in d_fasta.items()]
    _out_list = ['{} {}\n\n'.format(len(_tmp_list), max(_tmp_list))] + _out_list
    with open(_paml, 'w') as _f_out:
        _f_out.write(''.join(_out_list))


def unroot_tree(tree_):
    # check if rooted tree
    leaf_num = len([_.name for _ in tree_ .clade.clades if _.name is not None])
    if leaf_num > 1:
        return
    else:
        for _idx, _clade in enumerate(tree_.clade.clades):
            if _clade.name is None:
                tree_.collapse(tree_.clade.clades[_idx])
                break
        unroot_tree(tree_)
