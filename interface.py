# -*- coding: utf-8 -*-
# @Time : 2021/7/15 12:55
# @Author : Zhongyi Hua
# @FileName: interface.py
# @Usage: Interface for OrthoFinder results
# @Note:
# @E-mail: njbxhzy@hotmail.com

import os
import subprocess
import pandas as pd
from operator import itemgetter
from functools import partial
from multiprocessing import Pool
from paml import PerPAML, convertfasta2pml
from gadget import read_fasta, del_dir
from tempfile import mktemp


class FromOrthoFinder:
    def __init__(self, cds_fasta, single_list, orthogroups, tree, wd):
        self.fasta = read_fasta(cds_fasta)
        self.single_list = self.__tidy2(single_list)
        self.orthogroups = self.__tidy1(orthogroups)
        self.workingdir = wd
        self.tree = tree

    @staticmethod
    def __tidy1(orthofinder_tsv):
        _tmp_df = pd.read_table(orthofinder_tsv)
        _tmp_dict = {}
        for _ in _tmp_df.to_dict(orient='records'):
            _group = _.pop('Orthogroup')
            _tmp_dict[_group]={}
            for _key, _value in _.items():
                _tmp_dict[_group][_value] =  _key
        return _tmp_dict

    @staticmethod
    def __tidy2(single_file):
        with open(single_file) as f_in:
            single_list = f_in.read().split()
        return single_list

    def extract_and_rename(self):
        for _family in self.single_list:
            _family_ids = tuple(self.orthogroups[_family].keys())
            _family_seqs = itemgetter(*_family_ids)(self.fasta)
            _family_outs = []
            for _oldid, _seq in zip(_family_ids, _family_seqs):
                _id = self.orthogroups[_family][_oldid]
                _family_outs.append('>' + _id + '\n' + _seq + '\n')
            with open(os.path.join(self.workingdir, _family + '.fasta'), 'w') as f_out:
                f_out.write(''.join(_family_outs))

def _calculate_p(_msa, tree, cml_ins, icode):
    _wd = mktemp()
    os.mkdir(_wd)
    _family_name = os.path.split(_msa)[1].split('_')[0]
    convertfasta2pml(_msa, os.path.join(_wd, 'align.pml'))
    paml_ins = PerPAML(os.path.join(_wd, 'align.pml'), tree, _wd, icode)
    p_value = paml_ins.calculate_p_value(cml_ins)
    # output
    out_file = open(os.path.join(os.path.split(_msa)[0], _family_name + '.txt'), 'w')
    print(_family_name, str(p_value), sep="\t", file=out_file)
    out_file.close()
    del_dir(_wd)

def _align(in_out, macse_bin, icode):
    subprocess.run(['java', '-jar', macse_bin, '-prog', 'alignSequences', '-gc_def', str(icode),
                    '-seq', in_out[0], '-out_NT', in_out[1]],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.STDOUT)

def ortho_calculate(macse_bin, cml_ins, cds_fasta, single_list, orthogroups, tree, threads, icode):
    _tmp_dir = './BatchPAML_Results'
    os.mkdir(_tmp_dir)
    print('Tidy OrthoFinder result start')
    ortho_ins = FromOrthoFinder(cds_fasta, single_list, orthogroups, tree, _tmp_dir)
    ortho_ins.extract_and_rename()
    print('Tidy OrthoFinder result done')
    # align
    print('Align start')
    func1 = partial(_align, macse_bin=macse_bin, icode=icode)
    _in_out = [(os.path.join(_tmp_dir, _ + '.fasta'), os.path.join(_tmp_dir, _ + '_aligned.fasta')) for _ in ortho_ins.single_list]
    __pool1 = Pool(threads)
    __pool1.map(func1, _in_out)
    __pool1.close()
    __pool1.join()
    print('Align done')
    
    # calculate p
    print('Calculate start')
    aligned_fastas = [os.path.join(_tmp_dir, file) for file in os.listdir(_tmp_dir) if file.endswith('_aligned.fasta')]
    __pool2 = Pool(threads)
    for fasta in aligned_fastas:
        __pool2.apply_async(_calculate_p, (fasta, tree, cml_ins, icode))

    __pool2.close()
    __pool2.join()

    positive_selective = []
    for o_p in [os.path.join(_tmp_dir, file) for file in os.listdir(_tmp_dir) if file.endswith('.txt')]:
        with open(o_p, 'r') as fh:
            positive_selective.append(fh.read().strip().split())

    #del_dir(_tmp_dir)
    print('Calculate done')
    return positive_selective
