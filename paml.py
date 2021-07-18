# -*- coding: utf-8 -*-
# @Time : 2021/7/13 10:51
# @Author : Zhongyi Hua
# @FileName: paml.py
# @Usage: 
# @Note:
# @E-mail: njbxhzy@hotmail.com

import os
import subprocess
from tempfile import mktemp
from Bio.Phylo.PAML import codeml
from gadget import del_dir, convertfasta2pml
from multiprocessing import Pool
from functools import partial

__location__ = os.path.abspath(os.path.join(__file__, '..'))


def write_ctl(ctl_path, ctl_ins):
    """
    Write codeml config file for running
    :param ctl_path: config file path
    :param ctl_ins: a Bio.Phylo.PAML.codeml.Codeml instance (from codeml.Codeml())
    :return:
    """
    with open(ctl_path, "w") as ctl_handle:
        ctl_handle.write("seqfile = %s\n" % ctl_ins.alignment)
        ctl_handle.write("outfile = %s\n" % ctl_ins.out_file)
        ctl_handle.write("treefile = %s\n" % ctl_ins.tree)
        for option in ctl_ins._options.items():
            if option[1] is None:
                continue
            if option[0] == "NSsites":
                NSsites = " ".join(str(site) for site in option[1])
                ctl_handle.write("%s = %s\n" % (option[0], NSsites))
            else:
                ctl_handle.write("%s = %s\n" % (option[0], option[1]))


class CodemlExecutor:
    def __init__(self, _executable, _chi2):
        self.executable = _executable
        self.chi2 = _chi2

    def run_codeml(self, _ctl):
        """
        run codeml bypass the environment path
        :param _ctl:
        :return:
        """
        subprocess.run([self.executable, _ctl], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    def run_chi2(self, _chi2_in):
        """
        run chi2 bypass the environment path
        :param _chi2_in: 2*|lnF1 - lnF2|
        :return: p-value
        """
        _tmp_r = subprocess.run([self.chi2, '2', str(_chi2_in)], capture_output=True, text=True)
        _text = _tmp_r.stdout
        return float(_text.split()[-1])


class PerPAML():
    def __init__(self, _alignment_path, _tree, _wd, _icode):
        """
        :param _alignment_path:
        :param _tree:
        :param _wd: Working Directory. Using temporary directory.
        :param _outfile:
        """
        self.alignment = _alignment_path
        self.tree = _tree
        self.wd = _wd
        self.icode = _icode

    def runAmodel(self, _cmlexcutor):
        """
        :param _cmlexcutor: a CodemlExecutor instance
        """
        cml = codeml.Codeml()
        cml.read_ctl_file(os.path.join(__location__, 'configs', 'modelA.ctl'))
        cml.set_options(icode=self.icode-1)
        cml.tree = self.tree
        cml.alignment = self.alignment
        cml.working_dir = self.wd
        cml.out_file = os.path.join(self.wd, 'modelA.out')
        write_ctl(os.path.join(self.wd, 'modelA.ctl'), cml)
        _cmlexcutor.run_codeml(os.path.join(self.wd, 'modelA.ctl'))
        results = codeml.read(os.path.join(self.wd, 'modelA.out'))
        return results['NSsites'][2]['lnL']

    def runAnullmodel(self, _cmlexcutor):
        """
        :param _cmlexcutor: a CodemlExecutor instance
        """
        cml = codeml.Codeml()
        cml.read_ctl_file(os.path.join(__location__, 'configs', 'modelAnull.ctl'))
        cml.set_options(icode=self.icode-1)
        cml.tree = self.tree
        cml.alignment = self.alignment
        cml.working_dir = self.wd
        cml.out_file = os.path.join(self.wd, 'modelAnull.out')
        write_ctl(os.path.join(self.wd, 'modelAnull.ctl'), cml)
        _cmlexcutor.run_codeml(os.path.join(self.wd, 'modelAnull.ctl'))
        results = codeml.read(os.path.join(self.wd, 'modelAnull.out'))
        return results['NSsites'][2]['lnL']

    def calculate_p_value(self, _cmlexcutor):
        chi2_in = str(2*abs(self.runAnullmodel(_cmlexcutor)-self.runAmodel(_cmlexcutor)))
        return _cmlexcutor.run_chi2(chi2_in)

def _calculate_p(_item, _tmp_dir, _cmlexecutor, _tree, _icode):
    _tmp_dir2 = mktemp()
    os.mkdir(_tmp_dir2)
    convertfasta2pml(_item[1], os.path.join(_tmp_dir2, 'align.pml'))
    paml_ins = PerPAML(os.path.join(_tmp_dir2, 'align.pml'), _tree, _tmp_dir2, _icode)
    p_value = paml_ins.calculate_p_value(_cmlexecutor)
    # output
    out_file = open(os.path.join(_tmp_dir, _item[0] + '.txt'), 'w')
    print(_item[0], str(p_value), sep="\t", file=out_file)
    out_file.close()
    del_dir(_tmp_dir2)
        

def batch_calculate(meta, tree, cmlexecutor, threads, icode):
    """

    :param _meta: two columns with no head: family_name alignment_path
    :param _tree:
    :return:
    """
    # PAML
    with open(meta) as _f_meta:
        _item_list = [_.split() for _ in _f_meta.read().split('\n')]
    positive_selective = []
    _tmp_dir = './BatchPAML_Results'
    os.mkdir(_tmp_dir)
    __pool1 = Pool(threads)
    for _item in _item_list:
        __pool1.apply_async(_calculate_p, (_item, _tmp_dir, cmlexecutor, tree, icode))
    __pool1.close()
    __pool1.join()

    positive_selective = []
    for o_p in [os.path.join(_tmp_dir, file) for file in os.listdir(_tmp_dir) if file.endswith('.txt')]:
        with open(o_p, 'r') as fh:
            positive_selective.append(fh.read().strip().split())
    
    del_dir(_tmp_dir)
    
    return positive_selective
