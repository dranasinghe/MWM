import os
import logging
import numpy as np

class Gaussian(object):
    """
    parse and create Gaussian files.
    """
    def __init__(self, logfile=None, config_file=None, config_script=None, out_dir=None, quiet=False, raise_error=True):
        """
        Args:
        :param logfile: path to Gaussian log file
        :param config_file: path to file template to create Gaussian input files
        :param config_script: a string template to create Gaussian input files
        :param out_dir: path where Gaussian inout files are going to written
        :param quiet: controls printing
        :param raise_error: control error messages
        """
        self.logfile = logfile
        self.errors = []
        with open(self.logfile) as f:
            self.log = f.read().splitlines()
            for i, line in enumerate(self.log):
                if 'Error termination' in line:
                    self.errors.append(self.log[i])
        if self.errors:
            err_msg = f'Error in {self.logfile}: {self.errors}'
            if raise_error:
                raise Exception(err_msg)
            if not quiet:
                print(err_msg)

        self.out_dir = out_dir
        if self.out_dir is not None:
            os.makedirs(self.out_dir, exist_ok=True)

        if config_file is not None:
            with open(config_file) as f:
                self.config = [line.strip() for line in f]
        else:
            self.config = None

        if config_script is not None:
            self.config_script = config_script
        else:
            self.config_script = """%chk=check.chk
%mem=1000mb
%NProcShared=4
#P guess=mix apfd/defsvp scf=xqc
 
title
 
{charge} {multiplicity}
{xyz}


"""

    def get_energy(self, first=True):
        """We would like to parse and put diffent level of theory to dictionary"""
        if first:
            iterable = self.log
        else:
            iterable = reversed(self.log)
        E0 = {}
        e_elect, e0_composite, scaled_zpe = None, None, None
        for line in iterable:
            if 'SCF Done:' in line:
                e_elect = float(line.split()[4])
                elect_energy_source = 'SCF'
                E0['SCF'] = e_elect
            elif ' E2(' in line and ' E(' in line:
                e_elect = float(line.split()[-1].replace('D', 'E'))
                elect_energy_source = 'doublehybrd or MP2'
                E0['doublehybrd'] = e_elect
            elif 'MP2 =' in line:
                e_elect = float(line.split()[-1].replace('D', 'E'))
                elect_energy_source = 'MP2'
                E0['MP2'] = e_elect
            elif 'BD(T)= ' in line:
                e_elect = float(line.split()[1].replace('D', 'E'))
                elect_energy_source = 'BD(T)'
                E0['BD(T)'] = e_elect
            elif ' E(CORR)= ' in line:
                e_elect = float(line.split()[3])
                elect_energy_source = 'CCSD'
                E0['CCSD'] = e_elect
            elif 'CCSD(T)= ' in line:
                e_elect = float(line.split()[1].replace('D', 'E'))
                elect_energy_source = 'CCSD(T)'
                E0['CCSD(T)'] = e_elect
            elif 'CBS-QB3 (0 K)' in line:
                e0_composite = float(line.split()[3])
            elif 'G3(0 K)' in line:
                e0_composite = float(line.split()[2])
            elif 'E(ZPE)' in line:
                scaled_zpe = float(line.split()[1])

        if e0_composite is not None:
            logging.debug("Using the composite energy from the gaussian output file")
            if scaled_zpe is None:
                raise logging.error('Unable to find zero-point energy in Gaussian.py log file.')
            E0['composite'] = e0_composite - scaled_zpe
        elif e_elect is not None:
            logging.debug("Using the {0} energy from the gaussian output file".format(elect_energy_source))
            return E0
        else:
            raise logging.error('Unable to find energy in Gaussian.py log file.')

    def get_geometry(self, first=False):
        if first:
            iterable = range(len(self.log))
        else:
            iterable = reversed(range(len(self.log)))
        for i in iterable:
            line = self.log[i]
            if 'Input orientation' in line:
                symbols, coords = [], []
                for line in self.log[(i + 5):]:
                    if '------------------------------------------------------------------' in line:
                        return symbols, np.array(coords)
                    # print(line)
                    data = line.split()
                    symbols.append(data[1])
                    coords.append([float(c) for c in data[3:]])

        else:
            raise Exception(f'Geometry not found in {self.logfile}')

    def get_frequencies(self):
        freq = []
        freqs = []
        for i, line in enumerate(self.log):
            if 'Frequencies --' in line:
                freq.extend([float(f) for f in line.split()[2:]])
            elif 'Thermochemistry' in line:
                freqs.append(freq)
                freq = []

        return freqs

    def get_zpe(self):
        for line in reversed(self.log):
            if 'Zero-point correction' in line:
                return float(line.split()[-2])  # in Hartree
        else:
            raise Exception(f'ZPE not found in {self.logfile}')

    def get_negative_frequency(self):
        freqs = self.get_frequencies()
        negative_freqs = []
        for freq in freqs:
            negative_freqs.append([f for f in freq if f < 0])
        if len(negative_freqs) > 0:
            return negative_freqs
        else:
            raise Exception(f'negative frequancy not found in {self.logfile}')

    def get_charge_and_multiplicity(self):
        for i, line in enumerate(self.log):
            if 'Charge =' in line:
                charge = line.split()[2]
            if 'Multiplicity =' in line:
                multiplicity = line.split()[-1]
                return int(charge), int(multiplicity)
        else:
            raise Exception(f'Multiplicity not found in {self.logfile}')

    def make_input_file(self, name, mol, charge=0, multiplicity=1, cid=0, comment=None):
        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        coords = mol.GetConformer(cid).GetPositions()
        self.make_input_file_from_xyz(name, symbols, coords, charge=charge, multiplicity=multiplicity, comment=comment)

    def make_input_file_from_xyz(self, name, symbols, coords, charge=0, multiplicity=1, comment=None):
        xyz_str = ''
        for s, c in zip(symbols, coords):
            xyz_str = xyz_str + f'{s}  {c[0]: .10f}  {c[1]: .10f}  {c[2]: .10f}\n'

        script = self.config_script
        with open(os.path.join(self.out_dir, name + '.gjf'), 'w') as f:
            f.write(script.format(charge=charge, multiplicity=multiplicity, xyz=xyz_str))
