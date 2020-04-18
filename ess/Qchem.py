import os
import logging

class QChem(object):
    """
    parse and  create Qchem files
    """
    def __init__(self, logfile=None, config_file=None, config_script=None, out_dir=None, quiet=False, raise_error=True):
        """

        :param logfile: path to Qchem log file
        :param config_file: path to file template to create Qchem input files
        :param config_script: a string template to create Qchem input files
        :param out_dir: path where Qchem inout files are going to written
        :param quiet: controls printing
        :param raise_error: control error messages
        """
        self.logfile = logfile
        self.errors = []
        if self.logfile is None:
            self.log = None
        else:
            with open(self.logfile) as f:
                self.log = f.read().splitlines()
                for i, line in enumerate(self.log):
                    if 'Q-Chem fatal error' in line:
                        self.errors.append(self.log[i + 2])
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

    def make_input_file(self, name, mol, charge=0, multiplicity=1, cid=0, comment=None):
        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        coords = mol.GetConformer(cid).GetPositions()
        self.make_input_file_from_xyz(name, symbols, coords, charge=charge, multiplicity=multiplicity, comment=comment)

    def make_input_file_from_xyz(self, name, symbols, coords, charge=0, multiplicity=1, comment=None):
        xyz_str = ''
        for s, c in zip(symbols, coords):
            xyz_str = xyz_str + f'{s}  {c[0]: .10f}  {c[1]: .10f}  {c[2]: .10f}\n'

        script = self.config_script
        with open(os.path.join(self.out_dir, name + '.in'), 'w') as f:
            f.write(script.format(charge=charge, multiplicity=multiplicity, xyz=xyz_str))

    def get_energy(self, first=False):
        if first:
            iterable = self.log
        else:
            iterable = reversed(self.log)
        for line in iterable:
            if 'total energy' in line:  # Double hybrid methods
                return float(line.split()[-2])
            elif 'energy in the final basis set' in line:  # Other DFT methods
                return float(line.split()[-1])
        else:
            raise Exception(f'Energy not found in {self.logfile}')

    def get_geometry(self, first=False):
        if first:
            iterable = range(len(self.log))
        else:
            iterable = reversed(range(len(self.log)))
        for i in iterable:
            line = self.log[i]
            if 'Standard Nuclear Orientation' in line:
                symbols, coords = [], []
                for line in self.log[(i + 3):]:
                    if '----------' not in line:
                        data = line.split()
                        symbols.append(data[1])
                        coords.append([float(c) for c in data[2:]])
                    else:
                        return symbols, np.array(coords)
        else:
            raise Exception(f'Geometry not found in {self.logfile}')

    def get_frequencies(self):
        freqs = []
        for line in reversed(self.log):
            if 'Frequency' in line:
                freqs.extend([float(f) for f in reversed(line.split()[1:])])
            elif 'VIBRATIONAL ANALYSIS' in line:
                freqs.reverse()
                return np.array(freqs)
        else:
            raise Exception(f'Frequencies not found in {self.logfile}')

    def get_zpe(self):
        for line in reversed(self.log):
            if 'Zero point vibrational energy' in line:
                return float(line.split()[-2]) / 627.5095  # Convert to Hartree
        else:
            raise Exception(f'ZPE not found in {self.logfile}')

    def get_charge_and_multiplicity(self):
        for i, line in enumerate(self.log):
            if '$molecule' in line:
                charge, multiplicity = self.log[i + 1].strip().split()
                return int(charge), int(multiplicity)
        else:
            raise Exception(f'Multiplicity not found in {self.logfile}')

    def get_title(self):
        for i, line in enumerate(self.log):
            if 'User input:' in line:
                # name, smiles = self.log[i+2].strip().split()
                return self.log[i + 2]
        else:
            raise Exception(f'title not found in {self.logfile}')