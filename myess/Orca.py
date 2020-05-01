import os

class Orca(object):
    """
    parse and create Orca files
    """
    def __init__(self, logfile=None, config_file=None, config_script=None, out_dir=None, quiet=False, raise_error=True):
        """

        :param logfile: path to Orca log file
        :param config_file: path to file template to create Orca input files
        :param config_script: a string template to create Orca input files
        :param out_dir: path where Orca inout files are going to written
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
                    if 'ORCA finished by error termination' in line:
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
        else:
            self.config_script = """!UHF DLPNO-CCSD(T) def2-TZVP def2-TZVP/C TightSCF normalPNO
%maxcore 10000
#
%pal
nproc {nproc}
end
#how to do thermo
*xyz {charge} {multiplicity}
{xyz}*
"""

    def make_input_file(self, name, mol, charge=0, multiplicity=1, cid=0, comment=None):
        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        coords = mol.GetConformer(cid).GetPositions()
        self.make_input_file_from_xyz(name, symbols, coords, charge=charge, multiplicity=multiplicity, comment=comment)

    def make_input_file_from_xyz(self, name, symbols, coords, charge=0, multiplicity=1, comment=None):
        #        config = self.config[:]

        # setting number of proc depending on size of molecule
        if len(symbols) < 3:
            nproc = 1
        elif len(symbols) < 8:
            nproc = 8
        elif len(symbols) >= 4 and len(symbols) < 15 and 'H' in symbols:
            nproc = 16
        elif len(symbols) >= 15 :
            nproc = 20
        else:
            nproc = 4

        xyz_str = ''
        for s, c in zip(symbols, coords):
            xyz_str = xyz_str + f'{s}  {c[0]: .10f}  {c[1]: .10f}  {c[2]: .10f}\n'

        script = self.config_script
        with open(os.path.join(self.out_dir, name + '.inp'), 'w') as f:
            # for line in config:
            # f.write(line + '\n')
            f.write(script.format(nproc=nproc, charge=charge, multiplicity=multiplicity, xyz=xyz_str))

    def get_energy(self, first=False):
        if first:
            iterable = self.log
        else:
            iterable = reversed(self.log)
        for line in iterable:
            if 'FINAL SINGLE POINT ENERGY' in line:  # for all methods in Orca.py
                return float(line.split()[-1])
            # elif 'energy in the final basis set' in line:  # Other DFT methods
            #    return float(line.split()[-1])
        else:
            raise Exception(f'Energy not found in {self.logfile}')

    def get_geometry(self, first=False):
        if first:
            iterable = range(len(self.log))
        else:
            iterable = reversed(range(len(self.log)))
        for i in iterable:
            line = self.log[i]
            if 'CARTESIAN COORDINATES (ANGSTROEM)' in line:
                symbols, coords = [], []
                for line in self.log[(i + 2):]:
                    if not line.strip():
                        return symbols, np.array(coords)
                    if '---------------------------------' not in line:
                        # print(line)
                        data = line.split()
                        symbols.append(data[0])
                        coords.append([float(c) for c in data[1:]])

        else:
            raise Exception(f'Geometry not found in {self.logfile}')

    def get_frequencies(self):
        freqs = []
        for i, line in enumerate(self.log):
            if ' Mode    freq' in line:
                for line in self.log[(i + 2):]:
                    if not line.strip():
                        return np.array(freqs)
                    freqs.extend([float(line.split()[1])])

        else:
            raise Exception(f'Frequencies not found in {self.logfile}')

    def get_zpe(self):
        for line in reversed(self.log):
            if 'Zero point energy' in line:
                return float(line.split()[-4])  # in Hartree
        else:
            raise Exception(f'ZPE not found in {self.logfile}')

    def get_charge_and_multiplicity(self):
        for i, line in enumerate(self.log):
            if 'Total Charge' in line:
                charge = line.split()[-1]
            if 'Multiplicity' in line:
                multiplicity = line.split()[-1]
                return int(charge), int(multiplicity)
        else:
            raise Exception(f'Multiplicity not found in {self.logfile}')

    def get_mayer_pop(self):
        pop = []
        for i, line in enumerate(reversed(self.log)):
            if ' Mayer\'s free valence' in line:
                for line in self.log[(len(self.log) - i + 2):]:
                    if not line.strip():
                        # return np.array(pop)
                        return self.pop_array_to_dic(pop)
                    pop.extend([line.split()])

    def pop_array_to_dic(self, pop):
        print('NA   - Mulliken gross atomic population')
        print('ZA   - Total nuclear charge')
        print('QA   - Mulliken gross atomic charge')
        print('VA   - Mayer\'s total valence')
        print('BVA  - Mayer\'s bonded valence')
        print('FA   - Mayer\'s free valence')
        pop_dic = {}
        atom = []
        element = []
        NA = []
        ZA = []
        QA = []
        VA = []
        BVA = []
        FA = []
        for i, l in enumerate(pop):
            atom.append(int(l[0]))
            element.append(l[1])
            NA.append(float(l[2]))
            ZA.append(float(l[3]))
            QA.append(float(l[4]))
            VA.append(float(l[5]))
            BVA.append(float(l[6]))
            FA.append(float(l[7]))

        pop_dic['atom'] = atom
        pop_dic['element'] = element
        pop_dic['NA'] = NA
        pop_dic['ZA'] = ZA
        pop_dic['QA'] = QA
        pop_dic['VA'] = VA
        pop_dic['BVA'] = BVA
        pop_dic['FA'] = FA

        pop = pd.DataFrame.from_dict(pop_dic)
        pop.to_csv(str(self.logfile) + '.csv')

        return pop_dic

    def get_mayer_bond(self):
        bond = []
        for i, line in enumerate(reversed(self.log)):
            if ' Mayer bond orders larger than' in line:
                for line in self.log[(len(self.log) - i):]:
                    if not line.strip():
                        # return np.array(bond)
                        return self.rearrange_bond(bond)
                    bond.extend([line.replace('B(', '').replace(')', '').replace(',', '').replace(':', '').split()])

    def rearrange_bond(self, bond):
        new_bond = []
        for l in bond:
            if len(l) == 3:
                new_bond.append([[l[0], l[1]], float(l[2])])
            if len(l) == 6:
                new_bond.append([[l[0], l[1]], float(l[2])])
                new_bond.append([[l[3], l[4]], float(l[5])])
            if len(l) == 9:
                new_bond.append([[l[0], l[1]], float(l[2])])
                new_bond.append([[l[3], l[4]], float(l[5])])
                new_bond.append([[l[6], l[7]], float(l[8])])
        return new_bond

