import re
from Queue import Queue

class Molecule(object):
    def __init__(self, name=None):
        self.name = name
        self.atoms = {}
    def readitp(self, itp_file):
        """ load molecule information from given itp_file """
        with open(itp_file) as f:
            mode = None
            for each_line in f:
                each_line = self._remove_comments(each_line)
                if each_line.strip():
                    m = re.search(r'\[\s*([^\s]+)\s*\]', each_line)
                    if m:
                        section = m.group(1)
                        if section == 'atoms':
                            mode = 'atoms'
                        elif section == 'bonds':
                            mode = 'bonds'
                        else:                    
                            mode = None
                    elif mode:
                        tokens = each_line.strip().split()
                        if mode == 'atoms':
                            number = int(tokens[0])
                            type = tokens[1]
                            resid = int(tokens[2])
                            resname = tokens[3]
                            name = tokens[4]
                            atom = Atom(number, type, resid, resname,\
                                        name)
                            try:
                                if tokens[6]:
                                    charge = float(tokens[6])
                                    atom.charge = charge
                                if tokens[7]:
                                    mass = float(tokens[7])
                                    atom.mass = mass
                            except IndexError:
                                pass
                            if type == 'Sh':
                                self.atoms[number-1].shell = atom
                            else:
                                self.atoms[number] = atom
                        if mode == 'bonds':
                            atomi = int(tokens[0])
                            atomj = int(tokens[1])
                            self.atoms[atomi].bonds.append(atomj)
                            self.atoms[atomj].bonds.append(atomi)
    def find_atoms_at_distance(self, index_i, distance): 
        """ returns a list of atom indices of those separated specified distance
            away from given atom of index_i. distance=1 for nearest neighbors"""
        atomi = self.atoms[index_i]
        visited = set()
        def bfs(atomi, rem_dist=distance, visited=None):
            if not visited:
                visited = set()
            if rem_dist == 1:
                return list(atomi.bonds)
            else:
                retval = set()
                q = Queue()
                q.put((atomi, rem_dist))
                while not q.empty():
                    atomi, rem_dist = q.get()
                    if rem_dist < 1:
                        break
                    elif rem_dist == 1:
                        for neighbor_index in atomi.bonds:
                            if not neighbor_index in visited:
                                retval.add(neighbor_index)
                    else:
                        rem_dist -= 1
                        for neighbor_index in atomi.bonds:
                            visited.add(neighbor_index)
                            atomj = self.atoms[neighbor_index]
                            q.put((atomj, rem_dist))
                return list(retval)
        return bfs(atomi)
    @staticmethod
    def _remove_comments(line):
        return re.sub(r';[^\n]*', r'', line)
 


class Atom(object):
    def __init__(self, number=None, type=None, resid=None, resname=None,\
                 name=None, charge=0.0, mass=None, shell=None):
        self.number = number
        self.type = type
        self.resid = resid
        self.resname = resname
        self.name = name
        self.charge = charge
        self.mass = mass
        self.shell = shell
        self.bonds = []
    def is_shell(self):
        return self.type == 'Sh'
 
