

class Molecule(object):
    def __init__(self, name=None):
        self.name = name
        self.atoms = []
    def readitp(self, itp_file):
        """ load molecule information from given itp_file """




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
 
