import argparse
from molecule import *
import re
import os

# generates tables for custom 1-4 bonds which effectively take into account
# scaled down nonbonded 1-4 interactions.
# Given original index, itp, and table files, this code finds appropriate 1-4
# pairs, read non-bonded interactions between those pairs, and create
# tabulated bonds for these pairs using scale factor and corresponding 
# non-bonded potentials.

new_itp_file = 'bmim_new.itp'
index_file = 'index.ndx'
itp_file = 'bmim.itp'
table_preffix = 'table'
scale_factor = 0.5

E_CONVERSION = 138.935

def _check_nrexcl(itp_file):
    with open(itp_file) as f:
        contents = Molecule._remove_comments(f.read())
        m = re.search(r'\[\s*moleculetype\s*\]\s*\w+\s+(\d+)', contents)
        if m:
            if m.group(1) != '3':
                raise ValueError('ERROR! nrexcl is not set to 3!!!')
        else:
            raise ValueError('ERROR! cannot read nrexcl from ' + itp_file )


def build14bonds(itp_file, table_preffix, index_file, scale_factor, new_itp_file):

    warning_msg = '#####################################################################\n'+\
                  '#                                                                   #\n'+\
                  '#  Build14bonds code generates tables for custom 1-4 bonds          #\n'+\
                  '# which effectively takes into account scaled down nonbonded 1-4    #\n'+\
                  '# interactions. Given original index, itp, and table files,         #\n'+\
                  '# this code finds appropriate 1-4 pairs, read non-bonded            #\n'+\
                  '# interactions between those pairs, and creates tabulated bonds     #\n'+\
                  '# for these pairs using scale factors and appropriate non-bonded    #\n'+\
                  '# potentials.                                                       #\n'+\
                  '#                                                                   #\n'+\
                  '# to see description, type : python build14bonds.py --help          #\n'+\
                  '#                                                                   #\n'+\
                  '# WARNING! : 1) ENSURE THAT ALL C6 AND C12 VALUES ASSIGNED TO ATOM  #\n'+\
                  '#               TYPES USING TABLULATED NON-BONDED POTENTIALS ARE    #\n'+\
                  '#               SET TO 1.000.                                       #\n'+\
                  '#            2) CHARGES SHOULD BE EXPLICITLY DEFINED IN THE         #\n'+\
                  '#               PROVIDED ITP FILE.                                  #\n'+\
                  '#            3) BASE ITP FILE SHOULD EXCLUDE 1-4 INTERACTIONS BY    #\n'+\
                  '#               SETTING "nrexcl : 3" AND EXPLICIT EXCLUSIONS.       #\n'+\
                  '#                                                                   #\n'+\
                  '####################################################################\n'
    print(warning_msg)

    _check_nrexcl(itp_file)

    molecule = Molecule()
    molecule.readitp(itp_file)

    # assign appropriate interaction type for each atom type
    # using index file
    interaction_types = {}
    with open(index_file) as f:
        section = None
        for each_line in f:
            m = re.search(r'\[\s*([^\s]+)\s*\]', each_line)
            if m:
                section = m.group(1)
            else:
                tokens = each_line.strip().split()
                for each_token in tokens:
                    interaction_types[int(each_token)] = section

    # create list of 1-4 nonshell pairs(pairs) and 1-4 shell pairs(pairs_wshell)
    # note that shell pairs don't have VDW interaction
    pairs = []
    pairs_wshell = []
    atoms = molecule.atoms
    for index_i in atoms:
        lst = molecule.find_atoms_at_distance(index_i, 3)
        for index_j in lst:
            if index_i < index_j:
                pairs.append( (index_i, index_j) )
                if atoms[index_j].shell:
                    index_j_shell = atoms[index_j].shell.number
                    pairs_wshell.append( (index_i, index_j_shell) )
                if atoms[index_i].shell:
                    index_i_shell = atoms[index_i].shell.number
                    pairs_wshell.append( (index_i_shell, index_j) )
                    if atoms[index_j].shell:
                        index_j_shell = atoms[index_j].shell.number
                        pairs_wshell.append( (index_i_shell, index_j_shell) )

    # load copmlete pair information(potential is uniquely determined for given
    # interaction type of each atom and their charges)
    # e.g. pair_info dictionary maps pair of atom indices i and j to unique pair 
    # information : i's type, j's type, and q_i*q_j.      
    pair_info = {}
    pair_wshell_info = {}
    for each_pair in pairs:
        index_i = each_pair[0]
        index_j = each_pair[1]
        q_i = atoms[index_i].charge
        q_j = atoms[index_j].charge
        q_ij = q_i*q_j
        type_i = interaction_types[index_i]
        type_j = interaction_types[index_j]
        pair_info[each_pair] = (type_i, type_j, q_ij)
    for each_pair in pairs_wshell:
        index_i = each_pair[0]
        index_j = each_pair[1]
        atom_i = molecule.find_atom(index_i)
        atom_j = molecule.find_atom(index_j)
        q_i = atom_i.charge
        q_j = atom_j.charge
        q_ij = q_i*q_j
        type_i = interaction_types[index_i]
        type_j = interaction_types[index_j]
        pair_wshell_info[each_pair] = (type_i, type_j, q_ij)

    
    # developed_bonds_index dictionary maps unique pair information to corresponding
    # custom bond index. If no bond index is assigned yet, create a custom bond with
    # new bond index and assign the index to developed_bonds_index.
    developed_bonds_index = {}
    developed_bonds_index_wshell = {}
    for each_pair in pairs:
        developed_bonds_index[pair_info[each_pair]] = None
    for each_pair in pairs_wshell:
        developed_bonds_index_wshell[pair_wshell_info[each_pair]] = None

    bond_index = 0
    for each_pair in pairs:
        if not developed_bonds_index[pair_info[each_pair]]:
            type_i = pair_info[each_pair][0]
            type_j = pair_info[each_pair][1]
            q_ij = pair_info[each_pair][2]
            build_tableb(table_preffix, bond_index, type_i, type_j, q_ij, False)
            developed_bonds_index[pair_info[each_pair]] = bond_index
            bond_index += 1 

    for each_pair in pairs_wshell:
        if not developed_bonds_index_wshell[pair_wshell_info[each_pair]]:
            type_i = pair_wshell_info[each_pair][0]
            type_j = pair_wshell_info[each_pair][1]
            q_ij = pair_wshell_info[each_pair][2]
            build_tableb(table_preffix, bond_index, type_i, type_j, q_ij, True)
            developed_bonds_index_wshell[pair_wshell_info[each_pair]] = bond_index
            bond_index += 1 


    # add list of 1-4 pairs in [ bonds ] section and write them to new_itp_file
    with open(itp_file) as fin:  
        with open(new_itp_file, 'w') as fout:
            header = '; This itp file is created using build14bonds.py code with\n'+\
                     '; base itp file : ' + itp_file + '\n'+\
                     '; 1-4 custom bonds are created in [ bonds ] section\n'
            fout.write(header)
            section = None
            for each_line in fin:
                m = re.search(r'\[\s*([^\s]+)\s*\]', each_line)
                if m:
                    section = m.group(1)
                    if section == 'bonds':
                        bonds14 = '; custom bonds account for scaled non-bonded 1-4 interaction\n'
                        bonds14 += ';%4s %5s %5s %5s %5s\n'%('i', 'j', 'func', 'bond_idx',\
                                  'scale_factor')
                        for each_pair in pairs:
                            i = each_pair[0]
                            j = each_pair[1]
                            bond_idx = developed_bonds_index[pair_info[each_pair]]
                            bonds14 += '%5d %5d %5d %5d %5.2f\n'\
                                       %(i, j, 9, bond_idx, scale_factor)
                        for each_pair in pairs_wshell:
                            i = each_pair[0]
                            j = each_pair[1]
                            bond_idx = developed_bonds_index_wshell[pair_wshell_info[each_pair]]
                            bonds14 += '%5d %5d %5d %5d %5.2f\n'\
                                       %(i, j, 9, bond_idx, scale_factor)
                        bonds14 += '; normal harmonic bonds below..\n'
                        each_line += bonds14 
                fout.write(each_line)                     
                       
 
    
def build_tableb(input_table, bond_index, type_i, type_j, q_ij, shell_type=False):
    input_file = input_table + '_' + type_i + '_' + type_j +'.xvg'
    if not os.path.isfile('./' + input_file):
        input_file = input_table + '_' + type_j + '_' + type_i + '.xvg'
    if not os.path.isfile('./' + input_file):
        raise ValueError('No matching table file found!')
    output_file = input_table + '_b' + str(bond_index) + '.xvg'
    with open(input_file) as fin:
        with open(output_file, 'w') as fout: 
            for each_line in fin:
                if not '#' in each_line:
                    tokens = each_line.strip().split()
                    distance = '%11.10e'%(float(tokens[0]))
                    if not shell_type:
                        new_energy = '%11.10e'%((q_ij*E_CONVERSION*float(tokens[1]) +\
                                         float(tokens[3]) + float(tokens[5]) ))
                        new_force = '%11.10e'%((q_ij*E_CONVERSION*float(tokens[2]) +\
                                         float(tokens[4]) + float(tokens[6])) )
                    else:
                        new_energy = '%11.10e'%(q_ij*E_CONVERSION*float(tokens[1]))
                        new_force = '%11.10e'%(q_ij*E_CONVERSION*float(tokens[2]))
                    new_line = '%s %s %s\n'%(distance, new_energy, new_force)
                    fout.write(new_line)

 

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--itp', '-i', \
                        help='original itp file containing molecule information')
    parser.add_argument('--table', '-t',\
                        help='prefix of table files used for non-bonded interactions')
    parser.add_argument('--index', '-n',\
                        help='gromacs ndx file containing particle indices assigned to'+\
                             ' particle types')
    parser.add_argument('--scale', '-s', type=float,\
                        help='scale down factor for 1-4 interactions.')
    parser.add_argument('-output', '-o',\
                        help='file name for new itp file including 1-4 custom bonds')
    arg = parser.parse_args()
    build14bonds(arg.itp, arg.table, arg.index, arg.scale, arg.output)
