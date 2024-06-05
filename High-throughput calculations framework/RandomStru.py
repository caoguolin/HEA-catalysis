# here put the import lib
import ase
import random
import numpy as np
from ase.io import read,write
from pymatgen import Element
import os

ele = ['Ag','Ir','Pd','Pt','Ru']   
atom_num = [8,8,8,8,7]             
template_path = 'temp/HEA.vasp'  

'''
Automated modeling of HEAs
three steps：
1.Analyze the atomic positions in the original structrual file
2.Randomly shuffle the order of atoms and obtain a set of index values
3.Replace atomic types based on index values and output HEA structure files
'''

def RandomStru(order,ele,atom_num,template_path,sample_path):
    '''
    Args:
        order         : Random seed
        ele           : element type
        atom_num      : atom number
        template_path : teplate file path
        sample_path   : output file path
    Return:
        None
    '''
    ## Obtain atomic sequences based on random seeds
    seed = order * 11
    atoms_num = sum(atom_num)
    ele_num = [Element[x].number for x in ele]
    random.seed(seed)
    atoms = list(range(0,atoms_num))
    first_index = sorted(random.sample(atoms,atom_num[0]))

    for i in range(len(first_index)):
        atoms.remove(first_index[i])

    ## Obtain original file structure information
    stru = read(template_path)
    position = stru.get_positions()
    lattice = stru.cell.array
    number = stru.numbers

    ## Filter and replace atoms
    first_pos,other_pos = [],[]
    for i in range(len(position)):
        if i in first_index:
            first_pos.append(position[i])
        else:
            other_pos.append(position[i])
    
    ## Replace atoms
    random.shuffle(other_pos)
    newpos = np.array(first_pos + other_pos)
    atom_index = [[ele_num[0]]*atom_num[0],[ele_num[1]]*atom_num[1],[ele_num[2]]*atom_num[2],[ele_num[3]]*atom_num[3],[ele_num[4]]*atom_num[4]]
    new_num = [x for xx in atom_index for x in xx]
    newpos[:,2] += 7

    stru.positions = newpos
    stru.numbers = new_num

    ## Output HEA structure files
    if not os.path.exists(sample_path):
        os.mkdir(sample_path)
        slab = os.path.join(sample_path,'slab/')
        os.mkdir(slab)
        step1 = os.path.join(slab,'step1')
        step2 = os.path.join(slab,'step2')
        os.mkdir(step1)
        os.mkdir(step2)
        # ads = os.path.join(sample_path)
        write(os.path.join(step1,'HEA.vasp'),stru)
    # write('temp/HEA.vasp',stru)

    return

## test
for i in range(0,1):
    order = i+1
    sample_path = 'temp/87777order{}'.format(order)
    RandomStru(order,ele,atom_num,template_path,sample_path)
