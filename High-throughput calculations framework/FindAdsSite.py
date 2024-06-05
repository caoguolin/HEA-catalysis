# here put the import lib
import os
import random
from ase.io import read,write
from ase import Atoms
import numpy as np



'''
Traversing adsorption sites in HEA
Three steps：
1.Obtain the coordinates of adsorption sites
2.Add adsorbed molecules and add additional information such as relaxation
3.Output calculation files
'''

class FindSite(object):
    def __init__(self,poscar,path):
        '''
        Args:
            poscar : HEA stucture file path
            path   : Output file path
        '''
        self.poscar = poscar
        self.path = os.path.abspath(path)
        if not os.path.exists(self.path):
            os.mkdir(self.path)
    
    ## Obtain structural infromation of the two layers of HEA
    def GetTop(self):
        stru = read(self.poscar)
        pos = stru.get_positions()
        lat = stru.cell.array

        sortz = sorted([x[2] for x in pos],reverse=True)
        meanz = np.mean(sortz)
        for i in range(len(sortz)):
            if sortz[i] - meanz < 10:
                z_max = sortz[i]
                break
        top_pos,second_top_pos,third_top_pos = [],[],[]
        for i in range(len(pos)):
            if z_max - pos[i][2] > -1 and z_max - pos[i][2] < 1:
                top_pos.append(pos[i])
            elif z_max - pos[i][2] > 1.5 and z_max - pos[i][2] < 3:
                second_top_pos.append(pos[i])
            elif z_max - pos[i][2] > 4 and z_max - pos[i][2] < 5.5:
                third_top_pos.append(pos[i])

        return [top_pos,second_top_pos,third_top_pos]

    ## Obtain the position information of the top adsorption site
    def TopSite(self):
        top_pos = self.GetTop()[0]
        top_site = np.array(top_pos) + [0,0,1.8]

        return top_site

    ## Obtain the position information of the adsorption sites at the bridge section
    def BridgeSite(self):
        top_pos = self.GetTop()[0]
        index = []
        for i in range(len(top_pos)):
            for j in range(i+1,len(top_pos)):
                distance = ((top_pos[i][0] - top_pos[j][0])**2 + (top_pos[i][1] - top_pos[j][1])**2)**(1/2)
                if distance < 3:
                    index.append([i,j])

        bridget_site = [(top_pos[x[0]] + top_pos[x[1]])/2 + [0,0,1.5] for x in index]

        return bridget_site

    ## Obtain location information of hollow adsorption sites
    def HollowSite(self):
        second_top_pos = self.GetTop()[1]
        hollow_site = np.array(second_top_pos) + [0,0,3.8]

        return hollow_site
    
    def fccHollowSite(self):
        third_top_pos = self.GetTop()[2]
        fcchollow_site = np.array(third_top_pos) + [0,0,6.2]

        return fcchollow_site

    ## Add relaxation and other additional information to the file
    def AddTF(self,poscarpath,slab):
        structure = read(poscarpath)
        pos = structure.get_positions()
        z = [x[2] for x in pos]
        # Oz = sorted()
        zz = min(z)
        zzz = max(z)
        TF = []
        for i in range(len(z)):
            if z[i] == zzz:
                TF.append('F F T')
            elif z[i] - zz < slab[2]*slab[0] or z[i] -zz > (slab[1]+2)*slab[0]:
                TF.append('F F F')
            else:
                TF.append('T T T')
        
        with open(poscarpath,'r') as fi:
            flist = fi.readlines()
            coord = flist[8:]
            lattice_list = flist[:8]
            lattice_list.insert(7,'Selective dynamics\n')
            coord_flist = []
            for i in range(len(TF)):
                coord_flist.append(coord[i].replace('\n','     {} \n'.format(TF[i])))
            newflist = lattice_list + coord_flist
        with open(poscarpath,'w') as fi:
            for i in range(len(newflist)):
                fi.writelines(newflist[i])

    ## Create and output calculation files
    def MakeFile(self,site,slab):
        if site == 'fixtopO':
            pos_site = self.TopSite()
        elif site == 'fixbridgeO':
            pos_site = self.BridgeSite()
        elif site == 'fixhollowO':
            pos_site = self.HollowSite()
        elif site == 'fixfcchollowO':
            pos_site = self.fccHollowSite()
        
        site_path = os.path.join(self.path,'{}/'.format(site))
        if not os.path.exists(site_path):
            os.mkdir(site_path)        
            for i in range(len(pos_site)):
                os.mkdir(os.path.join(site_path,'{}/'.format(i)))
                fpath = os.path.join(site_path,'{}/step1/'.format(i))
                os.mkdir(fpath)
                poscar_path = os.path.join(fpath,'POSCAR')
                stru = read(self.poscar)
                atom = Atoms('O',[pos_site[i]])
                stru.extend(atom.copy())
                write(poscar_path,stru)
                self.AddTF(poscar_path,slab)


## test
slab = [2,4,2]
path = 'cal'
for i in range(10,20):
    poscar = '{}/87777order{}/slab/step2/CONTCAR'.format(path,i+1)
    zpath = '{}/87777order{}/ads/'.format(path,i+1)
    t = FindSite(poscar,zpath)
    t.MakeFile('fixtopO',slab)
    t.MakeFile('fixhollowO',slab)
    t.MakeFile('fixbridgeO',slab)
    t.MakeFile('fixfcchollowO',slab)

