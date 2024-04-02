#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   SubTask.py
@Time    :   2021/11/16 21:03:30
@Author  :   Cao Guolin 
@Version :   1.0
@Contact :   kafeicgl@sina.com
@License :   (C)Copyright 2020-2021, Liugroup-DFT
@Desc    :   None
'''

# here put the import lib

import os
import time


def SubTask(caltype,path):
    path = os.path.abspath(path)
    initpath = os.getcwd()
    steps = ['step1','step2']
    wlog = os.path.join(initpath,'log')   ##
    for i in range(len(steps)):
        step = steps[i]
        # filepath = os.path.abspath('cal/filepath/{}/{}/'.format(caltype,step))
        filepath = os.path.abspath('commonfile/{}/{}/'.format(caltype,step))
        works = []
        for root,dirs,files in os.walk(path):
            works.append(dirs)
        works_num = works[0]
        # tasks_name = []
        # for j in range(len(works_num)):
        #     tasks_name.append('87777order{}'.format(j+1))

        if caltype == 'slab':
            for j in range(len(works_num)):
                task_name = works_num[j]
                workpath = os.path.abspath('{}/{}/slab/{}/'.format(path,task_name,step))
                if not os.path.exists(workpath):
                    os.mkdir(workpath)
                fpath = os.path.join(workpath,'OUTCAR')
                if not os.path.exists(fpath):
                    # os.system('cd {} ; cp INCAR POTCAR KPOINTS {} ; cd {}'.format(filepath,workpath,initpath))
                    os.system('cd {} ; cp INCAR POTCAR KPOINTS {} ; cp ../../vasp2.pbs {} ; cd {}'.format(filepath,workpath,workpath,initpath))
                    if step == 'step1':
                        # os.system('cd {} ; mv HEA.vasp POSCAR ; mpirun -np 20 /public/software/vasp/vasp_std > log ; touch end.txt ; cd {}'.format(workpath,initpath))
                        os.system('cd {} ; mv HEA.vasp POSCAR ; qsub vasp2.pbs ; cd {}'.format(workpath,initpath))
                    elif step == 'step2':
                        # os.system('cd {} ; cd ../step1/ ; cp CONTCAR ../step2/ ; cd ../step2/ ; mv CONTCAR POSCAR ; mpirun -np 20 /public/software/vasp/vasp_std > log ; touch end.txt ; cd {}'.format(workpath,initpath))
                        os.system('cd {} ; cd ../step1/ ; cp CONTCAR ../step2/ ; cd ../step2/ ; mv CONTCAR POSCAR ; qsub vasp2.pbs ; cd {}'.format(workpath,initpath))
                    time.sleep(10)
                    completed = False
                    while completed == False:
                        with open(fpath,'r') as fi:
                            flist = fi.read()
                        if 'Total CPU time' in flist:
                            completed = True
                            ff = open(wlog,'a')
                            print('{} {} completed'.format(task_name,step),file=ff)
                            ff.close()
                        else:
                            time.sleep(10)
                        # if os.path.exists(fpath):
                        #     completed = True
                        # else:
                        #     time.sleep(10)

        elif caltype == 'adso':
        # elif caltype == 'adsoh':
            for j in range(len(works_num)):
                task_name = works_num[j]
                # sites = ['hollow','bridge','top']
                sites = ['fixtopO']
                # sites = ['fixtopO','fixbridgeO','fixhollowO','fixfcchollowO']
                # sites = ['hollowOH','bridgeOH','topOH']
                for k in range(len(sites)):
                    site = sites[k]
                    sitepath = os.path.abspath('{}/{}/ads/{}/'.format(path,task_name,site))
                    tmp = []
                    for root,dirs,files in os.walk(sitepath):
                        tmp.append(dirs)
                    site_num = tmp[0]
                    for l in range(len(site_num)):
                        workpath = os.path.abspath('{}/{}/ads/{}/{}/{}/'.format(path,task_name,site,site_num[l],step))
                        if not os.path.exists(workpath):
                            os.mkdir(workpath)
                        fpath = os.path.join(workpath,'OUTCAR')
                        if not os.path.exists(fpath):
                            # os.system('cd {} ; cp INCAR POTCAR KPOINTS {} ; cd {}'.format(filepath,workpath,initpath))
                            os.system('cd {} ; cp INCAR POTCAR KPOINTS {} ; cp ../../vasp.pbs {} ; cd {}'.format(filepath,workpath,workpath,initpath))     ##
                            if step == 'step1':
                                # os.system('cd {} ; mpirun -np 20 /public/software/vasp/vasp_std > log ; touch end.txt ; cd {}'.format(workpath,initpath))
                                os.system('cd {} ; qsub vasp.pbs ; cd {}'.format(workpath,initpath))
                            elif step == 'step2':
                                # os.system('cd {} ; cd ../step1/ ; cp CONTCAR ../step2/ ; cd ../step2/ ; mv CONTCAR POSCAR ; mpirun -np 20 /public/software/vasp/vasp_std > log ; touch end.txt ; cd {}'.format(workpath,initpath))
                                os.system('cd {} ; cd ../step1/ ; cp CONTCAR ../step2/ ; cd ../step2/ ; mv CONTCAR POSCAR ; qsub vasp.pbs ; cd {}'.format(workpath,initpath))
                            time.sleep(10)
                            completed = False
                            while completed == False:
                                with open(fpath,'r') as fi:
                                    flist = fi.read()
                                if 'Total CPU time' in flist:
                                    completed = True
                                    ff = open(wlog,'a')
                                    print('{} {} {} {} have completed'.format(task_name,site,site_num[l],step),file=ff)
                                    ff.close()
                                else:
                                    time.sleep(10)
                                # if os.path.exists(fpath):
                                #     completed = True
                                # else:
                                #     time.sleep(10)
   
    return

path = 'cal2/'    ##
caltype = 'adso'
t = SubTask(caltype,path)
