# here put the import lib
import ase
import random
import numpy as np
import os 
import math


## Parameter
Pure_Fillings = {'Ag':0.978,'Ir':0.819,'Pd':0.913,'Pt':0.875,'Ru':0.802,'Au':0.962,'Rh':0.826,'Os':0.776,'Cu':0.967}
Pure_antiFillings = {'Ag':0.132,'Ir':0.181,'Pd':0.087,'Pt':0.125,'Ru':0.198,'Au':0.038,'Rh':0.174,'Os':0.224,'Cu':0.033}
Pure_Xs = {'Ag':1.93,'Ir':2.20,'Pd':2.20,'Pt':2.28,'Ru':2.2,'Au':2.54,'Rh':2.28,'Os':2.2,'Cu':1.9}

##设计100*100的HEA slab，进行元素填充占据

##根据index转换为坐标
def IndextoCoord(indexs):
    coords = []
    for i in range(len(indexs)):
        cur_index = indexs[i]
        if (cur_index + 1)%100 == 0:
            cur_x = (cur_index + 1)//100 -1
            cur_y = 99
        else:
            cur_x = (cur_index + 1)//100
            cur_y = (cur_index + 1)%100 - 1
        coords.append([cur_x,cur_y])
    coords = np.array(coords)

    return coords

##判断相邻度
def NeighborDegree(coords):
    closed_atom_num = 0
    for i in range(len(coords)):
        cur_coord = coords[i]
        if cur_coord[0]%2 == 0:
            neighbor_coord = cur_coord + [[0,-1],[0,1],[-1,0],[-1,-1],[1,0],[1,-1]]
        else:
            neighbor_coord = cur_coord + [[0,-1],[0,1],[-1,0],[-1,1],[1,0],[1,1]]
        neighbor_num = 0
        for j in range(i+1,len(coords)):
            point_coord = coords[j]
            if point_coord.tolist() in neighbor_coord.tolist():
                neighbor_num += 1
            if j == 150:
                break
        if neighbor_num >= 4:
            closed_atom_num += 1
    
    return closed_atom_num



# ##判断结构是否合理
def GenerateSlab(ratios,seed):
    all_index = list(range(0,10000))
    whether_reasonable = 0
    random.seed(seed)
        ##先取第一种元素
    ele_first_index = sorted(random.sample(list(range(0,10000)),int(ratios[0]/sum(ratios)*10000)))
    # if NeighborDegree(IndextoCoord(ele_first_index)) == 0:
    #     whether_reasonable += 1

    ##取第二种元素
    all_index_first = list(set(all_index)-set(ele_first_index))
    ele_sec_index = sorted(random.sample(all_index_first,int(ratios[1]/sum(ratios)*10000)))
    # if NeighborDegree(IndextoCoord(ele_sec_index)) == 0:
    #     whether_reasonable += 1
    
    ##取第三种元素
    all_index_sec = list(set(all_index_first)-set(ele_sec_index))
    ele_third_index = sorted(random.sample(all_index_sec,int(ratios[2]/sum(ratios)*10000)))
    # if NeighborDegree(IndextoCoord(ele_third_index)) == 0:
    #     whether_reasonable += 1

    ##取第四种元素
    all_index_third = list(set(all_index_sec)-set(ele_third_index))
    ele_fourth_index = sorted(random.sample(all_index_third,int(ratios[3]/sum(ratios)*10000)))
    # if NeighborDegree(IndextoCoord(ele_fourth_index)) == 0:
    #     whether_reasonable += 1

    # all_index_fourth = list(set(all_index_third)-set(ele_fourth_index))
    # ele_fifth_index = sorted(random.sample(all_index_fourth,int(ratios[4]/sum(ratios)*10000)))

    # all_index_fifth = list(set(all_index_fourth)-set(ele_fifth_index))
    # ele_sixth_index = sorted(random.sample(all_index_fifth,int(ratios[5]/sum(ratios)*10000)))
    
    # all_index_sixth = list(set(all_index_fifth)-set(ele_sixth_index))
    # ele_seventh_index = sorted(random.sample(all_index_sixth,int(ratios[6]/sum(ratios)*10000)))

    # ele_eighth_index = sorted(list(set(all_index_sixth)-set(ele_seventh_index)))

    ##剩下的为第五种元素
    ele_fifth_index = sorted(list(set(all_index_third)-set(ele_fourth_index)))
    # if NeighborDegree(IndextoCoord(ele_fifth_index)) == 0:
    #     whether_reasonable += 1
        
        # seed = seed + 1 
    ##10W次循环来检查结构合理性,运行时间过长
    return [ele_first_index,ele_sec_index,ele_third_index,ele_fourth_index,ele_fifth_index],seed,whether_reasonable
    # return [ele_first_index,ele_sec_index,ele_third_index,ele_fourth_index,ele_fifth_index,ele_sixth_index,ele_seventh_index,ele_eighth_index],seed


##定义top位点和fcc位点的局域环境坐标组合

##对于top位点，包含在位1，表面近邻6，次表面近邻3
##当在位x为奇数行时(计算机语言的奇数):
##[x,y]  x为行，y为列
##[x,y-1];[x,y+1];[x-1,y];[x-1,y+1];[x+1,y];[x+1,y+1]
##[x-1,y];[x,y-1];[x,y]   ##次
##而当在位x为偶数行时:
##[x,y]
##[x,y-1];[x,y+1];[x-1,y-1];[x-1,y];[x+1,y-1];[x+1,y]
##[x-1,y-1];[x,y-1];[x,y]   ##次

##而对于fcc位点，包含在位3，表面近邻9，次表面近邻6
##fcc为正三角
##[x,y],上顶点
##当x为奇数：
##[x,y];[x+1,y];[x+1,y+1]
##[x-1,y];[x-1,y+1];[x,y-1];[x,y+1];[x+1,y-1];[x+1,y+2];[x+2,y-1];[x+2,y];[x+2,y+1]
##[x-1,y];[x,y-1];[x,y];[x+1,y-1];[x+1,y];[x+1,y+1]  ##次
##当x为偶数：
##[x,y];[x+1,y-1];[x+1,y]
##[x-1,y-1];[x-1,y];[x,y-1];[x,y+1];[x+1,y-2];[x+1,y+1];[x+2,y-1];[x+2,y];[x+2,y+1]
##[x-1,y-1];[x,y-1];[x,y];[x+1,y-2];[x+1,y-1];[x+1,y]   ##次


def GetLocalSiteCoord(index_site_coord,ads_type):
    ##根据顶点索引来得到局域环境的所有原子的坐标
    x,y = index_site_coord[0],index_site_coord[1]
    if ads_type == 'top':
        if x%2 == 0:
            sites = [[x,y]]
            slab_neighbor = [[x,y-1],[x,y+1],[x-1,y-1],[x-1,y],[x+1,y-1],[x+1,y]]
            subslab_neighbor = [[x-1,y-1],[x,y-1],[x,y]]
        else:
            sites = [[x,y]]
            slab_neighbor = [[x,y-1],[x,y+1],[x-1,y],[x-1,y+1],[x+1,y],[x+1,y+1]]
            subslab_neighbor = [[x-1,y],[x,y-1],[x,y]]
    elif ads_type == 'fcc':
        if x%2 == 0:
            sites = [[x,y],[x+1,y-1],[x+1,y]]
            slab_neighbor = [[x-1,y-1],[x-1,y],[x,y-1],[x,y+1],[x+1,y-2],[x+1,y+1],[x+2,y-1],[x+2,y],[x+2,y+1]]
            subslab_neighbor = [[x-1,y-1],[x,y-1],[x,y],[x+1,y-2],[x+1,y-1],[x+1,y]]
        else:
            sites = [[x,y],[x+1,y],[x+1,y+1]]
            slab_neighbor = [[x-1,y],[x-1,y+1],[x,y-1],[x,y+1],[x+1,y-1],[x+1,y+2],[x+2,y-1],[x+2,y],[x+2,y+1]]
            subslab_neighbor = [[x-1,y],[x,y-1],[x,y],[x+1,y-1],[x+1,y],[x+1,y+1]]
    
    return sites,slab_neighbor,subslab_neighbor

def GetIndexSiteCoord(ads_type):
    ##获得表面所有的索引顶点坐标
    index_sites = []
    if ads_type == 'top':
        for i in range(1,99):
            for j in range(1,99):
                index_sites.append([i,j])
    elif ads_type == 'fcc':
        for i in range(2,98):
            for j in range(2,98):
                index_sites.append([i,j])

    return index_sites

def CoordtoEle(slab,coord,eles):
    ##根据原子序列坐标来得知其元素信息
    x,y = coord[0],coord[1]
    index = x*100 + y
    if index in slab[0]:
        ele = eles[0]
    elif index in slab[1]:
        ele = eles[1]
    elif index in slab[2]:
        ele = eles[2]
    elif index in slab[3]:
        ele = eles[3]
    elif index in slab[4]:
        ele = eles[4]
    ##
    elif index in slab[5]:
        ele = eles[5]
    elif index in slab[6]:
        ele = eles[6]
    elif index in slab[7]:
        ele = eles[7]

    return ele

def GetLocalInfor(slab1,slab2,index_site,ads_type,eles):
    ##slab1为表面，slab2为次表面
    ##根据顶点求得某一局域的元素信息和索引信息，索引信息用于后续的筛选
    if ads_type == 'top':
        local_coords = GetLocalSiteCoord(index_site,'top')
        point_index = [local_coords[0][0][0]*100 + local_coords[0][0][1]]
        point_envir = [CoordtoEle(slab1,local_coords[0][0],eles)]
        slab_neighbor_index = [local_coords[1][0][0]*100 + local_coords[1][0][1],
                                local_coords[1][1][0]*100 + local_coords[1][1][1],
                                local_coords[1][2][0]*100 + local_coords[1][2][1],
                                local_coords[1][3][0]*100 + local_coords[1][3][1],
                                local_coords[1][4][0]*100 + local_coords[1][4][1],
                                local_coords[1][5][0]*100 + local_coords[1][5][1]]
        slab_neighbor_envir = [CoordtoEle(slab1,local_coords[1][0],eles),
                                CoordtoEle(slab1,local_coords[1][1],eles),
                                CoordtoEle(slab1,local_coords[1][2],eles),
                                CoordtoEle(slab1,local_coords[1][3],eles),
                                CoordtoEle(slab1,local_coords[1][4],eles),
                                CoordtoEle(slab1,local_coords[1][5],eles)]
        subslab_neighbor_index = [local_coords[2][0][0]*100 + local_coords[2][0][1],
                                local_coords[2][1][0]*100 + local_coords[2][1][1],
                                local_coords[2][2][0]*100 + local_coords[2][2][1]]
        subslab_neighbor_envir = [CoordtoEle(slab2,local_coords[2][0],eles),
                                    CoordtoEle(slab2,local_coords[2][1],eles),
                                    CoordtoEle(slab2,local_coords[2][2],eles)]
    elif ads_type == 'fcc':
        local_coords = GetLocalSiteCoord(index_site,'fcc')
        point_index = [local_coords[0][0][0]*100 + local_coords[0][0][1],
                        local_coords[0][1][0]*100 + local_coords[0][1][1],
                        local_coords[0][2][0]*100 + local_coords[0][2][1]]
        point_envir = [CoordtoEle(slab1,local_coords[0][0],eles),
                        CoordtoEle(slab1,local_coords[0][1],eles),
                        CoordtoEle(slab1,local_coords[0][2],eles)]
        slab_neighbor_index = [local_coords[1][0][0]*100 + local_coords[1][0][1],
                                local_coords[1][1][0]*100 + local_coords[1][1][1],
                                local_coords[1][2][0]*100 + local_coords[1][2][1],
                                local_coords[1][3][0]*100 + local_coords[1][3][1],
                                local_coords[1][4][0]*100 + local_coords[1][4][1],
                                local_coords[1][5][0]*100 + local_coords[1][5][1],
                                local_coords[1][6][0]*100 + local_coords[1][6][1],
                                local_coords[1][7][0]*100 + local_coords[1][7][1],
                                local_coords[1][8][0]*100 + local_coords[1][8][1]]
        slab_neighbor_envir = [CoordtoEle(slab1,local_coords[1][0],eles),
                                CoordtoEle(slab1,local_coords[1][1],eles),
                                CoordtoEle(slab1,local_coords[1][2],eles),
                                CoordtoEle(slab1,local_coords[1][3],eles),
                                CoordtoEle(slab1,local_coords[1][4],eles),
                                CoordtoEle(slab1,local_coords[1][5],eles),
                                CoordtoEle(slab1,local_coords[1][6],eles),
                                CoordtoEle(slab1,local_coords[1][7],eles),
                                CoordtoEle(slab1,local_coords[1][8],eles)]
        subslab_neighbor_index = [local_coords[2][0][0]*100 + local_coords[2][0][1],
                                local_coords[2][1][0]*100 + local_coords[2][1][1],
                                local_coords[2][2][0]*100 + local_coords[2][2][1],
                                local_coords[2][3][0]*100 + local_coords[2][3][1],
                                local_coords[2][4][0]*100 + local_coords[2][4][1],
                                local_coords[2][5][0]*100 + local_coords[2][5][1]]
        subslab_neighbor_envir = [CoordtoEle(slab2,local_coords[2][0],eles),
                                    CoordtoEle(slab2,local_coords[2][1],eles),
                                    CoordtoEle(slab2,local_coords[2][2],eles),
                                    CoordtoEle(slab2,local_coords[2][3],eles),
                                    CoordtoEle(slab2,local_coords[2][4],eles),
                                    CoordtoEle(slab2,local_coords[2][5],eles)]

    return [point_index,point_envir],[slab_neighbor_index,slab_neighbor_envir],[subslab_neighbor_index,subslab_neighbor_envir]

def ModeltoEn(point_envir,slab_envir,subslab_envir,ads,ads_type):
    ##描述符模型预测吸附能
    if ads == 'OH' and ads_type == 'top':
        point_xs = [Pure_Xs[point_envir[0]]]
        if point_envir[0] == 'Ag':
            point_fill = [Pure_Fillings[point_envir[0]]*1.035]
        elif point_envir[0] == 'Au':
            point_fill = [Pure_Fillings[point_envir[0]]*0.915]
        else:
            point_fill = [Pure_Fillings[point_envir[0]]]
        slab_xs = [Pure_Xs[x] for x in slab_envir]
        subslab_xs = [Pure_Xs[x] for x in subslab_envir]
        X = sum(point_xs) + (sum(slab_xs) + sum(subslab_xs))/9 + (5*point_fill[0])
        En = 1.62*X - 13.240
    elif ads == 'O' and ads_type == 'fcc':
        point_xs = [Pure_Xs[x] for x in point_envir]
        point_fill = [Pure_Fillings[x] for x in point_envir]
        point_antifill = [Pure_antiFillings[x] for x in point_envir]
        slab_xs = [Pure_Xs[x] for x in slab_envir]
        subslab_xs = [Pure_Xs[x] for x in subslab_envir]
        W = [point_antifill[0]/sum(point_antifill),point_antifill[1]/sum(point_antifill),point_antifill[2]/sum(point_antifill)]
        X = sum(point_xs)/3 + (sum(slab_xs)+sum(subslab_xs))/15 + (5*(point_fill[0]*W[0]+point_fill[1]*W[1]+point_fill[2]*W[2]))
        En = 3.45*X - 28.754

    return En

def GetSurfaceEn(slab1,slab2,ads,eles):
    ##获取一个表面的所有吸附能
    Ens = []
    Indexs = []
    if ads == 'O':
        indexsites = GetIndexSiteCoord('fcc')
        for i in range(len(indexsites)):
            cur_indexsite = indexsites[i]
            point_infors,slab_infors,subslab_infors = GetLocalInfor(slab1,slab2,cur_indexsite,'fcc',eles)
            point_envir,slab_envir,subslab_envir = point_infors[1],slab_infors[1],subslab_infors[1]
            Ens.append(ModeltoEn(point_envir,slab_envir,subslab_envir,'O','fcc'))
            Indexs.append(point_infors[0])
    elif ads == 'OH':
        indexsites = GetIndexSiteCoord('top')
        for i in range(len(indexsites)):
            cur_indexsite = indexsites[i]
            point_infors,slab_infors,subslab_infors = GetLocalInfor(slab1,slab2,cur_indexsite,'top',eles)
            point_envir,slab_envir,subslab_envir = point_infors[1],slab_infors[1],subslab_infors[1]
            Ens.append(ModeltoEn(point_envir,slab_envir,subslab_envir,'OH','top'))
            Indexs.append(point_infors[0])

    return Ens,Indexs

def SortFilter(O_ens,O_indexs,OH_ens,OH_indexs):
    ##对O和OH进行排序和筛选
    O_ens,OH_ens = np.array(O_ens),np.array(OH_ens)
    Ens_O,Ens_OH = np.sort(O_ens),np.sort(OH_ens)
    sorted_EnO_indexs,sorted_EnOH_indexs = np.argsort(O_ens),np.argsort(OH_ens)
    real_O_ens,real_OH_ens,O_filter_indexs,OH_filter_indexs = [],[],[],[]
    real_O_ens.append(Ens_O[0])
    OH_filter_indexs = OH_filter_indexs + O_indexs[sorted_EnO_indexs[0]]
    insert_O,insert_OH = 1,0
    while insert_O < len(Ens_O)-1 and insert_OH < len(Ens_OH)-1:
        while len(set(set(OH_indexs[sorted_EnOH_indexs[insert_OH]])&set(OH_filter_indexs))) > 0:
            if insert_OH < len(Ens_OH)-1:
                insert_OH += 1
            else:
                insert_OH += 1
                break
        if insert_OH < len(Ens_OH)-1:
            real_OH_ens.append(Ens_OH[insert_OH])
            O_filter_indexs = O_filter_indexs + OH_indexs[sorted_EnOH_indexs[insert_OH]]
            insert_OH += 1
        if insert_OH < len(Ens_OH)-1 and insert_O < len(Ens_O)-1:
            while len(set(set(O_indexs[sorted_EnO_indexs[insert_O]])&set(O_filter_indexs))) > 0:
                if insert_O < len(Ens_O)-1:
                    insert_O += 1
                else:
                    insert_O += 1
                    break
            if insert_O < len(Ens_O)-1:
                real_O_ens.append(Ens_O[insert_O])
                OH_filter_indexs = OH_filter_indexs + O_indexs[sorted_EnO_indexs[insert_O]]
                insert_O += 1
    if insert_O == len(Ens_O)-1 and insert_OH < len(Ens_OH)-1:
        while insert_OH < len(Ens_OH)-1:
            while len(set(set(OH_indexs[sorted_EnOH_indexs[insert_OH]])&set(OH_filter_indexs))) > 0:
                if insert_OH < len(Ens_OH)-1:
                    insert_OH += 1
                else:
                    insert_OH += 1
                    break
            if insert_OH < len(Ens_OH)-1:
                real_OH_ens.append(Ens_OH[insert_OH])
                insert_OH += 1
    elif insert_OH == len(Ens_OH)-1 and insert_O < len(Ens_O)-1:
        while insert_O < len(Ens_O)-1:
            while len(set(set(O_indexs[sorted_EnO_indexs[insert_O]])&set(O_filter_indexs))) > 0:
                if insert_O < len(Ens_O)-1:
                    insert_O += 1
                else:
                    insert_O += 1
                    break
            if insert_O < len(Ens_O)-1:
                real_O_ens.append(Ens_O[insert_O])
                insert_O += 1
    
    return real_O_ens,real_OH_ens

def GetSurfaceCurrent(O_ens,OH_ens,Pt_O_Eads=1.439,Pt_OH_Eads=0.522,U=0.82):
    ##求表面的电流密度以衡量其活性
    all_j = []
    for i in range(len(O_ens)):
        cur_jk = -math.e**((U-0.86-abs(O_ens[i]-(Pt_O_Eads+0.2)))*(1.6e-19)/((1.38e-23)*(300)))
        cur_j = (-cur_jk)/(cur_jk-1)
        all_j.append(cur_j) 
    for i in range(len(OH_ens)):
        cur_jk = -math.e**((U-0.86-abs(OH_ens[i]-(Pt_OH_Eads+0.1)))*(1.6e-19)/((1.38e-23)*(300)))
        cur_j = (-cur_jk)/(cur_jk-1) 
        all_j.append(cur_j)

    all_j
    res = sum(all_j)/10000
    return res

def GenerateComponent(ratio_step):
    coms = []
    fold = int(100/ratio_step)
    for i in range(0,fold+1):
        com_A = i*ratio_step
        for j in range(0,fold+1):
            com_B = j*ratio_step
            if com_A + com_B < 101:
                com_C = 100-com_A-com_B
                coms.append([com_A,com_B,com_C])
    
    return coms

##Main
##设定表面生成条件
# eles = ['Ag','Ir','Pd','Pt','Ru']
# eles_ratio = [20,20,20,20,20]


##生成表面
# slab1 = GenerateSlab(eles_ratio,100)[0]
# slab2 = GenerateSlab(eles_ratio,150)[0]

##生成表面所有吸附能信息
# O_Ens = GetSurfaceEn(slab1,slab2,'O')
# OH_Ens =  GetSurfaceEn(slab1,slab2,'OH')

##对表面吸附情况进行筛选
# filtered_Ens = SortFilter(O_Ens[0],O_Ens[1],OH_Ens[0],OH_Ens[1])

##求表面电流密度
# J = GetSurfaceCurrent(filtered_Ens[0],filtered_Ens[1])

#生成系列数据并保存本地
# eles = ['Ag','Ir','Pd','Pt','Ru']

# ele_list = [['Pd','Ir','Ag','Pt','Ru'],
#             ['Pt','Ru','Ag','Ir','Pd'],
#             ['Au','Ir','Pd','Pt','Ru'],
#             ['Au','Pd','Ir','Pt','Ru'],
#             ['Au','Pt','Ir','Pd','Ru'],
#             ['Au','Ru','Ir','Pd','Pt'],
#             ['Cu','Ir','Pd','Pt','Ru'],
#             ['Cu','Pd','Ir','Pt','Ru'],
#             ['Cu','Pt','Ir','Pd','Ru'],
#             ['Cu','Ru','Ir','Pd','Pt'],
#             ['Pt','Rh','Ag','Pd','Ru'],
#             ['Pd','Rh','Ag','Pt','Ru']]

# ele_list = [['Au','Ag','Pd','Pt','Ir','Ru','Rh','Os']]
ele_list = [['Pd','Ru','Ag','Ir','Pt']]

for k in range(len(ele_list)):
    eles = ele_list[k]
    coms = GenerateComponent(5)
    # coms = [[0,0,100],[0,100,0],[100,0,0]]
    current_densitys = []
    Ens_num = []
    for i in range(len(coms)):
        # ele_ratio_1,ele_ratio_3,ele_ratio_6 = int(coms[i][0]/2),int(coms[i][1]/3),int(coms[i][2]/3)
        # ele_ratio_4,ele_ratio_7 = int(coms[i][1]/3),int(coms[i][2]/3)
        # ele_ratio_2,ele_ratio_5,ele_ratio_8 = coms[i][0] - ele_ratio_1,coms[i][1]-ele_ratio_3-ele_ratio_4,coms[i][2]-ele_ratio_6-ele_ratio_7
        ele_ratio_1 = coms[i][0]
        ele_ratio_2 = coms[i][1]
        ele_ratio_3,ele_ratio_4 = int(coms[i][2]/3),int(coms[i][2]/3)
        ele_ratio_5 = coms[i][2] - ele_ratio_3 - ele_ratio_4
        # ele_ratio_2,ele_ratio_4 = int(coms[i][1]/2),int(coms[i][2]/2)
        # ele_ratio_3,ele_ratio_5 = coms[i][1]-ele_ratio_2,coms[i][2]-ele_ratio_4
        eles_ratio = [ele_ratio_1,ele_ratio_2,ele_ratio_3,ele_ratio_4,ele_ratio_5]
        # eles_ratio = [ele_ratio_1,ele_ratio_2,ele_ratio_3,ele_ratio_4,ele_ratio_5,ele_ratio_6,ele_ratio_7,ele_ratio_8]
        slab1 = GenerateSlab(eles_ratio,100)[0]
        slab2 = GenerateSlab(eles_ratio,150)[0]
        O_Ens = GetSurfaceEn(slab1,slab2,'O',eles)
        OH_Ens =  GetSurfaceEn(slab1,slab2,'OH',eles)
        filtered_Ens = SortFilter(O_Ens[0],O_Ens[1],OH_Ens[0],OH_Ens[1])
        En_num = len(filtered_Ens[0]) + len(filtered_Ens[1])
        J = GetSurfaceCurrent(filtered_Ens[0],filtered_Ens[1])

        current_densitys.append(J)
        Ens_num.append(En_num)
        
    # np.save('expernpy/coms/coms{}{}.npy'.format(eles[0],eles[1]),np.array(coms))
    np.save('expernpy/coms/comsPdRu.npy',np.array(coms))
    # np.save('expernpy/current/current{}{}.npy'.format(eles[0],eles[1]),np.array(current_densitys))
    np.save('expernpy/current/currentPdRu.npy',np.array(current_densitys))
    # np.save('npy/en_num/ennum{}{}.npy'.format(eles[0],eles[1]),np.array(Ens_num))


for k in range(len(ele_list)):
    eles = ele_list[k]
    coms = GenerateComponent(2)
    # coms = [[0,0,100],[0,100,0],[100,0,0]]
    current_densitys = []
    Ens_num = []
    for i in range(len(coms)):
        # ele_ratio_1,ele_ratio_3,ele_ratio_6 = int(coms[i][0]/2),int(coms[i][1]/3),int(coms[i][2]/3)
        # ele_ratio_4,ele_ratio_7 = int(coms[i][1]/3),int(coms[i][2]/3)
        # ele_ratio_2,ele_ratio_5,ele_ratio_8 = coms[i][0] - ele_ratio_1,coms[i][1]-ele_ratio_3-ele_ratio_4,coms[i][2]-ele_ratio_6-ele_ratio_7
        ele_ratio_1 = coms[i][0]
        ele_ratio_2 = coms[i][1]
        ele_ratio_3,ele_ratio_4 = int(coms[i][2]/3),int(coms[i][2]/3)
        ele_ratio_5 = coms[i][2] - ele_ratio_3 - ele_ratio_4
        # ele_ratio_2,ele_ratio_4 = int(coms[i][1]/2),int(coms[i][2]/2)
        # ele_ratio_3,ele_ratio_5 = coms[i][1]-ele_ratio_2,coms[i][2]-ele_ratio_4
        eles_ratio = [ele_ratio_1,ele_ratio_2,ele_ratio_3,ele_ratio_4,ele_ratio_5]
        # eles_ratio = [ele_ratio_1,ele_ratio_2,ele_ratio_3,ele_ratio_4,ele_ratio_5,ele_ratio_6,ele_ratio_7,ele_ratio_8]
        slab1 = GenerateSlab(eles_ratio,100)[0]
        slab2 = GenerateSlab(eles_ratio,150)[0]
        O_Ens = GetSurfaceEn(slab1,slab2,'O',eles)
        OH_Ens =  GetSurfaceEn(slab1,slab2,'OH',eles)
        filtered_Ens = SortFilter(O_Ens[0],O_Ens[1],OH_Ens[0],OH_Ens[1])
        En_num = len(filtered_Ens[0]) + len(filtered_Ens[1])
        J = GetSurfaceCurrent(filtered_Ens[0],filtered_Ens[1])

        current_densitys.append(J)
        Ens_num.append(En_num)
        
    # np.save('expernpy/coms/coms{}{}.npy'.format(eles[0],eles[1]),np.array(coms))
    np.save('expernpy/coms/comsPdRu2.npy',np.array(coms))
    # np.save('expernpy/current/current{}{}.npy'.format(eles[0],eles[1]),np.array(current_densitys))
    np.save('expernpy/current/currentPdRu2.npy',np.array(current_densitys))
    # np.save('npy/en_num/ennum{}{}.npy'.format(eles[0],eles[1]),np.array(Ens_num))


print()