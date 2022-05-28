### 1. Import libraries
# -----------------------------------------------------
from osgeo import gdal
import numpy as np
from bisect import bisect_right
from math import *
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.factory import get_sampling, get_crossover, get_mutation
from pymoo.factory import get_termination
from pymoo.optimize import minimize
from pymoo.core.problem import Problem
import copy
import matplotlib.pyplot as plt
import matplotlib as mlt
from swmm5.swmm5tools import SWMM5Simulation
import re
import pandas as pd
from osgeo import gdal
import os
import pickle


### 2. Functions 
# -----------------------------------------------------

'''INP EDIT'''

def writeinp(opath, dpath, formedlist, pos1, pos2, titlelist):
    '''写入inp文件'''
    file = open(opath, "r", encoding = "UTF-8" )
    content = file.read() 
    file.close()
    pos_1 = content.find(pos1)
    pos_2 = content.find(pos2)


    temp = content[pos_2:]
    content = content[:pos_1]

    content += titlelist[0]
    content += titlelist[1]
    content += titlelist[2]
    for i in range(len(formedlist)):
        content += formedlist[i]
    content += '\n'

    content += temp
    file = open(dpath, "w" ) 
    file.write(content) 
    file.close() 

    return 'Finished!'
  
def countx(X,iden):
    countx = []
    for i in range(len(X)):
        #print(np.sum(X[i] == iden))
        countx.append(np.sum(X[i] == iden))
    return np.array(countx)
  
def extractinp(file, pos1, pos2, startline=3, endline=2):
    ppos1 = content.find(pos1)
    ppos2 = content.find(pos2)
    contentnew = content[ppos1:ppos2]
    conlist = []
    context = contentnew.split('\n')

    for i in range(startline,len(context)-endline):
        conlist.append(re.split(r"[ ]+", context[i]))
    data = pd.DataFrame(conlist)
    return data

def formatinp(charnum, contentdic):
    contentdicK = list(contentdic.keys())
    Formedcontent = []
    for i in range(len(contentdic[contentdicK[0]])):
        temp = ''
        for j in range(len(contentdicK)):
            temp += f'{contentdic[contentdicK[j]][i]}'.ljust(charnum[j])
        temp += '\n '
        Formedcontent.append(temp)
    return Formedcontent

def insertLID(pos1, pos2, Formedcontent, tINP, INP2):
    file = open(tINP, "r", encoding = "UTF-8" )
    cont = file.read() 
    pos = cont.find(pos1)
    ppos2 = cont.find(pos2)

    temp = cont[ppos2:]
    cont = cont[:pos]

    cont += pos1 + "\n"
    cont += ";;Subcatchment   LID Process      Number  Area       Width      InitSat    FromImp    ToPerv     RptFile                  DrainTo         \n"
    cont += ";;-------------- ---------------- ------- ---------- ---------- ---------- ---------- ---------- ------------------------ ----------------\n"
    for i in range(len(Formedcontent)):
        cont += Formedcontent[i]
    
    cont += temp

    file = open(INP2, "w" ) 
    file.write(cont) 
    file.close() 
    return 'Finished!'
  
'''D8 process'''

def standardized8(d8array):
    d8b = []
    co = 0
    co2 = 0
    for i in d8array:
        d8b.append([])
        co2 = 0
        for j in i:
            d8d = bisect_right(tempvar1, j)
            d8b[co].append(d8d)
            co2 += 1
        co += 1
    
    d8c = np.array(d8b)
    return d8c

def to_tif(arr, tran):

    row = arr.shape[0]  
    columns = arr.shape[1]  
    dim = 1  
    driver = gdal.GetDriverByName('GTiff')
    dst_ds = driver.Create(r'D:/TEST/0727SWMM/data/standardd8.tif', columns, row, dim, gdal.GDT_UInt16) 
    dst_ds.SetGeoTransform(tran)
    dst_ds.GetRasterBand(1).WriteArray(arr)
    dst_ds.FlushCache()
    dst_ds = None

def formalizeX(res, n):
    X1 = res.X[n]
    fron = []
    for i in range(len(d8b)):
        fron.append([])
        for j in range(len(d8b[i])):
            if d8b[i][j] == 0:
                fron[i].append(-1)
            else:
                fron[i].append(X1[d8b[i][j]-1])
    X = np.array(fron)
    return X
    
def oriarray(d8array):
    idtemp = 0
    d8e = []
    for i in range(d8array.shape[0]):
        d8e.append([])
        for j in range(d8array.shape[1]):
            idtemp += 1 
            d8e[i].append(idtemp)
    d8f = np.array(d8e)
    return d8f


def desarray(oriarray, sdd8array):
    x = oriarray.shape[0]
    y = oriarray.shape[1]
    # oriarray2 = np.c_[oriarray, [ 0 * y]]
    des = []
    for i in range(x):
        des.append([])
        for j in range(y):
            try:
                if sdd8array[i][j] == 1:
                    des[i].append(oriarray[i][j+1])
                elif sdd8array[i][j] == 2:
                    des[i].append(oriarray[i+1][j+1])
                elif sdd8array[i][j] == 3:
                    des[i].append(oriarray[i+1][j])
                elif sdd8array[i][j] == 4:
                    des[i].append(oriarray[i+1][j-1])
                elif sdd8array[i][j] == 5:
                    des[i].append(oriarray[i][j-1])
                elif sdd8array[i][j] == 6:
                    des[i].append(oriarray[i-1][j-1])
                elif sdd8array[i][j] == 7:
                    des[i].append(oriarray[i-1][j])
                elif sdd8array[i][j] == 8:
                    des[i].append(oriarray[i-1][j+1])
                else:
                    des[i].append(-999)
            except:
                des[i].append(-999)
        
    desout = np.array(des)
    return desout

def subcatch(ginf, oriarray,desarray):
    siteid = oriarray
    namedid = []
    test = ginf
    siteidlist = siteid.tolist()
    oripX = test[0] 
    oripY = test[3] 
    width = test[1] 
    subcatchment = dict() 
    
    for i in range(siteid.shape[0]):
        namedid.append([])
        for j in range(siteid.shape[1]):
            temptext = 'subcatchment'+ str(siteidlist[i][j])
            namedid[i].append(temptext)
            
    for i in range(len(namedid)):
        for j in range(len(namedid[i])):
            subcatchment[namedid[i][j]] = [(j*width,i*width*(-1)),((j+1)*width,i*width*(-1)),((j+1)*width,(i+1)*width*(-1)),(j*width,(i+1)*width*(-1))]

    for i in range(len(desout)):
        for j in range(len(desout[i])):
            if desarray[i][j] == -999:
                subcatchment[namedid[i][j]].append('outfall')
            else:
                subcatchment[namedid[i][j]].append('subcatchment'+str(desarray[i][j]))
    
    return subcatchment
  

def Compile(X):
    popid[0] = 0
    R = []
    C = []
    for i in range(len(X)):
        popid[0] += 1
        temp1, temp2 = SWMM4lid(X[i])
        R.append(temp1)
        C.append(temp2)       
    GEN[0] += 1
    
    return R,C

def SWMM4lid(X, prop, width, cos, indi=4):

    '''1. extract current data'''
    Subcatchment = extractinp(content,'[SUBCATCHMENTS]','[SUBAREAS]' )  
    Lidcontrols = extractinp(content,'[LID_CONTROLS]','[LID_USAGE]' )  

    '''2.1 LIDs data'''
    Prop = prop 
    Lidnames = list(filter(None,list(Lidcontrols[0].drop_duplicates().dropna())))
    LW = width
    LidC = cos 
    EMP = len(Lidnames)


    '''2.2 LIDs info '''
    Lidprocess = []  
    Lidcost = []
    for i in range(len(X)):
        if X[i] < EMP:
            Lidprocess.append(Lidnames[X[i]])
            Lidcost.append(LidC[X[i]])
    Subname = [Subcatchment[0][i] for i in range(len(X)) if X[i] < EMP]
    Subarea = [Subcatchment[3][i] for i in range(len(X)) if X[i] < EMP]
    Lidpropotion = [Prop[i] for i in X if i < EMP] 
    Lidwidth = [LW[i] for i in X if i < EMP] 
    LidInitSat = [0 for i in X if i < EMP] 
    LidFromImp = [100 for i in X if i < EMP]
    LidToPerv = [0 for i in X if i < EMP]
    LidUsageN = [17, 17, 8, 11, 11, 11, 11, 10]
    LidNumber = [1 for i in X if i < EMP]
    Lidarea = [round(Lidpropotion[i]*float(Subarea[i])*10000,2) for i in range(len(Subarea))]

    '''2.3 Form LID dictionary'''

    LIDU = dict()
    LIDU['Subcatchment'] = Subname
    LIDU['LID_Process'] = Lidprocess
    LIDU['Number'] = LidNumber
    LIDU['Area'] = Lidarea
    LIDU['Width'] = Lidwidth
    LIDU['InitSat'] = LidInitSat
    LIDU['FromImp'] = LidFromImp
    LIDU['ToPerv'] = LidToPerv
    NUM = len(Subname)

    '''3. Formalization'''

    formedLIDU = formatinp(LidUsageN, LIDU)

    '''4. WriteINP'''

    insertLID('[LID_USAGE]', '[OUTFALLS]', formedLIDU, testINP, testINP2)

    '''5. Simulation and hydrological indicators evaluation'''

    st = SWMM5Simulation(testINP2)
    r = max(list(st.Results('NODE','outfall', indi)))

    '''6. Cost evaluation'''
    
    c = 0
    for i in range(NUM):
        c += LIDU['Area'][i] * Lidcost[i]
    
    return r,c

            
### 3. Main Process
# -----------------------------------------------------    
# 3.1 Pre-treatment
# -----------------------------------------------------  
print('INPUT YOUR DATA PATH')
D8_raster = input('The path of D8 flow direction(tif)')
DEM_raster = input('The path of DEM data(tif)')
Slope_raster = input('The path of Slope data(tif)')
Template_path = input('Template_path')
Output_path = input('Output_path')



d8 = gdal.Open(D8_raster) 
tempvar1 = [2**i for i in range(9)]
dem = gdal.Open(DEM_raster) 
slop = gdal.Open(Slope_raster) 

rows = d8.RasterYSize
cols = d8.RasterXSize
test = d8.GetGeoTransform()

d8a = d8.ReadAsArray()
dema = dem.ReadAsArray()
slo = slop.ReadAsArray()

d8c = standardized8(d8a)  
d8f = oriarray(d8a)  
desout = desarray(d8f, d8c) 

sub = subcatch(test, d8f, desout)
d8bool = d8c.tolist()

for i in range(len(d8bool)):
    for j in range(len(d8bool[i])):
        if d8bool[i][j] == 0:
            sub.pop('subcatchment'+ str(d8f[i][j]))
for i in sub.keys():
    if sub[i][4] in sub.keys():
        pass
    else:
        sub[i][4] = 'outfall'

subname = list(sub.keys())
subtext = []
raingage = 'R20' 
area = round(test[1] * test[1]*0.0001,3)
Width = [300 for i in range(len(subname))]
Slope = slo.reshape(1444).tolist()
Slope = [round(i+1,3) for i in Slope]
CurbLen = 0

for i in range(len(subname)):
    subtext.append([subname[i], raingage, sub[subname[i]][4], area, 5, Width[i], Slope[i],3,CurbLen])
SUBT = []
for i in range(len(subtext)):
    SUBTEMP = ''
    SUBTEMP += f'{subtext[i][0]}'.ljust(17)
    SUBTEMP += f'{subtext[i][1]}'.ljust(17)
    SUBTEMP += f'{subtext[i][2]}'.ljust(17)
    SUBTEMP += f'{subtext[i][3]}'.ljust(9)
    SUBTEMP += f'{subtext[i][4]}'.ljust(9)
    SUBTEMP += f'{subtext[i][5]}'.ljust(9)
    SUBTEMP += f'{subtext[i][6]}'.ljust(9)
    SUBTEMP += f'{subtext[i][7]}'.ljust(25)
    
    SUBTEMP += '\n '
    SUBT.append(SUBTEMP)

subtitle = ["[SUBCATCHMENTS]\n", 
             ";;Name           Rain Gage        Outlet           Area     %Imperv  Width    %Slope   CurbLen  SnowPack        \n",
             ";;-------------- ---------------- ---------------- -------- -------- -------- -------- -------- ----------------\n"]
opath = Template_path
dpath = Output_path 
pos1 = '[SUBCATCHMENTS]'
pos2 = '[OUTFALLS]' 

writeinp(opath, dpath, SUBT, pos1, pos2, subtitle)

poly = []
for i in sub.keys():
    for j in range(4):
        templ = []
        templ.append(i)
        templ.append(round(sub[i][j][0],3))
        templ.append(round(sub[i][j][1],3))
        poly.append(templ)

charnum = [17,19,18]
Formedpoly = []
for i in range(len(poly)):
    temp = ''
    for j in range(len(poly[i])):
        temp += f'{poly[i][j]}'.ljust(charnum[j])
    temp += '\n '
    Formedpoly.append(temp)
        
potitle = [ "[Polygons]\n", 
             ";;Subcatchment   X-Coord            Y-Coord           \n",
             ";;-------------- ------------------ ------------------\n"]
opath = Template_path
dpath = Output_path
pos1 = '[Polygons]'
pos2 = '[SYMBOLS]'
writeinp(opath, dpath, Formedpoly, pos1, pos2, potitle)

# 3.2 set LID parameters
# -----------------------------------------------------  
INP = input('INP path after per-treated')
testINP2 = input('output INP path')
testINP = input('backup INP path')

'''global variables'''

GEN = [0]
popid = [0]


file = open(INP, "r", encoding = "UTF-8" )
content = file.read() 
Subcatchment = extractinp(content,'[SUBCATCHMENTS]','[SUBAREAS]' )
NUMBER = Subcatchment.shape[0]
Lidcontrols = extractinp(content,'[LID_CONTROLS]','[LID_USAGE]' ) 
Lidnames = list(filter(None,list(Lidcontrols[0].drop_duplicates().dropna())))

print("These are LID types in your file:" + Lidnames)
print("Please set their parameters")
prop = input("Input a List of each LID's propotion from 0-1, for example: [0, 0.2, 0.3, 1]")
width = input("Input a List of each LID's width, for example: [10, 2, 5, 10]")
cos = input("Input a List of each LID's  cost (yuan / square meter), for example: [300, 500, 300, 200]"
indi = input("Input numbers of hydrological indicators: 0 Depth of water above invert (ft or m); 1 Hydraulic head (ft or m); 2 Volume of stored + ponded water (ft3 or m3); 3 Lateral inflow (flow units); 4 Total inflow (lateral + upstream) (flow units); 5 Flow lost to flooding (flow units); 6 Concentration of TSS (mg/l)")
    
# 4 Multi-objective spatial optimization
# -----------------------------------------------------  
            
Popu = input("The population number")
Gene = input("The generation number")

class SWMMOPTIM(Problem):

    def __init__(self):
        super().__init__(n_var=NUMBER,
                         n_obj=2,
                         n_constr=0,
                         xl=np.array([0 for i in range(NUMBER)]),
                         xu=np.array([20 for i in range(NUMBER)]),
                         type_var=int)

    def _evaluate(self, x, out, *args, **kwargs):
        f1,f2 = Compile(x)

        out["F"] = np.column_stack([f1, f2])



problem = SWMMOPTIM()

algorithm = NSGA2(
    pop_size=Popu,
    n_offsprings=Popu,
    sampling=get_sampling("int_random"),
    crossover=get_crossover("int_sbx", prob=0.9, eta=15),
    mutation=get_mutation("int_pm", eta=20),
    eliminate_duplicates=True
)
termination = get_termination("n_gen", Gene)
            
Flag = input("Whether to start the optimization process?(Y/N)")

            
if Flag == "Y":
    res = minimize(problem,
               algorithm,
               termination,
               seed=1,
               save_history=True,
               verbose=True)


# 5 Visualization & Analysis
# -----------------------------------------------------  

            
d8a = d8.ReadAsArray().tolist()
d8b = []
counter = 0
for i in range(len(d8a)):
    d8b.append([])
    for j in range(len(d8a[i])):
        if d8a[i][j] == 0:
            d8b[i].append(0)
        else:
            counter += 1
            d8b[i].append(counter)            

