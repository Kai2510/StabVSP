# stabvsp.py
# StabVSP - Python script for evaluate stablity using OpenVSP's stab file.
# Version: v2.0
# Author: LI Zhikai
import numpy as np
from pathlib import Path
import pandas as pd
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
mpl.rcParams.update({'figure.dpi':150})
# path = mpl.rcParams['datapath']
from io import StringIO
import os
import warnings
warnings.filterwarnings('ignore')

# #theta, gamma in degrees
# thetaR=5 # [degree]
# gammaR=3 # [degree]
G=9.81 # [m/s^2]

def readmass(mass_in_name):
    MassProperties={}
    data_in=pd.read_csv(mass_in_name,delimiter='\s+');#这里分隔符是任意多空白字符，包括空格、制表符、换页符
    m_tot=data_in['m_tot'][0]
    CGx=data_in['CGx'][0]
    CGy=data_in['CGy'][0]
    CGz=data_in['CGz'][0]
    Iyy=data_in['Iyy'][0]
    Ixz=data_in['Ixz'][0]
    Izz=data_in['Izz'][0]
    Ixx=data_in['Ixx'][0]
    #Series to floats
    Izx=Ixz;
    MassProperties['m']=m_tot
    MassProperties['Ixx']=Ixx
    MassProperties['CGx']=CGx
    MassProperties['Iyy']=Iyy
    MassProperties['CGz']=CGz
    MassProperties['CGy']=CGy
    MassProperties['Izz']=Izz
    MassProperties['Ixz']=Ixz
    MassProperties['Izx']=Izx
    return MassProperties

class StabDataWC:  #Each Working Condition
    def __init__(self,Refs,Coefs):
        self.Refs=Refs
        self.Coefs=Coefs
        self.pDynCoefs = pd.DataFrame(data=None)
        self.qDynCoefs = pd.DataFrame(data=None)
        self.rDynCoefs = pd.DataFrame(data=None)
    def __add_pDynCoefs__(self,pDynCoefs):
        self.pDynCoefs=pDynCoefs
    def __add_qDynCoefs__(self,qDynCoefs):
        self.pDynCoefs=qDynCoefs
    def __add_rDynCoefs__(self,rDynCoefs):
        self.pDynCoefs=rDynCoefs
    def __LateralData__(self,MassProps,thetaRinDegrees,gammaRinDegrees):
        self.EigValuesLtrl,self.EigVectorsLtrl,self.A2 = CalLateral(self.Refs,self.Coefs,MassProps,thetaRinDegrees,gammaRinDegrees,self.pDynCoefs,self.rDynCoefs)
        # print(self.EigValues)
    def __LongitudinalData__(self,MassProps,thetaRinDegrees,gammaRinDegrees):
        self.EigValuesLgl,self.EigVectorsLgl,self.A1= CalLongitudinal(self.Refs,self.Coefs,MassProps,thetaRinDegrees,gammaRinDegrees,self.qDynCoefs)
    def __CalDynCoefs__(self):
        # 定义要处理的变量前缀（CL, CD, CS, CMl, CMm, CMn）
        prefixes = ["CL", "CD", "CS", "CMl", "CMm", "CMn"]
                # 遍历每个前缀，计算 alpha_dot 并添加到 DataFrame
        for prefix in prefixes:
            if f"{prefix}_(q+alpha_dot)" in self.qDynCoefs.index:
                alpha_dotted_value = self.qDynCoefs.loc[f"{prefix}_(q+alpha_dot)", "Value"] - self.Coefs.loc[prefix,'p']
                self.qDynCoefs.loc[f"{prefix}_alpha_dot", "Value"] = alpha_dotted_value
        for prefix in prefixes:
            if f"{prefix}_(r-beta_dot)" in self.rDynCoefs.index:
                beta_dotted_value = - self.rDynCoefs.loc[f"{prefix}_(r-beta_dot)", "Value"] + self.Coefs.loc[prefix,'r']
                self.rDynCoefs.loc[f"{prefix}_beta_dot", "Value"] = beta_dotted_value


class StabDataCon: #Each Configuration
    def __init__(self,MassProp,StabDataList,dirname):
        self.MassProp=MassProp
        self.StabDataList=StabDataList
        self.dirname=dirname


# def ReadCase_UNDONE(linescache):
#     varnames=[]
#     reflines=[]
#     coeflines=[];
#     for iline in range(len(linescache)):
#         if (linescache[iline].startswith('*')
#             or linescache[iline].startswith('#')
#             or linescache[iline].startswith('\n')):
#             continue
#         elif linescache[iline].startswith('Case'):
#             iline+=10;
#             # lines = lines[11:]
#             continue
#         elif linescache[iline].startswith('Coef'):
#             coeflines.append(linescache[iline])
#             for jline in range(1,17):
#                 if (linescache[iline+jline].startswith('#')
#                     or linescache[iline+jline].startswith('\n')):
#                     continue;
#                 else:
#                     coeflines.append(linescache[iline+jline])
#             iline+=jline
#             # print(iline)
#         else:
#             coeflines.append(linescache[iline])
#             for jline in range(1,13):
#                 if (linescache[iline+jline].startswith('#')
#                     or linescache[iline+jline].startswith('\n')):
#                     continue;
#                 else:
#                     # !!! UNDONE: Result line
#                     parts = linescache[iline].split()
#                     print(parts)
#                     if len(parts) >= 3:
#                         name = parts[0]
#                         value = float(parts[1])
#                         units = ' '.join(parts[2:])  # Units 可能包含空格（如 "Lunit^2"）
#                         reflines.append([name, value, units])
#                     # reflines.append(linescache[iline+jline])
#             iline+=jline
            
#     # print(RefVars)
#     # print(coeflines)
#     coefs= pd.read_csv(StringIO(''.join(coeflines)), sep="\s+",index_col="Coef")
#     RefVars = pd.DataFrame(RefVars, columns=['Name', 'Value', 'Units'])
#     RefVars = RefVars.set_index('Name') 
#     # print(coefs.info())
#     return RefVars, coefs
def ReadCase(linescache):
    varnames=[]
    RefVars={}
    coeflines=[];
    for iline in range(len(linescache)):
        if (linescache[iline].startswith('*')
            or linescache[iline].startswith('#')
            or linescache[iline].startswith('\n')):
            continue
        elif linescache[iline].startswith('Case'):
            iline+=10;
            # lines = lines[11:]
            continue
        elif linescache[iline].startswith('Coef'):
            coeflines.append(linescache[iline])
            for jline in range(1,17):
                if (linescache[iline+jline].startswith('#')
                    or linescache[iline+jline].startswith('\n')):
                    continue;
                else:
                    coeflines.append(linescache[iline+jline])
            iline+=jline
            # print(iline)
        else:
            parts = linescache[iline].split()
            # print(parts,'\n')
            # var_name = parts[0]
            # var_value = float(parts[1])
            # exec('parts[0]=float(parts[1])')
            varnames.append(parts[0]);
            RefVars[parts[0]]=parts[1];
    # print(RefVars)
    # print(coeflines)
    coefs= pd.read_csv(StringIO(''.join(coeflines)), sep="\s+",index_col="Coef")
    # print(coefs.info())
    return RefVars, coefs


# Read data from the file
def ReadStabFile(stbname):
    StabData=[];
    try:
        with open(stbname, 'r') as file:
            lines = file.readlines()
        ncase=0;
        # print(len(lines)//57)
        for iline in range(0,len(lines)):
            if lines[iline].startswith('*'):
                ncase+=1;
                cache1, cache2=ReadCase(lines[iline:iline+57])
                StabData.append(StabDataWC(cache1,cache2))
                iline+=50;
            else:
                # iline+=50;
                continue
    except FileNotFoundError:
        print(f"Warning: {stbname} Not exist. Skip the file.")
    except PermissionError:
        print(f"Warning: no permission to read {stbname}. Skip the file.")
    except Exception as e:
        print(f"Warning: Unknown error: {e}. Skip the file.")
    else:
        print(f"File {stbname} has been read successfully.")

    
    return StabData

# Read Dynamic Stability derivatives
def Read_Add_DynCase(file_path,StabData):
    try:
        with open(file_path, 'r') as file:
            ext = Path(file_path).suffix[1:]
            lines = file.readlines()
        ncase=0;
        # print(len(lines)//34)
        for iline in range(0,len(lines)):
            if lines[iline].startswith('*'):
                ncase+=1;
                cache_Refs, cache_Derivatives=ReadDynStabCase(lines[iline:iline+34])
                for stabcase in StabData:
                    if (abs(float(stabcase.Refs['AoA_']) - cache_Refs.loc['AoA_','Value']) <0.01 and 
                        abs(float(stabcase.Refs['Beta_']) - cache_Refs.loc['Beta_','Value']) <0.01 ):
                            match ext:
                                case "pstab":
                                    stabcase.pDynCoefs = cache_Derivatives
                                case "qstab":
                                    stabcase.qDynCoefs = cache_Derivatives
                                case "rstab":
                                    stabcase.rDynCoefs = cache_Derivatives
                            break;
                iline+=34;
            else:
                # iline+=50;
                continue
    except FileNotFoundError:
        print(f"Warning: {file_path} Not exist. Skip the file.")
    except PermissionError:
        print(f"Warning: no permission to read {file_path}. Skip the file.")
    except Exception as e:
        print(f"Warning: Unknown error: {e}. Skip the file {file_path}.")
    else:
        print(f"File {file_path} has been read successfully.")

    return #StabData

def ReadDynStabCase(lines):
    # # 找到第二部分开始的标志
    # start_index = 0
    # for i, line in enumerate(lines):
    #     if line.startswith('# Name \t\t Value'):
    #         start_index = i + 1
    #         break
    # # 提取第二部分数据
    # data = []
    # for line in lines[start_index:]:
    #     line = line.strip()
    #     # 跳过空行和注释行
    #     if not line or line.startswith('#'):
    #         continue
    #     # 分割名称和值
    #     parts = line.split()
    #     if len(parts) >= 2:
    #         name = parts[0]
    #         # 值可能是负数，所以需要小心处理
    #         # 从右向左找第一个数值
    #         for i in range(len(parts)-1, -1, -1):
    #             try:
    #                 value = float(parts[i])
    #                 name = ' '.join(parts[:i]).strip()
    #                 data.append([name, value])
    #                 break
    #             except ValueError:
    #                 continue
    # df = pd.DataFrame(data, columns=['Name', 'Value'])
    # df['Name'] = df['Name'].str.replace(' ', '') # Delete spaces
    # df = df.set_index('Name') 
    # return df
        # 初始化两个数据存储列表
    data_part1 = []  # 存储第一部分数据（Name, Value, Units）
    data_part2 = []  # 存储第二部分数据（Name, Value）

    # 先处理第一部分（带 Units 的数据）
    part1_start = False
    for line in lines:
        line = line.strip()
        if line.startswith('# Name') and 'Units' in line:
            part1_start = True
            continue  # 跳过标题行
        elif line.startswith('# Name') and 'Value' in line:
            part1_start = False  # 第一部分结束，进入第二部分
            break  # 跳出循环，准备处理第二部分

        if part1_start and line and not line.startswith('#'):
            parts = line.split()
            if len(parts) >= 3:
                name = parts[0]
                value = float(parts[1])
                units = ' '.join(parts[2:])  # Units 可能包含空格（如 "Lunit^2"）
                data_part1.append([name, value, units])

    # 处理第二部分（仅 Name 和 Value）
    part2_start = False
    for line in lines:
        line = line.strip()
        if line.startswith('# Name') and 'Value' in line and not 'Units' in line:
            part2_start = True
            continue  # 跳过标题行
        elif part2_start:
            if not line or line.startswith('#'):
                continue  # 跳过空行和注释
            parts = line.split()
            if len(parts) >= 2:
                name = parts[0]
                # 从右向左找第一个数值（可能含负号）
                for i in range(len(parts)-1, -1, -1):
                    try:
                        value = float(parts[i])
                        name = ' '.join(parts[:i]).strip()
                        data_part2.append([name, value])
                        break
                    except ValueError:
                        continue
    df_Refs = pd.DataFrame(data_part1, columns=['Name', 'Value', 'Units'])
    df_Derivatives= pd.DataFrame(data_part2, columns=['Name', 'Value'])
    df_Derivatives['Name'] = df_Derivatives['Name'].str.replace(' ', '') # Delete spaces
    df_Derivatives = df_Derivatives.set_index('Name') 
    df_Refs = df_Refs.set_index('Name') 
    return df_Refs,df_Derivatives

def CalLateral(Refs,Coefs,MassProps,thetaRinDegrees,gammaRinDegrees,pDynCoefs,rDynCoefs):
    thetaR=thetaRinDegrees      * math.pi/180
    gammaR=gammaRinDegrees      * math.pi/180
    alphaR= float(Refs['AoA_']) * math.pi/180
    VR=     float(Refs['Vinf_'])
    rhoR=   float(Refs['Rho_'])
    # rhoR=   1.225
    qR=     0.5*rhoR*VR**2
    cR=     float(Refs['Cref_'])
    bR=     float(Refs['Bref_'])
    SR=     float(Refs['Sref_'])
    # print(rhoR)

    Ixz=MassProps['Ixz']
    Ixx=MassProps['Ixx']
    Iyy=MassProps['Iyy']
    Izz=MassProps['Izz']
    m=  MassProps['m']  
    Izx=Ixz
    
    Cal_L_bar=lambda Li,Ni:(Li+(Ixz/Izz)*Ni)/(Ixx-Ixz**2/Izz)
    Cal_N_bar=lambda Li,Ni:(Ni+(Ixz/Ixx)*Li)/(Izz-Ixz**2/Ixx)
    Cyb=Coefs.loc['CFy','Beta']  ;Yb_bar=Cyb*qR*SR/m /VR;
    Cyp=Coefs.loc['CFy','p']   ;Yp_bar=Cyp*qR*SR/m /VR * bR /(2*VR)
    Cyr=Coefs.loc['CFy','r']   ;Yr_bar=Cyr*qR*SR/m /VR * bR /(2*VR)
    Cnb=Coefs.loc['CMn','Beta']   ;Nb=Cnb*qR*SR*bR;
    Clb=Coefs.loc['CMl','Beta']  ;Lb=Clb*qR*SR*bR;
    Lb_bar=Cal_L_bar(Lb,Nb);Nb_bar=Cal_N_bar(Lb,Nb);
    # print(Lb_bar,Nb_bar)
    Cnp=Coefs.loc['CMn','p']  ;Np=Cnp*qR*SR*bR* bR/(2*VR)
    Clp=Coefs.loc['CMl','p']  ;Lp=Clp*qR*SR*bR* bR/(2*VR)
    Lp_bar=Cal_L_bar(Lp,Np);Np_bar=Cal_N_bar(Lp,Np)
    Cnr=Coefs.loc['CMn','r'] ;Nr=Cnr*qR*SR*bR* bR/(2*VR)
    Clr=Coefs.loc['CMl','r'];Lr=Clr*qR*SR*bR* bR/(2*VR)
    Lr_bar=Cal_L_bar(Lr,Nr);Nr_bar=Cal_N_bar(Lr,Nr)


    A2=np.asmatrix([[Yb_bar, alphaR+Yp_bar, Yr_bar - 1, G*math.cos(thetaR)/VR],
            [Lb_bar, Lp_bar, Lr_bar, 0],
            [Nb_bar, Np_bar, Nr_bar, 0],
            [0, 1, math.tan(thetaR), 0]])
    # print(A2)
    # print('CGz',MassProps['CGz'],'AoA',Refs['AoA_'])
    # print('AoA',Refs['AoA_'],'Beta',Refs['Beta_'])
    # print('spiral',Lb_bar*Nr_bar- Lr_bar*Nb_bar)
    # print('dutch roll',Nb_bar/Lb_bar)

    eigValue, eigVector = np.linalg.eig(A2);
    return eigValue, eigVector, A2


def CalLongitudinal(Refs,Coefs,MassProps,thetaRinDegrees,gammaRinDegrees,qDynCoefs):
    thetaR=thetaRinDegrees      * math.pi/180
    gammaR=gammaRinDegrees      * math.pi/180
    alphaR= float(Refs['AoA_']) * math.pi/180
    VR=     float(Refs['Vinf_'])
    rhoR=   float(Refs['Rho_'])
    # rhoR=   1.225
    qR=     0.5*rhoR*VR**2
    cR=     float(Refs['Cref_'])
    bR=     float(Refs['Bref_'])
    SR=     float(Refs['Sref_'])
    C_LR=   Coefs.loc['CL','Total']
    C_DR=   Coefs.loc['CD','Total']
    C_mR=   Coefs.loc['CMm','Total']
    # T_R=    C_DR * qR * SR;
    T_R = 0 # Thrust_ref can be specified

    Ixz=MassProps['Ixz']
    Ixx=MassProps['Ixx']
    Iyy=MassProps['Iyy']
    Izz=MassProps['Izz']
    m=  MassProps['m']  
    Izx=Ixz
    
    C_DV=   Coefs.loc['CD','U']
    C_LV=   Coefs.loc['CL','U']
    X_V=    -(C_DV+2*C_DR)* qR *SR /m /VR; #Ignore T_V
    Z_V=    (C_LV+2*C_LR)*qR*SR /m /VR /VR #Ignore T_V
    C_La=   Coefs.loc['CL','Alpha'];Z_a=(C_DR+C_La)*qR*SR/m /VR;
    C_Lq=   Coefs.loc['CL','q'];Z_q=C_Lq*(cR/2 /VR) *qR*SR/m/VR;
    C_mV=   Coefs.loc['CMm','U'];M_V_bar=(C_mV+2*C_mR)*qR*SR*cR /VR /Iyy;
    C_mq=   Coefs.loc['CMm','q'];M_q_bar=C_mq*(cR/2 /VR)*qR*SR*cR/Iyy;
    C_Da=   Coefs.loc['CD','Alpha'];X_a=(-T_R*math.sin(alphaR)-C_Da*qR*SR)/m
    C_ma=   Coefs.loc['CMm','Alpha'];M_a_bar=C_ma*qR*SR*cR/Iyy;

    if(qDynCoefs.empty):
        C_La_dot = 0; Z_a_dot = 0;
        C_ma_dot = 0; M_a_dot_bar = 0;
    else:
        C_La_dot = qDynCoefs.loc['CL_(q+alpha_dot)','Value']-C_Lq
        Z_a_dot = C_La_dot*(cR/2/VR)*qR*SR/m /VR
        C_ma_dot = qDynCoefs.loc['CMm_(q+alpha_dot)','Value']-C_mq 
        M_a_dot_bar = C_ma_dot*(cR/2/VR)*qR*SR*cR/Iyy
        
    # Exact Coef Matrix in Body Coord.
    A1=np.asmatrix([[X_V,X_a+G*math.cos(gammaR), 0, -G*math.cos(gammaR)],
            [-Z_V/(1+Z_a_dot), -(Z_a-G*math.sin(gammaR)/VR)/(1+Z_a_dot), 
             (1-Z_q)/(1+Z_a_dot), -G*math.sin(gammaR)/VR/(1+Z_a_dot)],
            [M_V_bar-M_a_dot_bar*Z_V/(1+Z_a_dot), 
             M_a_bar-M_a_dot_bar*(Z_a-G*math.sin(gammaR)/VR)/(1+Z_a_dot), 
             M_q_bar+M_a_dot_bar*(1-Z_q)/(1+Z_a_dot), -M_a_dot_bar*G*math.sin(gammaR)/VR/(1+Z_a_dot)],
            [0, 0, 1, 0]])
    
    # Approx. Way: Ignore Z_a_dot
    # A1=np.asmatrix([[X_V,X_a+G*math.cos(gammaR), 0, -G*math.cos(gammaR)],
    #         [-Z_V, -(Z_a-G*math.sin(gammaR)/VR), 
    #          (1-Z_q), -G*math.sin(gammaR)/VR],
    #         [M_V_bar-M_a_dot_bar*Z_V, 
    #          M_a_bar-M_a_dot_bar*(Z_a-G*math.sin(gammaR)/VR), 
    #          M_q_bar+M_a_dot_bar*(1-Z_q), -M_a_dot_bar*G*math.sin(gammaR)/VR],
    #         [0, 0, 1, 0]])
    # print(A1)

    eigValue, eigVector = np.linalg.eig(A1);
    return eigValue, eigVector, A1


def cal_w_c(eig):
    eigR,eigI=eig.real,eig.imag;
    w=math.sqrt(eigR**2+eigI**2)
    c=-eigR/w
    return w,c

def Cal_Lat_Lgl(MassProperties,StabData,thetaR,gammaR):
    for stabcase in StabData:
        stabcase.__LateralData__(MassProperties,thetaR,gammaR)
        stabcase.__LongitudinalData__(MassProperties,thetaR,gammaR)

def Cal_Dyn_Coefs_Batch(StabData):
    for stabcase in StabData:
        stabcase.__CalDynCoefs__()

        #Plot function
def PlotStabScatter(fig,EigValues,
                    color,mkr
                    ):
    # for ite in EigValues:
    #     fig.scatter(ite.real, ite.imag,
    #                 c=color,
    #                 marker='x')
    plt.plot(EigValues.real, EigValues.imag,
                 marker=mkr, c=color,ms=5,linewidth=0)
    return