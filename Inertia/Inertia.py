import numpy as np
import pandas as pd

in_name="Mass_Inertia_imp.txt"#这里改输入文件名
out_name="Mass_Inertia_out.txt"#这里改输出文件名
data_in=pd.read_csv(in_name,delimiter='\s+');#这里分隔符是任意多空白字符，包括空格、制表符、换页符

print(data_in)
m_tot=data_in['m'].sum()
# CGx=((data_in.loc[0,'m']*data_in.loc[0,'CGx'])+(data_in.loc[1,'m']*data_in.loc[1,'CGx']))/m_tot
# tmp=0
# for i in range (len(data_in)):
#     tmp+=(data_in.loc[i,'m']*data_in.loc[i,'CGx'])
#     print(tmp)
CGx=sum((data_in.loc[i,'m']*data_in.loc[i,'CGx'])for i in range(len(data_in)))/m_tot
CGy=sum((data_in.loc[i,'m']*data_in.loc[i,'CGy'])for i in range(len(data_in)))/m_tot
CGz=sum((data_in.loc[i,'m']*data_in.loc[i,'CGz'])for i in range(len(data_in)))/m_tot
Iyy=sum( (data_in.loc[i,'Iyy']+data_in.loc[i,'m'] * ((data_in.loc[i,'CGx']-CGx)**2+(data_in.loc[i,'CGz']-CGz)**2)) for i in range(len(data_in)))
Izz=sum( (data_in.loc[i,'Izz']+data_in.loc[i,'m'] * ((data_in.loc[i,'CGx']-CGx)**2+(data_in.loc[i,'CGy']-CGy)**2)) for i in range(len(data_in)))
Ixx=sum( (data_in.loc[i,'Ixx']+data_in.loc[i,'m'] * ((data_in.loc[i,'CGy']-CGy)**2+(data_in.loc[i,'CGz']-CGz)**2)) for i in range(len(data_in)))
Ixz=sum( (data_in.loc[i,'Ixz']+data_in.loc[i,'m'] * (data_in.loc[i,'CGz']-CGz)*(data_in.loc[i,'CGx']-CGx)) for i in range(len(data_in)) )
print(m_tot,CGx,CGy,CGz,Ixx,Iyy,Izz,Ixz)
lines=['m_tot\tCGx\tCGy\tCGz\tIxx\tIyy\tIzz\tIxz\n','{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t'.format(m_tot,CGx,CGy,CGz,Ixx,Iyy,Izz,Ixz)]
with open(out_name, 'w') as file:
    file.writelines(lines)
