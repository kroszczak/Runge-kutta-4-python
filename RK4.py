import matplotlib.pyplot as plt
import numpy as np
from sympy import sympify, symbols
from yaml import safe_load
import os 
# from alive_progress import alive_bar
import math
import copy
fin_dict = dict()
funtions = dict()
print("reading config file...")
data = safe_load(open("model.yaml", 'r'))

print("preparing function definitions...")
for key, value in data['factors'].items():
    fin_dict[key] = value['starting_value']
    funtions[key] = eval(f"lambda {','.join([i for i in list(sorted(data['factors'].keys()))])}: {sympify(value['function']).subs({symbols(key): value for key, value in data['parameters'].items()})}")
    
def rk4_vec(dt, fin_dict):
    for j, factor in enumerate([0, dt/2, dt/2, dt]):
        for i in funtions.keys():rk_step_matrix[4][i]=fin_dict[i]+factor*rk_step_matrix[j-1].get(i,0)
        for key, val in funtions.items(): rk_step_matrix[j][key] = val(**rk_step_matrix[4])
    for i in funtions.keys(): fin_dict[i]=fin_dict[i]+dt/6*(rk_step_matrix[0][i]+2*rk_step_matrix[1][i]+2*rk_step_matrix[2][i]+rk_step_matrix[3][i])
    return {key:fin_dict[key] for key in fin_dict}

print("setting simulation variables...")
rk_step_matrix = [dict() for i in range(5)]
h =  data['options']['dt']
tmax = 3600*data['options']['hours']
N = int((3600*data['options']['hours'])/h)

wynik = [rk4_vec(h, fin_dict) for i in range(N)]
# wynik = []
# with alive_bar(N) as bar:
#     for i in range(N):
#         wynik.append(rk4_vec(h, fin_dict))
#         bar()
# PLOT
print("exporting.jpg file...")
[plt.plot(range(N), [column[key] for column in wynik], label=key) for key in funtions.keys()]
plt.xlabel('Czas (s)')
plt.ylabel('Liczba cząsteczek')
plt.title(f"Symulacja {data['options']['hours']}h")
plt.xticks(np.arange(0, N, step=int(math.ceil(N/10))), labels=[str(int(data['options']['hours']*i/10)) for i in range(10)])
plt.legend()
plt.grid(True)
plt.savefig('wykres.jpg')