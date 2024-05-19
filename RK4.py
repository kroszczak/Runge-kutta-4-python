import matplotlib.pyplot as plt
from RPNparser import infix_to_rpn
import numpy as np
import copy
from sympy import sympify, symbols
import yaml
from datetime import datetime, date
import math

# Wczytaj plik YAML
sub_dict = dict()
fin_dict = dict()
funtions = dict()
data = None
with open("model.yaml", 'r') as file:
    data = yaml.safe_load(file)
for key, value in data['parameters'].items():
    sub_dict[symbols(key)] = value
for key, value in data['factors'].items():
    fin_dict[key] = value['starting_value']
for key, value in data['factors'].items():
    args = ",".join([i for i in list(sorted(data['factors'].keys()))])
    funtions[key] = eval(f"lambda {args}: {sympify(value['function']).subs(sub_dict)}")
# RK SIMULATION
def pochodne(s,k):
    for key, val in funtions.items(): k[key] = val(**fin_dict)

def rk4_vec(hd):
    for j, factor in loop:
        prev = rk_step_matrix[j-1]
        for i in funtions.keys():
            ifactor = factor*prev.get(i,0)
            rk_step_matrix[4][i]=fin_dict[i]+ifactor
        pochodne(rk_step_matrix[4], rk_step_matrix[j])
    for i in funtions.keys(): fin_dict[i]=fin_dict[i]+hd*(rk_step_matrix[0][i]+2*rk_step_matrix[1][i]+2*rk_step_matrix[2][i]+rk_step_matrix[3][i])
    return 0

# MAIN LOOP
t = 0;
dt = data['options']['dt']
time = []
wynik = []
dt = dt/6
tmax = 3600*data['options']['hours']
N = (3600*data['options']['hours'])/dt
rk_step_matrix = [dict() for i in range(5)]
loop = list(zip(list(range(4)), [0, dt/2, dt/2, dt]))
pct = 0
for i in range(int(N)):
    if  (pct := math.floor(100*i/N)) > pct: print(pct)
    rk4_vec(dt) # to moj step
    wynik.append(copy.deepcopy(fin_dict))
    time.append(t:=t+dt)
    

# PLOT
input(wynik[-1])
[plt.plot(time, [column[key] for column in wynik], label=key) for i, key in enumerate(funtions.keys())]
plt.xlabel('Czas (s)')
plt.ylabel('Liczba cząsteczek')
plt.title(f"Symulacja {data['options']['hours']}h")
plt.xticks(np.arange(0, tmax, step=(lbs := int(tmax/10))), labels=[str(int(data['options']['hours']*i/10)) for i in range(10)])
plt.legend()
plt.grid(True)
plt.savefig('wykres.jpg')