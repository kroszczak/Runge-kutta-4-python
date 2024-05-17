import matplotlib.pyplot as plt
import numpy as np
import copy

# tablica zmiennych - dane początkowe
s = [ 250000,250000,250000,250000 ]
labels = ['p53', 'MDMcyt', 'MDMn', 'PTEN']

def pochodne(s,k):
    # Stałe
    p1 = 8.8
    p2 = 440
    # p2 = 440 * 0.02 # włączone siRNA
    p3 = 100
    # p3 = 0 # bez PTEN
    d1 = 1.375 * 10**-14
    d2 = 1.375 * pow(10,-4) * 0.1
    # d2 = 1.375 * pow(10,-4) # uszkodzenie DNA
    d3 = 3 * pow(10,-5)
    k1 = 1.925 * pow(10,-4)
    k2 = pow(10,5)
    k3 = 1.5 * pow(10,5)

    # wzory na kolejne zmienne
    k[0] = p1 - (d1 * s[0] * pow(s[2],2))
    k[1] = p2 * ( pow(s[0],4) / (pow(s[0],4) + pow(k2,4)) ) - ( k1 * ( pow(k3,2) / (pow(k3,2) + pow(s[3],2)) ) * s[1] ) - (d2 * s[1])
    k[2] = (k1 * pow(k3,2) * s[1]) / (pow(k3,2) + pow(s[3],2)) - (d2 * s[2])
    k[3] = p3 * (pow(s[0],4) / (pow(s[0],4) + pow(k2,4))) - (d3 * s[3])


def rk4_vec(t, dt, n, s, f):
    M = len(s)
    k1 = np.zeros(M)
    k2 = np.zeros(M)
    k3 = np.zeros(M)
    k4 = np.zeros(M)
    w = np.zeros(M)

    for i in range(n): w[i]=s[i] # do W zapisz aktualny stan danej zmiennej z tablicy S
    pochodne(w, k1) # podajemy tablibce W ORAZ k1 - k1 to jest tablica stanów k1 dla wszytskich elementów, i odpowiednio kolejne tez; 
    for i in range(n): w[i]=s[i]+dt/2*k1[i]
    pochodne(w,k2)
    for i in range(n): w[i]=s[i]+dt/2*k2[i]
    pochodne(w,k3)
    for i in range(n): w[i]=s[i]+dt*k3[i]
    pochodne(w,k4)
    for i in range(n): s[i]=s[i]+dt/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i])
    return 0

# ustawienia
n = 4 # ilość zmiennych
hours = 64
tmax = 3600*hours# czas symulacji
dt = 1
N = tmax/dt; # ilość kroktów czasowych

t = 0;
time = []
wynik = []
for i in range(int(N)):
    rk4_vec(t, dt, n, s, pochodne) # to moj step
    wynik.append(copy.deepcopy(s))
    t += dt
    time.append(t)

# PLOT
[plt.plot(time, [column[i] for column in wynik], label=labels[i]) for i in range(n)]
plt.xlabel('Czas (s)')
plt.ylabel('Liczba cząsteczek')
plt.title(f'Symulacja {hours}h')
plt.xticks(np.arange(0, tmax, step= (lbs := int(tmax/10))), labels=[str(int(hours*i/10)) for i in range(10)])
plt.legend()
plt.savefig('wykres.jpg')