import matplotlib.pyplot as plt
import numpy as np
import propagação_de_incertezas as P_I
from uncertainties import ufloat
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
from math import pi,sqrt,exp 
from uncertainties import ufloat
from pratica_5 import H_i,H_s,H_i_graf,H_s_graf,vazão_convertida_i,vazão_convertida_s
#H_i,H_s -> Conjunto de ufloats para operações
#H_i_graf, H_s_graf -> Valores nominais usados na plotagem dos dados
#vazão_convertida_i, vazão_convertida_s -> Vazões associadas aos pontos anteriores

def plot_curve_HxQ(H,Q):
    plt.figure(figsize=(10, 6))
    plt.plot(Q, H, marker='o', linestyle='--', color='b', label='H vs Q')
    # Adicionando títulos e rótulos
    plt.title('Gráfico Altura(H) e Vazão(Q)')
    plt.xlabel('Q(m^3/s)')
    plt.ylabel('H(m)')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_all_curves(vazão_s,vazão_i,H_s,H_i):
    plt.figure(figsize=(10, 6))
    plt.plot(vazão_i, H_i, marker='o', linestyle='--', color='b', label='Bomba Inferior')
    plt.plot(vazão_s, H_s, marker='h', linestyle='-', color='r', label='Bomba Superior')
    # Adicionando títulos e rótulos
    plt.title('Gráfico Altura(H) e Vazão(Q)')
    plt.xlabel('Q(m^3/h)')
    plt.ylabel('H(m)')
    plt.legend()
    plt.grid(True)
    plt.show()

def show_results(vazão_s,vazão_i,H_s,H_i):
    print(f'Condições atmosféricas: {pressao_lab} Pa')
    print(f'Vazões: {vazão_s}, {vazão_i} m^3/s')
    print(f'Alturas manométricas(H): {H_s}, {H_i} m')

def plot_series_curve(vazão_s_importado,vazão_i_importado,vazão_serie_teorico,vazão_serie_pratico,H_s_importado,H_i_importado,H_serie_teorico,H_serie_pratico):
    plt.figure(figsize=(10, 6))
    plt.plot(vazão_i_importado, H_i_importado, marker='o', linestyle='--', color='b', label='Bomba Inferior')
    plt.plot(vazão_s_importado, H_s_importado, marker='h', linestyle='-', color='r', label='Bomba Superior')
    plt.plot(vazão_serie_teorico, H_serie_teorico, marker='o', linestyle='--', color='green', label='Série Teórico')
    plt.plot(vazão_serie_pratico, H_serie_pratico, marker='h', linestyle='-', color='black', label='Série Prático')
    # Adicionando títulos e rótulos
    plt.title('Gráfico Altura(H) e Vazão(Q)')
    plt.xlabel('Q(m^3/h)')
    plt.ylabel('H(m)')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_parallel_curve(vazão_s_importado,vazão_i_importado,vazão_paralelo_teorico,vazão_paralelo_pratico,H_s_importado,H_i_importado,H_paralelo_teorico,H_paralelo_pratico):
    plt.figure(figsize=(10, 6))
    plt.plot(vazão_i_importado, H_i_importado, marker='o', linestyle='--', color='b', label='Bomba Inferior')
    plt.plot(vazão_s_importado, H_s_importado, marker='h', linestyle='-', color='r', label='Bomba Superior')
    plt.plot(vazão_paralelo_teorico, H_paralelo_teorico, marker='o', linestyle='--', color='green', label='Série Teórico')
    plt.plot(vazão_paralelo_pratico, H_paralelo_pratico, marker='h', linestyle='-', color='black', label='Série Prático')
    # Adicionando títulos e rótulos
    plt.title('Gráfico Altura(H) e Vazão(Q)')
    plt.xlabel('Q(m^3/h)')
    plt.ylabel('H(m)')
    plt.legend()
    plt.grid(True)
    plt.show()

#-------------INCERTEZAS DO PAQ-----------#
i_paq = 0.03/1000 #m
i_regua = 0.5/1000 #m
i_vazão = 500 #L/h
#---------------CONDIÇÕES AMBIENTE-------------#
temperatura_lab = ufloat(22, 0.5) #°C
densidade_mercurio = ufloat(13600,2720) #kg/m^3
h_baro = ufloat(0.696,0.001) #m
gravidade = ufloat(9.784,0.005) #m/s^2
pressao_lab = gravidade*h_baro*densidade_mercurio #Pa(ufloat)
densidade_agua_CP = CP.PropsSI("D", "T", temperatura_lab.nominal_value + 273.15, "P", pressao_lab.nominal_value, 'water')
densidade_agua = ufloat(densidade_agua_CP, 0)
densidade_ar_CP = CP.PropsSI("D", "T", temperatura_lab.nominal_value + 273.15, "P", pressao_lab.nominal_value, 'air') #kg/m^3
densidade_ar = ufloat(densidade_ar_CP,0)
viscosidade_ar_CP = CP.PropsSI("V", "T", temperatura_lab.nominal_value + 273.15, "P", pressao_lab.nominal_value, 'air')
viscosidade_ar = ufloat(viscosidade_ar_CP,0)
#-------------ALTURAS---------------#

#--------------ASSOCIAÇÃO EM SÉRIE------------#
pressão_recalque_serie = [5.6,5.3,5.1,4.75,2.9] #atm
vazão_serie = [3000,4000,5000,6500,10500] #L/h
vazão_serie_graf = [a/1000 for a in vazão_serie] # m^3/h
desnivel_serie = 'calcular'
diferença_de_nível_dos_man_serie = [48.6/100,48.2/100,47.1/100,45.8/100,35.7/100] #cm -> m

H_serie_teorico = [(a+b) for a,b in zip(H_i,H_s)]
H_serie_graf_teorico = [a.nominal_value for a in H_serie_teorico]
H_serie_pratico = [(((P*pressao_lab)/(densidade_agua*gravidade)) + alt) for P,alt in zip(pressão_recalque_serie,diferença_de_nível_dos_man_serie)]
H_serie_graf_pratico = [a.nominal_value for a in H_serie_pratico]
#------------ASSOCIAÇÃO EM PARALELO-----------#
vazão_paralelo = [10500,6500,5000,4000,3000] #L/h
medição_vacuometro = [-0.05,-0.05,-0.04,-0.03,-0.03] #atm
medição_man_bomba2 = [2.55,2.8,2.85,2.9,2.95] #atm
diferença_de_nível_dos_man_paralelos = [46/100,48.5/100,49.4/100,49.8/100,50.4/100] #cm -> m
lista_a = [2,2.5,3,4,7]

vazão_paralelo_pratico_graf = [2*a/1000 for a in vazão_paralelo] # m^3/h
vazão_paralelo_teorico_graf = [(a+b) for a,b in zip(vazão_convertida_s,vazão_convertida_i)]
H_paralelo_pratico = [((((P_r - P_s)*pressao_lab)/(densidade_agua*gravidade)) + alt - a) for P_r,P_s,alt,a in zip(medição_man_bomba2,medição_vacuometro,diferença_de_nível_dos_man_paralelos,lista_a[::-1])]
H_paralelo_graf_pratico = [a.nominal_value for a in H_paralelo_pratico]
H_paralelo_teorico_graf = H_i_graf
#-----------RESULTADOS-------------#
print(f'Resultados em paralelo:')
print(f'H teórico: {H_i}')
print(f'H prático: {H_paralelo_pratico[::-1]}')
print(f'Q teórico: {vazão_paralelo_teorico_graf}')
print(f'Q prático: {vazão_paralelo_pratico_graf[::-1]}')
print(f'Resultados em serie:')
print(f'H teórico: {H_serie_teorico}')
print(f'H prático: {H_serie_pratico}')
print(f'Q teórico: {vazão_serie_graf}')
print(f'Q prático: {vazão_serie_graf}')
plot_series_curve(vazão_convertida_s,vazão_convertida_i,vazão_serie_graf,vazão_serie_graf,H_s_graf,H_i_graf,H_serie_graf_teorico,H_serie_graf_pratico)
plot_parallel_curve(vazão_convertida_s,vazão_convertida_i,vazão_paralelo_teorico_graf,vazão_paralelo_pratico_graf,H_s_graf,H_i_graf,H_paralelo_teorico_graf,H_paralelo_graf_pratico)
#H_i,H_s -> Conjunto de ufloats para operações
#H_i_graf, H_s_graf -> Valores nominais usados na plotagem dos dados
#vazão_convertida_i, vazão_convertida_s -> Vazões associadas aos pontos anteriores
#[2.55,2.8,2.85,2.9,2.95]
#[1.45,2.4,2.55,2.65,2.8] -> Resolveu manipulado