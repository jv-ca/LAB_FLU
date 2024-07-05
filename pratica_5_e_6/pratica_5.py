import matplotlib.pyplot as plt
import numpy as np
import propagação_de_incertezas as P_I
from uncertainties import ufloat
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
from math import pi,sqrt,exp 
from uncertainties import ufloat

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


#-------------INCERTEZAS DO PAQ-----------#
i_paq = 0.03/1000 #m
i_regua = 0.5/1000 #m
i_vazão = 500/(3.6*10**6) #m^3/s
#---------------CONDIÇÕES AMBIENTE-------------#
temperatura_lab = ufloat(23, 0.5) #°C
densidade_mercurio = ufloat(13600,2720) #kg/m^3
h_baro = ufloat(0.699,0.001) #m
gravidade = ufloat(9.784,0.005) #m/s^2
pressao_lab = gravidade*h_baro*densidade_mercurio #Pa(ufloat)
densidade_agua_CP = CP.PropsSI("D", "T", temperatura_lab.nominal_value + 273.15, "P", pressao_lab.nominal_value, 'water')
densidade_agua = ufloat(densidade_agua_CP, 0)
densidade_ar_CP = CP.PropsSI("D", "T", temperatura_lab.nominal_value + 273.15, "P", pressao_lab.nominal_value, 'air') #kg/m^3
densidade_ar = ufloat(densidade_ar_CP,0)
viscosidade_ar_CP = CP.PropsSI("V", "T", temperatura_lab.nominal_value + 273.15, "P", pressao_lab.nominal_value, 'air')
viscosidade_ar = ufloat(viscosidade_ar_CP,0)
#----------MEDIÇÕES BOMBA SUPERIOR--------#
pressões_de_recalque_s = [2.7,2.6,2.5,2.35,1.9] #atm
pressões_de_sucção_s = [-0.05,-0.05,-0.05,-0.05,-0.05] #atm
vazão_volumétrica_s = [3000,4000,5000,6500,8750] #L/h
vazão_convertida_s = [round(v_v_s/(10**3),5) for v_v_s in vazão_volumétrica_s] #m^3/h
altura_vacuometro = 16.4 #cm, a partir do meio da tubulação
cota_manometro_bourbon = 1
cota_vacuometro = 1
diferença_de_nivel_dos_manometros_s = 36.3/100 #cm -> m, incerteza de tudo e mais um pouco

H_s = [(((P_r - P_s)*pressao_lab)/(densidade_agua*gravidade) + diferença_de_nivel_dos_manometros_s) for P_r,P_s in zip(pressões_de_recalque_s,pressões_de_sucção_s)]
H_s_graf = [a.nominal_value for a in H_s]
#----------MEDIÇÕES BOMBA INFERIOR----------#
pressões_de_recalque_i = [2.7,2.6,2.5,2.35,1.9] #atm
vazão_volumétrica_i = [3000,4000,5000,6500,8850] #L/h
vazão_convertida_i = [round(v_v_i/(10**3),5) for v_v_i in vazão_volumétrica_i] #m^3/h
centro_do_tubo_ate_o_centro_do_recalque = 51.4/100 #cm -> m
desnível_de_agua_ = [48.7/100,47.6/100,46.9/100,45.1/100,40.1/100] #cm -> m

H_i = [(((P*pressao_lab)/(densidade_agua*gravidade)) + alt) for P,alt in zip(pressões_de_recalque_i,desnível_de_agua_)]
H_i_graf = [a.nominal_value for a in H_i]

#----------RESULTADOS---------------#
show_results(vazão_convertida_s,vazão_convertida_i, H_s, H_i)
#plot_all_curves(vazão_convertida_s, vazão_convertida_i, H_s_graf, H_i_graf)