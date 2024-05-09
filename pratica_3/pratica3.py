import matplotlib.pyplot as plt
import numpy as np
import math
import propagação_de_incertezas as P_I
from uncertainties import ufloat
import CoolProp.CoolProp as CP

i_paq = 0.03 #mm
i_regua = 0.5 #mm

#-----------------DADOS LAB--------------#
temperatura_com_incerteza = ufloat(26,0.5) #°C
gravidade = ufloat(9.784,0.005) #m/s^2
h_baro = ufloat(0.694,0.001) #m
densidade_mercurio = ufloat(13600,2720) #kg/m^3
pressao_lab = gravidade*h_baro*densidade_mercurio #Pa

massa_especifica_agua_CP = CP.PropsSI("D", "T", temperatura_com_incerteza.nominal_value + 273.15, "P", pressao_lab.nominal_value, 'Water')
massa_especifica_agua = ufloat(massa_especifica_agua_CP, 0) #kg/m^3

#-----------------DADOS------------------#
diametro_interno_placa_de_orificio = ufloat(47.02/1000,i_paq/1000) #mm
diametro_tubo_marrom = ufloat(74.76/1000,i_paq/1000) #mm
distancia_man_antes_da_placa = ufloat(82.07/1000,i_paq/1000) #mm -> distância da parede do tubo até a placa de orificio
distancia_man_depois_da_placa = ufloat(45.73/1000,i_paq/1000) #mm -> ''
diametro_depois = ufloat(8.64/1000,i_paq/1000) #mm
diametro_antes = ufloat(4.73/1000,i_paq/1000) #mm
area_placa_de_orificio = ufloat(math.pi,0)*(diametro_interno_placa_de_orificio/2)**2

ponto_antes = distancia_man_antes_da_placa - (diametro_antes/2)
ponto_depois = distancia_man_depois_da_placa - (diametro_depois/2)
dh_15 = [9.07, 10.37] #mm
dh_30 = [39.07, 39.79] #mm
dh_45 = [94.66, 94.89] #mm
dh_60 = [171.5, 172] #mm -> medido na régua
#-----------------INCERTEZAS--------------------#
n_de_medidas = 2
confiança = 0.9545
k_15Hz = P_I.calcular_k(confiança,n_de_medidas)
k_30Hz = P_I.calcular_k(confiança,n_de_medidas)
k_45Hz = P_I.calcular_k(confiança,n_de_medidas)
k_60Hz = P_I.calcular_k(confiança,n_de_medidas)

dh_15_i = P_I.medida_com_incerteza(dh_15,k_15Hz,i_paq)/1000
dh_30_i = P_I.medida_com_incerteza(dh_30,k_30Hz,i_paq)/1000
dh_45_i = P_I.medida_com_incerteza(dh_45,k_45Hz,i_paq)/1000
dh_60_i = P_I.medida_com_incerteza(dh_60,k_60Hz,i_regua)/1000

#-------------DADOS DO AR------------#
viscosidade_ar_CP = CP.PropsSI("V", "T", temperatura_com_incerteza.nominal_value + 273.15, "P", pressao_lab.nominal_value, 'Air')
densidade_ar_CP = CP.PropsSI("D", "T", temperatura_com_incerteza.nominal_value + 273.15, "P", pressao_lab.nominal_value, 'Air') #kg/m^3
densidade_ar = ufloat(densidade_ar_CP,0)
viscosidade_ar = ufloat(viscosidade_ar_CP,0)

dP_15Hz = dh_15_i*massa_especifica_agua*gravidade
dP_30Hz = dh_30_i*massa_especifica_agua*gravidade
dP_45Hz = dh_45_i*massa_especifica_agua*gravidade
dP_60Hz = dh_60_i*massa_especifica_agua*gravidade

#-----------HIPOTESE------------#
R = 287 # J/kg*K
T_c = temperatura_com_incerteza.nominal_value + 273.15
k = 1.4
v_som = math.sqrt(k*R*T_c)

P_r = (pressao_lab.nominal_value/1000)/3400000 #usando a P_c do Nitrogênio pois 78% do ar atm é
#-------CALCULO DA VAZÃO------------#
beta = ufloat(round(diametro_interno_placa_de_orificio.nominal_value/diametro_tubo_marrom.nominal_value,4),
            round(diametro_interno_placa_de_orificio.std_dev/diametro_tubo_marrom.std_dev,4))
Re = ufloat(1000,0)
Re_prev = ufloat(0,0)
vazao_massica_prev = ufloat(0,0)
vazao_massica = ufloat(0,0)
L2 = 0.47
M2 = ((2*L2)/1-beta)
def calcular_C(Re, beta):
    L1 = 1
    L2 = 0.47
    M2 = ((2 * L2) / (1 - beta))
    A = (19000 * beta / Re) ** 0.8
    C = (0.5961 + (0.0261 * beta ** 2) - (0.216 * beta ** 8) + 0.000521 * (10 ** (6) * beta / Re) ** 0.7
         + (0.0188 + 0.0063 * A) * (beta ** 3.5) * (10 ** 6 / Re) ** 0.3
         + (0.043 + 0.080 * math.exp(-10 * L1) - 0.123 * math.exp(-7 * L1)) * (1 - 0.11 * A) * (
                     beta ** 4 / (1 - beta ** 4))
         - 0.031 * (M2 - 0.8 * M2 ** 1.1) * beta ** 1.3 + 0.011 * (0.75 - beta) * (
                    2.8 - diametro_tubo_marrom / 25.4))
    return C

n_iter = 0
while (Re.nominal_value - Re_prev.nominal_value) == 0:
    vazao_massica_prev = vazao_massica
    coeficiente = calcular_C(Re, beta)
    vazao_massica = ((coeficiente * area_placa_de_orificio) * (2 * densidade_ar * dP_15Hz) ** (1 / 2)) / (
                1 - beta ** 4) ** (1 / 2)
    Re_prev = Re
    Re = (4 * vazao_massica) / (viscosidade_ar * math.pi * diametro_interno_placa_de_orificio)
    n_iter += 1
    print('----------------------------')
    print(f'Vazão: {vazao_massica}')
    print(f'Re previo: {Re_prev}')
    print(f'Re: {Re}')
    print(f'Iterações: {n_iter}')
    print(f'---------------------------')


#-----------CONFERIR RESULTADOS----------#
print(f'Beta: {beta}')
print(f'Área da placa de orificio: {area_placa_de_orificio} m^2')
print(f'Alturas: {dh_15_i}, {dh_30_i}, {dh_45_i}, {dh_60_i} m')
print(f'Pressões: {dP_15Hz}, {dP_30Hz}, {dP_45Hz}, {dP_60Hz} Pa')
print(f'C: {calcular_C(Re,beta)}')