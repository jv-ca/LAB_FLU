import matplotlib.pyplot as plt
import numpy as np
import math
import propagação_de_incertezas as P_I
from uncertainties import ufloat
import CoolProp.CoolProp as CP

#---------------PONTOS DE MEDIÇÃO-------------#
pontos =[]
pontos_de_medição = [0.044,0.146,0.296,0.5,0.704,0.854,0.956] #%
diametro_do_tubo_pitot = 85.14 #mm
diametro_interno = 3.21 #mm
diametro_do_tubo_externo = 75.51 #mm -> 0.07551
R_marrom = (diametro_do_tubo_externo/2)/1000 #m
for ponto in pontos_de_medição:
    pontos.append(((diametro_do_tubo_externo*ponto)-(diametro_interno/2))/1000) #m
pontos = [round(num,6) for num in pontos] #m

#------------------DADOS DE PRESSÃO--------------#
dh_30 = [2.32,3.03,3.58,3.63,3.06,2.76,2.20] #cm, e a index representa o ponto de medição
dh_45 = [5.53,6.96,7.64,8.16,6.90,6.27,5.18] #cm
dh_60 = [9.54,11.61,12.45,13.95,12.28,11.25,8.93] #cm
"""
dh_30_i = [ufloat(ponto,i_paq) for ponto in dh_30]
dh_45_i = [ufloat(ponto,i_paq) for ponto in dh_45]
dh_60_i = [ufloat(ponto,i_paq) for ponto in dh_60]
"""
#-----------------INCERTEZAS--------------------#
n_de_medidas = 1
confiança = 0.9545
i_paq = 0.03 #mm
k_30Hz = P_I.calcular_k(confiança,n_de_medidas)
k_45Hz = P_I.calcular_k(confiança,n_de_medidas)
k_60Hz = P_I.calcular_k(confiança,n_de_medidas)

dh_30_i = [(P_I.medida_com_incerteza(h/1000,k_30Hz,i_paq/1000)) for h in dh_30]
dh_45_i = [(P_I.medida_com_incerteza(h/1000,k_45Hz,i_paq/1000)) for h in dh_45]
dh_60_i = [(P_I.medida_com_incerteza(h/1000,k_60Hz,i_paq/1000)) for h in dh_60]

dh_30_i = [ufloat(round(a1.nominal_value,4),round(a1.std_dev,4)) for a1 in dh_30_i] #m
dh_45_i = [ufloat(round(a2.nominal_value,4),round(a2.std_dev,4)) for a2 in dh_45_i] #m
dh_60_i = [ufloat(round(a2.nominal_value,4),round(a2.std_dev,4)) for a2 in dh_60_i] #m

#-------------DADOS LAB ------------------#
temperatura_com_incerteza = ufloat(26.5,0.5) #°C
gravidade = ufloat(9.784,0.005) #m/s^2
h_baro = ufloat(0.695,0.001) #m
densidade_mercurio = ufloat(13600,2720) #kg/m^3
pressao_lab = gravidade*h_baro*densidade_mercurio #Pa
massa_especifica_agua_CP = CP.PropsSI("D", "T", temperatura_com_incerteza.nominal_value + 273.15, "P", pressao_lab.nominal_value, 'Water')
massa_especifica_agua = ufloat(massa_especifica_agua_CP, 0) #kg/m^3

#dP = ro*g*dH
dP_30Hz_arrendondar = [(h1*massa_especifica_agua*gravidade) for h1 in dh_30_i]
dP_45Hz_arrendondar = [(h2*massa_especifica_agua*gravidade) for h2 in dh_45_i]
dP_60Hz_arrendondar = [(h3*massa_especifica_agua*gravidade) for h3 in dh_60_i]
dP_30Hz = [ufloat(round(pressao.nominal_value,4),round(pressao.std_dev,4)) for pressao in dP_30Hz_arrendondar] #Pa
dP_45Hz = [ufloat(round(pressao.nominal_value,4),round(pressao.std_dev,4)) for pressao in dP_45Hz_arrendondar] #Pa
dP_60Hz = [ufloat(round(pressao.nominal_value,4),round(pressao.std_dev,4)) for pressao in dP_60Hz_arrendondar] #Pa

viscosidade_ar_CP = CP.PropsSI("V", "T", temperatura_com_incerteza.nominal_value + 273.15, "P", pressao_lab.nominal_value, 'Air')
densidade_ar_CP = CP.PropsSI("D", "T", temperatura_com_incerteza.nominal_value + 273.15, "P", pressao_lab.nominal_value, 'Air') #kg/m^3
densidade_ar = ufloat(densidade_ar_CP,0)
viscosidade_ar = ufloat(viscosidade_ar_CP,0)

#-------------VELOCIDADES------------------#
# V = (2*(delta_P)/densidade_agua)**1/2
v_30Hz_arredondar = [(((2*dP)/densidade_ar)**(1/2)) for dP in dP_30Hz]
v_45Hz_arredondar = [(((2*dP)/densidade_ar)**(1/2)) for dP in dP_45Hz]
v_60Hz_arredondar = [(((2*dP)/densidade_ar)**(1/2)) for dP in dP_60Hz]
v_30Hz = [ufloat(round(v.nominal_value,4),round(v.std_dev,4)) for v in v_30Hz_arredondar]
v_45Hz = [ufloat(round(v.nominal_value,4),round(v.std_dev,4)) for v in v_45Hz_arredondar]
v_60Hz = [ufloat(round(v.nominal_value,4),round(v.std_dev,4)) for v in v_60Hz_arredondar]

#-----------HIPOTESE------------#
R = 287 # J/kg*K
T_c = temperatura_com_incerteza.nominal_value + 273.15
k = 1.4
v_som = math.sqrt(k*R*T_c)

P_r = (pressao_lab.nominal_value/1000)/3400000 #usando a P_c do Nitrogênio pois 78% do ar atm é
#-----------REYNOLDS------------#
#Re = densidade*V*diametro/viscosidade
Re_30Hz = (densidade_ar*(v_30Hz[3]/2)*(diametro_do_tubo_externo/1000))/viscosidade_ar
Re_45Hz = (densidade_ar*(v_45Hz[3]/2)*(diametro_do_tubo_externo/1000))/viscosidade_ar
Re_60Hz = (densidade_ar*(v_60Hz[3]/2)*(diametro_do_tubo_externo/1000))/viscosidade_ar
Re = [Re_30Hz,Re_45Hz,Re_60Hz]

#-----------FUNÇÕES DO PERFIL DE VELOCIDADE--------------#
conjunto_de_pontos_y = [(round(R_marrom-ponto,4)) for ponto in pontos[0:3:1]]
u_u_max_inverter = [(round(vel.nominal_value/v_60Hz[3].nominal_value,3)) for vel in v_60Hz[-1:3:-1]]
u_u_max = u_u_max_inverter[::-1]
r_R_marro = [round(r/R_marrom,3) for r in conjunto_de_pontos_y[::-1]]

log_u_u_max = [math.log(valor) for valor in u_u_max]
log_r_R_marro = [math.log(valor) for valor in conjunto_de_pontos_y]

slopes = [0.2293, 0.1889, 0.2180]
n_fit = [1/s for s in slopes]
n_eq = [(-1.7 + 1.8*(math.log(rr.nominal_value)))/3 for rr in Re]

#-----------CONFERIR RESULTADOS----------#
print(f'Pontos: {pontos} m')
print(f'Alturas 30Hz: {dh_30_i} m')
print(f'Pressões em 60Hz: {dP_60Hz} Pa')
print(f'Pressão ambiente: {pressao_lab/1000} kPa')
print(f'Velocidade para 30Hz: {v_30Hz} m/s')
print(f'Velocidade para 45Hz: {v_45Hz} m/s')
print(f'Velocidade para 60Hz: {v_60Hz} m/s')
print(f'Características do ar: {viscosidade_ar}(viscosidade), {densidade_ar}(densidade)')
print(f'Reynolds: {Re_30Hz}(30Hz), {Re_45Hz}(45Hz), {Re_60Hz}(60Hz)')
print(f'Velocidade do som: {v_som}')
print(f'Incompressibilidade: {round(v_30Hz[3].nominal_value/v_som,4)}(Ma_30Hz), {round(v_45Hz[3].nominal_value/v_som,4)}(Ma_45Hz), {round(v_60Hz[3].nominal_value/v_som,4)}(Ma_60Hz)')
print(f'u/u_max: {log_u_u_max}')
print(f'Y: {conjunto_de_pontos_y}')
print(f'r/R: {log_r_R_marro}')
print(f'n = {n_fit}(fit), {n_eq}(eq)')

''''
plt.figure(figsize=(8, 6))
plt.plot(u_u_max, r_R_marro, marker='o', linestyle='-')
plt.xlabel('u/u_max')
plt.ylabel('r/R')
plt.title('u/u_max vs. r/R')
plt.grid(True)


plt.figure(figsize=(8, 6))
plt.plot(log_r_R_marro,log_u_u_max, marker='o', linestyle='-')
plt.xlabel('log(r/R)')
plt.ylabel('log(u/u_max)')
plt.title('log(y/R) vs. log(u/u_max)')
plt.grid(True)

# Calculando o ajuste linear usando NumPy
slope, intercept = np.polyfit(log_r_R_marro, log_u_u_max, 1)
y_fit = slope * np.array(log_r_R_marro) + intercept
print(slope)
# Plotando o gráfico com os dados e a linha de ajuste
plt.figure(figsize=(8, 6))
plt.scatter(log_r_R_marro, log_u_u_max, label='Dados Originais')
plt.plot(log_r_R_marro, y_fit, color='red', label='Ajuste Linear')
plt.xlabel('log(r/R)')
plt.ylabel('log(u/u_max)')
plt.title('Ajuste Linear')
plt.legend()
plt.grid(True)
equation_text = f'Equação do Ajuste Linear:\n y = {slope:.2f}x + {intercept:.2f}'
plt.text(0.5, min(log_u_u_max) - 0.5, equation_text, fontsize=12, bbox=dict(facecolor='lightgray', alpha=0.5), horizontalalignment='right')
plt.show()
'''