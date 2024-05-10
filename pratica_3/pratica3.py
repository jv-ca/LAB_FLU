import matplotlib.pyplot as plt
import numpy as np
import propagação_de_incertezas as P_I
from uncertainties import ufloat
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
from math import pi,sqrt,exp #importa o número pi e a função raiz quadradada da biblioteca matemática
from uncertainties import ufloat

i_paq = 0.03 #mm
i_regua = 0.5 #mm

#-----------------------PRATICA--------------------------#
##DADOS LAB
temperatura_com_incerteza = ufloat(26,0.5) #°C
gravidade = ufloat(9.784,0.005) #m/s^2
h_baro = ufloat(0.694,0.001) #m
densidade_mercurio = ufloat(13600,2720) #kg/m^3
pressao_lab = gravidade*h_baro*densidade_mercurio #Pa

massa_especifica_agua_CP = CP.PropsSI("D", "T", temperatura_com_incerteza.nominal_value + 273.15, "P", pressao_lab.nominal_value, 'water')
massa_especifica_agua = ufloat(massa_especifica_agua_CP, 0) #kg/m^3

#-----------------DADOS------------------#
diametro_interno_placa_de_orificio = ufloat(47.02/1000,i_paq/1000) #mm
diametro_tubo_marrom = ufloat(74.76/1000,i_paq/1000) #mm
distancia_man_antes_da_placa = ufloat(82.07/1000,i_paq/1000) #mm -> distância da parede do tubo até a placa de orificio
distancia_man_depois_da_placa = ufloat(45.73/1000,i_paq/1000) #mm -> ''
diametro_antes = ufloat(8.64/1000,i_paq/1000) #mm
diametro_depois = ufloat(4.73/1000,i_paq/1000) #mm
area_placa_de_orificio = ufloat(pi,0)*(diametro_interno_placa_de_orificio/2)**2
area_antes = ufloat(pi,0)*(diametro_antes/2)**2
area_depois = ufloat(pi,0)*(diametro_depois/2)**2

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
dh_60_i = P_I.medida_com_incerteza(dh_60,k_60Hz,i_paq)/1000

dP_15Hz = dh_15_i*massa_especifica_agua*gravidade
dP_30Hz = dh_30_i*massa_especifica_agua*gravidade
dP_45Hz = dh_45_i*massa_especifica_agua*gravidade
dP_60Hz = dh_60_i*massa_especifica_agua*gravidade

#-------------DADOS DO AR------------#
viscosidade_ar_CP = CP.PropsSI("V", "T", temperatura_com_incerteza.nominal_value + 273.15, "P", pressao_lab.nominal_value, 'air')
densidade_ar_CP = CP.PropsSI("D", "T", temperatura_com_incerteza.nominal_value + 273.15, "P", pressao_lab.nominal_value, 'air') #kg/m^3
densidade_ar = ufloat(densidade_ar_CP,0)
viscosidade_ar = ufloat(viscosidade_ar_CP,0)

#----------------------HIPOTESE BERNOULI----------------#

v_15Hz_arredondar = (((2*dP_15Hz/densidade_ar*(1-(area_depois/area_antes)**2)))**(1/2))
v_30Hz_arredondar = (((2*dP_30Hz/densidade_ar*(1-(area_depois/area_antes)**2)))**(1/2))
v_45Hz_arredondar = (((2*dP_45Hz/densidade_ar*(1-(area_depois/area_antes)**2)))**(1/2))
v_60Hz_arredondar = (((2*dP_60Hz/densidade_ar*(1-(area_depois/area_antes)**2)))**(1/2))

v_15Hz = ufloat(round(v_15Hz_arredondar.nominal_value,4),round(v_15Hz_arredondar.std_dev,4))
v_30Hz = ufloat(round(v_30Hz_arredondar.nominal_value,4),round(v_30Hz_arredondar.std_dev,4))
v_45Hz = ufloat(round(v_45Hz_arredondar.nominal_value,4),round(v_45Hz_arredondar.std_dev,4))
v_60Hz = ufloat(round(v_60Hz_arredondar.nominal_value,4),round(v_60Hz_arredondar.std_dev,4))

R = 287 # J/kg*K
T_c = temperatura_com_incerteza.nominal_value + 273.15
k = 1.4
v_som = sqrt(k*R*T_c)
#------------------------VARIAVEIS DA ROTINA DE PYTHON-------------#
D_1= diametro_tubo_marrom #Diâmetro Interno da Linha de Ar em m
D_t= diametro_interno_placa_de_orificio #Diâmetro da Placa de Orifício em m
P_1= pressao_lab #Pressão a montante em Pa
T_ar= temperatura_com_incerteza #Temperatura do ar em K
h= dh_60_i #Deflexão do manômetro a H2O em m
T_amb= temperatura_com_incerteza #Temperatura ambiente em K
P_atm= pressao_lab #Pressão atmosférica em Pa
g= gravidade #gravidade em m/s^2
rho_água= massa_especifica_agua
DELTA_P= rho_água*g*dh_60_i #ALTERAR AQUI PARA TER DIFERENTES VAZÕES
beta= D_t/D_1
rho_ar= densidade_ar
mu_ar= viscosidade_ar
A_t= area_placa_de_orificio
L1 = 1
L2 = 0.47
M2 = ((2 * L2) / (1 - beta))

def f(vazão_real): #definição de uma função que satisfaz a equação da vazão mássica real
  Re=(4*vazão_real)/(densidade_ar*pi*D_1)
  A = (19000 * beta / Re) ** 0.8
  C = (0.5961 + (0.0261 * beta ** 2) - (0.216 * beta**8) + 0.000521 * (10 ** (6) * beta / Re) ** 0.7
         + (0.0188 + 0.0063 * A) * (beta ** 3.5) * (10 ** 6 / Re) ** 0.3
         + (0.043 + 0.080 * exp(-10 * L1) - 0.123 * exp(-7 * L1)) * (1 - 0.11 * A) * (beta ** 4 / (1 - beta ** 4))
         - 0.031 * (M2 - 0.8 * M2 ** 1.1) * beta ** 1.3 + 0.011 * (0.75 - beta) * (2.8 - diametro_tubo_marrom / 25.4))
  f=(vazão_real-C*A_t*((2*rho_ar*DELTA_P)**(0.5)))/((1-beta**4)**(0.5))
  return f

def df(vazão_real): #definição da derivada da função que satisfaz a equação da vazão mássica real
  h=0.00001
  df=(f(vazão_real+h)-f(vazão_real))/h #derivada numérica da função
  return df

#Newton-Raphson para cálculo da vazão:
def vazão_real():
  vazão_real=1 #estimativa inicial da vazao
  for i in range(50): #o comando "for" é utilizado para criar um loop. nesse caso, a função será calculada inicialmente com 1 termo, n=1, depois será testada a tolerância especificada. caso não seja atendida, serão utilizados dois termos, n=2, e assim sucessivamente até que o critério seja atendido.
  #observe que o comando acima permite até 50 iterações. caso não convirja com esse número, provavelmente algo está errado na programação do método.
    vazão_real_new=vazão_real-f(vazão_real)/df(vazão_real) #expressão matemática de recorrência do método de Newton-Raphson
    if abs(f(vazão_real))<0.0000001:break #define a condição de convergência. o loop será encerrado apenas quando o valor da função seja ínfimo, já que a vazão procurada é a raiz dessa função.
  return vazão_real_new

def Q(vazão_real): #definição da função vazão volumétrica
  Q= vazão_real()/densidade_ar
  return Q

print(f'Beta: {beta}')
print(f'Condições Ambientes: {round(pressao_lab.nominal_value,6)}+/-{round(pressao_lab.std_dev,4)} Pa, {temperatura_com_incerteza} °C')
print(f'Densidade da água: {massa_especifica_agua} kg/m^3')
print(f'Densidade do ar: {densidade_ar} kg/m^3')
print(f'Área da placa de orificio: {area_placa_de_orificio} m^2')
print(f'Alturas: {dh_15_i}, {dh_30_i}, {dh_45_i}, {dh_60_i} m')
print(f'Pressões: {dP_15Hz}, {dP_30Hz}, {dP_45Hz}, {dP_60Hz} Pa')
print(f'Velocidades: {v_15Hz}, {v_30Hz}, {v_45Hz}, {v_60Hz} m/s')
print(f'Velocidade do som: {round(v_som,3)} m/s')
print(f'Mach: {round(v_15Hz.nominal_value/v_som,4)}, {round(v_30Hz.nominal_value/v_som,4)}, {round(v_45Hz.nominal_value/v_som,4)}, {round(v_60Hz.nominal_value/v_som,4)}')
print(f'Vazão real = {vazão_real()} kg/s')
print(f'Vazão volumétrica = {Q(vazão_real)} m^3/s')
