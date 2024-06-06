import matplotlib.pyplot as plt
import numpy as np
import propagação_de_incertezas as P_I
from uncertainties import ufloat
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
from math import pi,sqrt,exp 
from uncertainties import ufloat

i_paq = 0.03/1000 #m
i_regua = 0.5/1000 #m
#---------------DADOS DO FABRICANTE-------------#
vazão_fabricante = 6/60 #m^3/s
H_fabricante = 0.5 #m
rot_fabricante = 3438 #rpm
potência_fabricante = 2*745.7 #W
#---------------CONDIÇÕES AMBIENTE-------------#
temperatura_lab = ufloat(24.5, 0.5) #°C
densidade_mercurio = ufloat(13600,2720) #kg/m^3
h_baro = ufloat(0.694,0.001) #m
gravidade = ufloat(9.784,0.005) #m/s^2
pressao_lab = gravidade*h_baro*densidade_mercurio #Pa(ufloat)
#----------------DADOS DA TUBULAÇÃO-------------#
diametro_interno_placa_de_orificio = ufloat(47.02/1000,i_paq) #m
diametro_tubo_marrom = ufloat(74.76/1000,i_paq) #m
freq_30 = ufloat(30,0)
freq_55 = ufloat(55,0)
rot_30 = freq_30*60 #1800 RPM
rot_55 = freq_55*60 #3300 RPM
#---------------MEDIÇÕES--------------#
d_altura_30_man_calculo_H = [74.12/1000,69.17/1000,64.2/1000,59.25/1000,54.28/1000] #medido em mm, convertido em m 
d_altura_30_man_calculo_Q = [9.75/1000,20.75/1000,33.69/1000,44.48/1000,54.48/1000] # medido em mm, convertido em m

d_altura_55_man_calculo_H = [254.5/1000,236/1000,211/1000,192/1000,177/1000] #medido em mm, convertido em m-> com a régua
d_altura_55_man_calculo_Q = [30.21/1000,65.87/1000,114.6/1000,156.5/1000,193/1000] #medido em mm, convertido em m -> 5° ponto com a régua
#--------------RENDIMENTO-------------#
densidade_agua_CP = CP.PropsSI("D", "T", temperatura_lab.nominal_value + 273.15, "P", pressao_lab.nominal_value, 'water')
densidade_agua = ufloat(densidade_agua_CP, 0)
peso_especifico_agua = densidade_agua*gravidade #ufloat
rend_maq = (peso_especifico_agua*vazão_fabricante*H_fabricante)/potência_fabricante #ufloat
#-------------RENDIMENTO POR RATEUX---------------#
rend_30 = 1-(1-rend_maq)*(rot_30/rot_fabricante)**(0.1)
rend_55 = 1-(1-rend_maq)*(rot_55/rot_fabricante)**(0.1)
#-------------ALTURAS MANOMETRICAS---------------#
densidade_ar_CP = CP.PropsSI("D", "T", temperatura_lab.nominal_value + 273.15, "P", pressao_lab.nominal_value, 'air') #kg/m^3
densidade_ar = ufloat(densidade_ar_CP,0)

dP_H_30 = [(densidade_agua*gravidade*ufloat(d,i_paq)) for d in d_altura_30_man_calculo_H]
H_30 = [(d/(densidade_ar*gravidade)) for d in dP_H_30]

dP_H_55 = [(densidade_agua*gravidade*ufloat(d,i_paq)) for d in d_altura_55_man_calculo_H]
H_55 = [(d/(densidade_ar*gravidade)) for d in dP_H_55]
#------------CALCULO DAS VAZÕES VOLUMÉTRICAS-----------#
A_placa_de_orificio = 0.25*pi*diametro_interno_placa_de_orificio**2 #ufloat
beta = diametro_interno_placa_de_orificio/diametro_tubo_marrom #ufloat
viscosidade_ar_CP = CP.PropsSI("V", "T", temperatura_lab.nominal_value + 273.15, "P", pressao_lab.nominal_value, 'air')
viscosidade_ar = ufloat(viscosidade_ar_CP,0)
#-----------VAZÃO PRA 1° FREQUÊNCIA-------------#
dP_Q_30 = [(densidade_agua*gravidade*ufloat(d,i_paq)) for d in d_altura_30_man_calculo_Q]
#-----------VAZÃO PRA 2° FREQUÊNCIA-------------#
dP_Q_55 = [(densidade_agua*gravidade*ufloat(d,i_paq)) for d in d_altura_55_man_calculo_Q]
#-----------ITERAÇÃO PRA DETERMINAR A VAZÃO-----------------------#
D_1=diametro_tubo_marrom #Diâmetro Interno da Linha de Ar em m
D_t=diametro_interno_placa_de_orificio #Diâmetro da Placa de Orifício em m
T_ar=temperatura_lab+273.15 #Temperatura do ar em K
T_amb=temperatura_lab+273.15 #Temperatura ambiente em K
P_atm=pressao_lab #Pressão atmosférica em Pa
g=gravidade #gravidade em m/s^2
DELTA_P=dP_Q_55[4] # -> trocar para a vazão de cada altura do manometro
beta=D_t/D_1
rho_ar=densidade_ar
mu_ar=viscosidade_ar
A_t=A_placa_de_orificio
#A PARTIR DAQUI, É PRECISO FAZER OS DADOS 1 POR 1

def f(vazão_real): #definição de uma função que satisfaz a equação da vazão mássica real
  Re=(4*vazão_real)/(mu_ar*pi*D_1)
  C=0.5959+0.0312*beta**2.1-0.184*beta**8+(91.71*beta**2.5)/(Re**0.75) #atenção: a operação de potenciação é programada em python como **, não como ^, assim como em diversas outras linguagens
  f=vazão_real-C*A_t*(2*rho_ar*(DELTA_P))**(1/2)/(1-beta**4)**(1/2)
  return f

def df(vazão_real): #definição da derivada da função que satisfaz a equação da vazão mássica real
  h=0.00001
  df=(f(vazão_real+h)-f(vazão_real))/h #derivada numérica da função
  return df

#Newton-Raphson para cálculo da vazão:
def vazão_real():
  vazão_real=1 #estimativa inicial da vazão
  for i in range(50): #o comando "for" é utilizado para criar um loop. nesse caso, a função será calculada inicialmente com 1 termo, n=1, depois será testada a tolerância especificada. caso não seja atendida, serão utilizados dois termos, n=2, e assim sucessivamente até que o critério seja atendido.
  #observe que o comando acima permite até 50 iterações. caso não convirja com esse número, provavelmente algo está errado na programação do método.
    vazão_real_new=vazão_real-f(vazão_real)/df(vazão_real) #expressão matemática de recorrência do método de Newton-Raphson
    if abs(f(vazão_real_new))<0.00000001:break #define a condição de convergência. o loop será encerrado apenas quando o valor da função seja ínfimo, já que a vazão procurada é a raiz dessa função.
  return vazão_real_new

def Q(vazão_real): #definição da função vazão volumétrica
  Q=vazão_real()/rho_ar
  return Q

#------------VAZÕES DEPOIS DA ITERAÇÃO------------#
Q_30 = [ufloat(0.015170,0.000033),ufloat(0.02213,0.00004),ufloat(0.02820,0.00005),ufloat(0.03240,0.00005),ufloat(0.03586,0.00006)] #m^3/s
Q_55 = [ufloat(0.02670,0.00004),ufloat(0.04263,0.00007),ufloat(0.05200,0.00008),ufloat(0.06077,0.00009),ufloat(0.06749,0.0001)] #m^3/s
#------------CALCULO DAS POTÊNCIAS---------------#
peso_especifico_ar = densidade_ar*gravidade
potência_30Hz = [((peso_especifico_ar*vazao*H/rend_30)) for vazao,H in zip(Q_30,H_30)]
potência_55Hz = [((peso_especifico_ar*vazao*H/rend_55)) for vazao,H in zip(Q_55,H_55)]
#------------CALCULO POR RATEAUX-----------------#
H_55_R = [(((rot_55/rot_30)**(2))*H) for H in H_30]
Q_55_R = [(((H_55_/H_30_)**(1/2))*vazao) for H_55_,H_30_,vazao in zip(H_55_R,H_30,Q_30)]
potência_55_R = [(((H_55_/H_30_)**(3/2))*pot) for H_55_,H_30_,pot in zip(H_55_R,H_30,potência_30Hz)]
#------------FORMATAÇÃO DAS LISTAS PARA OS GRAFICOS------------#
H_30_graf = [a.nominal_value for a in H_30]
Q_30_graf = [a.nominal_value for a in Q_30]
potência_30_graf = [a.nominal_value for a in potência_30Hz]

H_55_graf = [a.nominal_value for a in H_55]
Q_55_graf = [a.nominal_value for a in Q_55]
potência_55_graf = [a.nominal_value for a in potência_55Hz]

H_55_R_graf = [a.nominal_value for a in H_55_R]
Q_55_R_graf = [a.nominal_value for a in Q_55_R]
potência_55_R_graf = [a.nominal_value for a in potência_55_R]
#------------VERIFICAR RESULTADOS----------------#
print(f'Condições do Laboratório: {temperatura_lab} °C, {round(pressao_lab.nominal_value,2)}+/-{round(pressao_lab.std_dev,2)} Pa')
print(f'Rendimento da maquina: {rend_maq}')
print(f'Rendimento por Rateaux: {rend_30}(30Hz), {rend_55}(55Hz)')
print(f'Densidades: {round(densidade_ar.nominal_value,2)}(Ar), {round(densidade_agua.nominal_value,2)}(Água), {round(densidade_mercurio.nominal_value,2)}(Hg) kg/m^3')
print(f'Viscosidades: {viscosidade_ar}')
print(f'beta: {beta}')
print(f'vazão_real= {vazão_real()} kg/s')
print(f'vazão_volumétrica= {Q(vazão_real)} m^3/s')
print(f'H: {H_30}(30Hz), {H_55}(55Hz) m')
print(f'Q: {Q_30}(30Hz), {Q_55}(55Hz)')
print(f'Potência de 30Hz: {potência_30Hz} W')
print(f'Potência de 55Hz: {potência_55Hz} W')
print(f'H de 55 por Rateaux: {H_55_R}')
print(f'Q de 55 por Rateaux: {Q_55_R}')
print(f'Potência de 55 por Rateaux: {potência_55_R}')

plt.figure(figsize=(10, 6))
plt.plot(Q_55_R_graf, H_55_R_graf, marker='o', linestyle='--', color='b', label='H_55_R vs Q_55_R')
plt.plot(Q_30_graf, H_30_graf, marker='h', linestyle='-', color='r', label='H_30 vs Q_30')
plt.plot(Q_55_graf, H_55_graf, marker='h', linestyle='-', color='g', label='H_55 vs Q_55')
# Adicionando títulos e rótulos
plt.title('Gráfico Altura(H) e Vazão(Q)')
plt.xlabel('Q(m^3/s)')
plt.ylabel('H(m)')
plt.legend()
plt.grid(True)

plt.figure(figsize=(10, 6))
plt.plot(Q_55_R_graf, potência_55_R_graf, marker='o', linestyle='--', color='b', label='P_55_R vs Q_55_R')
plt.plot(Q_30_graf, potência_30_graf, marker='h', linestyle='-', color='r', label='P_30 vs Q_30')
plt.plot(Q_55_graf, potência_55_graf, marker='h', linestyle='-', color='g', label='P_55 vs Q_55')
# Adicionando títulos e rótulos
plt.title('Gráfico Potência(N) e Vazão(Q)')
plt.xlabel('Q(m^3/s)')
plt.ylabel('N(W)')
plt.legend()
plt.grid(True)
# Mostrando o gráfico
plt.show()