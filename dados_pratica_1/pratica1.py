import propagação_de_incertezas as P_I
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
from uncertainties import ufloat
import math

#---------IMPORT DE DADOS---------#
n_de_medidas = 3
confiança = 0.9545

k_10Hz= P_I.calcular_k(confiança,n_de_medidas)
k_20Hz= P_I.calcular_k(confiança,n_de_medidas)
k_30Hz= P_I.calcular_k(confiança,n_de_medidas)
k_10Hz_inclinado = P_I.calcular_k(confiança,n_de_medidas)

dados_altura_10Hz = P_I.ler_dados("alturas10Hz.txt")
media_altura_10Hz = P_I.medida_com_incerteza(dados_altura_10Hz,k_10Hz)/100 #m

dados_altura_20Hz = P_I.ler_dados("alturas20Hz.txt")
media_altura_20Hz = P_I.medida_com_incerteza(dados_altura_20Hz,k_20Hz)/100 #m

dados_altura_30Hz = P_I.ler_dados("alturas30Hz.txt")
media_altura_30Hz = P_I.medida_com_incerteza(dados_altura_30Hz,k_30Hz)/100 #m

dados_altura_10Hz_inclinado = P_I.ler_dados("altura_inclinado_10Hz.txt")
media_altura_10Hz_inclinado = P_I.medida_com_incerteza(dados_altura_10Hz_inclinado,k_10Hz_inclinado)/100

#----------------VARIÁVEIS GERAIS------------------#
temperatura_lab_com_incerteza = ufloat(26.5,0.5) #°C
gravidade = ufloat(9.784,0.005) #m/s^2
h_baro = ufloat(0.693, 0.001) #m
densidade_mercurio = ufloat(13600,680) #kg/m^3
pressao_lab = gravidade*h_baro*densidade_mercurio #Pa -> ro_merc*g*h_baro
pressao_inmet = ufloat(91820,0) #Pa

#---------------ÂNGULO------------------#
h_angulo = ufloat(0.012,0.001) #m
l_angulo = ufloat(0.163,0.001) #m
argumento_arcoseno = h_angulo/l_angulo #ufloat
angulo_nominal_rad = math.asin(argumento_arcoseno.nominal_value) #rad
angulo_nominal = math.degrees(angulo_nominal_rad) #graus

incerteza_argumento = argumento_arcoseno.std_dev
derivada_arcsin = 1 / math.sqrt(1 - angulo_nominal_rad**2) #graus -> incerteza para a função arcsin
angulo_incerteza = derivada_arcsin*incerteza_argumento
angulo = ufloat(angulo_nominal,angulo_incerteza) #graus

print(angulo,angulo_nominal,angulo_incerteza)

#----------------CALCULO DE PRESSÃO VERTICAL-----------------#
massa_especifica_agua_CP = CP.PropsSI("D", "T", temperatura_lab_com_incerteza.nominal_value + 273.15, "P", pressao_lab.nominal_value, 'Water')
massa_especifica_agua = ufloat(massa_especifica_agua_CP, 0)
P_ar_10Hz = pressao_lab + (massa_especifica_agua*gravidade*media_altura_10Hz)
P_ar_20Hz = pressao_lab + (massa_especifica_agua*gravidade*media_altura_20Hz)
P_ar_30Hz = pressao_lab + (massa_especifica_agua*gravidade*media_altura_30Hz)

#-------------------CALCULO DE PRESSAO INCLINADO--------------#
altura_inclinado = media_altura_10Hz_inclinado*math.sin(math.radians(angulo.nominal_value))
massa_especifica_alcool = ufloat(786,39.3) #kg/m^3
P_ar_10Hz_inclinado = pressao_lab + (gravidade*massa_especifica_alcool*altura_inclinado)

print("Pressão atmosférica: {:.4f} Pa| Densidade: {} kg/m^3".format(pressao_lab,massa_especifica_agua))
print("Altura em 10Hz: {} m| Pressão em 10Hz: {:.4f} kPa".format(media_altura_10Hz,P_ar_10Hz/1000))
print("Altura em 20Hz: {} m| Pressão em 20Hz: {:.4f} kPa".format(media_altura_20Hz, P_ar_20Hz/1000))
print("Altura em 30Hz: {} m| Pressão em 30Hz: {:.4f} kPa".format(media_altura_30Hz, P_ar_30Hz/1000))
print("-------------------------------------------------------------------------")
print("Comprimento em 10 Hz: {} m| Pressão em 10Hz: {:.4f} kPa".format(media_altura_10Hz_inclinado, P_ar_10Hz_inclinado/1000))

#---------------PLOTAGEM------------------#
frequencias = [10, 20, 30]
pressoes_kPa = [P_ar_10Hz.nominal_value / 1000, P_ar_20Hz.nominal_value / 1000, P_ar_30Hz.nominal_value / 1000]
alturas_vertical = [media_altura_10Hz.nominal_value,media_altura_20Hz.nominal_value,media_altura_30Hz.nominal_value]
# Plotagem
plt.figure(figsize=(10, 8))  # Define o tamanho da figura
plt.plot(frequencias, pressoes_kPa, marker='H')  # Plota a pressão em relação à frequência

# Adiciona rótulos e título
plt.xlabel('Frequência (Hz)')
plt.ylabel('Pressão (kPa)')
plt.title('Pressão(kPa) vs Frequência(Hz)')

plt.figure(figsize=(10, 8))  # Define o tamanho da figura
plt.plot(alturas_vertical, pressoes_kPa, marker='H')  # Plota a pressão em relação à frequência

plt.xlabel('Altura (m)')
plt.ylabel('Pressão (kPa)')
plt.title('Pressão(kPa) vs Altura(m)')

# Mostra o gráfico
plt.show()
