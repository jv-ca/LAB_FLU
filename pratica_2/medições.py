import numpy as np
#------------PONTOS DE MEDIÇÃO---------------#
pontos_de_medição = [0.044,0.146,0.296,0.5,0.704,0.854,0.956]
diametro_do_tubo_pitot = 85.14 #mm
diametro_interno = 3.21 #mm
diametro_do_tubo_externo = 75.51 #mm
distancia = []
for ponto in pontos_de_medição:
    print("Ponto {}: {:.5} mm".format(pontos_de_medição.index(ponto)+1,(diametro_do_tubo_externo*ponto)-(diametro_interno/2)))

temperatura_da_sala = 25 #°C
#d_interno = 85.14 mm
#-----------PRESSÕES----------#

#--------------INCERTEZAS----------#
i_regua = 1 #m
i_paquimetro = 0.03 #mm
i_temperatura = 1 #°C

#---------------DADOS-----------#
dh_30 = [2.32,3.03,3.58,3.63,3.06,2.76,2.20] #cm, e a index representa o ponto de medição
dh_45 = [5.53,6.96,7.64,8.16,6.90,6.27,5.18] #cm
dh_60 = [9.54,11.61,12.45,13.95,12.28,11.25,8.93] #cm

print(diametro_do_tubo_pitot/2,diametro_do_tubo_externo/2)