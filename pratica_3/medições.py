import numpy as np

diametro_interno_placa_de_orificio = 47.02 #mm
diametro_tubo_marrom = 74.76 #mm
distancia_man_antes_da_placa = 82.07 #mm -> distância da parede do tubo até a placa de orificio
distancia_man_depois_da_placa = 45.73 #mm -> ''
diametro_depois = 8.64 #mm
diametro_antes = 4.73 #m

ponto_antes = distancia_man_antes_da_placa - (diametro_antes/2)
ponto_depois = distancia_man_depois_da_placa - (diametro_depois/2)
i_paq = 0.03 #mm
dh_15Hz = [9.07, 10.37] #mm
dh_30Hz = [39.07, 39.79] #mm
dh_45Hz = [94.66, 94.89] #mm
dh_60Hz = [171.5, 172] #mm -> medido na régua