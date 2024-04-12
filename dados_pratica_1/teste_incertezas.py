from uncertainties import ufloat
from scipy.stats import t
import numpy as np

# Valores medidos e incerteza
valores_medidos = [10.1, 10.2, 10.3, 10.2, 10.0]
incerteza = 0.1

# Calcular o número de medidas
n_medidas = len(valores_medidos)

# Calcular a média e desvio padrão dos valores medidos
media = np.mean(valores_medidos)
desvio_padrao = np.std(valores_medidos)

# Graus de liberdade para a distribuição t de Student
graus_liberdade = n_medidas - 1

# Calcular k para uma confiabilidade de p = 0.9545
confiabilidade = 0.9545
probabilidade_t = 1 - (1 - confiabilidade) / 2
k = t.ppf(probabilidade_t, graus_liberdade)

# Calcular a incerteza expandida
incerteza_expandida = k * (desvio_padrao / np.sqrt(n_medidas))

# Resultados
print(f"Valores medidos: {valores_medidos}")
print(f"Incerteza: {incerteza}")
print(f"Média: {media}")
print(f"Desvio padrão: {desvio_padrao}")
print(f"Graus de liberdade: {graus_liberdade}")
print(f"Valor de k para p = {confiabilidade}: {k}")
print(f"Incerteza expandida: {incerteza_expandida}")

