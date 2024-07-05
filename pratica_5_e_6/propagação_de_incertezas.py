import numpy as np
import uncertainties as UN
from scipy.stats import t
dados = []
#---------------------LEITURA DO TXT-------------------------#
def ler_dados(nome_arquivo):
    #abrir o arquivo txt no modo leitura
    arquivo = open(nome_arquivo, 'r') 
    #ler as linhas e armazenar cada linha em um elemento de uma lista, considerando todos os caracteres
    linhas = arquivo.readlines() 

    #editar as linhas lidas em outra lista, dessa vez com formatação
    linhas_processadas = [linha.strip() for linha in linhas]

    #separar somente os números
    separador = [parte.split(": ")[1] for parte in linhas_processadas]

    #converter cada termo para float
    dados = [float(parte) for parte in separador]
    return dados

def calcular_k(nivel_confianca, graus_liberdade):
    # Calcular a probabilidade de distribuição t de Student
    probabilidade_t = 1 - nivel_confianca
    
    # Calcular o fator de abrangência k (inverso da probabilidade t)
    k = t.ppf(1 - (probabilidade_t / 2), graus_liberdade)
    return k

# Definindo os valores medidos com suas incertezas
def medida_com_incerteza(dados,k,i):
    incerteza_do_instrumento = i
    U = incerteza_do_instrumento*k
    dados_com_incerteza = UN.ufloat(np.mean(dados), U) 
    return dados_com_incerteza


