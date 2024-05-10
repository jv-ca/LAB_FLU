from CoolProp.CoolProp import PropsSI #importa as propriedades termodinâmicas no Sistema Internacional
from math import pi,sqrt #importa o número pi e a função raiz quadradada da biblioteca matemática

D_1=150E-3 #Diâmetro Interno da Linha de Ar em m
D_t=100E-3 #Diâmetro da Placa de Orifício em m
P_1=600E3 #Pressão a montante em Pa
T_ar=25+273.15 #Temperatura do ar em K
h=750E-3 #Deflexão do manômetro a H2O em m
T_amb=23+273.15 #Temperatura ambiente em K
P_atm=92E3 #Pressão atmosférica em Pa
g=9.81 #gravidade em m/s^2

rho_água=PropsSI("D","P",P_atm,"T",T_amb,'water')
DELTA_P=rho_água*g*h
beta=D_t/D_1
rho_ar=PropsSI("D","P",P_1,"T",T_ar,'air')
mu_ar=PropsSI("V","P",P_1,"T",T_ar,'air')
A_t=0.25*pi*D_t**2

def f(vazão_real): #definição de uma função que satisfaz a equação da vazão mássica real
  Re=(4*vazão_real)/(mu_ar*pi*D_1)
  C=0.5959+0.0312*beta**2.1-0.184*beta**8+(91.71*beta**2.5)/(Re**0.75) #atenção: a operação de potenciação é programada em python como **, não como ^, assim como em diversas outras linguagens
  f=vazão_real-C*A_t*sqrt(2*rho_ar*(DELTA_P))/sqrt(1-beta**4)
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
    if abs(f(vazão_real))<0.0000001:break #define a condição de convergência. o loop será encerrado apenas quando o valor da função seja ínfimo, já que a vazão procurada é a raiz dessa função.
  return vazão_real_new

def Q(vazão_real): #definição da função vazão volumétrica
  Q=vazão_real()/rho_ar
  return Q

print("vazão_real=","%.3f" % vazão_real(),"kg/s")
print("vazão_volumétrica=","%.3f" % Q(vazão_real),"m^3/s")
print(f'Densidade do ar: {rho_ar}')
