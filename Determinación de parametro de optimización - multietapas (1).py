import numpy as np
import matplotlib.pyplot as plt

g_0 = 9.81
n = int(input("Inserte el numero de etapas: ")) #Numero de etapas
ISP = [] #Se crean las listas correspondientes al impulso específico y a la razon de masa estructural de cada etapa
EPS_pre = [] 
F_star = []
M_0 = []
M_P = []
MU_pre = []
GAMMA = []

dv = float(input("Inserte el delta de velocidad total (en m/s): ")) #Se llenan las listas mencionadas con anterioridad

for i in range(1,n+1):
    isp = float(input(f"Inserte el valor del impulso especifico para la etapa %d: " %i))
    ISP.append(isp*g_0)
    f_star = float(input("Ingrese la fuerza de propulsión total de la etapa %d (en kN): " %i))
    F_star.append(f_star)
    eps = float(input("Inserte la masa estructural de la etapa %d (en kg): " %i))
    EPS_pre.append(eps)
    m_0 = float(input("Inserte la masa inicial de la etapa %d (en kg): " %i))
    M_0.append(m_0)
    m_p = float(input("Ingrese la masa de propelente de la etapa %d (en kg): " %i))
    M_P.append(m_p)
    mu = float(input("Ingrese la masa final de la etapa %d (en kg): " %i))
    MU_pre.append(mu)
    gamma = float(input("Ingrese el ángulo de trayectoria de vuelo (en grados) de la etapa %d: " %i))
    GAMMA.append(gamma*np.pi/180)
    i+=1

EPS = []
MU = []

for i in range(n):
    razon_masa_total = MU_pre[i]/M_0[i]
    MU.append(razon_masa_total)
    razon_masa_es = EPS_pre[i]/(EPS_pre[i]+M_P[i])
    EPS.append(razon_masa_es)
    i+=1

class Cohete: #Se define la clase cohete, la cual se construye unicamente con los atributos pedidos al usuario con anterioridad.
    def __init__(self,n,ISP,EPS,dv,F_star,M_0,MU,GAMMA):
        self.etapas = n #Número de etapas
        self.isp = ISP #Valores de impulso especifico
        self.eps = EPS #Valores de razón de masa estructural
        self.dv = dv #Valor de delta v requerido
        self.fuerza_prop = F_star
        self.M_0 = M_0
        self.MU = MU
        self.GAMMA = GAMMA
        
    def delta_v(self,x): #Se define la funcion de optimizacion a la cual se le quieren encontrar las raices (los valores de optimizacion)
        delta = 0
        for i in range(n):
            delta = delta + self.isp[i]*np.log(np.abs((self.isp[i]-x)/(self.eps[i]*self.isp[i])))
            i += 1
        return dv-delta
    
    def integral_gravedad(self):
        t_combustion = []
        sin_gamma = []
        valor_ints = []
        for i in range(self.etapas):
            delta_t = self.M_0[i]*(1-self.MU[i])*g_0*self.isp[i]/(self.fuerza_prop[i]*1e+3)
            t_combustion.append(delta_t)
            s_gamma = np.sin(self.GAMMA[i])
            sin_gamma.append(s_gamma)
            i+=1
        for i in range(self.etapas):
            val_ints = sin_gamma[i]*t_combustion[i]*g_0
            valor_ints.append(val_ints)
        return valor_ints
    
    def delta_v_gravedad(self,x):
        delta = 0
        for i in range(n):
            delta = delta + self.isp[i]*np.log(np.abs((self.isp[i]-x)/(self.eps[i]*self.isp[i])))-self.integral_gravedad()[i]
            i += 1
        return dv-delta
    
    def d_delta_v(self,x): #Se obtiene la derivada numerica de la funcion anterior
        h = 1e-6
        return ((self.delta_v(x+h)-self.delta_v(x))/h)
    
    def d_delta_v_gravedad(self,x):
        h = 1e-6
        return (self.delta_v_gravedad(x+h)-self.delta_v_gravedad(x))/h
    
    def Newton_Raphson_sin_gravedad(self,x_0): #Se implementa un metodo de Newton-Raphson para la determinacion de las raices.
        n_i = 100
        datos=[]
        iteraciones = []
        for i in range(n_i):
            x = x_0 - self.delta_v(x_0)/self.d_delta_v(x_0)
            x_0 = x
        plt.plot(iteraciones,datos)
        plt.xlabel("Iteración")
        plt.ylabel("Convergencia")
        return x
    
    def Newton_Raphson_con_gravedad(self,x_0):
        n_i = 100
        datos=[]
        iteraciones = []
        for i in range(n_i):
            x = x_0 - self.delta_v_gravedad(x_0)/self.d_delta_v_gravedad(x_0)
            x_0 = x
            iteraciones.append(i)
            datos.append(x)
        plt.plot(iteraciones,datos)
        plt.xlabel("Iteración")
        plt.ylabel("Convergencia")
        return x
    
    def lambda_op(self,x_0):
        varphi_g_not = self.Newton_Raphson_sin_gravedad(x_0)
        varphi_g_aye = self.Newton_Raphson_con_gravedad(x_0)
        
        lambda_op_g_not = []
        lambda_op_g_aye = []
        
        for i in range(self.etapas):
            lambda_g_not = varphi_g_not*self.eps[i]/(self.isp[i]*g_0*(1-self.eps[i])-varphi_g_not)
            lambda_op_g_not.append(lambda_g_not)
            lambda_g_aye = varphi_g_aye*self.eps[i]/(self.isp[i]*g_0*(1-self.eps[i])-varphi_g_aye)
            lambda_op_g_aye.append(lambda_g_aye)
        
        def_g_not = 1
        def_g_aye = 1
        
        for i in range(self.etapas):
            def_g_not = def_g_not * lambda_op_g_not[i]/(1+lambda_op_g_not[i])
            def_g_aye = def_g_aye * lambda_op_g_aye[i]/(1+lambda_op_g_aye[i])
        
        return def_g_not, def_g_aye
    
x_0 = 4000 #Se plantea como posible valor inicial aquel recomendado por Danby(2003).
cohete_1 = Cohete(n,ISP,EPS,dv,F_star,M_0,MU,GAMMA) #Se crea la instancia de cohete
varphi = cohete_1.Newton_Raphson_sin_gravedad(x_0) #Se obtiene el valor de optimizacion
varphi_prime = cohete_1.Newton_Raphson_con_gravedad(x_0)

lambda_op_g_not = cohete_1.lambda_op(x_0)[0]
lambda_op_g_aye = cohete_1.lambda_op(x_0)[1]

print("Parámetro principal de optimización sin gravedad: ",varphi, "y con gravedad: ",varphi_prime,". Los valores de razón de masa óptimos son, para el caso sin gravedad: ",lambda_op_g_not,"y para el caso con gravedad: ",lambda_op_g_aye) #Tarea completa.
