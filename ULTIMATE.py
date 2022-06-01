# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 23:35:37 2021

@author: Abel Amado Gonzalez Bernad
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter1d

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'


def probability(R, r, n_try, n_bins, iterations):
    
    def generateC(R,r):
        p_2 = (R-r)**2
        mal = True
        while mal == True:
            xr_c, yr_c = np.random.uniform(-(R-r), R-r), np.random.uniform(-(R-r), R-r)
            if (xr_c)**2 + (yr_c)**2 < p_2:
                mal = False
        return xr_c, yr_c

    def generatePto(R):
        r_p = np.random.uniform(-R, R)
        return r_p
    
    
    def comprobacion(xr_c, yr_c, r_p, r):
        dentro = False
        if (r_p - xr_c)**2 + (yr_c)**2 < r**2:
               dentro = True
        return dentro
    
    probr = []
    #miss = 0
    no_deberia = 0
    for i in range(int(n_try)):
        
        xr_c, yr_c = generateC(R, r)
        r_p = generatePto(R)
        r_p1=r_p
        salir = False
        cont = 0
        R1=R
        r1=r
        #p1=p #es el chi
        
        while salir == False:
                                   
            if comprobacion(xr_c, yr_c, r_p1, r1):
                
                #Nueva distacia al centro
                r_p1 = np.sqrt((r_p1 - xr_c)**2 + (yr_c)**2)
                
                #Nuevo circulo grande
                R1=r1
                
                #Nuevo circulo pequeño
                r1 = 1/k * R1
                
                #Nuevo chi
                #p1 = R1-r1
                if r_p1 >R1:
                    print('No deberia')
                #Generar nuevo circulo, de acuerdo a la nueva disposicion
                xr_c, yr_c = generateC(R1, r1)
               
                #Sumamos 1 al contador que ha sobrevivido a la 1a iteracion
                cont = cont+1
                

                #print(r1)
            else:
                if abs(r_p)<99 and abs(r_p)>-99: #NO necesario, solo comprobacion
                    no_deberia = no_deberia+1
                salir = True
                
            if cont == iterations and salir == False:
                probr.append(r_p)
                salir = True
    
    n, bins, patches = plt.hist(probr, bins = n_bins, weights =  [n_bins*100/n_try] * len(probr), range=(-R, R), 
                                label = r'Hist. '+ str(iterations)+'º cierre')
    
    plt.xlabel(r"$r$ / u.a.", size=15)
    plt.ylabel(r"$P(r)$", size=15)
    
    
    plt.xlim(-230, 380)
    plt.hlines(100, -300, 300, colors='r',linestyle= 'dotted')
    plt.show()
    return n, bins
    #return(probr, weights)



def Iteracion(R, r, p, N, Mitera, Nitera):
    
    """
    Return:
        Probabilidades: array de longitud 2 con [0] = primera curva (primer cierre) & [1]= segunda curva (2 cierres)
    """
    """
    if R== abs(p):
        return np.array([0, 0])
    """
    def ProbPto(R, r, p):
        
        if p > R:
            print('Fallo en ProbPto')
            return False
        
        def P_ (R, r, p):
            
            """
            Meterle argumentos: R->Radio grande, r-> radio pequeño, p -> pto (valor) ó ptos (array)
            en los que se quiere saber el valor de la función de prbabilidad.
            
            SOLO VALIDO PARA R=2r
            """
            
            p = abs(p)
            P_ = (2*r**2*np.arctan(( (r-p/2)/(r+p/2) )**(1/2))-p*(r**2-(p**2)/4)**(1/2) + 
                  (np.pi*r**2)/2-r**2*np.arctan( p/(2*(r**2-(p**2)/4)**(1/2))))*100/(np.pi*r**2);
            return P_
        
        def P(R,r, p):
            """
            Meterle argumentos: R->Radio grande, r-> radio pequeño, p -> pto (valor) ó ptos (array)
            en los que se quiere saber el valor de la función de prbabilidad.
            
            SOLO VALIDO PARA R=!2r en la region de corte de circulos (es decir, la 'tapa' no se contempla)
            """
            chi=(R-r)
            b = ((chi)**2 -(p-r)**2)/(2*p) 
            c=(p-r)
            #a=(chi-p+r)
            P = (np.sqrt(b*(2*r-b))*(2*(r**2)*((np.arctan(np.sqrt(b/(2*r-b))))/(np.sqrt(b*(2*r-b)))) 
                 -r +b ) +(chi**2)*np.pi/2-(c+b)*np.sqrt(chi**2 -(c+b)**2)
                 -(chi**2)*np.arctan(( c+b )/(np.sqrt(chi**2 -(c+b)**2))))*100/(np.pi*chi**2)
            return P
            
        p = np.abs(p)
        chi=(R-r)
        if p == R:
            return 0
        
        if R == 2*r:
            P_ = P_(R, r, p)
            return P_
        
        elif R< 2*r:
            chi=(R-r)
            if p<=R-2*chi:
                return 100
            else:
                return P(R, r, p)
        else: #R > 2*r
            if p<=R-2*r:
                return (r/chi)**2 *100
            else:
                return P(R, r, p)    
       
    def IzLim1(x):
        return np.sqrt((r1)**2-(x-r1)**2)
        
    def IzLim2(x):
            
        return -np.sqrt((r1)**2-(x-r1)**2)
        
    def DerLim1(x):
        c = r_p1 -r1
        chi = R1-r1
        return np.sqrt(chi**2-(x+c)**2)
        
    def DerLim2(x):
        c = r_p1 -r1
        chi = R1-r1
        return -np.sqrt(chi**2-(x+c)**2) 
        
    def ProbPtoMAT(y, x, R1, r1 ):
        #Las defino aqui para que la funcion sea matematica, solo depender de la posicion:
        """
        Debe estar estar en el entorno en el que R1 y r1 son pequeñas
        """
        p = np.sqrt((x-(R1))**2+y**2) #donde pone R1 (ò k*r1) es  para poner ahi el radio, r, anterior
        #print('Susodicha distancia', p)
        chi=(R1-r1)
        if p == R1:
            return 0
        
        def P_ (R, r, p):
            """
            Meterle argumentos: R->Radio grande, r-> radio pequeño, p -> pto (valor) ó ptos (array)
            en los que se quiere saber el valor de la función de prbabilidad.
            
            SOLO VALIDO PARA R=2r
            """
            p = abs(p)
            P_ = (2*r**2*np.arctan(( (r-p/2)/(r+p/2) )**(1/2))-p*(r**2-(p**2)/4)**(1/2) + 
                  (np.pi*r**2)/2-r**2*np.arctan( p/(2*(r**2-(p**2)/4)**(1/2))))*100/(np.pi*r**2);
            return P_
        
        def P(R,r, p):
            """
            Meterle argumentos: R->Radio grande, r-> radio pequeño, p -> pto (valor) ó ptos (array)
            en los que se quiere saber el valor de la función de prbabilidad.
            
            SOLO VALIDO PARA R=!2r en la region de corte de circulos (es decir, la 'tapa' no se contempla)
            """
            chi=(R-r)
            b = ((chi)**2 -(p-r)**2)/(2*p) 
            c=(p-r)
            #a=(chi-p+r)
            
            P = (np.sqrt(b*(2*r-b))*(2*(r**2)*((np.arctan(np.sqrt(b/(2*r-b))))/(np.sqrt(b*(2*r-b)))) 
                 -r +b ) +(chi**2)*np.pi/2-(c+b)*np.sqrt(chi**2 -(c+b)**2)
                 -(chi**2)*np.arctan(( c+b )/(np.sqrt(chi**2 -(c+b)**2))))*100/(np.pi*chi**2)
            return P
        
        if R1 == 2*r1:
            P_ = P_(R1, r1, p)
            return P_
        
        elif R1< 2*r1:
            chi=(R1-r1)
            if p<=R1-2*chi:
                return 100
            else:
                return P(R1, r1, p)
            
        else: #R > 2*r
            if p<=R1-2*r1:
                return (r1/chi)**2 *100
            else:
                return P(R1, r1, p)
    """
    def GeneratePTO(x1, x2, f1, f2):
        radio = abs(x2-x1) #el chi
        salir = False
        while salir == False:
            x, y = np.random.uniform(x1, x2), np.random.uniform(-radio, +radio)
            if y> f1(x) or y<f2(x):
                salir = True
        return y, x
    """
    def GeneratePTO(x1, x2):
        chi = (x2-x1)/2
        pto_c =x1+chi
        mal = True
        while mal == True:
            x, y = np.random.uniform(x1, x2), np.random.uniform(-chi, chi)
            if (x-pto_c)**2 + (y)**2 < chi**2:
                mal = False 
        return x, y

    def GeneratePTO2(beta, alpha, R, r, p):
        """
        Genera puntos en la superficie de interseccion de circulos
        """
        def IzLim1(x):
            return np.sqrt((r)**2-(x-r)**2)
        
        def IzLim2(x):
            return -np.sqrt((r)**2-(x-r)**2)
        
        def DerLim1(x):
            c = p -r
            chi = R-r
            if chi**2-(x+c)**2<0:
                print(p)
            return np.sqrt(chi**2-(x+c)**2)
        
        def DerLim2(x):
            c = p -r
            chi = R-r
            if chi**2-(x+c)**2<0:
                print(p)
            return -np.sqrt(chi**2-(x+c)**2)   
        
        salir = False
        chi = R-r
        while salir == False:
            x, y = np.random.uniform(0, alpha), np.random.uniform(-chi, chi) #Esta y_min & y_max como limites solo valen en R<2R o R=2r 
            
            if x < beta and y < IzLim1(x) and y > IzLim2(x) or x>beta and y < DerLim1(x) and y > DerLim2(x):
                salir = True
        return x, y
    
    def GeneratePTO3(erre):
        """
        Para cuando el r esta totalmente contenido en el chi 
        Solo se da esta situacion en ciertos distancias si R>2r
        """
        mal = True
        while mal == True:
            x, y = np.random.uniform(-erre, erre), np.random.uniform(-erre, erre)
            if x**2 + y**2 < erre**2:
                mal = False 
        return x, y
    
    R1, r1 = R, r
    k = R/r
    chi1=R1-r1
    r_p1 = abs(p)
    
    if Mitera < Nitera:
        print('Mitera < Nitera')
        return 9999999999999999999999
    
    if abs(p) == R:
        return 0
    
    if Mitera == 1:
        return ProbPto(R, r, p)
        
    if abs(r_p1) <= r1-chi1: #Si el ciruclo CHI esta totalmente contenido en circulo r
        
        chi1 = R1-r1
        alpha1 = chi1 - r_p1 + r1
        alphita = r1-r_p1-chi1
 
        #Toca integrar: lo haremos montecarlo
        VALORES = np.array([])
        
        for vez in range(N):
            x, y = GeneratePTO(alphita, alpha1)
            #plt.plot(x, y, 'o')
            p1 = np.sqrt((x-r)**2 + y**2) #•nueva distancia entre ubicacion y punto aleatorio (CREO)
            if Mitera>Nitera+1:
                #print(Mitera, Nitera)
                VALORES = np.append(VALORES, Iteracion(R1/k, r1/k, p1, N, Mitera, Nitera+1)*ProbPtoMAT(y,x, r1, r1/k)/100)
            else:
                VALORES = np.append(VALORES, ProbPtoMAT(y,x, r1, r1/k))
                
        MEDIA = np.mean(VALORES) 
        #print('La media es', MEDIA, 'para el pto', r_p1)
        if Nitera ==1:
            return MEDIA*ProbPto(R, r, p)/100
        else:
            return MEDIA
        
    elif abs(r_p1) +r1 < chi1 and R1>2*r1: #Chi y r no intersectan: solo posible si R>2r
        #print('Hola')
        VALORES = np.array([])
        
        for vez in range(N):
            x, y = GeneratePTO3(r1)
            #plt.plot(x, y, 'o')
            #plt.plot(x, y, 'o')
            
            p1 = np.sqrt(x**2 + y**2) #nueva distancia entre ubicacion y punto aleatorio (CREO)
            if Mitera>Nitera+1:
                
                #print(Mitera, Nitera)
                VALORES = np.append(VALORES, Iteracion(R1/k, r1/k, p1, N, Mitera, Nitera+1)*ProbPto(R1/k, r1/k, p1)/100)
            else:
                VALORES = np.append(VALORES, ProbPto(R1/k, r1/k, p1))
                    
        MEDIA = np.mean(VALORES) 
        #print('La media es', MEDIA, 'para el pto', r_p1)
        if Nitera ==1:
            return MEDIA*ProbPto(R, r, p)/100
        else:
            return MEDIA     
    
    else:
        
        chi1= R1-r1
        beta1 = (chi1**2 -(r_p1-r1)**2)/(2*r_p1)
        alpha1 = chi1 - r_p1 + r1
        
        #El area nueva, interseccion ciruclos (ya no hace falta, que no hacemos la integral)
        #area = ProbPto(R1, r1, r_p1)*(np.pi*chi1**2)/100 
        #print('Area: ', area)
        
        VALORES = np.array([])
        
        #print('R1 = ', R1, 'y r1 = ', r1)
        
        #Puntos o evaluaciones de la funcion para hacer la 'integral'
        for vez in range(N):
            x, y = GeneratePTO2(beta1, alpha1, R1, r1, r_p1)
            #plt.plot(x, y, 'o')
            #print(x, y)
            p1 = np.sqrt((x-r)**2 + y**2) #•nueva distancia entre ubicacion y punto aleatorio (CREO)
            if Mitera>Nitera+1:
                #print(Mitera, Nitera)
                VALORES = np.append(VALORES, Iteracion(R1/k, r1/k, p1, N, Mitera, Nitera+1)*ProbPtoMAT(y,x, r1, r1/k)/100)
            else:
                VALORES = np.append(VALORES, ProbPtoMAT(y,x, r1, r1/k))
        MEDIA = np.mean(VALORES) 
        #print('La media es', MEDIA, 'para el pto', r_p1)
        if Nitera ==1:
            return MEDIA*ProbPto(R, r, p)/100
        else:
            return MEDIA
                
plt.clf()                

R, r = 200, 130
k =R/r
#p=90
#Potente= para que salga un poco bien aproximado
N = 2000 #PARA CALCULOS DE 1 Y 2 CIERRES. Potente: 3000
N_ = 100 #PARA CALCULOS DEL 3 CIERRE Potente: 300
N__ = 10#Potente: 100 
Mitera1 = 1 #Cantidad de iteraciones M  de maximo
Mitera2 = 2 #Cantidad de iteraciones M  de maximo
Mitera3 = 3 #Cantidad de iteraciones M  de maximo
Mitera4 = 4 
Nitera = 1 #Valor inicial de la iteraacion, es decir, 1, luego llegara a: itera, debe ser 1 siempre
n_bins = 100
#p=R-r #Radio interior de la corona circular

tries=1e5

n1, bin1 = probability(R, r, tries, n_bins, 1)

n2,  bin2 = probability(R, r, tries, n_bins, 2)

n3, bin3 = probability(R, r, tries, n_bins, 3)

n4, bin4 = probability(R, r, tries, n_bins, 4)


x = np.linspace(-R,R, 101)
x_ = np.linspace(-R,R, 30)


P1 = np.array([])
P2 = np.array([])
P3 = np.array([])
P4 = np.array([])


for i in x:
    P1 = np.append(P1, Iteracion(R, r, i, N, Mitera1, Nitera))
    P2 = np.append(P2, Iteracion(R, r, i, N, Mitera2, Nitera))


for j in x_:
    P3 = np.append(P3, Iteracion(R, r, j, N_, Mitera3, Nitera))
    P4 = np.append(P4, Iteracion(R, r, j, N__, Mitera4, Nitera))



#plt.ylim(0, 10)
plt.plot(x, P1, linewidth=2, color='darkblue', label='Calc. 1º cierre')#darkblue
plt.plot(x, gaussian_filter1d(P2, sigma=1), linewidth=2, color='orangered', label='Calc. 2º cierre')
plt.plot(x_, gaussian_filter1d(P3, sigma=1), linewidth=2, color='darkgreen', label='Calc. 3º cierre')
plt.plot(x_, gaussian_filter1d(P4, sigma=1), linewidth=2, color='red', label='Calc. 4º cierre')

plt.title(r'$R = 200$ u.a., $r = $'+str(r)+' u.a.:   $R=$'+str(np.round(k, 2))+'$r$')
plt.legend(loc='best')
plt.tight_layout()
#plt.savefig('imagen.png', dpi=500)

plt.show()


"""
Tarda 1 minuto en hacerse si no tocas nada
"""


















