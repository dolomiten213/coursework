
import matplotlib.pyplot as plt
import numpy as np

from math import cos, sin, sinh, log, exp, sqrt, pi

#from numpy.lib.shape_base import apply_along_axis

def printf(a):
    print("{0:.5f}".format(a))

class Point:
    x = 0
    z = 0

#============================================================================

Z = 0
k = 6
a = 2

switcher = 4

if(switcher == 0):
   
    z1 = 0
    z2 = 2
    
    def f(z):
        if z > 0 and z <= 2:
            return 0.5*(sin(pi*(z-0.5)) + 1)
        else:
            return 0

    def derivative(z):
        if z > 0 and z <= 2:
            return 0.5*pi*cos(pi*(z-0.5))
        else:
            return 0
    
    x = np.linspace(0, 2, 1001)
    y = 0.5*(np.sin(pi*(x-0.5)) + 1)
    plt.plot(x, y, color = '#4a28ea')
                      
elif(switcher == 1):

    z1 = 0
    z2 = 2
    aa = 0.64*pi 
    bb = -0.16
    
    def f(z):
        if z > 0 and z <= 1:
            return 0.5*(-2*(z**3)+3*(z**2))
        elif z > 1 and z <= 2:
            return 0.5*(2*((z-2)**3)+3*((z-2)**2))
        else:
            return 0

    def derivative(z):
        if z > 0 and z <= 1:
            return 0.5*(-6*(z**2)+6*z)
        elif z > 1 and z <= 2:
            return 0.5*(6*((z-2)**2)+6*(z-2))
        else:
            return 0

    x = np.linspace(0, 1, 1001)
    y = 0.5*(-2*(x**3)+3*(x**2))
    plt.plot(x, y, color = '#247719')
    x = np.linspace(1, 2, 1001)
    y = 0.5*(2*((x-2)**3)+3*((x-2)**2))
    plt.plot(x, y, color = '#247719')

elif(switcher == 2):

    z1 = 0
    z2 = 3
    
    aa = 0.56*pi 
    bb = -0.06

    def f(z):
        if z > 0 and z <= 1:
            return 0.5*(-2*(z**3)+3*(z**2))
        elif z > 1 and z <= 3:
            return 0.5*(0.25*((z-3)**3)+0.75*((z-3)**2))
        else:
            return 0

    def derivative(z):
        if z > 0 and z <= 1:
            return 0.5*(-6*(z**2)+6*z)
        elif z > 1 and z <= 3:
            return 0.5*(0.75*((z-2)**2)+1.5*(z-2))
        else:
            return 0

    x = np.linspace(0, 1, 1001)
    y = 0.5*(-2*(x**3)+3*(x**2))
    plt.plot(x, y, color = '#247719')
    x = np.linspace(1, 3, 1001)
    y = 0.5*(0.25*((x-3)**3)+0.75*((x-3)**2))
    plt.plot(x, y, color = '#247719') 

elif(switcher == 3):
    z1 = 0
    z2 = 4

    aa = 0.96*pi 
    bb = -0.04

    def f(z):
        if z > 0 and z <= 4:
            return 0.5*(sin(pi*(z-0.5)) + 1)
        else:
            return 0

    def derivative(z):
        if z > 0 and z <= 4:
            return 0.5*pi*cos(pi*(z-0.5))
        else:
            return 0

    x = np.linspace(0, 4, 1001)
    y = 0.5*(np.sin(pi*(x-0.5)) + 1)
    plt.plot(x, y, color = '#247719')

elif(switcher == 4):
    z1 = 0
    z2 = 6

    aa = 0.96*pi 
    bb = 0.04

    def f(z):
        if z > 0 and z <= 2:
            return 0.5*(sin(pi*(z-0.5)) + 1)
        if z > 2 and z <= 4:
            return 0
        if z > 4 and z <= 6:
            return 0.5*(sin(pi*(z-0.5)) + 1)
        else:
            return 0

    def derivative(z):
        if z > 0 and z <= 2:
            return 0.5*pi*cos(pi*(z-0.5))
        if z > 4 and z <= 6:
            return 0.5*pi*cos(pi*(z-0.5))
        else:
            return 0

    x = np.linspace(0, 2, 1001)
    y = 0.5*(np.sin(pi*(x-0.5)) + 1)
    plt.plot(x, y, color = '#247719')
    
    x = np.linspace(2, 4, 3)
    y = np.zeros(3)
    plt.plot(x, y, color = '#247719')
    
    x = np.linspace(4, 6, 1001)
    y = 0.5*(np.sin(pi*(x-0.5)) + 1)
    plt.plot(x, y, color = '#247719')

x = np.linspace(-2, 10, 3)
y = np.full(3, a)
plt.plot(x, y, color = 'black')

x = np.linspace(-2, 10, 3)
y = np.full(3, 0)
plt.plot(x, y, color = 'black')

plt.show()

def J(z):
    return sqrt(1+(derivative(z))**2)

#============================================================================

def fi(m, x):
    return sqrt(2/a)*sin(pi*m*x/a)

def pm(m):
    x = ((pi*m/a)**2)-(k**2)
    if x > 0:
        return complex(sqrt(x), 0)
    else:
        return complex(0, sqrt(-x))

def gamma(m):
    x = k**2-(pi*m/a)**2
    if x > 0:
        return complex(sqrt(x), 0)
    else:
        return complex(0, sqrt(-x))

#============================================================================

def G(M, P): 
    sum = 0
    for m in range(1, 51):        
        x = complex(0,0)
        b = ((pi*m/a)**2)-(k**2)
        if b < 0:
            x = complex(2*pi/a, 0) / pm(m) 
            x = x * complex(cos(pm(m).imag*abs(M.z-P.z)), -sin(pm(m).imag*abs(M.z-P.z)))
            x = x + complex((-2/m)*exp((-pi*m/a)*abs(M.z-P.z)), 0)
            x = x * sin((pi*m/a)*(M.x))
            x = x * sin((pi*m/a)*(P.x))

        else:
            x = (2*pi)/(a*pm(m).real) * exp(-pm(m).real*abs(M.z-P.z))        
            x = x + (-2/m)*exp((-pi*m/a)*abs(M.z-P.z))
            x = x * sin((pi*m/a)*(M.x))
            x = x * sin((pi*m/a)*(P.x))
        sum += x
        
    #sum = sum * sin((pi*m/a)*(M.x))
    #sum = sum * sin((pi*m/a)*(P.x))
    return sum

def psea(M, P):
    
    dividend = 0
    dividend = (sin((pi/(2*a))*(M.x+P.x)))**2
    dividend = dividend + (sinh((pi/(2*a))*(M.z-P.z)))**2

    divider = 0
    divider = (sin((pi/(2*a))*(M.x-P.x)))**2
    divider = divider + (sinh((pi/(2*a))*(M.z-P.z)))**2

    return dividend/divider

def K(i, j, Z):
    
    zm = (Z[i] + Z[i+1])/2
    zp = (Z[j] + Z[j+1])/2

    L = z2 - z1
    M = Point()
    M.z = zm
    M.x = f(zm)

    P = Point()
    P.z = zp
    P.x = f(zp)

    sum = 0
    sum = G(M,P)
    sum = sum * J(zp)
    sum = sum * (Z[j+1] - Z[j])
    #print(sum)

    if i == j:
         
        psi = 0
        psi = (Z[j+1]-zm)*(log(abs(Z[j+1] - zm))-1) - (Z[j]-zm)*(log(abs(Z[j] - zm))-1)
        psi = psi * J(zp)
        #print(psi)
        #print (sum-psi)
        return sum - psi
    
    else:
        res = 0
        res = psea(M, P) * ((abs(zp - zm))**2)
        res = 0.5*log(res)
        res = res * J(zp)
        res = res * (Z[j+1] - Z[j])

        psi = 0
        psi = (Z[j+1]-zm)*(log(abs(Z[j+1] - zm))-1) - (Z[j]-zm)*(log(abs(Z[j] - zm))-1)
        psi = psi * J(zp)

        return sum + res - psi

def R(i, Z):
    zm = (Z[i] + Z[i+1])/2
    return -(fi(1, f(zm)) * complex(cos(gamma(1).real*zm), sin(gamma(1).real*zm)))                                          

#============================================================================

def fia(m, x, a):
    return sqrt(2/a)*sin(pi*m*x/a)

def pma(m, a):
    x = ((pi*m/a)**2)-((pi+0.01)**2)
    if x > 0:
        return complex(sqrt(x), 0)
    else:
        return complex(0, sqrt(-x))

def gammaa(m, a):
    x = (pi+0.01)**2-(pi*m/a)**2
    if x > 0:
        return complex(sqrt(x), 0)
    else:
        return complex(0, sqrt(-x))

#def Rn(Ans, N, n, a):
    Z = np.linspace(z1, z2, N+1)
    sum = 0
    for i in range(N):
        zp = (Z[i] + Z[i+1])/2
        b = 0
        b = Ans[i] * fia(n, zp, a) * complex(cos(gammaa(n, a).real*zp), sin(gammaa(n, a).real*zp)) * J(zp)
        b = b * (Z[i+1]-Z[i])
        sum += b
    sum = sum / (2*pma(n, a))
    print("R",n," = ",abs(sum)**2)
    return sum

#def Tn(Ans, N, n, a):
    Z = np.linspace(z1, z2, N+1)
    sum = 0
    for i in range(N):
        zp = (Z[i] + Z[i+1])/2
        b = Ans[i] * fia(n, zp, a) * complex(cos(gammaa(n, a).real*zp), -sin(gammaa(n, a).real*zp)) * J(zp)
        b = b * (Z[i+1]-Z[i])
        sum += b
    sum = sum / (2*pma(n, a))
    if (n == 1): sum = sum + 1
    print("T",n," = ",abs(sum)**2)
    return sum

def rn(x, a, b):
    return a/x + b

#============================================================================

def buildGraph(N):
    Z = np.linspace(z1, z2, N+1)
    #print(Z)
    Kernels = np.zeros((N, N), dtype = np.complex)
    Right = np.zeros(N , dtype = np.complex)

    for i in range(N):
        s = ""
        for j in range(N):
            Kernels[i][j] = K(i,j,Z)/(2*pi)
            temp = Kernels[i][j]

            re = temp.real
            im = temp.imag
            
            if (abs(temp) < 0.01):
                s += "\033[37m {:>4.1f}{:>4.1f}i".format(re, im)     
            elif (abs(temp) < 0.1):
                s += "\033[36m {:>4.1f}{:>4.1f}i".format(re, im) 
            else:
                s += "\033[34m {:>4.1f}{:>4.1f}i".format(re, im)

        Right[i] = R(i, Z)
        temp = Right[i]
        re = temp.real
        im = temp.imag
            
        if (abs(temp) < 0.01):
            s += "  \033[37m {:>4.1f}{:>4.1f}i".format(re, im)     
        elif (abs(temp) < 0.1):
            s += "  \033[36m {:>4.1f}{:>4.1f}i".format(re, im) 
        else:
            s += "  \033[34m {:>4.1f}{:>4.1f}i".format(re, im)     
    
        #print(s)
        
    Ans = np.linalg.solve(Kernels, Right)

    x = np.linspace((Z[0]+Z[1])/2, (Z[N-1]+Z[N])/2, N)
    mylist = []

    for i in range(N):
        mylist.append(abs(Ans[i]))

    y = np.array(mylist)

    plt.plot(x, y)

    print(max(y))
    
    #R3 = gamma(1).real*abs(Rn(Ans,N,1))**2+gamma(2).real*abs(Rn(Ans,N,2))**2+gamma(3).real*abs(Rn(Ans,N,3))**2
    #T3 = gamma(1).real*abs(Tn(Ans,N,1))**2+gamma(2).real*abs(Tn(Ans,N,2))**2+gamma(3).real*abs(Tn(Ans,N,3))**2

    #print(R3+T3)
    #print(gamma(1))
    #print(gamma(2))
    #print(gamma(3))

def buildRTGraph():
    N = 50
    Z = np.linspace(z1, z2, N+1)
    Kernels = np.zeros((N, N), dtype = np.complex)
    Right = np.zeros(N , dtype = np.complex)

    x = np.linspace(3.14, 6.28, 20)
    mylist1 = []
    mylist2 = []

    for c in range (20):
        #for i in range(N):
            #s = ""
            #for j in range(N):
                #Kernels[i][j] = K(i,j,Z)/(2*pi)
                #temp = Kernels[i][j]
                #re = temp.real
                #im = temp.imag
                
                #if (abs(temp) < 0.01):
                #    s += "\033[37m {:>4.1f}{:>4.1f}i".format(re, im)     
                #elif (abs(temp) < 0.1):
                #    s += "\033[36m {:>4.1f}{:>4.1f}i".format(re, im) 
                #else:
                #    s += "\033[34m {:>4.1f}{:>4.1f}i".format(re, im) 
            #Right[i] = R(i, Z)
            #temp = Right[i]
            #re = temp.real
            #im = temp.imag
                
            #if (abs(temp) < 0.01):
            #    s += "  \033[37m {:>4.1f}{:>4.1f}i".format(re, im)     
            #elif (abs(temp) < 0.1):
            #    s += "  \033[36m {:>4.1f}{:>4.1f}i".format(re, im) 
            #else:
            #    s += "  \033[34m {:>4.1f}{:>4.1f}i".format(re, im)          
        #Ans = np.linalg.solve(Kernels, Right)

        r1 = rn(pi+c*(3.14/20), aa, bb)
        mylist1.append(r1)
        #r1 = abs(Rn(Ans,N,1,1+c*(3.14/20)))
        mylist2.append(1-(r1**2))
        #mylist1.append(r1)

    y = np.array(mylist1)
    plt.plot(x, y)
    y = np.array(mylist2)
    plt.plot(x, y)
    
    #mylist1 = []
    #mylist2 = []

    #for i in range(50):
    #    print(abs(Rn(Ans,N,2,a+i*(6.28/50),k)), abs(Tn(Ans,N,2,a+i*(6.28/50),k)))
    #    mylist1.append(abs(Rn(Ans,N,2,a+i*(6.28/50),k)))
    #    mylist2.append(abs(Tn(Ans,N,2,a+i*(6.28/50),k)))

    #y = np.array(mylist1)
    #plt.plot(x, y)
    #y = np.array(mylist2)
    #plt.plot(x, y)

    #mylist1 = []
    #mylist2 = []

    #for i in range(50):
    #    print(abs(Rn(Ans,N,1,a,k+i*(3.14/50))), abs(Tn(Ans,N,1,a,k+i*(3.14/50))))
    #    mylist1.append(abs(Rn(Ans,N,1,a,k+i*(3.14/50))))
    #    mylist2.append(abs(Tn(Ans,N,1,a,k+i*(3.14/50))))

    #y = np.array(mylist1)
    #plt.plot(x, y)
    #y = np.array(mylist2)
    #plt.plot(x, y)

    #mylist1 = []
    #mylist2 = []

    #for i in range(50):
    #    print(abs(Rn(Ans,N,1,a+i*(3.14/100),k+i*(3.14/100))), abs(Tn(Ans,N,1,a+i*(3.14/100),k+i*(3.14/100))))
    #    mylist1.append(abs(Rn(Ans,N,1,a+i*(3.14/100),k+i*(3.14/100))))
    #    mylist2.append(abs(Tn(Ans,N,1,a+i*(3.14/100),k+i*(3.14/100))))

    #y = np.array(mylist1)
    #plt.plot(x, y)
    #y = np.array(mylist2)
    #plt.plot(x, y)

    #R3 = gamma(1).real*abs(Rn(y,N,1,a))**2+gamma(2).real*abs(Rn(y,N,2,a))**2+gamma(3).real*abs(Rn(y,N,3,a))**2
    #T3 = gamma(1).real*abs(Tn(y,N,1,a))**2+gamma(2).real*abs(Tn(y,N,2,a))**2+gamma(3).real*abs(Tn(y,N,3,a))**2

    #print(R3+T3)
    #print(gamma(1))

def buildGGraph(N):
    Z = np.linspace(z1, z2, N+1)
    #print(Z)
    Kernels = np.zeros((N, N), dtype = np.complex)
    Right = np.zeros(N , dtype = np.complex)

    for i in range(N):
        s = ""
        sum = 0
        for j in range(N):

            zm = (Z[i] + Z[i+1])/2
            zp = (Z[j] + Z[j+1])/2

            L = z2 - z1
            M = Point()
            M.z = zm
            M.x = f(zm)

            P = Point()
            P.z = zp
            P.x = f(zp)

            temp = G(M,P)

            re = temp.real
            im = temp.imag
            
            if (abs(temp) < 0.01):
                s += "\033[37m {:>4.1f}{:>4.1f}i".format(re, im)     
            elif (abs(temp) < 0.1):
                s += "\033[36m {:>4.1f}{:>4.1f}i".format(re, im) 
            else:
                s += "\033[34m {:>4.1f}{:>4.1f}i".format(re, im)
            sum += temp

        temp = sum
        re = temp.real
        im = temp.imag
            
        if (abs(temp) < 0.01):
            s += "  \033[37m {:>4.1f}{:>4.1f}i".format(re, im)     
        elif (abs(temp) < 0.1):
            s += "  \033[36m {:>4.1f}{:>4.1f}i".format(re, im) 
        else:
            s += "  \033[34m {:>4.1f}{:>4.1f}i".format(re, im)

        #print(s)     
    
buildGraph(100)
buildGraph(200)

plt.show()

#buildRTGraph()

#plt.show()

#buildRTGraph(50)

#plt.show()


