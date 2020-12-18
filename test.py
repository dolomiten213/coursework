
import matplotlib.pyplot as plt
import numpy as np

from math import cos, sin, sinh, log, exp, sqrt, pi

def printf(a):
    print("{0:.5f}".format(a))

class Point:
    x = 0
    z = 0

#============================================================================

sumCount = 50
N = 45
Z = 0

switcher = 2

if(switcher == 0):
   
    k = 6
    a = 2
    z1 = 0
    z2 = 4
    Z = np.linspace(z1, z2, N+1)
    
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
    plt.plot(x, y)
    x = np.linspace(2, 4, 2)
    y = np.zeros(2)
    plt.plot(x, y)
                      
elif(switcher == 1):
   
    k = 6
    a = 2
    z1 = 0
    z2 = 2
    Z = np.linspace(z1, z2, N+1)
    
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
    plt.plot(x, y)
    x = np.linspace(1, 2, 1001)
    y = 0.5*(2*((x-2)**3)+3*((x-2)**2))
    plt.plot(x, y)

elif(switcher == 2):

    k = 6
    a = 2
    z1 = 0
    z2 = 3
    Z = np.linspace(z1, z2, N+1)
    
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
    plt.plot(x, y)
    x = np.linspace(1, 3, 1001)
    y = 0.5*(0.25*((x-3)**3)+0.75*((x-3)**2))
    plt.plot(x, y)  

elif(switcher == 3):
    
    k = 6
    a = 2
    z1 = 0
    z2 = 4
    Z = np.linspace(z1, z2, N+1)

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
    plt.plot(x, y)

def J(z):
    return sqrt(1+(derivative(z))**2)

#============================================================================

def fi(a, m, x):
    return sqrt(2/a)*sin(pi*m*x/a)

def pm(m):
    x = ((pi*m/a)**2)-(k**2)
    if x > 0:
        return complex(sqrt(x), 0)
    else:
        return complex(0, sqrt(-x))

def gamma(m):
    b = k**2-(pi*m/a)**2
    if(b > 0):
        return sqrt(b)
    else:
        return sqrt(b) 

#============================================================================

def G(M, P): 
    sum = 0
    for m in range(1, sumCount + 1):        
        x = complex(0,0)
        b = ((pi*m/a)**2)-(k**2)
       
        if b < 0:
            x = x + (complex(2*pi/a, 0) / pm(m) * complex(cos(pm(m).imag*abs(M.z-P.z)), -sin(pm(m).imag*abs(M.z-P.z))))
            x = x + complex((-2/m)*exp((-pi*m/a)*abs(M.z-P.z)), 0)
            x = x * complex(sin((pi*m/a)*M.x), 0)
            x = x * complex(sin((pi*m/a)*P.x), 0)
        else:
            x = x + (2*pi)/(a*pm(m).real) * exp(-pm(m).real*abs(M.z-P.z))        
            x = x + (-2/m)*exp((-pi*m/a)*abs(M.z-P.z))
            x = x * sin((pi*m/a)*M.x)
            x = x * sin((pi*m/a)*P.x)
        
        #print(x)
        #print("\033[37m {:>4.5f}".format(abs(x)))
        sum += x
        #print(abs(sum))
        
    #print("==============") 
    return sum

def psea(M, P):
    
    dividend = 0
    dividend = (sin((pi/(2*a))*(M.x+P.x)))**2
    dividend = dividend + (sinh((pi/(2*a))*(M.z-P.z)))**2

    divider = 0
    divider = (sin((pi/(2*a))*(M.x-P.x)))**2
    divider = divider + (sinh((pi/(2*a))*(M.z-P.z)))**2

    #print("psea")
    #printf(dividend/divider)

    return dividend/divider

#def lnF(M,P):
#   return 0.5*log(psea(M, P))

#def g(M, P): 
#    return G(M, P) + lnF(M, P)

def K(i, j):
    
    zm = (Z[i] + Z[i+1])/2
    zp = (Z[j] + Z[j+1])/2

    L = z2 - z1
    M = Point()
    M.z = zm
    M.x = f(zm)

    P = Point()
    P.z = zp
    P.x = f(zp)

    sum = G(M, P)
    sum = sum * J(zp)
    sum = sum * (Z[j+1] - Z[j])
    
    if i == j:
         
        psi = 0
        psi = (Z[j+1]-zm)*(log(abs(Z[j+1] - zm))-1) - (Z[j]-zm)*(log(abs(Z[j] - zm))-1)
        psi = psi * J(zp)

        return sum - psi
    
    elif abs(zm - zp) <= L/2:
        res = 0
        res = psea(M, P) * ((abs(zp - zm))**2)
        res = 0.5*log(res)
        res = res * J(zp)
        res = res * (Z[j+1] - Z[j])

        psi = 0
        psi = (Z[j+1]-zm)*(log(abs(Z[j+1] - zm))-1) - (Z[j]-zm)*(log(abs(Z[j] - zm))-1)
        psi = psi * J(zp)

        return sum + res - psi

    elif abs(zm - zp) > L/2:

        res = 0
        res = psea(M, P) * ((L - abs(zp - zm))**2)
        res = 0.5*log(res)
        res = res * J(zp)
        res = res * (Z[j+1] - Z[j])



        psi = 0
        psi = (L - abs(Z[j+1]-zm))*(log(L - abs(Z[j+1] - zm)) - 1)
        psi = psi - (L - abs(Z[j]-zm))*(log(L - abs(Z[j] - zm)) - 1)
        psi = psi * J(zp)
        
        if Z[j] > zm:
            return sum + res + psi
        else:           
            return sum + res - psi

def R(i):
    zm = (Z[i] + Z[i+1])/2
    return -complex(fi(a, 1, f(zm)) * cos(gamma(1)*zm), fi(a, 1, f(zm)) * sin(gamma(1)*zm))                                          

#============================================================================
Kernels = np.zeros((N, N), dtype = np.complex)
Right = np.zeros(N , dtype = np.complex)
#KernelsRe = np.zeros((N, N), dtype = np.float)
#KernelsIm = np.zeros((N, N), dtype = np.float)

for i in range(N):
    s = ""
    for j in range(N):
        Kernels[i][j] = K(i,j)/(2*pi)
        temp = Kernels[i][j]
        re = temp.real
        im = temp.imag
        
        if (abs(temp) < 0.01):
            s += "\033[37m {:>4.1f}{:>4.1f}i".format(re, im)     
        elif (abs(temp) < 0.1):
            s += "\033[36m {:>4.1f}{:>4.1f}i".format(re, im) 
        else:
            s += "\033[34m {:>4.1f}{:>4.1f}i".format(re, im) 
    Right[i] = R(i)
    temp = Right[i]
    re = temp.real
    im = temp.imag
        
    if (abs(temp) < 0.01):
        s += "  \033[37m {:>4.1f}{:>4.1f}i".format(re, im)     
    elif (abs(temp) < 0.1):
        s += "  \033[36m {:>4.1f}{:>4.1f}i".format(re, im) 
    else:
        s += "  \033[34m {:>4.1f}{:>4.1f}i".format(re, im)     
       
    print(s)

#print(np.linalg.solve(Kernels, Right))
Ans = np.linalg.solve(Kernels, Right)

#for i in range(N):
#    l = abs(Ans[i])
#    x = np.linspace(Z[i], Z[i+1], 10)
#    y = np.full(10, l)
#    plt.plot(x, y) 
#plt.show()

#for i in range(N):
#    x = np.linspace(Z[i], Z[i+1], 10)
#    y = np.full(10, Ans[i].imag)
#    plt.plot(x, y)
#    y = np.full(10, Ans[i].real)
#    plt.plot(x, y)
#plt.show()
x = np.linspace((Z[0]+Z[1])/2, (Z[N-1]+Z[N])/2, N)
mylist = []

for i in range(N):
    mylist.append(abs(Ans[i]))

y = np.array(mylist)

plt.plot(x, y)
plt.show()



KernelsRe = np.zeros((N, N), dtype = np.float)
KernelsIm = np.zeros((N, N), dtype = np.float)


#M = Point()
#M.x = 0.1
#M.z = 0.15
#P = Point()
#P.x = 0.12
#P.z = 0.16

#print(G(M,P))

#for i in range(N):
    #s = ""
    #for j in range(N):
        
        #KernelsRe[i][j] = K(i, j).real
        #KernelsIm[i][j] = K(i, j).imag
        #if (sqrt(KernelsRe[i][j]**2 + KernelsIm[i][j]**2) < 1):
        #    s += "\033[37m {:>4.1f}".format(KernelsRe[i][j])
        #    s += " +i"
        #    s += "\033[37m {:>4.1f}".format(KernelsIm[i][j])
        #else:
        #    s += "\033[36m {:>4.1f}".format(KernelsRe[i][j])
        #    s += " +i"
        #    s += "\033[36m {:>4.1f}".format(KernelsIm[i][j])
        
        
        #elif (abs(Kernels[i][j]) < 0.1):
        #   s += "\033[36m {:>4.1f}".format(Kernels[i][j])
        #else:
        #   s += "\033[34m {:>4.1f}".format(Kernels[i][j])

    #print(s)

#plt.show()

#a = 3+4j                     
#print(abs(a))