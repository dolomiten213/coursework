
import matplotlib.pyplot as plt
import numpy as np

from math import cos, sin, sinh, log, exp, sqrt, pi

def printf(a):
    print("{0:.5f}".format(a))

class Point:
    x = 0
    z = 0

#============================================================================

k = 6
a = 2
z1 = 0
z2 = 4
N = 70

Z = np.linspace(z1, z2, N+1)

sumCount = 20

#============================================================================

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

def J(z):
    return sqrt(1+(derivative(z))**2)

#============================================================================

def fi(a, m, x):
    return sqrt(2/a)*sin(pi*m*x/a)

def pm(m):
    x = (k**2)-((pi*m/a)**2)
    if x > 0:
        return sqrt(x)
    else:
        return sqrt(-x)

#============================================================================

def G(M, P): 
    sum = 0
    for m in range(1, sumCount + 1):
        x = 0

        x = x + (2*pi)/(a*pm(m))*exp(-pm(m)*abs(M.z-P.z))
        x = x + (-2/m)*exp((-pi*m/a)*abs(M.z-P.z))
        x = x * sin((pi*m/a)*M.x)
        x = x * sin((pi*m/a)*P.x)

        sum += x
        
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

#============================================================================

Kernels = np.zeros((N, N), dtype = np.float)

for i in range(N):
    s = ""
    for j in range(N):
        
        Kernels[i][j] = K(i, j)
        if (Kernels[i][j] < 0.01):
            s += "\033[37m {:>4.1f}".format(Kernels[i][j])
        elif (Kernels[i][j] < 0.1):
            s += "\033[36m {:>4.1f}".format(Kernels[i][j])
        else:
            s += "\033[34m {:>4.1f}".format(Kernels[i][j])

    print(s)


#x = np.linspace(0, 1, 1000)    # Create a list of evenly-spaced numbers over the range
#plt.plot(Z, f(Z))              # Plot the sine of each x point
#plt.show()                     # Display the plot