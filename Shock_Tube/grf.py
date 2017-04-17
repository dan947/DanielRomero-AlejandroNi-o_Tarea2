import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')

datos=np.genfromtxt('data.dat')
XD=datos[:,0]
RHO=datos[:,1]
PRESION=datos[:,2]
VELOCIDAD=datos[:,3]

#Solucion analitica:
L = 1.0
g = 1.4
dt = 0.001

rho1 = 1.0
rho2 = 0.125
u1 = 0.0
u2 = 0.0
p1 = 1.0
p2 = 0.1

Al = 2./((g+1)*rho1)
Ar = 2./((g+1)*rho2)
Bl = ((g-1)/(g+1))*p1
Br = ((g-1)/(g+1))*p2
al = (g*p1/rho1)**0.5
ar = (g*p2/rho2)**0.5

def F(p):
    return F_shock(p) + F_rare(p) + u2 - u1

def F_shock(p):
    return (p-p2)*(Ar/(p+Br))**0.5

def F_rare(p):
    return (2.0*al/(g-1))*((p/p1)**(0.5*(g-1)/g)-1)

def US(p):
    return 0.5*(u1 + u2 + F_shock(p) - F_rare(p))

def Fprime(p):
    return (al/(g*p1))*((p1/p)**((g+1)/(2.*g))) + ((Ar/(p+Br))**0.5)*((p2-p)/(2.*(p+Br)) + 1.)

def NR(p):
    pa = p
    pn = pa - F(pa)/Fprime(pa)
    CHA = abs(pn - pa)/(0.5*(pn + pa))
    while(CHA < 10E-6):
        pa = pn
        pn = pa - F(pa)/Fprime(pa)
        CHA = abs(pn - pa)/(0.5*(pn + pa))

    return pn

def Solucion():

    p0 = 0.315

    Ps = NR(p0)
    Us = US(Ps)
    RhoSL = rho1*(Ps/p1)**(1./g)
    aSL = al*(Ps/p1)**((g-1)/(2.*g))
    SHL = u1 - al
    STL = Us - aSL
    RhoSR = rho2*(((Ps/p2)+((g-1.)/(g+1.)))/(((g-1.)/(g+1.))*(Ps/p2) + 1.))
    SR = u2 + ar*((((g+1)/(2.*g))*(Ps/p2)) + ((g-1.)/(2.*g)))**(0.5)

    X = np.linspace(0,1,100)
    P = np.zeros(100)
    Rho = np.zeros(100)
    U = np.zeros(100)

    t=0
    run = True
    
    while(run == True):
        t = t + dt
        for i in range(1,100):
            p = 0
            r = 0
            u = 0 
            x = X[i] - 0.5
            xt = x/t

            if(xt <= SHL):
                p = p1
                r = rho1
                u = u1
            elif(SHL < xt <= STL):
                p = p1*((2./(g+1.)) + (g-1.)/(al*(g+1))*(u1 - xt))**((2.*g)/(g-1.))
                r = rho1*((2./(g+1.)) + (g-1.)/(al*(g+1))*(u1 - xt))**(2./(g-1.))
                u = (2./(g+1.))*(al + ((g-1)*u1)/2. + xt)
            elif(STL < xt <= Us):
                p = Ps
                r = RhoSL
                u = Us
            elif(Us < xt <= SR):
                p = Ps
                r = RhoSR
                u = Us
            else:
                p = p2
                r = rho2
                u = u2
            
            P[i] = p
            Rho[i] = r
            U[i] = u

        if(abs(U[90]-U[89]) > 0):
            run = False
     
    return P,Rho,U,X


PA, RA, UA, XA = Solucion()

plt.subplot(311)
plt.scatter(XD,PRESION,color='yellow',label='Computacional')
plt.plot(XA,PA,color='black',label='Analitica')
plt.title("Presion")
plt.legend()

plt.subplot(312)
plt.scatter(XD,RHO,color='green',label='Computacional')
plt.plot(XA,RA,color='black',label='Analitica')
plt.title("Densidad")
plt.legend()

plt.subplot(313)
plt.scatter(XD,VELOCIDAD,color='pink',label='Computacional')
plt.plot(XA,UA,color='black',label='Analitica')
plt.title("Velocidad")
plt.legend()

plt.savefig("Graficas.pdf")
plt.show()
    





