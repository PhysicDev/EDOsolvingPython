import numpy as np
#import seaborn as sb
import matplotlib.pyplot as plt

# a,b,c :les coeficients de la matrice de runge kutta
# x0,t0 : les conditions initiales
# f : fonction dependant de x et de t (equation differentielle)
# h : pas de temps
def RungeKutta(a,b,c,x0,t0,f,h):
  p=np.array([x0 for i in range(len(a))])
  x=x0
  t=t0
  for i in range(len(a)):
    t=t0+h*c[i]
    x=x0+h*sum(np.array([a[i,j]*p[j] for j in range(len(a))]))
    p[i]=f(x,t)
  return(x0+h*sum(np.array([b[j]*p[j] for j in range(len(a))])),t0+h);


# a,b,c :les coeficients de la matrice de runge kutta
# x0,t0 : les conditions initiales
# f : fonction dependant de x et de t (equation differentielle)
# h : pas de temps
# param : nombre d'iteration/ temps de simulation
# cond : vaut false si on veut simuler sur n iteraion et True si on veut simuler sur un certain temps
def RungeKuttaSim(a,b,c,x0,t0,f,h,param=100,cond=False):
  n=0
  if(cond):
    n=(int)(param/h)
  else:
    n=param
  t=[0]*n
  x=[0]*n
  x[0]=x0
  t[0]=t0
  for i in range(1,n):
    if(i%500==0):
      print("progress :",i/n*100,"%")
    (x[i],t[i])=RungeKutta(a,b,c,x[i-1],t[i-1],f,h)
  return (x,t)

# x0,t0 : les conditions initiales
# f : fonction dependant de x et de t (equation differentielle)
# derf : la jacobienne de f
# h : pas de temps
# solveIter : le nombre d'itération de la méthode de Newton à effectuer
def EulerImplicite(x0,t0,f,derf,h,solveIter=5):
  #schema euler implicite
  x=x0
  t=t0+h
  for i in range(solveIter):
    x=x-np.dot(np.linalg.inv(h*derf(x,t)-np.eye(len(x))),x0+h*f(x,t)-x)
  return(x,t)
  
# x0,t0 : les conditions initiales
# f : fonction dependant de x et de t (equation differentielle)
# derf : la jacobienne de f
# h : pas de temps
# solveIter : le nombre d'itération de la méthode de Newton à effectuer
# param : nombre d'iteration/ temps de simulation
# cond : vaut false si on veut simuler sur n iteraion et True si on veut simuler sur un certain temps
def EulerImpliciteSim(x0,t0,f,derf,h,solveIter,param=100,cond=False):
  n=0
  if(cond):
    n=(int)(param/h)
  else:
    n=param
  t=[0]*n
  x=[0]*n
  x[0]=x0
  t[0]=t0
  for i in range(1,n):
    if(i%500==0):
      print("progress :",i/n*100,"%")
    (x[i],t[i])=EulerImplicite(x[i-1],t[i-1],f,derf,h,solveIter)
  return(x,t)

#RK-4 coefficient
a=np.array([[0,0,0,0],
            [1/2,0,0,0],
            [0,1/2,0,0],
            [0,0,1,0]])
b=np.array([1/6,2/6,2/6,1/6])
c=np.array([0,1/2,1/2,1])

#sys1
#"""
alpha=1
beta=1/20
delta=1/200
gamma=1/2


"""
#sys2
alpha=3
beta=1/40
delta=2
gamma=1/10
epsilon=1/20
zeta=1/40
eta=1/20
theta=1
iota=1/20
kappa=1/20
#"""
#custom
hunt=0.05
C=0.5
debut=10
fin=30
A=0.25



def sys1(x,t):
  return np.array([x[0]*(alpha-beta*x[1])
                   ,x[1]*(delta*x[0]-gamma)])

def sys2(x,t):
  return np.array([x[0]*(alpha-beta*x[0]-gamma*x[1]),
                   x[1]*(delta-epsilon*x[1]-zeta*x[0]-eta*x[2]),
                   x[2]*(theta*x[1]-iota*x[2]-kappa)])

def derSys1(x,t):
  return np.array([[(alpha-beta*x[1]),-x[0]*beta],
                  [x[1]*delta,(delta*x[0]-gamma)]])

def derSys2(x,t):
  return np.array([[(alpha-2*beta*x[0]-gamma*x[1]),-x[0]*gamma,0],
                   [-x[1]*zeta,(delta-2*epsilon*x[1]-zeta*x[0]-eta*x[2]),-x[1]*eta],
                   [0,-x[2]*theta,(theta*x[1]-2*iota*x[2]-kappa)]])

def sysCustom(x,t):
  return np.array([x[0]*(alpha-beta*x[1])
                   ,x[1]*(delta*x[0]-gamma-t*hunt)])

def guerreMondiale(x,t):
    return np.array([x[0]*(alpha-beta*x[1]-chasse(t))
                    ,x[1]*(delta*x[0]-gamma)])
def chasse(t):
    out=C
    if(t>debut and t<fin):
        out=out-A*(np.sin(np.pi*((t-debut)/(fin-debut)))**(2))
    return out

#test
def exp(x,t):
  return x



#(x,t)=RungeKuttaSim(a,b,c,np.array([1,2,2]),0,sys2,0.05,20,True)


"""
(x,t)=EulerImpliciteSim(np.array([80,30]),0,sys1,derSys1,0.01,5,50,True)
plt.plot(t,[i[0] for i in x])
plt.plot(t,[i[1] for i in x])
#plt.plot(t,[i[2] for i in x])
plt.title("Système 2 Euler Implicite")
#"""#"""
(x1,t1)=RungeKuttaSim(a,b,c,np.array([100,10]),0,guerreMondiale,0.001,60,True)

plt.plot(t1,[i[0] for i in x1])
plt.plot(t1,[i[1] for i in x1])
plt.plot(t1,[90*chasse(i) for i in t1])
plt.plot(t1,[8*i[0]/i[1] for i in x1])
#plt.plot(t1,[i[2] for i in x1])
#"""
plt.title("chasse sur les proies")
plt.legend(["proie","predateur","chasse (pas à l'échelle)","rapport proie/prédateur (pas à l'échelle)"])
"""
(x1,t1)=RungeKuttaSim(a,b,c,np.array([100,20]),0,sysCustom,0.001,60,True)
plt.xlabel('Prédateur')
plt.ylabel('Proie')
plt.plot([i[1] for i in x1],[i[0] for i in x1])
#"""
#plt.title("Système 1 Euler Implicite")
"""
for j in range(20):
    (x1,t1)=RungeKuttaSim(a,b,c,np.array([100,20+j]),0,sys1,0.01,50,True)
    plt.plot([i[1] for i in x1],[i[0] for i in x1],color="#5999e0")
    plt.text(x1[1][1], x1[1][0], str(20+j))
    plt.title("Système 1")


#plt.legend(["Euler Implicite","RungeKutta"])
#"""
#plt.plot(t,[i[2] for i in x])
plt.grid(which='major', linestyle='-')
plt.grid(which='minor', linestyle='--',linewidth=0.3)
plt.minorticks_on()

plt.gca().set_ylim(bottom=0)

plt.show()

