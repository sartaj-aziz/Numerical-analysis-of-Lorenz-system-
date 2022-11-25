import numpy as np
from matplotlib import pyplot as plt

rhoe= 28
miu= 10
beta= 3
x0= -1
y0= -30
z0= 10
X0=[x0,y0,z0]

Nt=1000
tmax=10
t=np.linspace(0,tmax,Nt)

def derivative(X,t,rhoe,miu,beta):
    x,y,z = X
    dotx = miu*(y-x)
    doty = (x*(rhoe-z))-y
    dotz = x*y-beta*z
    return np.array([dotx,doty,dotz])
def euler(func,X0,t,rhoe,miu,beta):
    dt=t[1]-t[0]
    nt=len(t)
    X=np.zeros([nt,len(X0)])
    X[0]=X0
    for i in range(nt-1):
        X[i+1]=X[i]+func(X[i],t[i],rhoe,miu,beta)*dt
    return  X
Xe=euler(derivative,X0,t,rhoe,miu,beta)


plt.figure()
plt.title("Lorentz equation by euler method")
plt.plot(t,Xe[:,0],'xb',label='rate of convection')
plt.plot(t,Xe[:,1],'+r',label='horizontal temperature variation')
plt.plot(t,Xe[:,2],'r-',label='vertical temperature variation')
plt.grid()
plt.show()

ax=plt.axes(projection='3d')
x_data=Xe[:,0]
y_data=Xe[:,1]
z_data=Xe[:,2]
ax.plot(x_data,y_data,z_data)
plt.show()

