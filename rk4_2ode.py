import numpy as np
import matplotlib.pyplot as plt

class RK4:
    def __init__(self,N,del_t,f1,f2):
        self.N=N # Number of iterations
        self.del_t=del_t # discrete time step
        self.f1=f1 # function f(x,x_dot,t)
        self.f2=f2 # function g(x,x_dot,t)

    def kl_values(self,x,v,t):
        """
        
        """
        k_0=self.del_t*self.f1(x,v,t)
        l_0=self.del_t*self.f2(x,v,t)
        k_1=self.del_t*self.f1(x+0.5*k_0,v+0.5*l_0,t+0.5*self.del_t)
        l_1=self.del_t*self.f2(x+0.5*k_0,v+0.5*l_0,t+0.5*self.del_t)
        k_2=self.del_t*self.f1(x+0.5*k_1,v+0.5*l_1,t+0.5*self.del_t)
        l_2=self.del_t*self.f2(x+0.5*k_1,v+0.5*l_1,t+0.5*self.del_t)
        k_3=self.del_t*self.f1(x+k_2,v+l_2,t+self.del_t)
        l_3=self.del_t*self.f2(x+k_2,v+l_2,t+self.del_t)

        return (k_0+2.*k_1+2.*k_2+k_3)/6 , (l_0+2.*l_1+2.*l_2+l_3)/6

    def solution(self,x0,v0):
        t=np.empty(self.N) # time space 
        x=np.empty(self.N) # Solution space for position,x
        v=np.empty(self.N) # dx/dt = v, first derivative of position is velocity 
        x[0]=x0
        v[0]=v0
        t[0]=0.

        for i in range(1,self.N):
            k,l=self.kl_values(x[i-1],v[i-1],t[i-1])
            t[i]=t[i-1]+self.del_t
            v[i]=v[i-1]+l 
            x[i]=x[i-1]+k 

        return x,v,t 


if __name__=='__main__':

    def g(x,x_dot,t): # g(x,x_dot,time) = -x or dv/dt = -x, The oscillator equation without any damping
        return -x   

    def f(x,x_dot,t): # dx/dt = v
        return x_dot

    n=500 # Number of iterations
    dt=0.1 # discrete time step

    ob=RK4(n,dt,f,g)
    x_init=0.0 # Initial position of oscillator
    v_init=8.0 # Initial velocity of oscillator
    x,x_dot,time=ob.solution(x_init,v_init)

    fig=plt.figure(figsize=(12,4))
    
    plot_one=fig.add_subplot(111)
    plot_one.set_xlim(-5,55)
    plot_one.set_ylim(-10,10)
    result_one, = plot_one.plot([],[],linestyle='-.',color='purple')
    result_one.set_data(time,x)
    plt.show()