import matplotlib.pyplot as plt
from rk4_2ode import RK4

def g(x,x_dot,t):
    return 2.*x_dot*(1-x*x) - x 

def f(x,x_dot,t):
    return x_dot

n=1000 # Number of iterations
dt=0.1 # discrete time step

ob=RK4(n,dt,f,g)
x_init=2.0 # Initital position of oscillator
v_init=0.0 # Initial velocity of oscillator

x,x_dot,time=ob.solution(x_init,v_init)

# Visualizing the results for damped oscillator

fig = plt.figure(figsize=(10,5))

plot_one=fig.add_subplot(111)
plot_one.set_xlim(-5,105)
plot_one.set_ylim(-10,10)

result_one, = plot_one.plot([],[],linestyle='-.',color='black')
result_one.set_data(time,x)

plt.show()