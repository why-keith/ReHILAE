from pylab import *
from scipy.integrate import odeint

def dy_dx(y,x):
    return x - y

x = arange(0.,10.,0.001)
y = odeint(dy_dx,1.0,x) # the 2nd argument encodes y(0.) = 1.0
plot(x,y)
xlabel('x')
ylabel('y')
title('dy_dx = x - y')
show()
