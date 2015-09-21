from scitools.std import *
from os import system
import sympy as sym
import matplotlib.pyplot as mpl

q,w,u,f,l,sx,st = sym.symbols('q w u f l sx st')

def u_exact(x,t):	# Exact solution 
	return cos(pi*x/L)*cos(t)

def qfunction(x):
	return 1 + (x-L/2)**4

def V(x):	# Initial velocity
	return 0
	
def I(x):	# Initial condition
	return u_exact(x,0)

def WriteToFile(x):
	data = open("data.txt", "a")
	for value in x:
		data.write("{0:15g}".format(float(value)))
	data.write("\n")
	data.close()
	return

def Error(u,x,t,dt):
    E = sqrt(dt*sum((u_exact(x,t)-u)**2))
    return E


def solver(f, Nx=50, T=1, Vis=False):
    Beta=0.9
    x = linspace(0,L,Nx+1)
    dx = x[1]
    q = qfunction(x)
    dt = float(dx*Beta/sqrt(max(q)))
    Nt = int(round(T/dt))
    t = linspace(0, Nt*dt, Nt+1)

    C2 = (dt/dx)**2
    dt2 = dt**2
    u = zeros(Nx+1); u_1 = zeros(Nx+1); u_2 = zeros(Nx+1)
    
    f = f.subs(l,L)
    f = sym.lambdify((sx,st),f)   # Changes f to a lambda function dependent of x and t.
    
    Ix = range(0, Nx+1)
    It = range(0, Nt+1)
    
    
    
    # First initial condition
    u_1 = I(x)

    for i in Ix[1:-1]:
        u[i] = u_1[i] + dt*V(x[i]) + \
        0.5*C2*(0.5*(q[i] + q[i+1])*(u_1[i+1] - u_1[i]) - \
                0.5*(q[i] + q[i-1])*(u_1[i] - u_1[i-1])) + \
        0.5*dt2*f(x[i], t[0])

    i = Ix[0]

    # Set boundary values (x=0: i-1 -> i+1 since u[i-1]=u[i+1]
    # when du/dn = 0, on x=L: i+1 -> i-1 since u[i+1]=u[i-1])
    ip1 = i+1
    u[i] = u_1[i] + dt*V(x[i]) + \
    0.5*C2*(q[i] + q[ip1])*(u_1[ip1] - u_1[i]) + 0.5*dt2*f(x[i], t[0])

    i = Ix[-1]
    im1 = i-1
    u[i] = u_1[i] + dt*V(x[i]) + \
    0.5*C2*(q[i] + q[im1])*(u_1[im1] - u_1[i]) + 0.5*dt2*f(x[i], t[0])
    
    WriteToFile(u)
    u_2, u_1, u = u_1, u, u_2
    
    for n in It[1:-1]:
        
        for i in Ix[1:-1]:
                u[i] = - u_2[i] + 2*u_1[i] + \
                    C2*(0.5*(q[i] + q[i+1])*(u_1[i+1] - u_1[i])  - \
                        0.5*(q[i] + q[i-1])*(u_1[i] - u_1[i-1])) + \
                dt2*f(x[i], t[n])

        # Setting boundary values using equation (57)
        # x=0: i-1 -> i+1 since u[i-1]=u[i+1] when du/dn=0
        i = Ix[0]
        ip1 = i+1
        u[i] = 2*u_1[i] - u_2[i] + \
        C2*(q[i] + q[ip1])*(u_1[ip1] - u_1[i]) + dt2*f(x[i], t[n])
        
        # x=L: i+1 -> i-1 since u[i+1]=u[i-1] when du/dn=0
        i = Ix[-1]
        im1 = i-1
        u[i] = - u_2[i] + 2*u_1[i] + \
        C2*(q[i] + q[im1])*(u_1[i-1] - u_1[i]) + dt2*f(x[i], t[n])
    
        WriteToFile(u)
        u_2, u_1, u = u_1, u, u_2
    
    u = []
    datau = open("data.txt", "r")
    for line in datau:
         #if counter%skip_frames == 2: # plots every 10th frame
         u.append(map(float, line.split()))   # Converts a list of strings to list of floats.
    datau.close()
    system("rm data.txt")

    if Vis:
        visualize(u,x,t)
    
    return u, x, t, dt



def visualize(u,x,t):
     skip_frames=1
     filename = "banimation"

     n = 0
     
     # Make a first plot
     mpl.ion()
     y1 = u_exact(x,t[n])
     y2 = u[n]
     lines = mpl.plot(x, y1, "r-", x, y2, "b-")
     mpl.axis([0, 1, -1, 1])
     mpl.xlabel('x')
     mpl.ylabel('u')
     mpl.legend(["Exact","Numerical"])
     mpl.grid("on")
     mpl.draw()
     mpl.savefig('tmp_%05d.png' % n)
     mpl.title('time: {0:.3f}'.format(t[n]))
     # Show the movie, and make hardcopies of frames simultaneously

     for n in range(1,len(u)):
        y1 = u_exact(x,t[n])
        y2 = u[n]
        for line, y in zip(lines, [y1, y2]):
            line.set_ydata(y)
        mpl.title('time: {0:.3f}'.format(t[n]))
        mpl.draw()
        mpl.savefig('tmp_%04d.png' % n)
     #movie('tmp_*.png', encoder='convert', fps=25, output_file=filename+'.gif')
     system("rm tmp_*.png")
     return
    
    
    
if __name__ == "__main__":

    q = 1+(sx-l/2)**4
    #q = 1+sym.cos(pi*sx/l)    
    u = sym.cos(sym.pi*sx/l)*sym.cos(st)
    f = sym.diff(u, st, st) - sym.diff(q*sym.diff(u,sx),sx)

    L=1.0
    
    solver(f,Nx=50,T=4*pi, Vis=True)

    
    
    m = 6
    Nx = 5
    E_list = zeros(m)
    dt_list = zeros(m)
    Nx_list=zeros(m)
    for i in range(m):
        print i
        u,x,t,dt = solver(f,Nx,T=10*pi)
        E_list[i] = Error(u,x,t[-1],dt) 
        dt_list[i] = dt
        Nx_list[i] = Nx
        Nx = Nx*2
    
    
    #r = log(E_list[:-1]/E_list[1:]) / log(dt_list[:-1]/dt_list[1:])
    r = log(E_list[:-1]/E_list[1:]) / log(Nx_list[:-1]/Nx_list[1:])
    figure()
    plot(log(dt_list), log(E_list))
    xlabel("log(dt)")
    ylabel("log(Error)")    
    figure()
    plot(log(Nx_list), log(E_list))
    xlabel("log(Nx)")
    ylabel("log(Error)")
    print "Convergence rate: ", r
    #raw_input("hit return to exit")














