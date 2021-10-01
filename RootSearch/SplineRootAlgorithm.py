import math
from matplotlib.pyplot import legend
import numpy as np
from numpy.core.fromnumeric import size
# Created by Maximo Xavier DeLeon
# This code utilizes a natural spline algorithm to approximate the roots of an 'unknown' function 

test_function = lambda xi:  np.sin(xi)**4

def main():
    # PARAMETER DECLARATION
    VISUAL:bool = True # toggle for matplotlib plot
    m = 5000
    global test_function

    x = np.arange(0,20,.001) #x values
    f = test_function(x) # function to approximate

    # APPROXIMATE FUNCTIION WITH SPLINE
    spline = Spline() # create instance of spline class
    x_spline, f_spline = spline.fit(m=m,xi=x,fi=f,a_0=0,a_n=0) # perform spline approximation of the function
    print('Step size h: {}'.format((x_spline[m-1] + x_spline[0]) / m))
    
    x_spline, f_spline = np.delete(x_spline,0) , np.delete(f_spline,0) # delete the first values because they're 0,0. Basically I get under the hood bursts of stupidity and realize that after the fact

    
    #poly_array = spline.fitFormula(m=m,xi=x,fi=f,a_0=0,a_n=0)
    
    #terms = [alpha,beta,gamma,eta,xi[k],xi[k+1]]

    #formula = lambda poly_data,v:  poly_data[0]*(v-poly_data[4])*(v-poly_data[4])*(v-poly_data[4]) + poly_data[1]*(v-poly_data[5])*(v-poly_data[5])*(v-poly_data[5]) + poly_data[2]*(v-poly_data[4]) + poly_data[3]*(v-poly_data[5])


    #print('polynomial array size'.format(len(poly_array)))
    '''f_spline2 = np.zeros([len(x_spline)])

    for i in range(len(1,x_spline)):
        f_spline2[i] = formula(poly_array[i],f_spline[i])
        print(f_spline[i])'''
    
    # FIND ROOTS DYLAN METHOD
    roots, rootBins = findRoots(x_spline,f_spline,step=1,f=test_function)

    # FIND ROOT MAX METHOD
    error_tolerance = 1e-6
    ''' for root in roots:
        secant_roots = secant(n=100,dell=error_tolerance,x=root,f=f_spline)'''

    # VISUALIZE RESULTS IN MATPLOTLIB
    # if visual is true then matplotlib will be imported and then plots will be made
    if VISUAL:
        import matplotlib.pyplot as plt
        from matplotlib.pyplot import cm
        plt.figure(figsize=(8,4))
        plt.axhline(0,c='black',alpha=0.5)
        plt.plot(x_spline,f_spline,'k-',label='spline function approximation')
        plt.plot(x_spline, test_function(x_spline))
        n = len(roots)
        color = iter(cm.rainbow(np.linspace(0, 1, n)))
        for xc in roots:
            c = next(color)
            plt.axvline(x=xc, linestyle='--',c=c,label='root approximation x~{}'.format(round(xc,4)))
        #plt.plot(x_spline2,f_spline2)
        plt.scatter(x,f,label='function points') # scatter plot for the raw data points
        plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=True, ncol=4)
        plt.xlabel('x')
        plt.ylabel('f')
        plt.title('Spline Root Function Approximation')
        plt.grid()
        plt.show() # show the plot


# Credit to Dylan Bloodsworth for providing this root search idea
# for more 'intense' fucntions the mid point root approximation has a lower accuracy
# the bins are returned to perform a secant method approach since secant will allow for us to control the accuracy of each approximation
def findRoots(xi,fi,step:int,f):
    roots = []
    rootBins = []
    # iterate through x and f values while not indexing anything outside of the array because then something would break
    for i in range(0,len(xi)-step):
        f_lower = fi[i] # lower bound
        f_upper = fi[i+step] # upper bound
        # check to see if the bounds have a change of sign
        if np.sign(f_lower) * np.sign(f_upper) == -1:
            midpoint_guess = (xi[i] + xi[i+step])/2
            print('guess root (x): {} test function value: {}'.format(midpoint_guess, f(midpoint_guess)))
            roots.append((xi[i] + xi[i+step])/2)
            rootBins.append([f_lower,f_upper])
    return roots, rootBins

def secant(n:int,dell:float,x:float,dx:float,f):
        k:int = 0
        x1:float = x+dx
        while abs(dx) > dell and k<n:
            d:float = f(x1) - f(x)
            x2:float = x1-f(x1)*(x1-x)/d
            x = x1
            x1 = x2
            dx = x1-x
            k+=1
        if k==n:
            print('Convergence not found after {} iterations'.format(n))
            return 'NaN'
        else:
            return x1

class Spline:
    def __init__(self): # empty init
        pass
    # Reads the data file and outputs two numpy arrays
    def readFile(self,datafile='xy.data'):
        f = open(datafile, "r")
        non_empty_lines =[line.strip("\n") for line in f if line != "\n" and len(line) != 0]
        if non_empty_lines[0] == 'x,y': # ignore the column names in the file and carry on 
            non_empty_lines.pop(0)
        n = len(non_empty_lines) # set the n automatically
        # create zero spline approximation arrays
        xi = np.zeros([n])
        fi = np.zeros([n])
        # read in datapoints xi and fi
        # iterate trhough the data, split the text up, and make the arrays for spline interpolator
        i = 0
        for line in non_empty_lines:
            line = line.split(',')
            xi[i], fi[i] = float(line[0]), float(line[1])
            i+=1
        # return the data
        return xi, fi

    # Creates a spline fit to the data 
    def fit(self,m=100,xi=np.zeros([10]),fi=np.zeros([10]),a_0=0,a_n=0):
        n = len(xi) - 1
        p2 = np.zeros([n+1]) # second derivative of cubics
        p2 = self.cubicSpline(x=xi,f=fi) # coefficients of each polynomial
        p2[0] = a_0
        p2[n] = a_n

        # Find the approximation of the function
        h = (xi[n] - xi[0])/m # distance between points
        x = xi[0] # starting x
        

        a = np.zeros([m]) # approximation
        x1 = np.zeros([m]) # fitted x array

        for i in range(1,m):
            x+=h # step x by h
            # Find the interval where x resides
            k = 0 # k = 0 for counter start
            dx = x-xi[0]
            while dx>0: # iterate until point is found
                k+=1
                dx = x-xi[k]
            k-=1
            
            # Find the value of the function f(x)
            dx = xi[k+1]-xi[k] # steps
            alpha = p2[k+1]/(6*dx) # coef 1
            beta = -p2[k]/(6*dx) # coef 2 
            gamma = fi[k+1]/dx - dx*p2[k+1]/6 # coef 3 
            eta = dx*p2[k]/6 - fi[k]/dx # coef 4 

            # approximated cubic polynomial
            f = alpha*(x-xi[k])*(x-xi[k])*(x-xi[k]) + beta*(x-xi[k+1])*(x-xi[k+1])*(x-xi[k+1]) + gamma*(x-xi[k]) + eta*(x-xi[k+1])
            #print(x,f)
            # append the fucntion value for each value of x that is a step of the cubic spline algorithm
            a[i] = f
            x1[i] = x
        return x1,a

    def fitFormula(self,m=100,xi=np.zeros([10]),fi=np.zeros([10]),a_0=0,a_n=0):
        n = len(xi) - 1
        p2 = np.zeros([n+1]) # second derivative of cubics
        p2 = self.cubicSpline(x=xi,f=fi) # coefficients of each polynomial
        p2[0] = a_0
        p2[n] = a_n

        # Find the approximation of the function
        h = (xi[n] - xi[0])/m # distance between points
        x = xi[0] # starting x
        

        a = np.zeros([m]) # approximation
        x1 = np.zeros([m]) # fitted x array
        pw = np.array([m], dtype=object)
        
        for i in range(1,m):
            x+=h # step x by h 
            # Find the interval where x resides
            k = 0 # k = 0 for counter start
            dx = x-xi[0]
            while dx>0: # iterate until point is found
                k+=1
                dx = x-xi[k]
            k-=1
            
            # Find the value of the function f(x)
            dx = xi[k+1]-xi[k] # steps
            alpha = p2[k+1]/(6*dx) # coef 1
            beta = -p2[k]/(6*dx) # coef 2 
            gamma = fi[k+1]/dx - dx*p2[k+1]/6 # coef 3 
            eta = dx*p2[k]/6 - fi[k]/dx # coef 4 

            # approximated cubic polynomial
            #f = alpha*(x-xi[k])*(x-xi[k])*(x-xi[k]) + beta*(x-xi[k+1])*(x-xi[k+1])*(x-xi[k+1]) + gamma*(x-xi[k]) + eta*(x-xi[k+1])
            formula = lambda v:  alpha*(v-xi[k])*(v-xi[k])*(v-xi[k]) + beta*(v-xi[k+1])*(v-xi[k+1])*(v-xi[k+1]) + gamma*(v-xi[k]) + eta*(v-xi[k+1])

            terms = [alpha,beta,gamma,eta,xi[k],xi[k+1]]


            #print(x,f)
            # append the fucntion value for each value of x that is a step of the cubic spline algorithm

            pw[i] = terms

        return pw
    
    # Method to perform the cubic spline approximation
    def cubicSpline(self,x,f):
        n = len(x)-1
        p = np.zeros([n+1])
        d = np.zeros([n-1])
        b = np.zeros([n-1])
        c = np.zeros([n-1])
        g = np.zeros([n])
        h = np.zeros([n])

        # Assign the intervals and function differences

        for i in range(n):
            h[i] = x[i+1]-x[i] 
            g[i] = f[i+1]-f[i]
        
        # Evaluate the coefficient matrix elements

        for i in range(0,n-1):
            d[i] = 2*(h[i+1]+h[i])
            b[i] = 6*(g[i+1]/h[i+1]-g[i]/h[i])
            c[i] = h[i+1]


        # Obtain second order derivatives
        g = self.tridiagonalLinearEq(d, c, c, b)
        for i in range(1,n):
            p[i] = g[i-1]
        
        return p # return the 2nd order derivatives for each polynomial segment

    # performs the LU decomposition to solve for the spline's p values
    def tridiagonalLinearEq(self,d, e, c, b):
        m = len(b)
        w = np.zeros([m])
        y = np.zeros([m])
        z = np.zeros([m])
        v = np.zeros([m-1])
        t = np.zeros([m-1])

        # Evaluate the elements in the LU decomposition
        w[0] = d[0]
        v[0] = c[0]
        t[0] = e[0]/w[0]
        
        for i in range(1,m-1):
            w[i]  = d[i]-v[i-1]*t[i-1]
            v[i]  = c[i]
            t[i]  = e[i]/w[i]

        w[m-1]  = d[m-1]-v[m-2]*t[m-2]


        # Forward substitution to get y
        y[0] = b[0]/w[0]
        for i in range(1,m):
            y[i] = (b[i] - v[i-1]*y[i-1])/w[i]
        
        # Backward substitution to obtain z
        z[m-1] = y[m-1]
        i=m-2
        while i >=0: # potential source of error ----------------
            z[i] = y[i]-t[i]*z[i+1]
            i-=1
        
        return z


main()