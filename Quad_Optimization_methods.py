import math
import numpy as np
from sympy import*
def function(t):
    #F=t/math.log10(t)
    #F = 0.65 - (0.75/(1 + pow(t,2))) - 0.65*t*(math.atan(pow(t,-1)))
    F = pow(t,5) - 5*pow(t,3) - 20*t + 5
    return F

def simplex(x1,x2,x3,alpha,beta,gamma,epsilon):
    fs =[]
    xs=[]
    xs.append(x1)
    xs.append(x2)
    xs.append(x3)
    fs.append(multi_function(x1))
    fs.append(multi_function(x2))
    fs.append(multi_function(x3))
    xh = fs.index(max(fs))
    xl = fs.index(min(fs))
    fs[xh] = -9999
    fs[xl] = -9999
    xm = fs.index(max(fs))
    fs = []
    fs.append(multi_function(x1))
    fs.append(multi_function(x2))
    fs.append(multi_function(x3))
    xlp = xs[xl]
    xhp = xs[xh]
    xmp = xs[xm]

    daviation = epsilon
    # reflection
    while daviation >= epsilon:
        fs = []
        fs.append(multi_function(xlp))
        fs.append(multi_function(xhp))
        fs.append(multi_function(xmp))
        xs = []
        xs.append(xlp)
        xs.append(xhp)
        xs.append(xmp)
        xh = fs.index(max(fs))
        xl = fs.index(min(fs))
        fs[xh] = -9999
        fs[xl] = -9999
        xm = fs.index(max(fs))
        xlp = xs[xl]
        xhp = xs[xh]
        xmp = xs[xm]
        x_bar = (xlp+xmp)/2
        xrp = x_bar + alpha*(x_bar-xhp)
        f_r = multi_function(xrp)
        if f_r < multi_function(xlp):
            xep = x_bar + gamma*(xrp - x_bar)
            #if multi_function(xep) < f_r:
            xhp =xep
            xs[xh] = xhp
        elif f_r > multi_function(xlp) and f_r < multi_function(xhp):

            xcp = beta*xhp + (1-beta)*(x_bar)
            #if multi_function(xcp) < f_r and multi_function(xcp) < multi_function(xlp):
            xhp =xcp
            xs[xh] = xhp
        elif f_r > multi_function(xhp):
            xcp = beta * xhp + (1 - beta) * (x_bar)
            #if multi_function(xcp) < multi_function(xhp):
            xhp = xcp
            xs[xh] = xhp
        term1 = multi_function(xlp) - multi_function(x_bar)
        term2 = multi_function(xhp) - multi_function(x_bar)
        term3 = multi_function(xmp) - multi_function(x_bar)
        daviation = pow((1/(len(x1)+1))*(pow(term1,2) + pow(term2,2) + pow(term3,2)),0.5)
        print("simplex centroid and daviation:",x_bar,daviation)
    return x_bar, daviation
def multi_function(xi):
    #F =  xi[0] - xi[1] + 2*pow(xi[0],2) + 2*xi[0]*xi[1] + pow(xi[1],2) #pow(xi[0],4) - 2*xi[1]*pow(xi[0],2) + pow(xi[0],2) + pow(xi[1],2) - 2*xi[0] +1
    F = pow((xi[0] + 2*xi[1] - 7),2) + pow((2*xi[0] + xi[1] - 5),2)
    #F = pow(xi[0],4) - 2*pow(xi[0],2)*xi[1] + pow(xi[0],2) + pow(xi[1],2) - 2*xi[0] +1
    #F = 2*pow(xi[1] - pow(xi[0],2),2) + pow(1-xi[0],2)
    return F
def P_x(a0,a1,a2,n):
    P = a0 + a1*n + a2*(pow(n,2))
    return P

def unrestricted_fixed(lamda,step):
    lamdas = []
    lamdas.append(lamda + step)
    lamdas.append(lamdas[0] + step)
    lamdas.append(lamdas[1] + step)
    while function(lamdas[2]) < function(lamdas[1]):
        lamdas.append(lamdas[2] +step)
        del lamdas[0]
    avg = (lamdas[2] + lamdas[0])/2

    return avg

def unrestricted_accel(lamda,step):
    lamdas = []
    steps = []
    lamdas.append(lamda + step)
    steps.append(step)
    lamdas.append(lamdas[0] + 2*step)
    steps.append(2*step)
    lamdas.append(lamdas[1] + 4*step)
    steps.append(4*step)
    while function(lamdas[2]) < function(lamdas[1]):
        steps.append(1.5*steps[2])
        del steps[0]
        lamdas.append(lamdas[2] +steps[2])
        del lamdas[0]

    avg = (lamdas[2] + lamdas[0])/2

    return avg

def exhaustive_search(n,step):
    lamdas = []
    lamdas.append(0 + step)
    lamdas.append(lamdas[0] + step)
    lamdas.append(lamdas[1] + step)
    for i in range(0,n,1):
        lamdas.append(lamdas[2] + step)
        del lamdas[0]
        if function(lamdas[2]) > function(lamdas[1]):
            break
    avg = (lamdas[2] + lamdas[0]) / 2

    return avg
def Newton_method(df,d2f,lamda,epsilon):
    while abs(df(lamda)) > epsilon:
        lamda = abs(lamda - (df(lamda) / d2f(lamda)))
    return lamda
def seacant_method(df,epsilon,A,t0):
    A = A + t0
    while df(A) < 0:
        B = A+t0
        if df(B) > 0:
           break
        else:
            A = B
            t0 = 2 * A
    while abs(df(B)) > epsilon:
        B = A - ((df(A)*(B-A))/(df(B)-df(A)))
    return B

def Quasi_Newton(x1,del_x,epsilon):
    convergence = epsilon
    while convergence >= epsilon:
        f1 = function(x1)
        f_plus = function(x1+del_x)
        f_neg = function(x1 - del_x)
        book=f_plus+ f_neg
        x2 = x1 - ((del_x*(f_plus - f_neg))/(2*(f_plus - 2*f1 + f_neg)))
        convergence = abs((function(x2 + del_x) - function(x2 - del_x))/(2*del_x))
        x1 = x2
    return x2


def cubic_interpolation(f,df,acc,x1,x2):
    while true:
        b1 = f(x1)
        b2 = df(x1)
        b3 = f(x2)
        b4 = df(x2)
        B = np.array([[b1], [b2], [b3], [b4]])
        A = np.array([[1, x1, pow(x1, 2), pow(x1, 3)],[ 0, 1, 2*x1, 3*pow(x1, 2)], [1, x2, pow(x2, 2), pow(x2, 3)] , [0, 1, 2*x2, 3*pow(x2, 2)]])
        X = np.linalg.solve(A, B)
        Z=(3*(b1-b3))/(x2-x1) + b2+b4
        a1 = (1/pow((x1-x2),2))*(pow(x2,2)*b2 + pow(x1,2)*b4 + 2*Z*x1*x2)#X[1,0]
        a2 = -(1/(pow((x1-x2),2)))*((x1+x2)*Z + x2*b2 + x1*b4)#X[2,0]
        a3 = (1/(3*pow((x1-x2),2)))*(2*Z + b2 + b4)#X[3,0]

        x_star = (-a2 + pow((pow(a2,2)-3*a1*a3),0.5))/(3*a3)#float(x1 + ((-1*a1)/(a2 + sqrt((keera)))))
        check = abs(df(x_star))
        if check < acc:
            break
        else:
            x1 = x_star

    return x_star


    #print()

def quad_interpolation(x1,t_not,acc):
    fa = function(x1)
    temp = function(x1)
    f1 = function(t_not)
    while f1<temp:
        fb = f1
        x2 = t_not
        t_not = 2*t_not
        f1 = function(t_not)
    fc = f1
    x3 = t_not
    lamda_star = (4*fb - 3*fa - fc)/(4*fb - 2*fc - 2*fa)
    a0 = (function(x1)*x2*x3*(x3-x2) + function(x2)*x3*x1*(x1-x3) + function(x3)*x1*x2*(x2-x1))/((x1-x2)*(x2-x3)*(x3-x1))
    a1 = (function(x1)*(pow(x1,2)-pow(x3,2)) + function(x2)*(pow(x3,2)-pow(x1,2)) + function(x3)*(pow(x1,2)-pow(x2,2)))/((x1-x2)*(x2-x3)*(x3-x1))
    a2 = -1*(function(x1) * (x2 - x3) + function(x2)*(x3 - x1) + function(x3)*(x1- x2)) / ((x1 - x2) * (x2 - x3) * (x3 - x1))
    x_star = -1*(a1/(2*a2))#-1*(a1/(2*a2))
    error = acc
    x_prev = 9999
    while error >= acc:
        a0 = (function(x1) * x2 * x3 * (x3 - x2) + function(x2) * x3 * x1 * (x1 - x3) + function(x3) * x1 * x2 * (x2 - x1)) / ((x1 - x2) * (x2 - x3) * (x3 - x1))
        a1 = (function(x1) * (pow(x1, 2) - pow(x3, 2)) + function(x2) * (pow(x3, 2) - pow(x1, 2)) + function(x3) * (pow(x1, 2) - pow(x2, 2))) / ((x1 - x2) * (x2 - x3) * (x3 - x1))
        a2 = -1 * (function(x1) * (x2 - x3) + function(x2) * (x3 - x1) + function(x3) * (x1 - x2)) / ((x1 - x2) * (x2 - x3) * (x3 - x1))
        x_star = -1*(a1/(2*a2))#(function(x1)*(pow(x2,2) - pow(x3,2)) + function(x2)*(pow(x3,2) - pow(x1,2)) + function(x3)*(pow(x1,2) - pow(x2,2)))/(2*(function(x1)*(x2 - x3) + function(x2)*(x3 - x1) + function(x3)*(x1 - x2)))
        error = abs((x_star - x_prev)/x_prev)
        x_prev = x_star
        if x_star < x2 and function(x_star) < function(x2):
            x1 = x1
            x3 = x2
            x2 = x_star

        elif x_star < x2 and function(x_star) > function(x2):
            x1 = x_star
            x2 = x2
            x3 = x3
        elif x_star > x2 and function(x_star) < function(x2):
            x1 = x2
            x2 = x_star
            x3 =x3
        elif x_star > x2 and function(x_star) > function(x2):
            x3 = x_star
    return x_star

def feb_seq(n):
    feb = []
    feb.append(1)
    feb.append(1)
    for i in range(2,n,1):
        feb.append(feb[i-1] + feb[i - 2])
    return feb

def interval_halving(a,b,acc,max_min):
    L_not = abs(b - a)
    Ln = acc * L_not
    while Ln>acc:
        L_not = abs(b-a)
        Ln = acc*L_not
        quater = L_not/4
        x1 =  a + quater
        x0 = x1 + quater
        x2 = x0 + quater
        fx1 = function(x1)
        fx2 = function(x2)
        fx0 = function(x0)
        if max_min == 0:
            if fx2 > fx0 > fx1:
                b = x0
                a = a
            elif fx2 < fx0 < fx1:
                a = x0
                b = b
            elif fx1>fx0 and fx2>fx0:
                b =x2
                a =x1
        elif max_min == 1:
            if fx2 > fx0 > fx1:
                b = x0
                a = a
            elif fx2 < fx0 < fx1:
                b = x2
                a = a
            elif fx1<fx0 and fx2<fx0:
                b = x2
                a = x1

    return a,b

def dichotomous_search(a,b,acc,sigma,max_min):
    L_not = abs(b - a)
    Ln = acc * L_not
    while Ln>acc:
        L_not = abs(b-a)
        Ln = acc*L_not
        half = a + L_not/2
        x1 = half - sigma/2
        x2 = half + sigma/2
        fx1 = function(x1)
        fx2 = function(x2)
        if max_min == 0:
            if fx1 > fx2:
                a = x1
                b = b
            else:
                b = x2
                a = a
        elif max_min == 1:
            if fx1 < fx2:
                a = x1
                b = b
            else:
                b = x2
                a = a
    return a,b

def fibonacci_search(a,b,epsilon,feb,max_min):
    val = (pow(epsilon*2,-1))
    for i in range(0,len(feb),1):
        if feb[i]>=val:
            index = i
            break

    L_not = abs(b - a)
    i = 0

    while i != index-1:
        i = i+1
        L_star = ((feb[index-(i+1)])/(feb[index]))*L_not
        x1 = a + L_star
        x2 = b - L_star
        fx1 = function(x1)
        fx2 = function(x2)

        if max_min == 0:
            if fx1 > fx2:
                a = x1
                b = b
            else:
                b = x2
                a = a
        elif max_min == 1:
            if fx1 < fx2:
                a = x1
                b = b
            else:
                b = x2
                a = a
    return a,b

def golden_section(a,b,epsilon,max_min):
    L = epsilon
    L_not = abs(b - a)
    i =0
    while L>=epsilon:
        i = i+1
        L_star = (1/(pow(1.618,i+1)))*L_not
        L = (1/(pow(1.618,i)))*L_not
        x1 = a + L_star
        x2 = b - L_star
        fx1 = function(x1)
        fx2 = function(x2)

        if max_min == 0:
            if fx1 > fx2:
                a = x1
                b = b
            else:
                b = x2
                a = a
        elif max_min == 1:
            if fx1 < fx2:
                a = x1
                b = b
            else:
                b = x2
                a = a
    return a,b

x = symbols('x')
f = pow(x,5) - 5*pow(x,3) - 20*x +5#0.65 - (0.75/(1 + pow(x,2))) - 0.65*x*(atan(pow(x,-1)))#  x/log(x)
df = diff(f)
d2f = diff(df)
print(df)
f = lambdify(x, f)
df = lambdify(x, df)
d2f = lambdify(x, d2f)



feb = feb_seq(10)#create fibonacci sequence for 10 iterations

a,b = golden_section(0.1,3,0.1,0)#(a,b,accuracy,min or max)
print("Golden section Interval: ",a,b)

a,b = fibonacci_search(0.2,3,0.05,feb,0)#(a,b,accuracy,febonacci_sequence,min or max)
print("fibonacci Interval: ",a,b)

a,b = dichotomous_search(0.2,3,0.005,0.001,0)#(a,b,accuracy,delta,min or max)
print("dichotomous Interval: ",a,b)

a,b = interval_halving(0.2,3,0.005,0)#(a,b,accuracy,min or max)
print("Interval halving Interval: ",a,b)

a= quad_interpolation(0.1,0.5,0.01)#x1,x2,x3,acc
print("Quad interpolation Approx sol: ",a)

x_star = cubic_interpolation(f,df,0.1,0,3.2)#f,df,acc,x1,x2
print("Cubic interpolation Approx sol: ",x_star)

x_star = Quasi_Newton(0.1,0.01,0.1)#x1,del_x,epsilon
print("Quasi Newton Approx sol: ",x_star)

x_star = seacant_method(df,0.01,0,0.1)#df,epsilon,A,t0
print("seacant Approx sol: ",x_star)

x_star = Newton_method(df,d2f,0.2,0.01)#df,d2f,lamda,epsilon
print("Newton method Approx sol: ",x_star)

x_star = unrestricted_fixed(0.2,0.1)#lamda,step
print("unrestricted fixed step Approx sol: ",x_star)

x_star = exhaustive_search(20,0.1)#n,step
print("exhaustive search Approx sol: ",x_star)

x_star = unrestricted_accel(0.2,0.01)#lamda,step
print("unrestricted accel Approx sol: ",x_star)

'''
x1 = np.array([-2, -2])
x2 = np.array([-3, 0])
x3 = np.array([-1, -1])
opt_sol,daviation = simplex(x1,x2,x3,1,0.5,2,0.2)#x1,x2,x3,alpha,beta,gamma,epsilon
print("simplex centroid and daviation:",opt_sol,daviation)'''