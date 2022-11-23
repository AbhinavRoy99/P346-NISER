#!/usr/bin/env python
# coding: utf-8

# In[1]:


import math
import toolsar
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

def printhash():
    print('\n ##### \n')


def crossmat(a,b):
    c = []
    if len(a[0])==len(b):
        for i in range(0,len(a)):
            temp=[]
            for j in range(0,len(b[0])):
                s = 0
                for k in range(0,len(a[0])):
                    s += a[i][k]*b[k][j]
                temp.append(s)
            c.append(temp)
    else:
        print('Not possible.')

    return c

def matrixtranspose(l):
    r = len(l)
    c = len(l[0])
    lt = []
    for j in range(c):
        row = []
        for i in range(r):
            row.append(l[i][j])
        lt.append(row)
    return lt

def rnum(r):
    
    a = 1103515245
    c = 12345
    m = 32768

    x = ((a*r+c)%m)
        
    return x

def rnumlt1(r):
    
    a = 1103515245
    c = 12345
    m = 32768

    x = ((a*r+c)%m)/m
        
    return x

def rnumbwAB(r,A,B):
    
    a = 1103515245
    c = 12345
    m = 32768

    x = ((a*r+c)%m)/m
    
    y= A + (B-A)*x
        
    return y
    
# In[ ]:
def pivotmat(M):
    m = len(M)                                                                                                                                                                                           
    id_mat = [[float(i ==j) for i in range(m)] for j in range(m)]
    
    for j in range(m):
        row = max(range(j, m), key=lambda i: abs(M[i][j]))
        if j != row:                                                                                                                                                                          
            id_mat[j], id_mat[row] = id_mat[row], id_mat[j]

    return id_mat

def LUdecomposition(A,P):
    n = len(A)                                                                                                                                                                                                                 
    L = [[0] * n for i in range(n)]
    U = [[0] * n for i in range(n)]                                                                                                                                                     
    PA = toolsar.crossmat(P,A)                                                                                                                                                                                                                     
    for j in range(n):                                                                                                                                                                                                  
        L[j][j] = 1                                                                                                                                                                                      
        for i in range(j+1):
            s1 = sum(U[k][j] * L[i][k] for k in range(i))
            U[i][j] = PA[i][j] - s1
        
        for i in range(j, n):
            s2 = sum(U[k][j] * L[i][k] for k in range(j))
            L[i][j] = (PA[i][j] - s2) / U[j][j]
    
    print('\n The lower triangle matrix is: \n')
    for row in L:
        print(row)
    
    print('\n The upper triangle matrix is: \n')
    for row in U:
        print(row)
    return L, U

# In[ ]:

def gaussjordan(l,k):
    x=len(l)
    y=len(l[0])
    
    
    for ok in range(0,x):
        

        
        for i in range(ok+1, x):
            if l[i][ok]>l[ok][ok]:
                for j in range(ok,x):
                    l[ok][j],l[i][j]=l[i][j],l[ok][j]
                k[ok],k[i]=k[i],k[ok]
        
        t = l[ok][ok]
        
        if t != 0:
            k[ok][0] = k[ok][0]/t
            for a in range(ok,y):
                l[ok][a] = l[ok][a]/t
                
        for a in [ x for x in range(0,x) if x!=ok ]:
            lok=l[a][ok]
            k[a][0] = k[a][0] - (lok*k[ok][0])
            for b in range(0,y):
                l[a][b] = l[a][b] - (lok*l[ok][b])
    
    for a in range(0,x):
        for b in range(0,y):
            l[a][b] = round((l[a][b]),3)

    for a in range(0,x):
        k[a][0] = round((k[a][0]),3)
        
    print('\n System of linear equations after Gauss Jordan Elimination is: \n')
    for i in range(0,len(l)):
        print(l[i],'(','X',i+1,')',k[i])
    
    print('\n The answer matrix X is: \n',k)
    
    K1=[]
    for i in range(0,len(k)):
        K1.append(k[i][0])
    
    return K1


# In[ ]:

def printeq(l,k,x,key):
    print('\n',key,'system of linear equations is: \n')
    for i in range(0,len(l)):
        print(l[i],'(',x,i+1,')',k[i])

def checksym(l):
    def check(l):
        lt=toolsar.matrixtranspose(l)
        for i in range(0,len(l)):
            for j in range(0,len(l)):
                if lt[i][j] != l[j][i]:
                    return False
        return True
    if (check(l)):
        print('\n Yes the given matrix is symmetric.')
    else:
        print('\n No the given matrix is not symmetric.')

def CholeskyD(content):
    
    print('\n The matrix to be decomposed is: \n')
    for row in content:
        print(row)
    
    n=len(content)
        
    
    for j in range(0,n):
        for i in range(0,n):
            
            s=0
            
            if j == i:
                
                for k in range(0,j):
                    s = s + (content[i][k]**2)
                    
                content[i][i] = math.sqrt(content[i][i] - s)
                
                
                
            s=0
            
            if j < i:
                
                for k in range(0,j):
                    s = s + (content[i][k] * content[j][k])
                
                if(content[j][j] != 0):
                    content[i][j] = int(content[i][j] - s) / content[j][j]
                    
                
        
            
            
            if j > i:
                
                content[i][j] = 0.0
             
                
                
    print('\n The decomposed matrix is: \n')
    for row in content:
        print(row)   

    
    contenttranspose = toolsar.matrixtranspose(content)
    print('\n The decomposed matrix transpose is: \n')
    for row in contenttranspose:
        print(row)
    return content,contenttranspose
    
def forwardsub(l,k):
    k1 = []
    i=0
    for i in range(0,len(l)):
        k1.append(k[i][0])
        for j in range(0,i):
            k1[i]= k1[i]-(l[i][j]*k1[j])
        k1[i]= k1[i]/l[i][i]
    
    m, n = len(l), 1
    k2 = [[0] * n for i in range(m)]
    
    for i in range(0,len(l)):
        k2[i][0]=k1[i]
    
    return k2

def backwardsub(l,k):
    m, n = len(l), 1
    k1 = [[0] * n for i in range(m)]
    k2 = [[0] * n for i in range(m)]
    
    for i in range(0,len(l)):
        k2[i][0]=k[i][0]
    
    for i in range(len(l)-1,-1,-1):
        for j in range(i+1,len(l)):
            k2[i][0]= k2[i][0]-(l[i][j]*k1[j][0])
        
        k1[i][0]= k2[i][0]/l[i][i]
    
    return k1

# In[ ]:
def diagonaldom(l,k):
    lnew=[[0] * len(l) for i in range(len(l))]
    knew=[0] * len(l)
    for i in range(len(l)):
        for j in range(len(l)):
            lnew[i][j]=l[i][j]
            knew[i]=k[i]
    swapc=0
    for i in range(len(l)):
        summ=0
        for t in range(len(l)):
            summ = summ + abs(lnew[i][t])
        for j in range(len(l)):
            if abs(lnew[i][j])> (summ-l[i][j]):
                lnew[i],lnew[j] = lnew[j], lnew[i]
                knew[i], knew[j] = knew[j], knew[i]
    
    toolsar.printeq(lnew,knew,'X','The diaonally dominant')
    return lnew,knew

# In[ ]:

def arrclose(x,y,tol):
    count=0
    if len(x)==len(y):
        for i in range(0,len(x)):
            if (abs(x[i]-y[i]))/abs(y[i]) < tol:
                count =count+1
            else:
                return False
    if count==len(x):
        return True

# In[ ]:

def JacobiM(a,b,it,tol):
    print("\n Jacobi Calculations: \n")
    
    x=[1]*len(a)
    newx=[1]*len(a)
    
    
    for i in range(0,it):
        t=0
        for t in range(0,len(x)):
            x[t] = newx[t]
        
        newx = toolsar.jac(a,b,newx)
    
        print(i+1,".",newx)
    
        z = toolsar.arrclose(newx,x,tol)
        if z == True:
            print('\n The answer is: \n',newx,"\n after",i+1,"iterations.")
            break
        else:
            continue

def jac(a,b,x):
    for i in range(0, len(a)):
        d=b[i]
        for j in range(0, len(a)):  
            if (j != i):
                d = d - (a[i][j]*x[j])
        
        x[i] = d / a[i][i]
        
    return x
# In[ ]:
def gsm(a,x,b):
                        

    for j in range(0, len(a)):               
          
        d=b[j]
        
        for i in range(0, len(a)):     
            if(j != i):
                d = d - (a[j][i]*x[i])
        
        x[j] = d / a[j][j]
        
    return x

def GaussSeidel(a,x,b,it,tol):
    #print("\n Gauss Seidel Calculations: \n")
    
    y=[]
    
    for i in range(0,len(x)):
        y.append(x[i])
    
    for i in range(0,it):
        t=0
        for t in range(0,len(x)):
            y[t] = x[t]
    
        x = toolsar.gsm(a,x,b)
    
        #print(i+1,".",x)
    
        z = toolsar.arrclose(x,y,tol)
        if z == True:
            print('\n The answer is: \n',x,"\n after",i+1,"iterations.")
            break
        else:
            continue
            

# In[ ]:

def Bracket(a,b,func,e):
    
    co=0
    
    while (func(a)*func(b)>=0):
        if (abs(func(a)) < abs(func(b))):
            a = a - e*(b-a)
            co=co+1
            #print(a,b)
            
            
        if (abs(func(a)) > abs(func(b))):
            b = b + e*(b-a)
            co=co+1
            #print(a,b)
    
    print('Interval:',a,b)
    return (a,b)

            

def Bisect(a,b,e,d,func):
    print('BISECT')
    c=(a+b)/2
    t=0
    print(t,">>",'interval:',a,b,'start root c:',c)
    
        
    while abs(a-b)>e:
        while abs(func(a))>d:
            c=(a+b)/2
            
            
            if abs(func(c))<=d:
                print(t+1,">>",'c:',c)
                print("the solution is",'%10.6E' %c,"in",t+1,"iterations.")
                return c
            else:  
                if (func(a)*func(c))<0:
                    b=c
                    t=t+1
                    #print(t,">>",'c:',c)
                    
                else:
                    a=c
                    t=t+1
                    #print(t,">>",'c:',c)
# In[ ]:

def RegulaFalsi(a,b,e,d,fn):
    print('REGULA FALSI')
    c0=a
    c1=b
    t=0
    print(t,">>",'interval:',a,b,'start prev root c:',c0,'start next root c:',c1)
    
    while abs(c1-c0)>e:
        while abs(fn(a))>d:
            c1 = b-((b-a)*fn(b))/(fn(b)-fn(a))
            
            if abs(fn(c1))<=d:
                print(t+1,">>",'c:',c1)
                print("the solution is",'%10.6E' %c1,"in",t+1,"iterations.")
                return c1
            else:
                if (fn(a)*fn(c1))<0:
                    b=c1
                    c0=c1
                    t=t+1
                    #print(t,">>",'c:',c1)
                    
                else:
                    a=c1
                    c0=c1
                    t=t+1
                    #print(t,">>",'c:',c1)

# In[ ]:

def NewtonRaphson(x,e,d,fn,fd):
    print('NEWTON RAPHSON')
    x0=x-1
    x1=x
    t=0
    print(t,">>",'start prev root guess:',x0,'start new root guess:',x1)
    
    while abs(x1-x0)>e:
        while abs(fn(x1))>d:
            x1 = x1-(fn(x1)/fd(x1))
            
            if abs(fn(x1))<=d:
                print(t+1,">>",' new root:',x1)
                print("the solution for x is",'%10.6E' %x1,"in",t+1,"iterations.")
                return x1
            else:
                t=t+1
                #print(t,">>",'new root guess',x1)
                x0=x1
                
# In[ ]:

def LinearFit(X,Y):
    sigma = [1]*len(X)
    plt.scatter((X),Y)
    plt.show()  
    sum_x = 0
    sum_y = 0
    sum_x2 = 0
    sum_xy = 0
    sum_y2 = 0
    N = len(X)
    for i in range(0,len(X)):
        sum_x = sum_x + (X[i]/sigma[i])
        sum_x2 = sum_x2 + ((X[i]**2)/sigma[i])
        sum_y = sum_y + (Y[i]/sigma[i])
        sum_xy = sum_xy + ((X[i]*Y[i])/sigma[i])
        sum_y2 = sum_y2 + ((Y[i]**2)/sigma[i])
    a1 = (N*sum_xy-sum_x*sum_y)/(N*sum_x2-(sum_x**2))
    print('Slope:',a1)
    a2 = (sum_y - a1*sum_x)/(N)
    print('Y Intercept:',a2)
   
    delta_x = N*sum_x2-(sum_x**2)
    delta_y = N*sum_y2-(sum_y**2)
    r = ((N*sum_xy - sum_x*sum_y)**2)/(delta_x*delta_y)
    print('R^2:',r)
    return a1,a2,r

# In[ ]:

def Laguerre(W,x0,e):
    
    def f1(r,W):
        h=0
        for i in range(0,len(W)):
            h=h+(W[i]*(r**i))
        return h

    def f2(r,W):
        h=0
        for i in range(0,len(W)):
            h=h+(W[i]*i*(r**(i-1)))
        return h

    def f3(r,W):
        h=0
        for i in range(0,len(W)):
            h=h+(W[i]*i*(i-1)*(r**(i-2)))
        return h
    
    
    n=len(W)
    xk = x0
    while abs(f1(xk,W)) > e:
        G = f2(xk,W) / f1(xk,W)
        H = G ** 2 - (f3(xk,W) / f1(xk,W))
        root = math.sqrt((n - 1) * (n * H - (G ** 2)))
        d = max([G + root, G - root], key=abs)
        a = n / d
        xk -= a
    return xk

def SynDivision(W,root):
    WW = []
    for i in range(0,len(W)):
        WW.append(W[i])
    if len(W) != 1:
        WW[len(W)-2] = WW[len(W)-2] + WW[len(W)-1]*root
        for i in range(3,len(W)):
            WW[len(W)-i] = WW[len(W)-i] + root*WW[len(W)-i+1]
        WW.pop(0)
    else:
        print("cannot deflate")
    return WW

def LagRoots(W,x0,e):
    
    n0=len(W)
    
    def f11(r,A):
        h=0
        for i in range(0,len(A)):
            h=h+(A[i]*(r**i))
        return h

    def f22(r,A):
        h=0
        for i in range(0,len(A)):
            h=h+(A[i]*i*(r**(i-1)))
        return h

    def f33(r,A):
        h=0
        for i in range(0,len(A)):
            h=h+(A[i]*i*(i-1)*(r**(i-2)))
        return h
    
    WW = []
    RR = []
    for i in range(0,n0):
        WW.append(W[i])
    
    n00 = len(WW)
    
    for i in range(0,n0-2):
    
        r1 = toolsar.Laguerre(WW, x0, e)
        RR.append(r1)
        
        WW=toolsar.SynDivision(WW,r1)
        #print('M:',WW)
        
        n00=len(WW)
        #print('X:',n00,r1)
        
    if n00 == 2:
        RR.append(-1*(WW[0]/WW[1]))
        print('Calculations done. \n')
        print('Polynomial:')
        toolsar.printPoly(W)
        print('Roots: \n',RR)
        
    return RR

def printPoly(W):
    n=len(W)
    for i in range(0,n):
        print('+(',W[n-i-1],') x^',n-i-1, end=" ")
    print('= 0')
    print('\n')
    
#In[ ]:

def PolyFitLS(X,Y,n):
    def L(ds,n):
        L = [[0] * n for i in range(n)]
        for i in range(0,n):
            for j in range(0,n):
                for t in range(0,len(ds)):
                    L[i][j]=L[i][j]+(ds[t]**(i+j))
    
        return L
    
    def K(ds,ht,n):
        K = [[0] * 1 for i in range(n)]
        for i in range(0,n):
            for t in range(0,len(ds)):
                K[i][0] = K[i][0] + (ht[t]*(ds[t]**(i)))
        return K
    def ftrac(r):
        h=0
        for i in range(0,len(WW)):
            h=h+(WW[i]*(r**i))
        return h
    
    L1=L(X,n)

    K1=K(X,Y,n)

    toolsar.printeq(L1,K1,'A','The solved matrix form is')

    toolsar.printhash()
    
    WW=[]
    WW=toolsar.gaussjordan(L1,K1)
    print('\n')

    arr_f=[]

    for i in range(0,len(X)):
        arr_f.append(ftrac(X[i]))

    plt.plot(X,arr_f,'bo-')
    plt.plot(X,Y,'ro')
    
#In[ ]:

#def d2xdt2(t, x, v):
#    return -1*0.1*v - 0.51*0.51*x

#def dxdt(t, x, v):
#    return v

#x0_1, y0_1 = RKXY(d2xdt2, dxdt, 0, 1, 1, 50, 0.5)

def RKXY(d2ydx2, dydx, x0, y0, z0, xf, st, yes):
   
    x = [x0]
    y = [y0]
    z = [z0]      # dy/dx

    n = int((xf-x0)/st)     # no. of steps
    for i in range(n):
        x.append(x[i] + st)
        k1 = st * dydx(x[i], y[i], z[i])
        l1 = st * d2ydx2(x[i], y[i], z[i])
        k2 = st * dydx(x[i] + st/2, y[i] + k1/2, z[i] + l1/2)
        l2 = st * d2ydx2(x[i] + st/2, y[i] + k1/2, z[i] + l1/2)
        k3 = st * dydx(x[i] + st/2, y[i] + k2/2, z[i] + l2/2)
        l3 = st * d2ydx2(x[i] + st/2, y[i] + k2/2, z[i] + l2/2)
        k4 = st * dydx(x[i] + st, y[i] + k3, z[i] + l3)
        l4 = st * d2ydx2(x[i] + st, y[i] + k3, z[i] + l3)

        y.append(y[i] + (k1 + 2*k2 + 2*k3 + k4)/6)
        z.append(z[i] + (l1 + 2*l2 + 2*l3 + l4)/6)
        
    if str(yes)=='yes':
        plt.plot(x,y,'ro-',label='Final Plot')
        plt.legend()
        plt.show()
        print('Line here is not a fitting of the polynomial. Has been added to aid the eye to track the points.')

    return x, y, z

#def dxdt(x,y,z,t):
#    return 10*(y-x)

#def dydt(x,y,z,t):
#    return x*(28-z) - y

#def dzdt(x,y,z,t):
#    return x*y - 8*z/3

#arrx, arry, arrz, arrt = RKXYZ(dxdt, dydt, dzdt, 1, 0, 0, 0, 25, 0.1, 'yes')

def RKXYZ(dxdt, dydt, dzdt, x0, y0, z0, t0, tf, st):
   
    x = [x0]
    y = [y0]
    z = [z0]
    t = [t0] 

    n = int((tf-t0)/st)     # no. of steps
    for i in range(n):
        t.append(t[i] + st)
        k1x = st * dxdt(x[i], y[i], z[i],t[i])
        k1y = st * dydt(x[i], y[i], z[i],t[i])
        k1z = st * dzdt(x[i], y[i], z[i],t[i])
        
        k2x = st * dxdt(x[i] + k1x/2, y[i] + k1y/2, z[i] + k1z/2,t[i] + st/2)
        k2y = st * dydt(x[i] + k1x/2, y[i] + k1y/2, z[i] + k1z/2,t[i] + st/2)
        k2z = st * dzdt(x[i] + k1x/2, y[i] + k1y/2, z[i] + k1z/2,t[i] + st/2)
        
        k3x = st * dxdt(x[i] + k2x/2, y[i] + k2y/2, z[i] + k2z/2,t[i] + st/2)
        k3y = st * dydt(x[i] + k2x/2, y[i] + k2y/2, z[i] + k2z/2,t[i] + st/2)
        k3z = st * dzdt(x[i] + k2x/2, y[i] + k2y/2, z[i] + k2z/2,t[i] + st/2)
        
        k4x = st * dxdt(x[i] + k3x, y[i] + k3y, z[i] + k3z,t[i] + st)
        k4y = st * dydt(x[i] + k3x, y[i] + k3y, z[i] + k3z,t[i] + st)
        k4z = st * dzdt(x[i] + k3x, y[i] + k3y, z[i] + k3z,t[i] + st)
        

        x.append(x[i] + (k1x + 2*k2x + 2*k3x + k4x)/6)
        y.append(y[i] + (k1y + 2*k2y + 2*k3y + k4y)/6)
        z.append(z[i] + (k1z + 2*k2z + 2*k3z + k4z)/6)

    sphere = plt.axes(projection='3d')
    sphere.plot(x, y, z,'o-')
    
    return x, y, z, t

#def d2ydx2(x, y, z):
#    return 2*y

#def dydx(x, y, z):
#    return z

#x0_1, y0_1, z0_1 = Shoot(d2ydx2, dydx, 0, 1, 1.2, 0.9, -1.5, 0.1, 0.001, 'yes')

def ShootXY(d2ydx2, dydx, x0, xf, y0, yf, z0, dx, tol, yes):
    yff=yf+1
    r=10
    
    def LagI(zh,zl,yh,yl,yf):
        return zl + ((zh-zl)*(yf-yl)/(yh-yl))
    
    while abs(yff-yf) >= tol:
        
        zguess=z0
        
        x0_1, y0_1, z0_1 = toolsar.RKXY(d2ydx2, dydx, x0, y0, z0, xf, dx, 'no')
        plt.plot(x0_1,y0_1,'o',label='guess plot')
        
        n = len(y0_1)-1
        #print(n)
        
        zguess=toolsar.rnumlt1(zguess)
        x0_11, y0_11, z0_11 = toolsar.RKXY(d2ydx2, dydx, x0, y0, zguess, xf, dx, 'no')
        plt.plot(x0_11,y0_11,'o',label='guess plot')
        
        #print('\n',x0_1,'\n',y0_1,'\n',y0_11)
        
        if yf > y0_1[n]:
            #print('l0')
            if yf > y0_11[n]:
                zguess=toolsar.rnumlt1(zguess)
                x0_11, y0_11, z0_11 = toolsar.RKXY(d2ydx2, dydx, x0, y0, zguess, xf, dx, 'no')
                plt.plot(x0_11,y0_11,'o',label='guess plot')
                zl = z0_1[0]
                yl = y0_1[n]
                zh = z0_11[0]
                yh = y0_11[n]
                #print('l1')
            else:
                zh = z0_1[0]
                yh = y0_1[n]
                zl = z0_11[0]
                yl = y0_11[n]
                #print('l2')
                
        if yf < y0_1[n]:
            #print('h0')
            if yf < y0_11[n]:
                zguess=toolsar.rnumlt1(zguess)
                x0_11, y0_11, z0_11 = toolsar.RKXY(d2ydx2, dydx, x0, y0, zguess, xf, dx, 'no')
                plt.plot(x0_11,y0_11,'o',label='guess plot')
                zh = z0_1[0]
                yh = y0_1[n]
                zl = z0_11[0]
                yl = y0_11[n]
                #print('h1')
            else:
                zl = z0_1[0]
                yl = y0_1[n]
                zh = z0_11[0]
                yh = y0_11[n]
                #print('h2')
                
        z0 = LagI(zh,zl,yh,yl,yf)
        #print('L',z0)
        
        yff = y0_1[n]
    
    if str(yes)=='yes':
        plt.plot(x0_1,y0_1,'ro-',label="Final Plot")
        plt.legend()
        plt.show()
        print('Line here is not a fitting of the polynomial. Has been added to aid the eye to track the points.')
    return x0_1, y0_1, z0_1


#In[ ]:

#def ux0(x):
#    y = math.pi*x
#    return 20*abs(math.sin(y))

#T1,L1,T0 = ExpPDE(ux0, 2, 4, 20, 5000, 50)

def ExpPDE(ux0, lx, lt, nx, nt, nt1):
    ht = lt/nt
    hx = lx/nx
    #alpha = ht/(hx**2)
    alpha=0.2
    
    #print(ht,hx,alpha)
    
    Vt = [0]*(nx+1)
    V0 = [0]*(nx+1)
    V00 = [0]*(nx+1)
    LX = [0]*(nx+1)
    
    x1=0
    
    for i in range(0,nx+1):
        x1 = i*hx
        #print('x',x1)
        LX[i] = x1
        V0[i] = ux0(x1)
    
    #print(len(V0))
    for i in range(len(V0)):
        Vt[i] = V0[i]
        V00[i] = V0[i]
    
    for j in range(0,nt1):
        #print(j)
        for i in range(0,nx+1):
        
            if i == 0:
                Vt[i] = ((1-(2*alpha))*V00[i]) + (alpha*V00[i+1])
                #print('a',i,j)
            
            if i == nx:
                Vt[i] = ((1-(2*alpha))*V00[i]) + (alpha*V00[i-1])
                #print('b',Vt[i])
            
            else:
                Vt[i] = ((1-(2*alpha))*V00[i]) + (alpha*V00[i+1]) + (alpha*V00[i-1])
                #print('c',Vt[i])
        
        for i in range(len(V0)):
            V00[i] = Vt[i]
    
    plt.plot(LX,V0,'ro-',label='Time Step = 0')
    plt.plot(LX,Vt,'o-',label='Time Step Specified')
    plt.legend()
    plt.show()
    
    return Vt, LX, V0


#In[ ]:

#nA=500
#nB=0
#nC=0
#tA=25
#tB=35
#T=100
#dT=2
#NA,NB,NC,NT,NcountA,NcountB=DecayABC(nA,nB,nC,tA,tB,T,dT)


def DecayABC(nA,nB,nC,tA,tB,T,dT):
    lamA = math.log(2)/tA
    lamB = math.log(2)/tB
    
    Nmax=nA
    
    r = toolsar.rnum(T)
    t=0
    
    NA=[nA]
    NB=[0]
    NC=[0]
    NT=[0]
    NcountA=[0]
    NcountB=[0]
    i=0
    count=0
    while i in range(0,int(T/dT)) and nA >=0 and  nB<=Nmax and nB >=0 and nC<=Nmax: #defining conditions for the loop to stop after all decays possible
        countA=0
        countB=0
        for j in range(0,nA):#trying the decay for all the nuclei of A
            if r <= lamA*dT:
                nA = nA -1
                nB = nB +1
                countA += 1
        
            r=toolsar.rnumlt1(r)
        
        for j in range(0,nB):#trying the decay for all the nuclei of B
            if r <= lamB*dT:
                nB = nB -1
                nC = nC +1
                countB+=1
        
            r=toolsar.rnumlt1(r)
        
        NA.append(nA)
        NB.append(nB)
        NC.append(nC)
        t +=dT
        NT.append(t)
        NcountA.append(countA)
        NcountB.append(countB)
        i=i+1
    
    plt.plot(NT,NA,'bo-',label="A")
    plt.plot(NT,NB,'ro-',label="B")
    plt.plot(NT,NC,'go-',label="C")
    plt.legend()
    plt.show()
    #plt.plot(NT,NcountA,'o-')
    plt.hist(NcountA, bins = max(NcountA),rwidth=0.9)
    plt.show()
    #plt.plot(NT,NcountB,'o-')
    plt.hist(NcountB, bins = max(NcountB),rwidth=0.9)
    plt.show()
    
    return NA,NB,NC,NT,NcountA,NcountB

#In[ ]:

#import toolsar as ar
#import math

#A=[[1,-1,0], [-2,4,-2],[0,-1,2]]

#lambda1,v1=ar.EigenPI(A,20,'random')

def EigenPI(A, k_max, tol, guess):
    na = len(A)
    b_k=[]
    lambda1=0
    if str(guess)=='random':
        for i in range(0,na):
            R=[]
            r=toolsar.rnum(i)
            R.append(r)
            b_k.append(R)
    else:
        b_k = guess
    print('\n Guess vector:',b_k)
    
    e_new=1
    e_old=0
    i=0
    while (abs(e_new - e_old) >= tol) and (i in range(0,k_max)):
        e_old=lambda1
        b_k = toolsar.crossmat(A,b_k)
        b_k_n = 0
        for j in range(0,len(b_k)):
            b_k_n += (b_k[j][0])**2
        
        for j in range(0,len(b_k)):
            b_k[j][0] = b_k[j][0] / math.sqrt(b_k_n)
        
        i=i+1
        
        #finding the value of b_k*Ab_k
        x1=toolsar.matrixtranspose(b_k)
        x2=toolsar.crossmat(x1,A)
        x3=toolsar.crossmat(x2,b_k)
        #finding value of b_k*b_k
        v=toolsar.matrixtranspose(b_k)
        v1=toolsar.crossmat(v,b_k)
        #lambda1 = b_k*Ab_k / b_k*b_k
        lambda1=x3[0][0]/v1[0][0]
        
        e_new=lambda1
    
    
    
    print('\n Matrix A:',A,'\n \n Eignvector:',b_k,'\n \n Eigenvalue:',lambda1,'\n \n in', i,'iterations.')
    return b_k,lambda1