#!/usr/bin/env python
# coding: utf-8

# In[1]:


import math
import toolsar
import matplotlib.pyplot as plt

def printhash():
    print('\n ##### \n')

def crossmat(a,b):
    p_ab = [[0]*len(b[0]) for j in range(len(a))]
    
    for i in range(0,len(a)):
        for j in range(0,len(b[0])):
            for k in range(0,len(b)):
                p_ab[i][j] = p_ab[i][j] + a[i][k]*b[k][j]
    return p_ab

def rnum(a,c,m,r):

    x = ((a*r+c)%m)
        
    return x

def rnumlt1(a,c,m,r):

    x = ((a*r+c)%m)/m
        
    return x
    
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

def matrixtranspose(l):
    m,n = len(l),len(l)
    lt=[[0] * n for i in range(m)]
    
    for i in range(0,len(l)):
        for j in range(0,len(l)):
            lt[i][j]=l[j][i]
    
    return lt

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

def Bracket(a,b,func):
    
    while (func(a)*func(b)>=0):
        if (abs(func(a)) < abs(func(b))):
            a = a - 1.5*(b-a)
            print(a,b)
            
            
        if (abs(func(a)) > abs(func(b))):
            b = b + 1.5*(b-a)
            print(a,b)
    
    return (a,b)

            

def Bisect(a,b,e,d,func):
    c=0
    t=0
    print(t,">>",a,b,c)
    
        
    while abs(a-b)>e:
        while abs(func(a))>d:
            c=(a+b)/2
            
            
            if abs(func(c))<=d:
                print(t+1,">>",a,b,c)
                print("the solution is",'%10.4E' %c,"in",t+1,"iterations.")
                return c
            else:  
                if (func(a)*func(c))<0:
                    b=c
                    t=t+1
                    print(t,">>",a,b,c)
                    
                else:
                    a=c
                    t=t+1
                    print(t,">>",a,b,c)
# In[ ]:

def RegulaFalsi(a,b,e,d,fn):
    c0=0
    c1=1
    t=0
    print(t,">>",a,b,c0,c1)
    
    while abs(c1-c0)>e:
        while abs(fn(a))>d:
            c1 = b-((b-a)*fn(b))/(fn(b)-fn(a))
            
            if abs(fn(c1))<=d:
                print(t+1,">>",a,b,c0,c1)
                print("the solution is",c1,"in",t+1,"iterations.")
                return c1
            else:
                if (fn(a)*fn(c1))<0:
                    b=c1
                    c0=c1
                    t=t+1
                    print(t,">>",a,b,c0,c1)
                    
                else:
                    a=c1
                    c0=c1
                    t=t+1
                    print(t,">>",a,b,c0,c1)

# In[ ]:

def NewtonRaphson(x,e,d,fn,fd):
    x0=x-1
    x1=x
    t=0
    #print(t,">>",'prev',x0,'next',x1)
    
    while abs(x1-x0)>e:
        while abs(fn(x1))>d:
            x1 = x1-(fn(x1)/fd(x1))
            
            if abs(fn(x1))<=d:
                #print(t+1,">>",'prev',x0,'next',x1)
                print("the solution for x is",x1,"in",t+1,"iterations.")
                return x1
            else:
                t=t+1
                #print(t,">>",'prev',x0,'next',x1)
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