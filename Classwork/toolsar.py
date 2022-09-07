#!/usr/bin/env python
# coding: utf-8

# In[1]:


import math
import toolsar

def rnum(a,c,m,r):

    x = ((a*r+c)%m)
        
    return x

def rnumlt1(a,c,m,r):

    x = ((a*r+c)%m)/m
        
    return x
    



# In[ ]:

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
                
                if(content[j][j] > 0):
                    content[i][j] = int(content[i][j] - s) / content[j][j]
                    
                
        
            
            
            if j > i:
                
                content[i][j] = 0
             
                
                
                
    print('\n The decomposed matrix is: \n')
    for row in content:
        print(row)
    
    contenttranspose = [[content[j][i] for j in range(len(content))] for i in range(len(content[0]))]
    print('\n The decomposed matrix transpose is: \n')
    for row in contenttranspose:
        print(row)


# In[ ]:

def Jacobi(m):
    
    print('\n The matrix to be decomposed is: \n')
    for row in m:
        print(row)
    
    n=len(m)
    
    arrx=[]
    
    r=10
    
    for i in range(0,n):
        x= toolsar.rnumlt1(1103515245,12345,32768,r)
        r=x
        arrx.append(x)
        
    print(arrx)
    
    for i in range(0,n):
        
        sum_ax=0
        
            
        for t in range(0,n):
            for j in range(0,n):
                if j != t:
                    sum_ax = sum_ax + m[t][j]*arrx[t]
                
            arrx[t] = (1/m[t][t])((m[t][n])- sum_ax)
            

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

def gsm(a,x,b):
                        

    for j in range(0, len(a)):               
          
        d=b[j]
        
        for i in range(0, len(a)):     
            if(j != i):
                d = d - a[j][i] * x[i]
        
        x[j] = d / a[j][j]
        
    return x

def GaussSeidel(a,x,b,it,tol):
    
    y=[]
    
    for i in range(0,len(x)):
        y.append(x[i])
    
    for i in range(0,it):
        t=0
        for t in range(0,len(x)):
            y[t] = x[t]
    
        x = toolsar.gsm(a,x,b)
    
        print(x)
    
        z = toolsar.arrclose(x,y,tol)
        if z == True:
            print('The answer is:',x)
            break
        else:
            continue