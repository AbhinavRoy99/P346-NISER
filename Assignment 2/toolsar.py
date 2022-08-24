#!/usr/bin/env python
# coding: utf-8

# In[1]:


def rnum(content):

    if len(content) == 5:
    
        n=int(content[0])
        a=int(content[1])
        c=int(content[2])
        m=int(content[3])
        r=int(content[3])
        
        rnumarr=[]
        
        for i in range(0,n):
            r = (a*r+c)%m
            rnumarr.append(r)
            
        
        if n<= 100:
            print("The random number array is saved as: \n",rnumarr)
        else:
            print("The random number array is too large to print but has been saved as instructed.")
        
        return rnumarr
        
    else:
        print("Give all parameters please.")
        

def rnumlt1(content):

    if len(content) == 5:
    
        n=int(content[0])
        a=int(content[1])
        c=int(content[2])
        m=int(content[3])
        r=int(content[3])
        
        rnumarr=[]
        
        for i in range(0,n):
            r = ((a*r+c)%m)/m
            rnumarr.append(r)
            
        if n<= 100:
            print("The random number array is saved as: \n",rnumarr)
        else:
            print("The random number array is too large to print but has been saved as instructed.")
        
        return rnumarr
        
    else:
        print("Give all parameters please.")

