#!/usr/bin/env python
# coding: utf-8

# In[1]:

def rnum(a,c,m,r):

    x = ((a*r+c)%m)
        
    return x

def rnumlt1(a,c,m,r):

    x = ((a*r+c)%m)/m
        
    return x