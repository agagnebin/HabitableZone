#!/usr/bin/env python
# coding: utf-8

# In[33]:


import numpy as np
import matplotlib.pyplot as plt
import os


# In[2]:


pwd = os.getcwd()


# In[3]:


x = np.loadtxt("Final_Positions.txt", usecols = 0, skiprows=1)
y = np.loadtxt("Final_Positions.txt", usecols = 1, skiprows=1)


# In[5]:


specclass = input("Please input the same spectral class (capitalized): ")


# In[29]:


os.mkdir("Plots")
os.chdir("{0}/Plots".format(pwd))


# In[45]:


plt.style.use('dark_background')
if specclass == "G":    
    for i in range(1000):
        plt.figure()
        plt.plot(figsize = (8, 6), dpi = 80)
        plt.title("2D Orbit of an Earth-like Planet around a "+specclass+"-type star")
        plt.plot(x, y, color = "snow", linewidth = 0.3)
        plt.scatter(x[i], y[i], color = "aqua", s = 4, label = "Planet")
        plt.scatter(0, 0, color = "lightyellow", s = 100, label = "Star")
        plt.legend()
        plt.xlabel("Distance in x (AU)")
        plt.ylabel("Distance in y (AU)")
        plt.savefig("fig%04d.png" %i)
        
if specclass == "O":    
    for i in range(1000):
        plt.figure()
        plt.plot(figsize = (8, 6), dpi = 80)
        plt.title("2D Orbit of an Earth-like Planet around an "+specclass+"-type star")
        plt.plot(x, y, color = "snow", linewidth = 0.3)
        plt.scatter(x[i], y[i], color = "aqua", s = 4, label = "Planet")
        plt.scatter(0, 0, color = "cornflowerblue", s = 100, label = "Star")
        plt.legend()
        plt.xlabel("Distance in x (AU)")
        plt.ylabel("Distance in y (AU)")
        plt.savefig("fig%04d.png" %i)
        
if specclass == "B":    
    for i in range(1000):
        plt.figure()
        plt.plot(figsize = (8, 6), dpi = 80)
        plt.title("2D Orbit of an Earth-like Planet around a "+specclass+"-type star")
        plt.plot(x, y, color = "snow", linewidth = 0.3)
        plt.scatter(x[i], y[i], color = "aqua", s = 4, label = "Planet")
        plt.scatter(0, 0, color = "dodgerblue", s = 100, label = "Star")
        plt.legend()
        plt.xlabel("Distance in x (AU)")
        plt.ylabel("Distance in y (AU)")
        plt.savefig("fig%04d.png" %i)
        
if specclass == "A":    
    for i in range(1000):
        plt.figure()
        plt.plot(figsize = (8, 6), dpi = 80)
        plt.title("2D Orbit of an Earth-like Planet around an "+specclass+"-type star")
        plt.plot(x, y, color = "snow", linewidth = 0.3)
        plt.scatter(x[i], y[i], color = "aqua", s = 4, label = "Planet")
        plt.scatter(0, 0, color = "deepskyblue", s = 100, label = "Star")
        plt.legend()
        plt.xlabel("Distance in x (AU)")
        plt.ylabel("Distance in y (AU)")
        plt.savefig("fig%04d.png" %i)
        
if specclass == "F":    
    for i in range(1000):
        plt.figure()
        plt.plot(figsize = (8, 6), dpi = 80)
        plt.title("2D Orbit of an Earth-like Planet around a "+specclass+"-type star")
        plt.plot(x, y, color = "snow", linewidth = 0.3)
        plt.scatter(x[i], y[i], color = "aqua", s = 4, label = "Planet")
        plt.scatter(0, 0, color = "aliceblue", s = 100, label = "Star")
        plt.legend()
        plt.xlabel("Distance in x (AU)")
        plt.ylabel("Distance in y (AU)")
        plt.savefig("fig%04d.png" %i)
        
if specclass == "K":    
    for i in range(1000):
        plt.figure()
        plt.plot(figsize = (8, 6), dpi = 80)
        plt.title("2D Orbit of an Earth-like Planet around a "+specclass+"-type star")
        plt.plot(x, y, color = "snow", linewidth = 0.3)
        plt.scatter(x[i], y[i], color = "aqua", s = 4, label = "Planet")
        plt.scatter(0, 0, color = "coral", s = 100, label = "Star")
        plt.legend()
        plt.xlabel("Distance in x (AU)")
        plt.ylabel("Distance in y (AU)")
        plt.savefig("fig%04d.png" %i)
        
if specclass == "M":    
    for i in range(1000):
        plt.figure()
        plt.plot(figsize = (8, 6), dpi = 80)
        plt.title("2D Orbit of an Earth-like Planet around a "+specclass+"-type star")
        plt.plot(x, y, color = "snow", linewidth = 0.3)
        plt.scatter(x[i], y[i], color = "aqua", s = 4, label = "Planet")
        plt.scatter(0, 0, color = "red", s = 100, label = "Star")
        plt.legend()
        plt.xlabel("Distance in x (AU)")
        plt.ylabel("Distance in y (AU)")
        plt.savefig("fig%04d.png" %i)


# In[ ]:


os.chdir("{0}".format(pwd))

