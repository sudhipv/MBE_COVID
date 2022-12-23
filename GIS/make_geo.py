#!/usr/bin/env python
# coding: utf-8

# In[1]:


# dependencies
import os
import sys
import copy
import numpy as np
import pandas as pd
import geopandas as gpd
import pyvista as pv
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib as mpl
import seaborn as sns

sns.set()
sns.set_style("dark")
# sns.set_style("whitegrid")
sns.set_context("paper", font_scale=4, rc={"lines.linewidth": 4})
mpl.rcParams['lines.linewidth'] = 3



# # Load Ontario map with outer boundary

# In[3]:


ON_shp = gpd.read_file("./shp/ottawa/bd/boundary.shp")
print(ON_shp)



ON_geom = ON_shp['geometry']
# ON_geom.plot()
# ON_geom.head()

# plt.show()

ON_poly = ON_geom[0]
print(ON_geom)

x, y = ON_poly.exterior.coords.xy
x = np.array(x)
y = np.array(y)


# PC = np.loadtxt('./mesh/PC_coordinates.txt')

# exit(1)
# plt.plot(x,y,marker='o')
# plt.show()

# In[7]:


def create_geo_file(filename, x, y, mesh_size):
    # create .geo file for gmsh
    geof = open(filename+'.geo','w')
    geof.write("cl__1 = {};\n".format(mesh_size))
    geof.write("Mesh.MshFileVersion = 2.2;\n")

    # print points
    N = len(x)
    for i in range(N):
        geof.write('Point({}) = {{{}, {}, 0., cl__1}};\n'.format(i+1, x[i], y[i]))

    # print lines
    for i in range(N-1):
        geof.write('Line({}) = {{{},{}}};\n'.format(i+1, i+1, i+2))

    geof.write('Line({}) = {{{},{}}};\n'.format(N, N, 1))

    # define polygon surface
    geof.write('Line Loop(1) = {')
    for i in range(N-1):
        geof.write('{},'.format(i+1))
    geof.write('{}'.format(N))
    geof.write('};\n')
    geof.write('Plane Surface(1) = {1};\n')
    tag = '"' + "a" + '"'
    geof.write('Physical Surface({}) = {{-1}};\n'.format(tag))
    geof.close()



# get minimum of length of lines in polygon
def get_min_length(x,y):
    min_length = 1000.
    for i in range(len(x)-1):
        d = (x[i+1] - x[i])*(x[i+1] - x[i]) + (y[i+1] - y[i])*(y[i+1] - y[i])
        d = np.sqrt(d)
        if min_length > d:
            min_length = d

    return min_length



# create a array excluding last (x,y) which is equal to first in shapefile
# Last node and First Node are same
X = x[:-1]
Y = y[:-1]


min_length = get_min_length(X,Y)
print('Minimum length of line: {}'.format(min_length))



# create_geo_file('mesh/mesh_10h', X, Y, 10. * min_length)

create_geo_file('mesh/ottawa', X, Y, 1)



exit(1)
# In[11]:

# os.system("gmsh './mesh/mesh_10h'.geo -2")
# os.system("gmsh './mesh/mesh_10h'.geo -2 -o 'data/mesh/mesh_10h'.vtk")

# create_geo_file('data/mesh/mesh_15h', X, Y, 15. * min_length)
# get_ipython().system("gmsh 'data/mesh/mesh_15h'.geo -2")
# get_ipython().system("gmsh 'data/mesh/mesh_15h'.geo -2 -o 'data/mesh/mesh_15h'.vtk")


# # In[12]:


# create_geo_file('data/mesh/mesh_20h', X, Y, 20. * min_length)
# get_ipython().system("gmsh 'data/mesh/mesh_20h'.geo -2")
# get_ipython().system("gmsh 'data/mesh/mesh_20h'.geo -2 -o 'data/mesh/mesh_20h'.vtk")


# In[13]:


create_geo_file('data/mesh/mesh_5h', X, Y, 5. * min_length)
get_ipython().system("gmsh 'data/mesh/mesh_5h'.geo -2")
get_ipython().system("gmsh 'data/mesh/mesh_5h'.geo -2 -o 'data/mesh/mesh_5h'.vtk")


# In[14]:


# create_geo_file('data/mesh/mesh_4h', X, Y, 4. * min_length)
# get_ipython().system("gmsh 'data/mesh/mesh_4h'.geo -2")
# get_ipython().system("gmsh 'data/mesh/mesh_4h'.geo -2 -o 'data/mesh/mesh_4h'.vtk")


# # In[15]:


# create_geo_file('data/mesh/mesh_2h', X, Y, 2. * min_length)
# get_ipython().system("gmsh 'data/mesh/mesh_2h'.geo -2")
# get_ipython().system("gmsh 'data/mesh/mesh_2h'.geo -2 -o 'data/mesh/mesh_2h'.vtk")


# # In[16]:


create_geo_file('data/mesh/mesh_h', X, Y, 1. * min_length)
get_ipython().system("gmsh 'data/mesh/mesh_h'.geo -2")
get_ipython().system("gmsh 'data/mesh/mesh_h'.geo -2 -o 'data/mesh/mesh_h'.vtk")


# In[10]:


# create_geo_file('data/mesh/mesh_0.5h', X, Y, 0.5 * min_length)
# get_ipython().system("gmsh 'data/mesh/mesh_0.5h'.geo -2")
# get_ipython().system("gmsh 'data/mesh/mesh_0.5h'.geo -2 -o 'data/mesh/mesh_0.5h'.vtk")


# # find max and min range of map

# In[11]:


x_min = np.min(X)
x_max = np.max(X)

y_min = np.min(Y)
y_max = np.max(Y)

print('min: {0} max: {1}'.format([x_min, y_min], [x_max, y_max]))


# In[19]:


W = x_max - x_min
L = y_max - y_min

print('W: {}  L: {}  WL: {}'.format(W,L, W*L))


# In[20]:


texas_area_approx = 0.5*W*L
texas_area_approx


# In[22]:


texas_area_km = 696200. #km^2
texas_area_100km = texas_area_km / (100. * 100.)

print('area([km]^2): {0:<10} area([100km]^2): {1:<10} area approx: {2}'.format(texas_area_km,
                                        texas_area_100km, texas_area_approx))


# In[ ]:




