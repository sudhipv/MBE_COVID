
# Sudhi P V
# Civil and Environmental Engineering
# Carleton University

# INSPIRED FROM "https://github.com/prashjha/BayesForSEIRD"
######### SEIRD -One dimension #############


from dolfin import *
import numpy as np
import os as os
import math
import matplotlib.pyplot as plt


## For avoiding FEniCS internal node reordering optimized for speed
parameters['reorder_dofs_serial'] = False

SU = 0
EX = 1
IN = 2
RE = 3
DE = 4

#### Creating mesh
mesh = UnitIntervalMesh(2000)

coordinates = mesh.coordinates()

nv = mesh.num_vertices()

print('number of vertices', nv)

print(coordinates)

np.savetxt('./seird/nodes.txt',coordinates)

### Discretized Function space on the mesh

V = FunctionSpace(mesh, "Lagrange", 1)

Ve = FunctionSpace(mesh, "Lagrange", 3)


u_trial = TrialFunction(V)

v = TestFunction(V)

### Previous time step
u_0 = [Function(V) for i in range(5)]

u_Ve = [Function(Ve) for i in range(5)]


ue_V = [Function(V) for i in range(5)]


### Previous Iteration in Picard Loop
u_k = [Function(V) for i in range(5)]


u_k_1 = [Function(V) for i in range(5)]

## Current Iteration
u = [Function(V) for i in range(5)]


a = [None]*5
L = [None]*5


A = 0
beta_i = 0.01
beta_e = 0.01
nu_s = 5e-5
nu_e = 1e-3
nu_r = 5e-5
nu_i = 1e-10
gamma_r = 1/24
gamma_d = 1/160
gamma_e = 1/6
sigma=1/8


print(gamma_r)
print(gamma_d)
print(sigma)


T = 200
dt = 0.1

nt = round(T/dt)

print("Total timesteps is", nt)


### Initial Conditions

k = 10
c = 0.2

u_Truth = [Function(V) for i in range(5)]


u_Truth[0] = Expression(("500 + 25*sin(k*x[0] + c*t)"),t=0,k = k, c = c, degree = 3)

u_Truth[1] = Expression(("300 + 25*sin(k*x[0] + c*t)"),t=0,k = k, c = c, degree = 3)

u_Truth[2] = Expression(("200 + 25*sin(k*x[0] + c*t)"),t=0,k = k, c = c, degree = 3)

u_Truth[3] = Expression(("100 + 25*sin(k*x[0] + c*t)"),t=0,k = k, c = c, degree = 3)

u_Truth[4] = Expression(("80 + 25*sin(k*x[0] + c*t)"),t=0,k = k, c = c, degree = 3)


u_0[SU] = interpolate(u_Truth[0], V)
u_0[EX] = interpolate(u_Truth[1], V)
u_0[IN] = interpolate(u_Truth[2], V);
u_0[RE] = interpolate(u_Truth[3], V);
u_0[DE] = interpolate(u_Truth[4], V);


np.savetxt('./seird/sus_ini.txt', u_0[SU].vector().array())
np.savetxt('./seird/exp_ini.txt', u_0[EX].vector().array())

np.savetxt('./seird/inf_ini.txt', u_0[IN].vector().array())
np.savetxt('./seird/rec_ini.txt', u_0[RE].vector().array())
np.savetxt('./seird/dec_ini.txt', u_0[DE].vector().array())

############## Boundary Condition ####################

tolbc = 1e-14

def left(x, on_boundary):

    return on_boundary and abs(x[0]) < tolbc


def right(x, on_boundary):

    return on_boundary and abs(x[0]-1.00) < tolbc



for i in range(5):

    Gamma0 = DirichletBC(V, u_Truth[i], left)

    Gamma1 = DirichletBC(V, u_Truth[i], right)

    if(i == 0):
        bcsL = ([Gamma0])
        bcsR = ([Gamma1])

    else:
        bcsL.append(Gamma0)
        bcsR.append(Gamma1)


############ #########################################


u_force = [Function(V) for i in range(5)]

u_force[0] = Expression(("625*beta_e*(sin(0.2*t + 10*x[0]) + 12)*(sin(0.2*t + 10*x[0]) + 20) + 625*beta_i*(sin(0.2*t + 10*x[0]) + 8)*(sin(0.2*t + 10*x[0]) + 20) +\
 12500*nu_s*(25*sin(0.2*t + 10*x[0]) + 236)*sin(0.2*t + 10*x[0]) - 312500*nu_s*pow(cos(0.2*t + 10*x[0]),2) + 5.0*cos(0.2*t + 10*x[0])"),\
t=0, beta_i = beta_i, beta_e = beta_e, nu_s = nu_s,degree = 3)


u_force[1] = Expression(("-625*beta_e*(sin(0.2*t + 10*x[0]) + 12)*(sin(0.2*t + 10*x[0]) + 20) - 625*beta_i*(sin(0.2*t + 10*x[0]) + 8)*(sin(0.2*t + 10*x[0]) + 20) +\
 25*gamma_e*(sin(0.2*t + 10*x[0]) + 12) + 25*sigma*(sin(0.2*t + 10*x[0]) + 12) + 12500*nu_e*(25*sin(0.2*t + 10*x[0]) + 236)*sin(0.2*t + 10*x[0]) -\
  312500*nu_e*pow(cos(0.2*t + 10*x[0]),2) + 5.0*cos(0.2*t + 10*x[0])"),t = 0,beta_i = beta_i, beta_e = beta_e, nu_e = nu_e,sigma= sigma, gamma_e = gamma_e, degree=3)


u_force[2] = Expression(("25*gamma_d*(sin(0.2*t + 10*x[0]) + 8) + 25*gamma_r*(sin(0.2*t + 10*x[0]) + 8) - 25*sigma*(sin(0.2*t + 10*x[0]) + 12) + 12500*nu_i*(25*sin(0.2*t + 10*x[0]) +\
 236)*sin(0.2*t + 10*x[0]) - 312500*nu_i*pow(cos(0.2*t + 10*x[0]),2) + 5.0*cos(0.2*t + 10*x[0])"),t = 0,nu_i = nu_i,sigma= sigma, gamma_d = gamma_d,gamma_r = gamma_r, degree=3)


u_force[3] = Expression(("-25*gamma_e*(sin(0.2*t + 10*x[0]) + 12) - 25*gamma_r*(sin(0.2*t + 10*x[0]) + 8) + 12500*nu_r*(25*sin(0.2*t + 10*x[0]) + 236)*sin(0.2*t + 10*x[0]) -\
 312500*nu_r*pow(cos(0.2*t + 10*x[0]),2) + 5.0*cos(0.2*t + 10*x[0])"),t = 0,nu_r = nu_r,gamma_r = gamma_r,gamma_e = gamma_e, degree=3)


u_force[4] = Expression(("-25*gamma_d*(sin(0.2*t + 10*x[0]) + 8) + 5.0*cos(0.2*t + 10*x[0])"),t = 0,gamma_d = gamma_d, degree=3)

#####################################################


def weakform(u_k,u_0,u_force):

    F = [None]*5

    un_k = u_k[0] + u_k[1] + u_k[2] + u_k[3] + u_k[4]

    # Define weak form for phi_s evolution
    f_s = - (1.- A/un_k)*beta_i*u_trial*u_k[IN] - (1.- A/un_k)*beta_e*u_trial*u_k[EX] + u_force[0]


    F[SU]     = (u_trial-u_0[SU])*v*dx - dt*f_s*v*dx +\
                        dt*(un_k*nu_s*inner(grad(u_trial), grad(v)))*dx

    # Define the weak form for phi_e evolution
    f_e = - sigma*u_trial - gamma_e*u_trial + (1.- A/un_k)*beta_i*u_k[SU]*u_k[IN] + (1.- A/un_k)*beta_e*u_k[SU]*u_trial + u_force[1]
    F[EX]     = (u_trial-u_0[EX])*v*dx - dt*f_e*v*dx +\
                        dt*(un_k*nu_e*inner(grad(u_trial), grad(v)))*dx

    # Define the weak form for phi_i evolution
    f_i = - gamma_d*u_trial - gamma_r*u_trial + sigma*u_k[EX] + u_force[2]
    F[IN]     = (u_trial-u_0[IN])*v*dx - dt*f_i*v*dx +\
                        dt*(un_k*nu_i*inner(grad(u_trial), grad(v)))*dx

    # Define the weak form for phi_r evolution
    f_r = gamma_r*u_k[IN] + gamma_e*u_k[EX] + u_force[3]
    F[RE]     = (u_trial-u_0[RE])*v*dx - dt*f_r*v*dx +\
                        dt*(un_k*nu_r*inner(grad(u_trial), grad(v)))*dx

    # Define the weak form for phi_d evolution
    F[DE]     = (u_trial-u_0[DE])*v*dx - dt*gamma_d*u_k[IN]*v*dx - dt*u_force[4]*v*dx


    #### Bilinear and Linear terms assembly

    for i in range(5):
        a[i], L[i] = lhs(F[i]), rhs(F[i])

    return(a,L)


#### Solution #######


H = [None]*5
b = [None]*5


M = assemble(u_trial*v*dx)

z = assemble(Constant(1.0)*v*dx)


out_values = np.empty((nt, 5))
out_space = np.empty((nv,5))

out_truth = np.empty((nt, 5))
out_space_truth = np.empty((nv, 5))

out_error = np.empty((nt,1))
out_error_space = np.empty((nv,1))


u_k = u_0

t = dt

tolP = 1e-12
g = 0

### Time loop
for i in range(nt):

    print("time is", t)

    u_Truth[0].t = t
    u_Truth[1].t = t
    u_Truth[2].t = t
    u_Truth[3].t = t
    u_Truth[4].t = t

    u_force[0].t = t
    u_force[1].t = t
    u_force[2].t = t
    u_force[3].t = t
    u_force[4].t = t

#### Picard Loop
    for k in range(25):

        a,L = weakform(u_k,u_0,u_force)

        for j in range(5):
            H[j] = assemble(a[j])
            b[j] = assemble(L[j])

            bcs = [bcsL[j], bcsR[j]]

            for bc in bcs:
                bc.apply(H[j],b[j])

            solve(H[j], u[j].vector(), b[j],'lu')

        error = 0.0
        diff = 0.0
        for q in range(5):
            diff = u[q].vector().copy()
            diff = diff - u_k[q].vector()

            error += np.linalg.norm(diff, ord = 2)/np.linalg.norm(u_0[q].vector(), ord = 2)

        u_k = u


        print('error in iteration'+str(k), error)


        if(error < tolP):
            print('exit Picards loop after '+str(k)+'iterations')
            break

        if(k > 25):
            print("Exiting after 25 iterations")
            break

    # Finished Picard Loop

### Saving the new time step solution (picard final solution) as previous time step solution
    u_0 = u_k

    errl2 = 0
    errl2fun = 0
    for g in range(5):
        u_Ve[g] = interpolate(u_Truth[g], Ve)
        ue_V[g] = interpolate(u_Truth[g], V)

        diff2 = ue_V[g].vector().copy()
        diff2 = diff2 - u_k[g].vector()


        errl2 += np.linalg.norm(diff2, ord = 2)/np.linalg.norm(ue_V[g].vector(), ord = 2)

    # print('L2 Norm of error in Higher FE space:', errl2fun)
    print('L2 Norm of error:', errl2)


    for l in range (5):
        out_values[i, l] = u_k[l].vector().inner(z)      #### z is assembled test function
        out_truth[i,l]  = ue_V[l].vector().inner(z)


    if(i == 100):
        for h in range(5):
            out_space[:,h] = u_k[h].vector()      #### z is assembled test function
            out_space_truth[:,h]  = ue_V[h].vector()


            out_error_space[:,0] += (out_space_truth[:,h] - out_space[:,h])/max(out_space_truth[:,h])

        print(out_error_space[:,0])

    t += dt


################################################################

np.savetxt('./seird/numerical.txt',out_values)
np.savetxt('./seird/analytical.txt',out_truth)



np.savetxt('./seird/numerical_space.txt',out_space)
np.savetxt('./seird/analytical_space.txt',out_space_truth)


np.savetxt('./seird/errl2.txt',out_error)

np.savetxt('./seird/err_space.txt',out_error_space)








