#!/usr/bin/env python
# coding: utf-8

# In[1]:


from sympy import *
from sympy.solvers.solveset import linear_coeffs
import sympy.physics.mechanics as me
import numpy as np

init_printing()

def symb(x, y, z = ''):
    return symbols('{0}_{1}{2}'.format(x, y, z), type = float)

def apply2Eq(fun, eq):
    return Eq(fun(eq.lhs), fun(eq.rhs))


# Given the following equations:
# $$
# M \ddot{X} + D \dot{X} + K X = F \\
# \sigma_F = C_F \frac{T}{b d_1 m_n}
# $$
# The nondimensional groups are given by:

# In[2]:


X, x, m_n, = symbols('X x m_n')
tau, omega_0, t = symbols('tau omega_0, t')
omega, Omega = symbols('omega Omega')
zeta, d, k, m = symbols('zeta d k m')
p, P, F = symbols('p P F')
mu, sigma = symbols('mu sigma')
B, b, f = symbols('B b f')
rho, V, pi = symbols('rho V pi')


# In[3]:


pi1 = Eq(X, x/m_n)

pi2 = Eq(B, b/m_n)

pi3 = Eq(tau**2, (omega_0*t)**2)
pi3 = pi3.subs(omega_0**2, k/m)
pi3 = pi3.subs(m, rho*V)
pi3 = pi3.subs(V, pi*b*m_n**2)

pi4 = Eq(F, f/(k*m_n))

pi5 = Eq(4*zeta**2, d**2/(m*k))
pi5 = pi5.subs(m, rho*V)
pi5 = pi5.subs(V, pi*b*m_n**2)

pi6 = Eq(Sigma, (sigma*omega_0*b*m_n**2)/P)
pi7 = Eq(mu, (sigma*m_n**2)/f)

# pi6 = Eq(p, P/(F*m_n*Omega))

display(pi1, pi2, pi3, pi4, pi5, pi6, pi7)


# In[4]:


PI = [pi1, pi2, pi3, pi4, pi5, pi6, pi71, pi8] #, pi71]
PIrhs = list(map(lambda x: x.rhs, PI))
#list(map(, PIrhs))
#display(PIrhs[:2])


# In[5]:


log_simp = lambda x: expand_log(log(x), force = True)

PI = [pi1, pi2, pi3, pi4, pi5, pi6, pi71, pi8] #, pi71]
PI = list(map(lambda x: apply2Eq(log_simp,x), PI))
PI = list(map(lambda x: x.subs([(t, 1), # no way to scale time
                                (rho, 1), # same material
                                (log(2), 0), # get rid of this
                                (sigma, 1)]), # keep structural safety
                          #      (m_n, b)]), # no distortion for gears
                PI))

for q in PI:
    display(q)


# In[6]:


#rel1 = Eq(log(m_n), log(b))
#rel2 = Eq(log(x), log(b))
#rel3 = Eq(log(F), 4*log(b))
#PI2 = list(map(lambda x: x.subs([(rel1.lhs, rel1.rhs), (rel2.lhs, rel2.rhs)]),PI))
#for q in PI2:
#    display(q)
def coeff_matrix(equations, *symbols):
    # adapted from the source code of sympy.linear_eq_to_matrix(equations, *symbols)
    A , b = [], []
    for i, f in enumerate(equations):
        coeff_list = linear_coeffs(f, *symbols)
        b.append(-coeff_list.pop())
        A.append(coeff_list)

    A, b = map(Matrix, (A, b))

    return A, b


# In[7]:


PI3 = PI
#del PI3[0]
#del PI3[6]
vv = [m_n, b, k, d, x]
#vv = list(map(lambda x: (log(x), x),vv))

for v in vv:
    PI3 = list(map(lambda x: x.subs(log(v), v), PI3))
PI4 = list(map(lambda x: x.rhs - x.lhs, PI3))

A, bb = coeff_matrix(PI4, log(P),
                                  log(Omega),
                                  log(F),
                                  log(B),
                                  log(zeta),
                                  log(X),
                                  log(tau),
                                  log(f),
                                  log(omega),
                                  log(Sigma),
                                  log(p),
                                  m_n,
                                  b,
                                  k,
                                  d,
                                  x)

display(A, bb)
for q in PI4:
    display(q)

#list(map(lambda x: x.subs(),PI3))


# In[8]:


PI02 = Eq(PI3[0].lhs + PI3[1].lhs, PI3[0].rhs + PI3[1].rhs)
display(PI3[:4])


# In[9]:


type(A)


# In[10]:


get_ipython().system('jupyter nbconvert --to script DA.ipynb')


# In[11]:


printing.octave.octave_code(A)


# In[8]:


X[m_n]


# In[9]:


gamma, Gamma = IndexedBase('gamma Gamma')


# In[ ]:




