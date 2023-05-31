#!/usr/bin/env python
# coding: utf-8

# In[1]:


import diag_helpers as helpers
import itertools as it
from diag_optimization import weight
import numpy as np
import galois
GF = galois.GF(2)

import pstring
import tableau_operations as tabop



# In[2]:


def pstring_to_vector(mypstr):
    n = len(mypstr)
    res = np.zeros(2*n, dtype=int)
    for i in range(n):
        if mypstr[i] == "1":
            res[i] = 1
        if mypstr[i] == "2":
            res[i] = 1
            res[i + n] = 1
        if mypstr[i] == "3":
            res[i+n] = 1
    return GF(res)


# In[3]:


def pcommute(p1, p2):

    return p1.commute(p2)


# In[4]:


def commute_with_others(p, others):
    for other in others:
        if not p.commute(other):
            return False

    return True


# In[5]:


def insertpmatrices(n, qbits, pmatrices):

    pstr = ['0']*n

    for q in range(len(qbits)):
        pstr[qbits[q]] = str(pmatrices[q])

    return "".join(pstr)


# In[7]:


def candidate_weight_w(x,z,wght):

    n = x.shape[1]
    T = helpers.tableau_to_pstrings(z,x)

    qbitsList = it.combinations(range(n), r=wght)


    for qbits in qbitsList:

        pmatricesList = it.product(range(1,4), repeat=wght)

        for pmatrices in pmatricesList:

            p = pstring.pstring(insertpmatrices(n, qbits, pmatrices),1)

            if commute_with_others(p, T) and p.weight() != 0:

                return GF(p.to_vector())


# In[8]:


def search_for_vector(x,z):

    r, n = x.shape

    for wght in range(r//2 + 2):

        cand = candidate_weight_w(x,z,wght)

        if not cand is None:

            return cand

    print("we should never get here")

    return


# In[ ]:
