#!/usr/bin/env /usr/bin/python3.6

import numpy as np 
from sympy import *
import math
# from evaluator import *
def likelihood_given_mu(signal,bkg,bkg_err,events,event,mu,corr):
    corr_mat=np.flipud(corr)
    #########Getting covariance matirx from correlation matrix############

    std_dev_mat=np.diag(bkg_err)
    covar_mat_temp=np.matmul(std_dev_mat,corr_mat)
    covar_mat=np.matmul(covar_mat_temp,std_dev_mat)
    covar_mat=Matrix(covar_mat)
    covar_mat_inv=covar_mat**(-1)
    ###defining matrix
    theta_list=[]
    for i in range(len(bkg)):
        theta_list.append(Symbol('theta_'+str(i)))
    theta=Matrix(theta_list)
    #generating F[\theta]
    uu=symbols('\mu')
    kk=(Matrix(uu*signal+bkg)).T+theta
    kk=diag(*kk)
    # F[\theta]
    F_theta=kk*(covar_mat_inv*theta+Matrix([1 for i in range(len(bkg))]))-(Matrix(events)).T
    #Jacobian for F[\theta]
    JF_theta=F_theta.jacobian(theta)
    # initial guess for solution
    guess=Matrix([0 for i in range(len(bkg))])
    for ii in range(10):
         s1=JF_theta.subs(uu,mu)
         s2=s1.subs(zip(theta,guess))
         j1=F_theta.subs(uu,mu)
         j2=j1.subs(zip(theta,guess))
         nn=(s2**-1)*j2
         guess=guess-nn

    ###calculating likelihood
    oo=theta.subs(zip(theta,guess)) 
    pp=(0.5)*oo.T*covar_mat_inv*oo
    sum1=0
    for jj in range(len(bkg)):
            qq=(mu)*signal[jj]+bkg[jj]+guess[jj]
            rr=events[jj]*math.log(qq)-qq-math.log(math.factorial(event[jj]))
            sum1=sum1+rr

    

    return (sum1-pp[0])








