__author__ = 'shahar'

import numpy as np
from matplotlib import pyplot as plt

def tau(R1,R2,C1,C2):
    return (R1*R2*(C1+C2))/(R1+R2)

def A(R1,R2,C1,C2,alpha, V0):
    def a(Q1,t):
        return (Q1-((alpha*R1*C1)/(R1+R2))*((R2*(R1*C1-R2*C2)/(R1+R2))+t) -
                V0*(R1*C1)/(R1+R2))*np.exp(t/tau(R1,R2,C1,C2))
    return a

def B(R1,C1,C2,alpha):
    def b(Q1,t):
        return (Q1-alpha*R1*C1*C2)*np.exp(t/(R1*(C1+C2)))
    return b

def C(C1,C2,alpha, V0):
    def c(Q,t):
        return Q-(alpha*t+V0)*((C1*C2)/(C1+C2))
    return c

def Qth(Vth, C):
    return Vth*C

def Q1BothLow(C1,C2,alpha, V0):
    def q1BothLow(t,C):
        return (alpha*t+V0)*((C1*C2)/(C1+C2)) + C
    return q1BothLow

def Q2BothLow(C1,C2,alpha, V0):
    def q2BothLow(t,C):
        return (alpha*t + V0)*((C1*C2)/(C1+C2)) - C*(C2/C1)
    return q2BothLow

def Q1OneHigh(C1,C2,R1,alpha):
    def q1OneHigh(t,B):
        return B*np.exp(-t/(R1*(C1+C2)))+alpha*C1*R1*C2
    return q1OneHigh

def Q2OneHigh(C1,C2,R1,alpha, V0):
    def q2OneHigh(t,B):
        return (-C2/C1)*B*np.exp(-t/(R1*(C1+C2)))+alpha*C2*(t-R1*C2) + V0*C2
    return q2OneHigh


def Q1BothHigh(C1,C2,R1,R2,alpha, V0):
    def q1BothHigh(t,A):
        return A*np.exp(-t/tau(R1,R2,C1,C2))+((alpha*C1*R1)/(R1+R2))*(
            R2*(C1*R1-C2*R2)/(R1+R2)+t) + V0*(R1*C1)/(R1+R2)
    return q1BothHigh

def Q2BothHigh(C1,C2,R1,R2,alpha, V0):
    def q2BothHigh(t,A):
        return (-C2/C1)*A*np.exp(-t/tau(R1,R2,C1,C2))+((alpha*C2*R2)/(
            R1+R2))*((R1*(R2*C2-R1*C1)/(R1+R2))+t)+V0*(R2*C2)/(R1+R2)
    return q2BothHigh

def calcQ(tVec, C1, C2, R1, R2, alpha, Vth1, Vth2, Q1Init, Q2Init, V0):
    q11=Q1BothLow(C1,C2,alpha,V0)
    q12=Q1OneHigh(C1,C2,R1,alpha)
    q13=Q2OneHigh(C1,C2,R2,alpha, V0)
    q14=Q1BothHigh(C1,C2,R1,R2,alpha, V0)
    q21=Q2BothLow(C1,C2,alpha,V0)
    q22=Q2OneHigh(C1,C2,R1,alpha, V0)
    q23=Q1OneHigh(C1,C2,R2,alpha)
    q24=Q2BothHigh(C1,C2,R1,R2,alpha, V0)
    a=A(R1,R2,C1,C2,alpha, V0)
    b1=B(R1,C1,C2,alpha)
    b2=B(R2,C2,C1,alpha)
    c=C(C1,C2,alpha, V0)
    Qth1 = Qth(Vth1, C1)
    Qth2 = Qth(Vth2, C2)
    t0 = tVec[0]
    Q1 = [Q1Init]
    Q2 = [Q2Init]
    CurrA = a(Q1Init, t0)
    CurrB1 = b1(Q1Init, t0)
    CurrB2 = b2(Q2Init, t0)
    CurrC = c(Q1Init, t0)
    conduct1 = Q1Init >= Qth1
    conduct2 = Q2Init >= Qth2
    for t in tVec:
        if Q2[-1] < Qth2:
            if Q1[-1] < Qth1:
                if conduct1 or conduct2:
                    CurrC = c(Q1[-1],t)
                Q1.append(q11(t,CurrC))
                Q2.append(q21(t,CurrC))
                conduct1 = False
                conduct2 = False
            else:
                if not(conduct1 and not conduct2):
                    CurrB1 = b1(Q1[-1],t)
                Q1.append(q12(t,CurrB1))
                Q2.append(q22(t,CurrB1))
                conduct1 = True
                conduct2 = False
        else:
            if Q1[-1] < Qth1:
                if not(conduct2 and not conduct1):
                    CurrB2 = b2(Q2[-1],t)
                Q1.append(q13(t,CurrB2))
                Q2.append(q23(t,CurrB2))
                conduct1 = False
                conduct2 = True
            else:
                if not(conduct1 and  conduct2):
                    CurrA = a(Q1[-1],t)
                Q1.append(q14(t,CurrA))
                Q2.append(q24(t,CurrA))
                conduct1 = True
                conduct2 = True
    return Q1[1:], Q2[1:]

def calcI(Q1, R1, C1, V1th):
    Q1 = np.array(Q1)
    Q1th = V1th*C1
    Q1dot = Q1[1:]-Q1[:-1]
    Q1 = Q1[1:]
    I = -Q1dot
    I[Q1 >= Q1th] += Q1[Q1 >= Q1th] / (R1*C1)
    return I

def calcV(alpha, t):
    V = alpha*t
    tswitch = t[len(t)//2]
    Vswitch = V[len(V)//2]
    V[(len(V)//2):] = Vswitch - alpha*(t[(len(V)//2):] -tswitch)
    return V

if __name__ == "__main__":
    C1=1
    C2=1
    R1 = 1
    R2 = 1
    Vth1 = 5
    Vth2 = 10
    alpha = 1
    Q1Init = 0
    Q2Init = 0
    tVecUp = np.arange(0,30,0.001)
    tVecDown = np.arange(30,60,0.001)



    Q1Up,Q2Up = calcQ(tVecUp, C1, C2, R1, R2, alpha, Vth1, Vth2, Q1Init,
                      Q2Init, 0)
    Q1Init = Q1Up[-1]
    Q2Init = Q2Up[-1]
    Q1Down,Q2Down = calcQ(tVecUp, C1, C2, R1, R2, -alpha, Vth1, Vth2,
                      Q1Init, Q2Init, alpha*30)
    Q1 = np.array(Q1Up+ Q1Down)
    Q2 = np.array(Q2Up+Q2Down)
    tVec = np.hstack((tVecUp,tVecDown))
    plt.figure()
    plt.plot(tVec, Q1, '.b', tVec, Q2, '.r')

    I1 = calcI(Q1, R1, C1, Vth1)
    I2 = calcI(Q2, R2, C2, Vth2)
    V = calcV(alpha, tVec)
    plt.figure()
    plt.plot(V[1:len(V)//2], I1[:len(V)//2-1], '.b', V[len(V)//2+1:], I1[len(
        V)//2:],'.r')
    plt.figure()
    plt.plot(V[1:len(V)//2], I2[:len(V)//2-1], '.b', V[len(V)//2+1:], I2[len(
        V)//2:],'.r')
    plt.show()













