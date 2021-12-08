import sympy as sp
import numpy as np
from sympy import Matrix,shape,symbols, diff, sin, cos, tan, sinh, tanh, acos, atan, sqrt, limit, oo

import multiprocess as Multi

#for general purposes, zero array creation function
def ZeroArrayCreate(*integers):
    elemnum = np.prod(integers)
    zerolist = [0] * elemnum
    ZeroArray = sp.MutableDenseNDimArray(zerolist,integers)
    return ZeroArray
    pass

#for general purposes, zero matrix creation function
def ZeroMatrixCreate(*integers):
    ZeroMatrix = sp.zeros(integers[0],integers[1])
    return ZeroMatrix
    pass

#this part defines a function of a coordinate transform for a rank-2 tensor

def BaseTransform2DMat(xin,xout,inputMetric):
    outputMetric = ZeroMatrixCreate(len(xin),len(xout))
    for u in range(0,len(xin)):
        for v in range(0,len(xin)):
            for i in range(0,len(xin)):
                for j in range(0,len(xin)):
                    outputMetric[u,v] += diff(xin[i],xout[u])*diff(xin[j],xout[v])*inputMetric[i,j]

    return outputMetric
    pass


#this class defines various properties of relevant 2-rank tensors that we are interested in. The inverse matrix should be inputted manually.
class TwoRankTensor:
    def __init__(self, Matrix, basis):
        self.arr = Matrix
        self.basis = basis
        self.inv = Matrix**-1
        self.purebas = basis[basis[0]]

    def ComplexifyTo(self,target_basis):
        for n in range(0,len(self.basis[self.basis[0]])):
            self.arr = self.arr.subs(self.purebas[n], self.basis[target_basis[0]][n])
        self.inv = self.arr**-1

    def TransformTo(self,target_basis):
        if self.basis == target_basis:
            pass
        else:
            self.ComplexifyTo(target_basis)
            self.arr = BaseTransform2DMat(self.basis[target_basis[0]], target_basis[target_basis[0]], self.arr)
            self.basis = target_basis
            self.purebas = target_basis[target_basis[0]]
            self.inv = self.arr**-1

#a special tensor adder to hopefully prevent errors
def TensorSum(obj1, obj2, simplifyState = 'on'):
    if obj1.basis == obj2.basis:
        if strname == '+':
            TensorTotal = obj1.arr + obj2.arr
        elif strname == '-':
            TensorTotal = obj1.arr - obj2.arr
        if simplifyState == 'on':
            return TwoRankTensor(sp.simplify(TensorTotal),obj1.basis)
        elif simplifyState == 'off':
            return TwoRankTensor(TensorTotal,obj1.basis)
    else:
        print('Tensors are not of equal basis!')
    pass

#defining tensor functions
def ChristoffelSym(tensor, simplifyState = 'on'):           #Christoffel symbol function
    from sympy import diff
    repeater = 4
    paper = ZeroArrayCreate(repeater,repeater,repeater)
    for y in range(0,repeater):
        for b in range(0,repeater):
            for u in range(0,repeater):
                total = 0
                for a in range(0,repeater): #this is the dummy index
                    matderiv = diff(tensor.arr[a,b],tensor.purebas[u]) + diff(tensor.arr[a,u],tensor.purebas[b]) - diff(tensor.arr[b,u],tensor.purebas[a])
                    answer =  tensor.inv[a,y] *matderiv
                    print(f'Christoffel symbol for index {y,b,u} and dummy index {a} has been calculated with answer {answer}')
                    total += answer
                paper[y,b,u] = 0.5 * total
    if simplifyState == 'on':
        return sp.simplify(paper)
    elif simplifyState == 'off':
        return paper
    pass

def RiemannCurvature(tensor, CS, simplifyState = 'on'):     #Riemann Curvature tensor function
    from sympy import diff
    repeater = 4
    paper = ZeroArrayCreate(repeater,repeater,repeater,repeater)
    for a in range(0,repeater):
        for b in range(0,repeater):
            for u in range(0,repeater):
                for v in range(0,repeater):
                    T1 = diff(CS[a,b,v],tensor.purebas[u])
                    T2 = diff(CS[a,b,u],tensor.purebas[v])
                    SumLast = 0
                    for o in range(0,repeater): #this is the dummy index
                        T3 = CS[a,o,u] * CS[o,b,v]
                        T4 = CS[a,o,v] * CS[o,b,u]
                        sum = T3 - T4
                        SumLast += sum
                        print(f'Partial Riemann curvature tensor component for index {a,b,u,v} and dummy index {o} has been calculated with answer {SumLast}')

                    paper[a,b,u,v] = T1 - T2 + SumLast
    if simplifyState == 'on':
        return sp.simplify(paper)
    elif simplifyState == 'off':
        return paper
    pass

def RicciTensor(tensor, RCT, simplifyState = 'on'):
    Delta  = sp.eye(4)
    repeater = 4
    paper = ZeroArrayCreate(repeater,repeater)
    for a in range(0,repeater):
        for b in range(0,repeater):
            for y in range(0,repeater):
                answer = RCT[y,a,y,b]
                print(f'Ricci tensor component for index {a,b} and dummy index {y} has been calculated with answer {answer}')
                paper[a,b] += answer
    if simplifyState == 'on':
        return sp.simplify(paper)
    elif simplifyState == 'off':
        return paper
    pass

def RicciScalar(tensor, RT,  simplifyState = 'on'):
    paper = 0
    repeater = 4
    for a in range(0, repeater):
        for b in range(0, repeater):
            answer = tensor.inv[a,b] * RT[a,b]
            print(f'Ricci scalar value for dummy index {a,b} has been calculated with answer {answer}')
            paper += answer

    if simplifyState == 'on':
        return sp.simplify(paper)
    elif simplifyState == 'off':
        return paper
    pass
