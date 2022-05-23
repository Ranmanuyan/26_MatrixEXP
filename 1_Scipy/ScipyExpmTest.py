import numpy as np
from scipy.linalg import expm
from time import perf_counter



a = np.loadtxt('Mat169.txt')
a1= np.array(a)

s = (len(a1),1)
b = np.ones(s)

start = perf_counter()
ex = expm(450* a1)
cdf  = 1- np.matmul(ex,b)[0];
end = perf_counter()

print(end - start)