import numpy as np
import matplotlib.pyplot as plt
from statsmodels.distributions.empirical_distribution import ECDF
import math



# sum of two

# sample = np.random.exponential(scale=2, size=10000)+np.random.exponential(scale=1/0.7, size=10000)
# ecdf = ECDF(sample)

# x1 = np.arange(0, 7, 0.0001)
# F = [1-0.5 *(7* math.exp(-0.5 * x2)- 5* math.exp(-0.7 * x2)) for x2 in x1]


# plt.plot(x1, F)


# plt.plot(ecdf.x, ecdf.y)
# # h = plt.hist(np.random.exponential(scale=0.5, size=100000)+np.random.exponential(scale=0.7, size=100000), bins=200, density=True)
# plt.show()


# sum of three test

# lambda 0.5, 0.7, 0.2
# scale  2,  1/0.7, 5

# sample = np.random.exponential(scale=2, size=10000)+np.random.exponential(scale=1/0.7, size=10000)++np.random.exponential(scale=5, size=10000)
# ecdf = ECDF(sample)

# x1 = np.arange(0, 40, 0.0001)
# F = [1- (math.exp(-0.7 * x2) +  7/3 * math.exp(-0.2 * x2) -7/3 * math.exp(-0.5 * x2)) for x2 in x1]


# plt.plot(x1, F)


# plt.plot(ecdf.x, ecdf.y)
# # h = plt.hist(np.random.exponential(scale=0.5, size=100000)+np.random.exponential(scale=0.7, size=100000), bins=200, density=True)
# plt.show()




# dominance rule for M1 and M2
# sample1 = np.random.exponential(scale=2, size=10000)+np.random.exponential(scale=1/0.7, size=10000)++np.random.exponential(scale=5, size=10000)
# ecdf1 = ECDF(sample1)


# sample2 = np.random.exponential(scale=5, size=10000)+np.random.exponential(scale=2, size=10000)++np.random.exponential(scale=2, size=10000)
# ecdf2 = ECDF(sample2)

# plt.plot(ecdf1.x, ecdf1.y, color ='r')

# plt.plot(ecdf2.x, ecdf2.y, color ='b')


# plt.show()












sample1 = np.random.exponential(scale=1/0.48, size=10000)+ np.random.exponential(scale=1/0.27, size=10000)+ np.random.exponential(scale=1/0.52, size=10000)+ np.random.exponential(scale=1/0.17, size=10000)+ np.random.exponential(scale=1/0.09, size=10000)+ np.random.exponential(scale=1/0.07, size=10000)+np.random.exponential(scale=1/0.49, size=10000)+np.random.exponential(scale=1/0.03, size=10000)+np.random.exponential(scale=1/0.19, size=10000)+np.random.exponential(scale=1/0.18, size=10000)+np.random.exponential(scale=1/0.16, size=10000)

ecdf1 = ECDF(sample1)



plt.plot(ecdf1.x, ecdf1.y, color ='r')




plt.show()


