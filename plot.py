import numpy as np 
import matplotlib.pyplot as plt 

# with open("result_1000_b") as f:
# 	lines = f.readlines()

# lines = lines[1:]
# x = [float(line.split()[1]) for line in lines]
# y = [float(line.split()[3]) for line in lines]

# plt.plot(x,y, label = 'Babaioff 2018')

with open("UniformMatroids/result.txt") as f:
	lines = f.readlines()

lines = lines[1:]
x = [float(line.split()[1]) for line in lines]
y = [float(line.split()[2]) for line in lines]

plt.plot(x,y, label = 'Our algorithm')	

# with open("result_1000_k") as f:
# 	lines = f.readlines()

# lines = lines[1:]
# x = [float(line.split()[1]) for line in lines]
# y = [float(line.split()[2]) for line in lines]

# plt.plot(x,y, label = 'Kleinberg\'s algorithm')

leg = plt.legend()
plt.savefig('bruteforcecheck.png')
plt.show()		