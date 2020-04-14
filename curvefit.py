import numpy as np 
import pandas as pd 
import matplotlib

def get_features(filepath):

	df = pd.read_csv(filepath, delimiter = ' ')
	read = df.rename_axis('ID').values

	numrows = read.shape[0]
	numcols = read.shape[1]
	
	y1 = np.array(read[:,2:3])
	y2 = np.array(read[:,3:4])

	phi = np.zeros([numrows,8])

	phi[:,0:1] = np.ones([numrows,1])
	phi[:,1:2] = read[:,1:2]/read[:,0:1]
	phi[:,2:3] = read[:,0:1]/read[:,1:2]
	phi[:,3:4] = np.sqrt(read[:,1:2])
	phi[:,4:5] = 1/np.sqrt(read[:,1:2])	
	phi[:,5:6] = phi[:,4:5]**2
	phi[:,6:7] = phi[:,4:5]**3
	phi[:,7:8] = phi[:,4:5]**4

	return phi, y1, y2

def closed_soln(phi, y):
	return np.linalg.pinv(phi).dot(y)

def compute_RMSE(phi, w, y):
	numrows = phi.shape[0]
	error = 0

	for i in range(numrows):
		error += 1/len(phi)*((np.dot(phi[i],w)[0]-y[i])**2)

	error[0] = np.sqrt(error[0])
	return error[0]

def main():
	phi, y1, y2 = get_features('result_120')

	w1 = closed_soln(phi,y1)
	w2 = closed_soln(phi,y2)

	print(w1, compute_RMSE(phi,w1,y1))
	print(w2, compute_RMSE(phi,w2,y2))

main()
