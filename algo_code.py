import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
#Given l and time t, find out the support and resistance by finding out the time points, through which the lines should pass.    
#col is the input vector of time series. 
# The output is a vector m[0,1,2,3,4,5] where m[1] and m[2] are the time points so that the support line is obtained by joining
# (m[1],col[m[1]]) and (m[2], col[m[2]]). 
# m[0] gives the distance of CG from the support line. 
# Similarly m[3] gives the distance of CG from the resistance line and m[4],m[5] give the resistance line.

col=[9835.4, 9837.65, 9847, 9839.3, 9834.1, 9825.6, 9822.7, 9814.35, 9811.65, 9806.95, 9821.1, 9822.4, 9816.95, 9810.05, 9808.6, 9809.2, 9797.35, 9786.25, 9787.15, 9793.9, 9801.75, 9803.5, 9804.45, 9800, 9802.4, 9806.3, 9802.25, 9807.05, 9809, 9807.65, 9802.95, 9798.6, 9796.95, 9797.4, 9798.3, 9793.95, 9796.7, 9799.15, 9804.8, 9810.1, 9800.9, 9794.95, 9794.3, 9790.35, 9792.9, 9785, 9790.05, 9790.8, 9790.95, 9793.65, 9792.7, 9796.05, 9796.45, 9792.35, 9786.8, 9782.45, 9791.7, 9789.3, 9793.65, 9797.45, 9790, 9794.15, 9796.05, 9793.1, 9779.5, 9777.35, 9775, 9764.25, 9768.35, 9773.8, 9773.85, 9765.8, 9765.5, 9761.05, 9768.35, 9767.75, 9783.65, 9809.05, 9803.75, 9805.3, 9807.7, 9805.7, 9803.6, 9814.75, 9848.5, 9851.65, 9844.5, 9841.95, 9842.05, 9835.25, 9847.5, 9839.6, 9836.5, 9828.6, 9832.2, 9838.25, 9828.75, 9829.05, 9821.3, 9816.05]
t=80
l=65

#defining function SnR
def SnR(data,t,l):
	if l*0.5-(l/2)!=0:
	        l=l-1
	cond1=[[0 for x in range(l)] for y in range(l)]
	cond2=[[0 for x in range(l)] for y in range(l)]
	J=[[0 for x in range(l)] for y in range(l)]
	m=[0 for x in range(4)]
	x=list(data[t-l:t])
	g=[t-0.5*l,sum(x)/l]
	m[0]=10**8
	m[2]=10**8
	for i in range(0,l//2):
		for j in range(l//2,l):
	            a=np.abs((x[i]-x[j])*g[0]-(i-j)*g[1]+ (x[j]*i-x[i]*j))
	            b=np.sqrt((x[i]-x[j])**2+(i-j)**2)
	            J[i][j]=a/b
	            C= (x[i]-x[j])/(i-j)
	            for k in range (0,l):
	                if x[k]< x[j]+C*(k-j):
	                    cond1[i][j]+=1
	                elif x[k]== x[j]+C*(k-j):
	                    cond2[i][j]+=1
	            if cond1[i][j]==0:
	                    m[0]=min(m[0],J[i][j])
	                    m[1:2]=[i,j]
	            elif cond1[i][j]+cond2[i][j]==l:
	                m[3]=min(m[3],J[i][j])
	                if m[3]==J[i][j]:
	                        m[4:5]=[i,j]
		                    
	m[1]+=t-l
	m[2]+=t-l
	m[4]+=t-l
	m[5]+=t-l
	print("m = ", m)

		#to print the time series graph with line of support and resistance in continuation with above	code 
	print("Total no. of datapoints = ",len(data))
	#x=[]
	#for i in range (1,len(data)):
	#    x.append(i)
	#plt.scatter(x, data, label = 'Scatterplot',  color = 'k', s = 10)
	#plt.xlabel('x')
	#plt.ylabel('y')
	#plt.title('Interesting Graph')

	#plt.show()
	#Rounding off m to the nearest integer
	ar = np.array(m)
	ar = np.round(ar,decimals = 0)
	ar = ar.astype(np.int)   #Converts float into int
	print("m as rounded off to the nearest integer",ar)
	x1, y1 = [ar[1],ar[2]],[data[ar[1]],data[ar[2]]]
	x2, y2 = [ar[4],ar[5]],[col[ar[4]],col[ar[5]]]
	print("({},{}) and ({},{}) are points lying on the support, and ({},{}) and ({},{}) are points lying on the resistance.".format(x1[0],y1[0],x1[1],y1[1],x2[0],y2[0],x2[1],y2[1]))
	#plt.plot(x1,y1,x2,y2,marker = 'o')
	#plt.show()

#defining function to compute delta
def delta(data,t,l,x_sup, y_sup, x_res, y_res):
	x=list(data[t-l:t])
	# for calculating delta_plus (based on resistance equation)
	count1 = 0
	for i in range(t+1, len(data)):
		L_plus = y_res[1] + ((y_res[1]-y_res[0])*(i-x_res[1])/(x_res[1]-x_res[0]))
		if data[i] <= L_plus:
			count1 = count1 + 1
	if count1 == len(data) - l:
		count1 = min(count1,l)
	#for calculating delta_minus (based on support equation)
	count2 = 0
	for j in range(t+1, len(data)):
		L_minus = y_sup[1] + ((y_sup[1]-y_sup[0)*(j-x_sup[1])/(x_sup[1]-x_sup[0]))
		if data[i] >= L_minus:
			count2 = count2 + 1
	if count2 == len(data) - l:
		count2 = min(count2,l)

#defining function to determine the empirical distribution of delta_plus and delta_minus
def F_delta(delta_plus, delta_minus): #delta_plus and delta_minus are arrays containing delta values for t = l, l+1, ... , N
	A = sorted(delta_plus)
	B = sorted(delta_minus)
	F_delta_plus = [0 for x in range(max(delta_plus))]
	F_delta_minus = [0 for x in range(max(delta_plus))]
	#empirical distribution of delta_plus
	for k in range(max(delta_plus)):
		count = 0
		for i in range(len(A)):
			if A[i]<= k:
				count = count +1
			else:
				break
		F_delta_plus[k] = count/(len(delta_plus))
	
	#empirical distribution of delta_minus
	for k in range(max(delta_minus)):
		count = 0
		for i in range(len(B)):
			if B[i]<= k:
				count = count +1
			else:
				break
		F_delta_minus[k] = count/(len(delta_minus))
		
#calling the functions based on data			
	

