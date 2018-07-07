# code for converting the position file into a csv file
pylab
t1,xbh1, ybh1 = np.loadtxt("BH_diagnostics.ah1.gp", usecols=(1,2,3), unpack=True)
t2, xbh2, ybh2 = np.loadtxt("BH_diagnostics.ah2.gp", usecols=(1,2,3), unpack=True)


import csv
with open('positions-b4.csv', 'w') as f:
writer = csv.writer(f, delimiter=' ')
writer.writerows(zip(t1,xbh1,ybh1,t2,xbh2,ybh2))
 quit()

