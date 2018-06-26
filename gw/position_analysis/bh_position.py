import matplotlib
from StringIO import StringIO
xbh1, ybh1 = np.loadtxt("BH_diagnostics.ah1.gp", usecols=(2,3), unpack=True)
xbh2, ybh2 = np.loadtxt("BH_diagnostics.ah2.gp", usecols=(2,3), unpack=True)
axis (" equal ")
plot(xbh1, ybh1, linestyle="-", color="black") 
plot(xbh2, ybh2, linestyle="-.", color="red")

# write to file
>>> a=[1,2,3]
>>> b=[4,5,6]
>>> zip(a,b)
[(1, 4), (2, 5), (3, 6)]
>>> import csv
>>> with open('text.csv', 'w') as f:
...    writer = csv.writer(f, delimiter='\t')
...    writer.writerows(zip(a,b))
...
>>> quit()
$ cat text.csv
1       4
2       5
3       6
