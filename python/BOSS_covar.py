#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python


import matplotlib
import numpy
import matplotlib.pylab
from pylab import *
import numpy 
from matplotlib import rc

rc('font', family='sans serif', style='normal', variant='normal', stretch='normal', size='10.0')
lab_fontsize = 16
axes_fontsize = 14

F_arr = arange(0.575,0.772,0.004)
fsig8_arr = arange(0.25,0.60,0.007)
best_fsig8 = 0.4298
best_F = 0.6771
covar = numpy.zeros((2,2))
covar[0,0] = 0.004509868 
covar[1,1] = 0.001736087
covar[1,0] = 0.002435891
covar[0,1] = 0.002435891
covar = linalg.inv(covar)
diff = zeros(2)
prob = zeros((50,50))

for j in range(0,50):
    for i in range(0,50):
        fsig8 = fsig8_arr[i]
        F = F_arr[j]
        diff[0] = fsig8-best_fsig8
        diff[1] = F-best_F
        chi2 = dot(diff,dot(covar,diff))
        prob[j,i] = exp(-chi2/2.)


fig = plt.figure()
fig.subplots_adjust(bottom=0.2)
ax = fig.add_subplot(111)
plt.contour(fsig8_arr,F_arr,prob,[0.32,0.05])
ylabel(r'$F(z=0.57)$', fontsize=lab_fontsize)
xlabel(r'$f\sigma_8(z=0.57)$', fontsize=lab_fontsize)

savefig("BOSS_contour.pdf")
