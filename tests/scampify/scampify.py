import h5py
import scipy
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, "/Users/ajaruga/clones/libcloudphxx/build/bindings/python")
from libcloudphxx import blk_1m

# read in the data

#fname_cnst = 'const_dycoms.h5'
#fname = 'timestep0000007150.h5'
#dt   = 1
#dz   = 5
#nbin = 50


dt   = 1
dz   = 20
nbin = 10 
fname_cnst = 'const_eddy.h5'
fname = 'timestep0000009000.h5'
fname_w = 'dycoms_w.txt'

data_cnst = h5py.File(fname_cnst, 'r')
data      = h5py.File(fname, 'r')

print("Keys: %s" % data_cnst.keys())
print("Keys: %s" % data.keys())
#a_group_key = list(f.keys())[0]

#w   = np.array(data['w'])
w = np.loadtxt(fname_w) * 400

rc  = np.array(data['rc']) * 1000
rr  = np.array(data['rr']) * 1000
rv  = np.array(data['rv']) * 1000
th  = np.array(data['th'])
rho = np.array(data_cnst['G'][0])

#print w.shape
height  = np.arange(0, rc.shape[1], 1)

np.set_printoptions(threshold=np.inf)

w_10    = np.percentile(w,  10, axis=0)
w_40    = np.percentile(w,  40, axis=0)
w_50    = np.percentile(w,  50, axis=0)
w_60    = np.percentile(w,  60, axis=0)
w_90    = np.percentile(w,  90, axis=0)

rc_10   = np.percentile(rc, 10, axis=0)
rc_40   = np.percentile(rc, 40, axis=0)
rc_50   = np.percentile(rc, 50, axis=0)
rc_60   = np.percentile(rc, 60, axis=0)
rc_90   = np.percentile(rc, 90, axis=0)

rv_10   = np.percentile(rv, 10, axis=0)
rv_40   = np.percentile(rv, 40, axis=0)
rv_50   = np.percentile(rv, 50, axis=0)
rv_60   = np.percentile(rv, 60, axis=0)
rv_90   = np.percentile(rv, 90, axis=0)

rr_10   = np.percentile(rr, 10, axis=0)
rr_40   = np.percentile(rr, 40, axis=0)
rr_50   = np.percentile(rr, 50, axis=0)
rr_60   = np.percentile(rr, 60, axis=0)
rr_90   = np.percentile(rr, 90, axis=0)

th_10   = np.percentile(th, 10, axis=0)
th_40   = np.percentile(th, 40, axis=0)
th_50   = np.percentile(th, 50, axis=0)
th_60   = np.percentile(th, 60, axis=0)
th_90   = np.percentile(th, 90, axis=0)

plt.figure(1)

#plt.subplot(221)
#plt.ylim(0,1500)
##plt.axhline(y=545, color='orange')
##plt.axhline(y=855, color='orange')
#plt.plot(w_10, height*dz, c='lightgray')
#plt.plot(w_90, height*dz, c='lightgray')
#plt.fill_betweenx(height*dz, w_10, w_90, color='lightgray')
#plt.plot(w_40, height*dz, c='gray')
#plt.plot(w_60, height*dz, c='gray')
#plt.fill_betweenx(height*dz, w_40, w_60, color='gray')
#plt.plot(w_50, height*dz, color='black')
#plt.axvline(x=0, color='orange')
#plt.ylabel("height [m]")
#plt.xlabel("w [m/s]")
plt.subplot(221)
plt.ylim(0,1500)
#plt.axhline(y=545, color='orange')
#plt.axhline(y=855, color='orange')
plt.plot(w_10, height[0:-1]*dz, c='lightgray')
plt.plot(w_90, height[0:-1]*dz, c='lightgray')
plt.fill_betweenx(height[0:-1]*dz, w_10, w_90, color='lightgray')
plt.plot(w_40, height[0:-1]*dz, c='gray')
plt.plot(w_60, height[0:-1]*dz, c='gray')
plt.fill_betweenx(height[0:-1]*dz, w_40, w_60, color='gray')
plt.plot(w_50, height[0:-1]*dz, color='black')
plt.axvline(x=0, color='orange')
plt.ylabel("height [m]")
plt.xlabel("w [m/s]")

plt.subplot(222)
plt.ylim(0,1500)
#plt.axhline(y=545, color='orange')
#plt.axhline(y=855, color='orange')
plt.plot(rv_10, height*dz, c='lightgray')
plt.plot(rv_90, height*dz, c='lightgray')
plt.fill_betweenx(height*dz, rv_10, rv_90, color='lightgray', label=r"$P_{10} - P_{90}$")
plt.plot(rv_40, height*dz, c='gray')
plt.plot(rv_60, height*dz, c='gray')
plt.fill_betweenx(height*dz, rv_40, rv_60, color='gray', label=r"$P_{40} - P_{60}$")
plt.plot(rv_50, height*dz, color='black', label='median')
plt.xlabel("rv [g/kg]")
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    labelleft='off') # labels along the bottom edge are off
plt.legend(loc="lower left", frameon=False)

plt.subplot(223)
plt.ylim(0,1500)
#plt.axhline(y=545, color='orange')
#plt.axhline(y=855, color='orange')
plt.plot(rc_10, height*dz, c='lightgray')
plt.plot(rc_90, height*dz, c='lightgray')
plt.fill_betweenx(height*dz, rc_10, rc_90, color='lightgray')
plt.plot(rc_40, height*dz, c='gray')
plt.plot(rc_60, height*dz, c='gray')
plt.fill_betweenx(height*dz, rc_40, rc_60, color='gray')
plt.plot(rc_50, height*dz, color='black')
plt.xlabel("rc [g/kg]")
plt.ylabel("height [m]")

plt.subplot(224)
plt.ylim(0,1500)
#plt.axhline(y=545, color='orange')
#plt.axhline(y=855, color='orange')
plt.plot(rr_10 * 10, height*dz, c='lightgray')
plt.plot(rr_90 * 10, height*dz, c='lightgray')
plt.fill_betweenx(height*dz, rr_10 * 10, rr_90 * 10, color='lightgray')
plt.plot(rr_40 * 10, height*dz, c='gray')
plt.plot(rr_60 * 10, height*dz, c='gray')
plt.fill_betweenx(height*dz, rr_40 * 10, rr_60 * 10, color='gray')
plt.plot(rr_50 * 10, height*dz, color='black', label='median')
plt.xlabel("rr [mg/kg]")
plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    labelleft='off') # labels along the bottom edge are off

plt.tight_layout()
plt.show()

opts = blk_1m.opts_t()
print "cond =", opts.cond
print "cevp =", opts.cevp
print "revp =", opts.revp
print "conv =", opts.conv
print "accr =", opts.accr
print "sedi =", opts.sedi
print "r_c0 =", opts.r_c0
print "r_eps =", opts.r_eps

var_th,  var_rc,  var_rv,   var_rr    = np.zeros(height.shape[0]), np.zeros(height.shape[0]), np.zeros(height.shape[0]), np.zeros(height.shape[0])
mean_th, mean_rc, mean_rv,  mean_rr   = np.zeros(height.shape[0]), np.zeros(height.shape[0]), np.zeros(height.shape[0]), np.zeros(height.shape[0])
dot_rc,  dot_rr,  var_flux, mean_flux = np.zeros(height.shape[0]), np.zeros(height.shape[0]), np.zeros(height.shape[0]), np.zeros(height.shape[0])

for it in range(0, height.shape[0]):
    print "height = ", height[it]*5
    tmp = np.histogram(rv[:,it], bins=nbin, density=False)
    bin_len = tmp[1][1] - tmp[1][0]

    tmp_th = np.array(float(max(273.15, np.mean(th[:, it]))))
    tmp_rc = np.array(float(max(0,      np.mean(rc[:, it]))) / 1000)
    tmp_rr = np.array(float(max(0,      np.mean(rr[:, it]))) / 1000)
    tmp_rv = np.array(float(max(0,      np.mean(rv[:, it]))) / 1000)

    rc_old = tmp_rc.copy()
    rr_old = tmp_rr.copy()
    in_dot_rc = np.array(0.)
    in_dot_rr = np.array(0.)

    blk_1m.adj_cellwise(opts, np.array(float(rho[it])), tmp_th, tmp_rv, tmp_rc, tmp_rr, dt)

    mean_th[it] = tmp_th
    mean_rv[it] = tmp_rv
    mean_rc[it] = tmp_rc
    mean_rr[it] = tmp_rr

    blk_1m.rhs_cellwise(opts, in_dot_rc, in_dot_rr, rc_old, rr_old)
    mean_flux[it] = blk_1m.rhs_columnwise(opts, in_dot_rr, np.array(float(rho[it])), rr_old, dz)
    
    mean_rc[it] += in_dot_rc * dt
    mean_rr[it] += in_dot_rr * dt

    for bin_it in range(0, nbin):
        idx_tmp = np.where((rv[:,it] >= tmp[1][bin_it]) & (rv[:,it] < tmp[1][bin_it+1]))

        in_th  = np.array(float(max(273.15, np.mean(th[idx_tmp, it]))))
        in_rc  = np.array(float(max(0,      np.mean(rc[idx_tmp, it]))) / 1000)
        in_rr  = np.array(float(max(0,      np.mean(rr[idx_tmp, it]))) / 1000)
        in_rv  = np.array(float(max(0,      np.mean(rv[idx_tmp, it]))) / 1000)
        weight = tmp[0][bin_it]/float(np.sum(tmp[0]))

        blk_1m.adj_cellwise(opts, np.array(float(rho[it])), in_th, in_rv, in_rc, in_rr, dt)

        var_th[it] += weight * in_th
        var_rv[it] += weight * in_rv
        var_rc[it] += weight * in_rc
        var_rr[it] += weight * in_rr

    for bin_it in range(0, nbin):
        idx_tmp = np.where((rv[:,it] >= tmp[1][bin_it]) & (rv[:,it] < tmp[1][bin_it+1]))

        in_rc     = np.array(float(max(0,      np.mean(rc[idx_tmp, it]))) / 1000)
        in_rr     = np.array(float(max(0,      np.mean(rr[idx_tmp, it]))) / 1000)
        in_dot_rc = np.array(0.)
        in_dot_rr = np.array(0.)
        weight    = tmp[0][bin_it]/float(np.sum(tmp[0]))

        blk_1m.rhs_cellwise(opts, in_dot_rc, in_dot_rr, in_rc, in_rr)
        var_flux[it] = blk_1m.rhs_columnwise(opts, in_dot_rr, np.array(float(rho[it])), in_rr, dz)
 
        dot_rc[it] += weight * in_dot_rc
        dot_rr[it] += weight * in_dot_rr

    var_rc[it] += dot_rc[it] * dt
    var_rr[it] += dot_rr[it] * dt

from pylab import rcParams
rcParams['figure.figsize'] = 3.5, 10

plt.figure(2)

#plt.subplot(221)
##plt.ylim(0,1500)
#plt.plot(mean_th, height*dz, color='orange')
#plt.plot(var_th,  height*dz, color='red')
#plt.plot(th_50,   height*dz, color='black')
#plt.ylabel("height [m]")
#plt.xlabel("th [K]")

lw=2.5
plt.subplot(311)
plt.ylim(0,1500)
plt.yticks([0, 250, 500, 750, 1000, 1250, 1500])
plt.plot(mean_rv * 1000, height*dz, color='orange', lw=lw, label='mean')
plt.plot(var_rv  * 1000, height*dz, color='red',    lw=lw, label='integral')
plt.plot(rv_50         , height*dz, color='black',  lw=lw, label='LES')
plt.ylabel("height [m]")
plt.xlabel("rv [g/kg]")
plt.legend(loc="upper right", frameon=False)

plt.subplot(312)
plt.ylim(0,1500)
plt.yticks([0, 250, 500, 750, 1000, 1250, 1500])
plt.plot(mean_rc * 1000, height*dz, color='orange', lw=lw)
plt.plot(var_rc  * 1000, height*dz, color='red',    lw=lw)
plt.plot(rc_50         , height*dz, color='black',  lw=lw )
plt.ylabel("height [m]")
plt.xlabel("rc [g/kg]")

plt.subplot(313)
plt.ylim(0,1500)
plt.yticks([0, 250, 500, 750, 1000, 1250, 1500])
plt.plot(mean_rr * 1000 * 10, height*dz, color='orange', lw=lw)
plt.plot(var_rr  * 1000 * 10, height*dz, color='red',    lw=lw)
plt.plot(rr_50          * 10, height*dz, color='black',  lw=lw)
plt.ylabel("height [m]")
plt.xlabel("rr [mg/kg]")

plt.tight_layout()
plt.savefig("eddy_results.pdf")
#plt.show()
