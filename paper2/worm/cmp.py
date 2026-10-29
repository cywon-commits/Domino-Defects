import sys,numpy as np,subprocess
W,H,J,xt,yt=sys.argv[1:6]
ex=np.array(subprocess.run(['python3','brute.py',W,H,J,xt,yt],capture_output=True,text=True).stdout.split(),float)
subprocess.run(['./worm',W,H,J,xt,yt,'4e7','3','h.bin'],capture_output=True)
mc=np.fromfile('h.bin'); mc/=mc.sum()
print('W=%s H=%s J=%s tail=(%s,%s): max|mc-exact| = %.4f  (max prob %.3f)'%(W,H,J,xt,yt,np.abs(mc-ex).max(),ex.max()))
