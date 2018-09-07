from pylab import *
import ROOT
from ROOTtoPython import HistToHist
from ROOTtoPython import HistToLine



rf = ROOT.TFile("1M-SST_Albedo_DataBase.root")

hCamera = rf.Get("hAngleCamera")
hAlbedo = rf.Get("hAngleAlbedo")


fig, ax = subplots(1, 1)
hCamera.Scale(1./hCamera.Integral())
hAlbedo.Scale(1./hAlbedo.Integral())

_,x,ww = HistToHist(hCamera, ax, histtype="stepfilled", lw = 1, color = "blue")
#HistToHist(hAlbedo, ax, histtype="step", lw = 2, color = "black")
print x,ww

with file("lightguide_weights.txt",'w') as f:
    for i,w in enumerate(ww):
        f.write("{} {}\n".format(x[i],w))

ax.axvline(x=24, color="k", ls="--", lw=2)
ax.text(25, 0.030, "$\\"+"theta_{c} = 24^{\circ}$", fontsize=25, va="top", rotation=90, weight="bold")

ax.set_xlabel(r"$\Theta$ [$^{\circ}$]", fontsize=20)
ax.set_ylabel("Number of rays (norm.)", fontsize=20)
#ax.grid()
ax.set_xlim(0, 30)
ax.set_ylim(0, 0.035)

setp(ax.get_xticklabels(), fontsize = 15)
setp(ax.get_yticklabels(), fontsize = 15)


savefig("lightcamera.eps")

show()

