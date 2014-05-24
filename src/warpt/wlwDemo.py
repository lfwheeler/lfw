"""
Warping of well logs from Teapot Dome survey.
Author: Dave Hale, Colorado School of Mines
Version: 2013.12.28

Modified by: Loralee Wheeler, Colorado School of Mines
Date modified: 2014.2.17
"""

from java.awt import *
from java.io import *
from java.nio import *
from java.lang import *
from java.util import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.interp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mesh import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from tpt import *
from warpt import *

##########################################################################
wlw = WellLogWarping()
curve = "v" # velocity logs
logs = None


perr = 0.25 # power of norm used in alignment error calculation
ms = 250 # maximum shift
freplace = -2.0 # value to replace nulls for velocity
fclips = (2.0,6.0) 
cblabel = "Velocity (km/s)"

##########################################################################

def main(args):
  global logs,sz
  sz,logs = getLogs(curve)
  goShifts() # simultaneously aligns all logs
  #goWarping() # aligns pairs of logs
  #goErrors() # alignment errors for log pairs

"""
Performs simultaneous correlation of many logs. 
Plots logs before and after alignment.
"""
def goShifts():
  fs = logs
  nk,nl = len(fs[0]),len(fs)
  wlw.setPowError(perr) 
  wlw.setMaxShift(ms)
  s = wlw.findShifts(fs)
  gs = wlw.applyShifts(fs,s)
  s = mul(1000*sz.delta,s) # convert shifts to m
  fs = wlw.replaceNulls(fs,freplace)
  gs = wlw.replaceNulls(gs,freplace)

  # plots raw logs
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(700,380)
  sp.setVLabel("Depth (km)")
  sp.setHLabel("Log index")
  sp.addColorBar(cblabel)
  sp.plotPanel.setColorBarWidthMinimum(90)
  pv = sp.addPixels(sz,Sampling(nl,1,1),fs)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cjet)
  pv.setClips(fclips[0],fclips[1])
  
  # plots aligned logs
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(700,380)
  sp.setVLabel("Relative geologic time")
  sp.setHLabel("Log index")
  sp.addColorBar(cblabel)
  sp.plotPanel.setColorBarWidthMinimum(90)
  pv = sp.addPixels(sz,Sampling(nl,1,1),gs)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cjet)
  pv.setClips(fclips[0],fclips[1])


"""
Performs pairwise correlation between specified log pairs.
Plots log pairs before and after alignment.
"""
def goWarping():
  # define pairs of logs to compare below
  pairs = [(0,1),(1,2),(2,3),(3,4),(4,5)] 
  fs = logs
  wlw.setPowError(perr)
  wlw.setMaxShift(ms)
  for pair in pairs:
    jf,jg = pair[0],pair[1]
    fi,gj = fs[jf],fs[jg]
    e = wlw.computeErrors(fi,gj)
    nl,nk = len(e[0]),len(e)
    d = wlw.accumulateErrors(e)
    kl = wlw.findWarping(d)
    fk,gk = wlw.applyWarping(kl,fi,gj)
    fi = wlw.replaceNulls(fi,freplace)
    gj = wlw.replaceNulls(gj,freplace)
    fk = wlw.replaceNulls(fk,freplace)
    gk = wlw.replaceNulls(gk,freplace)
    ii,ff = removeZeros(fi)
    jj,gg = removeZeros(gj)
    ki,kf = removeZeros(fk)
    kj,kg = removeZeros(gk)
    for ij in range(len(ff)):
      if ff[ij] == -2.0:
        ff[ij] = 3.0
    for ij in range(len(gg)):
      if gg[ij] == -2.0:
        gg[ij] = 3.0
    for ij in range(len(kf)):
      if kf[ij] == -2.0:
        kf[ij] = 3.0
    for ij in range(len(kg)):
      if kg[ij] == -2.0:
        kg[ij] = 3.0
    sff = Sampling(len(ff),sz.getDelta(),sz.getFirst()+ii*sz.getDelta())
    sgg = Sampling(len(gg),sz.getDelta(),sz.getFirst()+jj*sz.getDelta())
    skf = Sampling(len(kf),0.5*sz.getDelta(),sz.getFirst()+ki*0.5*sz.getDelta())
    skg = Sampling(len(kg),0.5*sz.getDelta(),sz.getFirst()+kj*0.5*sz.getDelta())
    title = "("+str(jf)+","+str(jg)+")"
    if True:
      sp = SimplePlot()
      sp.setSize(750,500)
      sp.setTitle(title)
      sp.setHLabel("Depth (km)")
      sp.setVLabel(cblabel)
      sp.setVLimits(fclips[0],fclips[1])
      sp.setHLimits(0,1.9)
      pv = sp.addPoints(sff,ff)
      pv.setLineColor(Color.BLACK)
      pv = sp.addPoints(sgg,gg)
      pv.setLineColor(Color.RED)
    if True:
      sp = SimplePlot()
      sp.setSize(750,500)
      sp.setTitle(title)
      sp.setHLabel("Depth (km)")
      sp.setVLabel(cblabel)
      sp.setVLimits(fclips[0],fclips[1])
      sp.setHLimits(0,1.9)
      pv = sp.addPoints(skf,kf)
      pv.setLineColor(Color.BLACK)
      pv = sp.addPoints(skg,kg)
      pv.setLineColor(Color.RED)


"""
Computes alignment errors on the kl-coordinate system for specified log 
pairs and finds an optimal path.
Plots alignment error array with and without an optimal path.
"""
def goErrors():
  # define pairs of logs to compare below
  pairs = [(0,1),(1,2),(2,3),(3,4),(4,5)] 
  fs = logs
  wlw.setPowError(perr)
  wlw.setMaxShift(ms)
  for pair in pairs:
    ia,ib = pair[0],pair[1]
    e = wlw.computeErrors(fs[ia],fs[ib])
    wlw.interpolateOddErrors(e)
    nl,nk = len(e[0]),len(e)
    lmax = (nl-1)/2
    lmin = -lmax
    sl = Sampling(nl,1,lmin)
    sk = Sampling(nk,1,0)

    title = "("+str(ia)+","+str(ib)+")"
    # plots alignment errors 
    sp = SimplePlot()
    sp.setSize(750,500)
    sp.setTitle(title)
    sp.setHLabel("Depth index k")
    sp.setVLabel("Lag index l")
    sp.setHFormat("%5f")
    pv = sp.addPixels(sk,sl,transpose(e))
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setPercentiles(0,100)
    d = wlw.accumulateErrors(e)
    wlw.interpolateOddErrors(d)

    # plots alignment errors with optimal path
    sp = SimplePlot()
    sp.setSize(750,500)
    sp.setTitle(title)
    sp.setHLabel("Depth index k")
    sp.setVLabel("Lag index l")
    sp.setHFormat("%5f")
    pv = sp.addPixels(sk,sl,transpose(e))
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setPercentiles(0,100)
    d = wlw.accumulateErrors(e)
    wlw.interpolateOddErrors(d)
    kw,lw = wlw.findWarping(d)
    kw = wlw.toFloat(kw)
    lw = wlw.toFloat(lw)
    pv = sp.addPoints(kw,lw)
    pv.setLineColor(Color.WHITE)

    """
    # plot accumulated alignment errors
    sp = SimplePlot()
    sp.setSize(750,500)
    sp.setTitle(title)
    sp.setHLabel("Depth index k")
    sp.setVLabel("Lag index l")
    sp.setHFormat("%5f")
    pv = sp.addPixels(sk,sl,transpose(d))
    #pv.setColorModel(cjet)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setPercentiles(0,90)
    kw,lw = wlw.findWarping(d)
    kw = wlw.toFloat(kw)
    lw = wlw.toFloat(lw)
    pv = sp.addPoints(kw,lw)
    pv.setLineColor(Color.WHITE)
    """


#############################################################################
# utilities

def getLogs(type):
  fileName = type+"logs.txt"
  ifile = open(fileName,'r+')
  lines = ifile.readlines()
  nz = len(lines)
  ifile.close()
  ifile = open(fileName,'r+')
  c = 0
  start = False
  while (start != True):
    line = ifile.readline()
    c += 1
    if line == "":
      print 'End of file'
      break
    elif line[0] == '~':
      start = True

  line = ifile.readline()
  wdata = line.split('\t')
  nl = len(wdata)
  logs = zerofloat(nz-c,nl-1)
  depth = zerodouble(nz-c)
  depth[0] = float(wdata[0])
  for l in range(1,nl):
    logs[l-1][0] = float(wdata[l])

  i = 1
  start = False
  while (start != True):
    line = ifile.readline()
    if line == "":
      start = True
    else:
      wdata = line.split('\t')
      depth[i] = float(wdata[0])
      for l in range(1,nl):
        logs[l-1][i] = float(wdata[l])
      i += 1
  ifile.close()
  zs = Sampling(depth)
  sz = Sampling(zs.count,zs.delta*0.001,zs.first)
  return sz,logs

def removeZeros(f):
  n = len(f)
  i = 0
  while f[i] == -2.0:
    i += 1
  fs = zerofloat(n)
  c = 0
  for j in range(i,n):
    fs[c] = f[j]
    if (f[j] == -2.0 and j+2 >= n) or (f[j] == -2.0 and f[j+1] == -2.0
        and f[j+2] == -2.0):
      break
    c += 1
  #ft = zerofloat(c)
  ft = copy(c,fs)
  return i,ft

#############################################################################
# graphics

cjet = ColorMap.JET
alpha = fillfloat(1.0,256); alpha[0] = 0.0
ajet = ColorMap.setAlpha(cjet,alpha)

#############################################################################
# Run the function main on the Swing thread
import sys
from javax.swing import *
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
