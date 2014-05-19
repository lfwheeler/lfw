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
curve = "v" # well log type. v = velocity, d = density, etc
logs = None


perr = 0.25 # power of norm used in alignment error calculation
ms = 250 # maximum shift
freplace = 2.0 # value to replace nulls for velocity
fclips = (2.0,6.0) 
cblabel = "Velocity (km/s)"
if curve=="d": # settings for density
  freplace = 1.0
  fclips = (2.0,2.8)
  cblabel = "Density (g/cc)"
if curve=="p": # settings for porosity
  freplace = 0.0
  fclips = (0.0,0.45)
  cblabel = "Porosity"
if curve=="g": # settings for gamma ray
  freplace = 30.0
  fclips = (30.0,160.0)
  cblabel = "Gamma ray (API)"

##########################################################################

def main(args):
  global logs
  logs = getLogs("d",curve)
  goShifts() # simultaneously aligns all logs
  #goWarping() # aligns pairs of logs
  #goErrors() # alignment errors for log pairs
  #goSort() # plots well locations

"""
Performs simultaneous correlation of many logs. 
Plots logs before and after alignment.
"""
def goShifts():
  sz,fs = resample(logs,curve)
  # chose a subset of logs below if you do not desire the entire set
  fs = [fs[0],fs[4],fs[9],fs[14],fs[17],fs[20]] # deepest 6 velocity logs
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
  pairs = [(0,4),(4,9),(9,14),(14,17),(17,20)] # deepest 6 velocity logs
  sz,fs = resample(logs,curve)
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
    title = "("+str(jf)+","+str(jg)+")"
    sk = Sampling(2*sz.count-1,0.5*sz.delta,sz.first)
    if True:
      sp = SimplePlot()
      sp.setSize(750,500)
      sp.setTitle(title)
      sp.setHLabel("Depth (km)")
      sp.setVLabel(cblabel)
      sp.setVLimits(fclips[0],fclips[1])
      pv = sp.addPoints(sz,fi)
      pv.setLineColor(Color.BLACK)
      pv = sp.addPoints(sz,gj)
      pv.setLineColor(Color.RED)
    if True:
      sp = SimplePlot()
      sp.setSize(750,500)
      sp.setTitle(title)
      sp.setHLabel("Depth (km)")
      sp.setVLabel(cblabel)
      sp.setVLimits(fclips[0],fclips[1])
      pv = sp.addPoints(sk,fk)
      pv.setLineColor(Color.BLACK)
      pv = sp.addPoints(sk,gk)
      pv.setLineColor(Color.RED)


"""
Computes alignment errors on the kl-coordinate system for specified log 
pairs and finds an optimal path.
Plots alignment error array with and without an optimal path.
"""
def goErrors():
  # define pairs of logs to compare below
  pairs = [(0,4),(4,9),(9,14),(14,17),(17,20)]
  sz,f = resample(logs,curve)
  wlw.setPowError(perr)
  wlw.setMaxShift(ms)
  for pair in pairs:
    ia,ib = pair[0],pair[1]
    e = wlw.computeErrors(f[ia],f[ib])
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
    #pv.setColorModel(cjet)
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
    #pv.setColorModel(cjet)
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

"""
Plots locations of well logs on constant time slice.
"""
def goSort():
  x,y = [],[]
  for index,log in enumerate(logs):
    x.append(log.x2[0])
    y.append(log.x3[0])
  seis = readImage()
  nx = len(seis)
  ny = len(seis[0])
  dy = 0.025
  dx = 0.025
  sx = Sampling(nx,dx,0.0)
  sy = Sampling(ny,dx,0.0)
  ss = zerofloat(ny,nx)
  for ix in range(nx):
    for iy in range(ny):
      ss[ix][iy] = seis[ix][iy][500]
  sp = SimplePlot()
  sp.addPixels(sy,sx,ss)
  sp.setSize(700,380)
  sp.setHLabel("Crossline (km)")
  sp.setVLabel("Inline (km)")
  pv = sp.addPoints(x,y)
  pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  pv.setLineStyle(PointsView.Line.NONE)
  pv.setMarkSize(6)
  pv.setMarkColor(Color.RED)
  if curve=="v":
    x = [x[0],x[4],x[9],x[14],x[17],x[20]]
    y = [y[0],y[4],y[9],y[14],y[17],y[20]]
    pv = sp.addPoints(x,y)
    pv.setLineStyle(PointsView.Line.NONE)
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    pv.setMarkSize(10)
    pv.setMarkColor(Color.RED)


#############################################################################
# utilities

_tpDir = "/data/seis/tpd/"
_csmDir = _tpDir+"csm/"
_wellLogsDir = _csmDir+"welllogs/"
_seismicLogsDir = _csmDir+"seismicz/"
  
def getLogs(set,type):
  fileName = _wellLogsDir+"tpw"+set[0]+".dat"
  wdata = WellLog.Data.readBinary(fileName)
  logs = []
  logt = wdata.getLogsWith(type)
  for log in logt:
    logs.append(log)
  logs = sortLogs(logs)
  return logs

def sortLogs(logs):
  nlog = len(logs)
  xlog = zerodouble(nlog)
  ylog = zerodouble(nlog)
  for i,log in enumerate(logs):
    xlog[i] = log.x2[0]
    ylog[i] = log.x3[0]
  j = WellLogWarping.sortWells(xlog,ylog)
  logt = list(logs)
  for i,log in enumerate(logs):
    logt[i] = logs[j[i]]
  return logt

def resample(logs,curve):
  nlog = len(logs)
  zs = zerofloat(0,nlog)
  fs = zerofloat(0,nlog)
  for i,log in enumerate(logs):
    zs[i] = log.z
    fs[i] = log.getCurve(curve)
  zs = mul(zs,0.0003048) # ft to km
  sz = wlw.getDepthSampling(zs,fs)
  nz,dz,fz,lz = sz.count,sz.delta,sz.first,sz.last
  #print "resample before: nz =",nz," dz =",dz," fz =",fz
  dz = 0.001 # 1 m
  nz = 1+int((lz-fz)/dz)
  sz = Sampling(nz,dz,fz)
  #print "resample  after: nz =",nz," dz =",dz," fz =",fz
  fs = wlw.resampleLogs(sz,zs,fs)
  return sz,fs

def readLogSamples(set,type,smooth=0):
  """ 
  Reads log curves from the specified set that have the specified type.
  set: "s" for shallow, "d" for deep, or "a" for all
  type: "v" (velocity), "d" (density), "p" (porosity), or "g" (gamma)
  smooth: half-width of Gaussian smoothing filter
  Returns a tuple (f,x1,x2,x3) of lists of arrays of samples f(x1,x2,x3)
  """
  logs = getLogs(set,type)
  fl,x1l,x2l,x3l = [],[],[],[]
  for log in logs:
    if smooth: 
      log.smooth(smooth)
    samples = log.getSamples(type)
    if samples:
      f,x1,x2,x3 = samples
      fl.append(f)
      x1l.append(x1)
      x2l.append(x2)
      x3l.append(x3)
  return fl,x1l,x2l,x3l

def getWellIntersections(set,type,x1):
  fileName = _wellLogsDir+"tpw"+set[0]+".dat"
  wdata = WellLog.Data.readBinary(fileName)
  x2,x3 = wdata.getIntersections(type,x1)
  return x2,x3

def fgood(f):
  n = len(f)
  for i in range(n):
    if f[i]!=-999.2500:
      return i
def lgood(f):
  n = len(f)
  for i in range(n):
    if f[n-1-i]!=-999.2500:
      return n-1-i

def readImage():
  fileName = _seismicLogsDir+"tpsz.dat"
  n1,n2,n3 = 2762,357,161
  x = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName,ByteOrder.BIG_ENDIAN)
  ais.readFloats(x)
  ais.close()
  return x

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
