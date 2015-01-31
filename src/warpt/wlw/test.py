
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

from tp import *
from warpt import *

wlw = WellLogWarpingT()
curve = "p"
logs = None

def main(args):
  global logs
  logs = getLogs("d",curve)
  goTest()
  #goTestCG()

def goTest():
  sz,fs = resample(logs,curve)
  #fs = [fs[37],fs[38]]
  #fs = [fs[16],fs[18],fs[19],fs[28],fs[33],fs[34]]
  #fs = [fs[35],fs[37],fs[38]]
  #//fs = [fs[16],fs[18],fs[19],fs[28],fs[33],fs[34],fs[35],fs[37],fs[38]]
  fs = [fs[16],fs[18],fs[19],fs[28],fs[33],
        fs[34],fs[35],fs[37],fs[38],fs[39],
        fs[45],fs[50],fs[68]] # deepest 13 porosity logs
  nk,nl = len(fs[0]),len(fs)
  wlw.setPowError(0.25)
  wlw.setMaxShift(250)
  s = wlw.findShifts(fs)
  #for i in range(len(s)):
   # for j in range(len(s[0])):
    #  if (s[i][j] > 150 or s[i][j] < -150):
     #   print "l=",(i+1)," d=",j,"s=",s[i][j]
  #gs = wlw.applyShiftsZ(fs,s)
  gs = wlw.applyShifts(fs,s)
  s1 = mul(1000*sz.delta,s) # convert shifts to m
  freplace = 2.0
  fclips = (2.0,6.0)
  cblabel = "Velocity (km/s)"
  if curve=="d":
    freplace = 1.0
    fclips = (2.0,2.8)
    cblabel = "Density (g/cc)"
  if curve=="p":
    freplace = 0.0
    fclips = (0.0,0.45)
    cblabel = "Porosity"
  if curve=="g":
    freplace = 30.0
    fclips = (30.0,160.0)
    cblabel = "Gamma ray (API)"
  fs = wlw.replaceNulls(fs,freplace)
  gs = wlw.replaceNulls(gs,freplace)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(650,550)
  #sp.setSize(700,380) # for slides
  #sp.setSize(448,874) # for slides - each log
  #sp.setSize(617,376)
  sp.setVLabel("Depth (km)")
  sp.setHLabel("Log index")
  sp.addColorBar(cblabel)
  #sp.setVLimits(0.892,1.142)
  sp.plotPanel.setColorBarWidthMinimum(90)
  pv = sp.addPixels(sz,Sampling(nl,1,1),fs)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cjet)
  pv.setClips(fclips[0],fclips[1])
  #sp.setFontSizeForSlide(9.5/11.0,1.0,16.0/9.0)
  #sp.setFontSizeForPrint(8.0,255.6)
  #sp.paintToPng(720.0,3.55,"original.png")
  #sp.setVLimits(0.6,1.0) 
  #sp.paintToPng(720.0,6.51,"rawzoomshift.png")
  #sp.setFontSizeForPrint(8.0,234.5)
  #sp.paintToPng(720.0,3.25,"original.png")
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  #sp.setSize(650,508) # zoom & werror
  #sp.setSize(700,380) # for slides
  sp.setSize(650,550)
  #sp.setSize(617,376)
  #sp.setSize(650,550)
  sp.setVLabel("Relative geologic time")
  sp.setHLabel("Log index")
  #sp.setVLimits(1.25,1.4) # zoom for slides
  #sp.setVLimits(1.2,1.8) # zoom
  #sp.setVLimits(0.88,1.05) # werror
  sp.addColorBar(cblabel)
  sp.plotPanel.setColorBarWidthMinimum(90)
  pv = sp.addPixels(sz,Sampling(nl,1,1),gs)
  #pv.setClips(30,160) # zoom gamma
  #pv.setClips(2.8,6) # zoom velocity
  #pv.setClips(2.5,5) # zoom velocity error maybe
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cjet)
  pv.setClips(fclips[0],fclips[1])
  #pv.setClips(2.5,5.5)
  #sp.setFontSizeForSlide(9.5/11.0,1.0,16.0/9.0)
  #sp.setFontSizeForPrint(8.0,234.5)
  #sp.paintToPng(720.0,3.25,"warped.png")
  #sp.setFontSizeForPrint(8.0,222.0)
  #sp.paintToPng(720.0,3.08,"zoomv.png")
  #sp.paintToPng(720.0,3.08,"werror.png")
  sp1 = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp1.setSize(650,550)
  sp1.setVLabel("Depth (km)")
  sp1.setHLabel("Log index")
  sp1.addColorBar("Shifts (m)")
  sp1.plotPanel.setColorBarWidthMinimum(90)
  pv = sp1.addPixels(sz,Sampling(nl,1,1),s)
  #pv.setClips(-250,250)
  #pv.setClips(-225,175)
  #pv.setClips(-150,150)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cjet)
  #sp1.setFontSizeForPrint(8.0,156.3)
  #sp1.paintToPng(720.0,6.51,"shifts.png")

  
  """
  pp = PlotPanel(1,2,PlotPanel.Orientation.X1DOWN_X2RIGHT,PlotPanel.AxesPlacement.LEFT_BOTTOM)
  #pp = PlotPanel(1,2,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  pp.setSize(650,550)
  pp.setVLabel("Depth (km)")
  pp.setHLabel(0,"Log index")
  pp.setHLabel(1,"Log index")
  pp.addColorBar("Velocity (km/s)")
  pp.setColorBarWidthMinimum(90)
  pv = pp.addPixels(0,0,sz,Sampling(nl),fs)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cjet)
  pv.setClips(fclips[0],fclips[1])
  pv = pp.addPixels(0,1,sz,Sampling(nl),gs)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cjet)
  pv.setClips(fclips[0],fclips[1])
  pf = PlotFrame(pp)
  pf.setVisible(True)
  pf.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE)
  """

def goTestCG():
  l1 = [2,5,6,6,0,6,8]
  l2 = [3,0,2,6,6,5,6]
  l3 = [3,3,2,6,6,8,7]
  ll = [l1,l2,l3]
  nls = [[False,False,False,False,True,False,False],[False,True,False,False,False,False,False],[False,False,False,False,False,False,False]]

  s = wlw.testCG()
  sm = [0]*7
  for i in range(len(s)):
    for j in range(len(s[0])):
      if (i==0 and j!=4):
        sm[j] += s[i][j] 
      if (i==1 and j!=1):
        sm[j] += s[i][j] 
      if (i==2):
        sm[j] += s[i][j] 

  for j in range(len(s[0])):
    print sm[j]
  sp1 = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp1.setSize(650,550)
  sp1.setVLabel("Depth (km)")
  sp1.setHLabel("Log index")
  sp1.addColorBar("Values")
  sp1.plotPanel.setColorBarWidthMinimum(90)
  pv = sp1.addPixels(Sampling(7,1,1),Sampling(3,1,1),ll)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cjet)
  sp1 = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp1.setSize(650,550)
  sp1.setVLabel("Depth (km)")
  sp1.setHLabel("Log index")
  sp1.addColorBar("Shifts (m)")
  sp1.plotPanel.setColorBarWidthMinimum(90)
  pv = sp1.addPixels(Sampling(7,1,1),Sampling(3,1,1),s)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cjet)
  sp1 = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp1.setSize(650,550)
  sp1.setVLabel("Depth (km)")
  sp1.setHLabel("Log index")
  sp1.addColorBar("nulls (m)")
  sp1.plotPanel.setColorBarWidthMinimum(90)
  pv = sp1.addPixels(Sampling(7,1,1),Sampling(3,1,1),nls)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cjet)

  

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

def getXYLocation(set,type):
  logs = getLogs(set,type)
  nlog = len(logs)
  xlog = zerofloat(nlog)
  ylog = zerofloat(nlog)
  for i,log in enumerate(logs):
    xlog[i] = log.x2[0]
    ylog[i] = log.x3[0]
  return xlog,ylog

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

def writeFile(va):
  nz = len(va[0])
  nl = len(va)
  ofile = open('vlogs.txt','r+')
  for z in range(nz):
    ofile.write(str(z)+'\t'+str(va[0][z])+'\t'+str(va[1][z])+'\t'+
                str(va[2][z])+'\t'+str(va[3][z])+'\t'+str(va[4][z])+'\t'+
                str(va[5][z])+'\n')
  ofile.close()



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
