"""
Warping of well logs from Teapot Dome survey.
Author: Dave Hale, Colorado School of Mines
Version: 2013.12.28
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

from tp import *
from warpt import *

wlw = WellLogWarping()
curve = "v"
logs = None

def main(args):
  global logs
  logs = getLogs("d",curve)
  #goRGT()
  #goHorizon()
  goShifts()
  #goWarping()
  #goErrors()
  #goErrorsIJ()
  #goResample()
  #goSort()
  #goMesh()

def goRGT():
  sz,fs = resample(logs,curve)
  #fs = [fs[0],fs[4],fs[9],fs[14],fs[17],fs[20]] # deepest 6 velocity logs
  #fs = [fs[0],fs[4],fs[9],fs[11],fs[14],fs[17],fs[20]]
  #fs = [fs[ 1],fs[ 2],fs[ 3],fs[ 4],fs[ 7],
  #      fs[11],fs[21],fs[22],fs[33],fs[35],
  #      fs[43],fs[48],fs[50],fs[56],fs[66],
  #      fs[81],fs[88],fs[163]] # deepest 18 density logs
  nk,nl = len(fs[0]),len(fs)
  wlw.setPowError(0.25)
  wlw.setMaxShift(250)
  s = wlw.findShifts(fs)
  freplace = 2.0
  if curve=="d":
    freplace = 1.0
  fclips = (2.0,6.0)
  if curve=="d":
    fclips = (2.0,2.8)
  fs = wlw.replaceNulls(fs,freplace)
  tz = zerofloat(nk,nl)
  for il in range(nl):
    for ik in range(nk):
      #print ik, s[il][ik]
      tz[il][ik] = ik+s[il][ik] #FIX

  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(650,550)
  sp.setVLabel("Depth (km)")
  sp.setHLabel("Log index")
  sp.addColorBar("RGT")
  sp.plotPanel.setColorBarWidthMinimum(90)
  pv = sp.addPixels(sz,Sampling(nl),tz)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cjet)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(650,550)
  sp.setVLabel("Depth (km)")
  sp.setHLabel("Log index")
  sp.addColorBar("Velocity (km/s)")
  sp.plotPanel.setColorBarWidthMinimum(90)
  pv = sp.addPixels(sz,Sampling(nl),fs)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cjet)
  pv.setClips(fclips[0],fclips[1])

def goHorizon():
  sz,fs = resample(logs,curve)
  #fs = [fs[0],fs[4],fs[9],fs[14],fs[17],fs[20]] # deepest 6 velocity logs
  #fs = [fs[0],fs[4],fs[9],fs[11],fs[14],fs[17],fs[20]]
  #fs = [fs[ 1],fs[ 2],fs[ 3],fs[ 4],fs[ 7],
  #      fs[11],fs[21],fs[22],fs[33],fs[35],
  #      fs[43],fs[48],fs[50],fs[56],fs[66],
  #      fs[81],fs[88],fs[163]] # deepest 18 density logs
  nk,nl = len(fs[0]),len(fs)
  wlw.setPowError(0.25)
  wlw.setMaxShift(250)
  s = wlw.findShifts(fs)
  st = copy(s)
  wlw.invertShifts(st)
  """
  gs = wlw.applyShifts(fs,s)
  s = mul(1000*sz.delta,s) # convert shifts to m
  freplace = 2.0
  if curve=="d":
    freplace = 1.0
  fclips = (2.0,6.0)
  if curve=="d":
    fclips = (2.0,2.8)
  """

  tz = sz.delta*len(fs[0])/2
  zt = zerofloat(nl)
  ti = sz.indexOfNearest(tz)
  for i in range(nl):
    zt[i] = sz.getValue(int(ti - st[i][ti]))
    #zt[i] = ti - st[i][ti]
  x,y = getXYLocation("d",curve)
    

  """
  seis = readImage()
  nx = len(seis)
  ny = len(seis[0])
  nzz = len(seis[0][0])
  dz = 0.002
  dy = 0.025
  dx = 0.025
  sx = Sampling(nx,dx,0.0)
  sy = Sampling(ny,dx,0.0)
  szz = Sampling(nzz,dz,0.0)
  """
  ny = 50
  nx = 20
  dy = 10.0/ny
  dx = 5.0/nx
  sx = Sampling(nx,dx,0.0)
  sy = Sampling(ny,dy,0.0)

  #grd = SplinesGridder2(zt,y,x)
  #grd = DiscreteSibsonGridder2(zt,x,y)
  grd = SibsonGridder2(zt,x,y)
  zt = grd.grid(sy,sx)

  sp = SimplePlot()
  sp.setSize(700,380)
  sp.addColorBar("Depth (km)")
  pv = sp.addPixels(sy,sx,zt)
  pv.setColorModel(cjet)
  pv.setClips(0.87,1.00)
  ptv = sp.addPoints(x,y)
  ptv.setLineStyle(PointsView.Line.NONE)
  ptv.setMarkColor(Color.BLACK)
  ptv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  ptv.setMarkSize(12)
  cv = sp.addContours(sy,sx,zt)
  cv.setColorModel(ColorMap.GRAY)
  cv.setClips(0.87,1.00)
  #cv.setContours()
  """
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(650,550)
  sp.setVLabel("Depth (km)")
  sp.setHLabel("Log index")
  sp.addColorBar("Velocity (km/s)")
  sp.plotPanel.setColorBarWidthMinimum(90)
  pv = sp.addPixels(sz,Sampling(nl),fs)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cjet)
  pv.setClips(fclips[0],fclips[1])
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(650,550)
  sp.setVLabel("Depth (km)")
  sp.setHLabel("Log index")
  sp.addColorBar("Relative Geologic Time")
  sp.plotPanel.setColorBarWidthMinimum(90)
  pv = sp.addPixels(sz,Sampling(nl),gs)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cjet)
  pv.setClips(fclips[0],fclips[1])
  """

  """
  sinc = SincInterp()
  fx,fy = sx.first,sy.first
  horizon = zerofloat(sy.count,sx.count)
  for ix in range(nx):
    for iy in range(ny):
      zi= zt[ix][iy]
      xi = fx+ix*dx
      yi= fy+iy*dy
      horizon[ix][iy] = sinc.interpolate(szz,sy,sx,seis,zi,yi,xi)
  SimplePlot.asPixels(horizon)
  """


def goShifts():
  sz,fs = resample(logs,curve)
  #fs = [fs[0],fs[4]]
  #fs = [fs[ 2],fs[ 3],fs[ 6],fs[ 7],fs[ 8],fs[10],fs[12],
  #      fs[16],fs[18],fs[21],fs[22],fs[23],fs[24],
  #      fs[25],fs[26],fs[27],fs[28],fs[29]] # shortest velocity logs
  fs = [fs[0],fs[4],fs[9],fs[14],fs[17],fs[20]] # deepest 6 velocity logs
  #fs = [fs[0],fs[9]]
  #fs = [fs[0],fs[4],fs[9],fs[11],fs[14],fs[17],fs[20]]
  #fs = [fs[ 1],fs[ 2],fs[ 3],fs[ 4],fs[ 7],
  #      fs[11],fs[21],fs[22],fs[33],fs[35],
  #      fs[43],fs[48],fs[50],fs[56],fs[66],
  #      fs[81],fs[88],fs[163]] # deepest 18 density logs
  #fs = [fs[ 6],fs[16],fs[18],fs[31],fs[32],fs[33],fs[36],
  #      fs[38],fs[48],fs[55],fs[86]] # deepest 11 gamma logs
  #fs = [fs[16],fs[18],fs[19],fs[28],fs[33],
  #      fs[34],fs[35],fs[37],fs[38],fs[39],
  #      fs[45],fs[50],fs[68]] # deepest 13 porosity logs
  nk,nl = len(fs[0]),len(fs)
  """
  wlw.setPowError(0.25)
  wlw.setMaxShift(250)
  s = wlw.findShifts(fs)
  gs = wlw.applyShifts(fs,s)
  s = mul(1000*sz.delta,s) # convert shifts to m
  """
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
  """
  gs = wlw.replaceNulls(gs,freplace)
  """
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  #sp.setSize(650,550)
  sp.setSize(700,380) # for slides
  sp.setSize(448,874) # for slides - each log
  #sp.setSize(617,376)
  sp.setVLabel("Depth (km)")
  sp.setHLabel("Log index")
  sp.addColorBar(cblabel)
  sp.setVLimits(0.892,1.142)
  sp.plotPanel.setColorBarWidthMinimum(90)
  pv = sp.addPixels(sz,Sampling(nl,1,1),fs)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cjet)
  pv.setClips(fclips[0],fclips[1])
  sp.setFontSizeForSlide(9.5/11.0,1.0,16.0/9.0)
  #sp.setFontSizeForPrint(8.0,255.6)
  sp.paintToPng(720.0,3.55,"original.png")
  #sp.setFontSizeForPrint(8.0,234.5)
  #sp.paintToPng(720.0,3.25,"original.png")
  """
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  #sp.setSize(650,508) # zoom & werror
  sp.setSize(700,380) # for slides
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
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(cjet)
  #sp1.setFontSizeForPrint(8.0,156.3)
  #sp1.paintToPng(720.0,6.51,"shifts.png")
  """

  
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

def goWarping():
  #pairs = [(0,4),(4,9),(9,14),(14,17),(17,20)] # deepest 6 velocity logs
  pairs = [(0,4),(0,14)] #slides
  #pairs = [(0,4),(4,20),(9,14)] # paper
  #pairs = [(0,4),(4,9),(17,20)] 
  #pairs = [(4,0),(4,9),(4,14),(4,17),(4,20)] # log 4 is nearest to centroid
  #pairs = [(4,14)]
  sz,fs = resample(logs,curve)
  wlw.setPowError(0.25)
  #wlw.setPowError(2.00)
  wlw.setMaxShift(250)
  freplace = -2.0
  #freplace = 2.0
  vlabel = "Velocity (km/s)"
  if curve=="d":
    freplace = 1.0
    vlabel = "Density (g/cc)"
  if curve=="p":
    freplace = 0.0
    vlabel = "Porosity"
  if curve=="g":
    freplace = 30.0
    vlabel = "Gamma ray (API)"
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
    for ij in range(len(ff)):
      if ff[ij] == -2.0:
        ff[ij] = 3.0
    for ij in range(len(gg)):
      if gg[ij] == -2.0:
        gg[ij] = 3.0
    sff = Sampling(len(ff),sz.getDelta(),sz.getFirst() +ii*sz.getDelta())
    sgg = Sampling(len(gg),sz.getDelta(),sz.getFirst() +jj*sz.getDelta())
    title = "("+str(jf)+","+str(jg)+")"
    sk = Sampling(2*sz.count-1,0.5*sz.delta,sz.first)
    #pp = PlotPanel(2,1)
    pp = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
    #pp.setVLimits(0,2,6.5)
    #pp.setVLimits(1,2,6.5)
    #pp.setTitle(title)
    #pp.setVLabel("Relative geologic time")
    pp.setVLabel("Depth (km)")
    pp.setHLabel("Velocity (km/s)")
    #pp.setHLimits(3,4.5)
    pp.setHLimits(2,6.4)
    pp.setVLimits(0,1.9)
    #pp.setVLimits(0.892,1.142)
    #pp.setVLimits(1.005,1.025)
    #pp.setHLabel("Depth (km)")
    #pp.setVLabel(0,"Velocity (km/s)")
    #pp.setVLabel(1,"Velocity (km/s)")
    #pp.setVLimits(0,2,6.4)
    #pp.setVLimits(1,2,6.4)
    if True:
      """
      sp = SimplePlot()
      sp.setSize(750,500)
      sp.setTitle(title)
      sp.setHLabel("Depth (km)")
      sp.setVLabel(vlabel)
      sp.setVLimits(2,6.5)
      pv = sp.addPoints(sz,fi)
      pv.setLineColor(Color.BLACK)
      pv = sp.addPoints(sz,gj)
      pv.setLineColor(Color.RED)
      """
      """
      pv = pp.addPoints(0,0,sz,fi)
      pv.setLineColor(Color.BLACK)
      #pv = pp.addPoints(0,0,sz,gj)
      #pv.setLineColor(Color.RED)
      pv.setLineWidth(2.0)
      """
      pv = pp.addPoints(0,0,sff,ff) ## to start plotting at first depth
      pv.setLineColor(Color.BLACK)
      #pv = pp.addPoints(0,0,sgg,gg)
      #pv.setLineColor(Color.RED)
      pv.setLineWidth(2.0)
    """
    if True:
      sp = SimplePlot()
      sp.setSize(750,500)
      sp.setTitle(title)
      sp.setHLabel("Depth (km)")
      sp.setVLabel(vlabel)
      sp.setVLimits(2,6.5)
      pv = sp.addPoints(sk,fk)
      pv.setLineColor(Color.BLACK)
      pv = sp.addPoints(sk,gk)
      pv.setLineColor(Color.RED)

      #pv = pp.addPoints(0,0,sk,fk)
      #pv.setLineColor(Color.BLACK)
      pv = pp.addPoints(0,0,sk,gk)
      pv.setLineColor(Color.RED)
      pv.setLineWidth(2.0)
    """

    #pp = PlotPanel(2,1,PlotPanel.Orientation.X1DOWN_X2RIGHT,PlotPanel.AxesPlacement.LEFT_BOTTOM)
    #pp = PlotPanel(1,2,PlotPanel.Orientation.X1DOWN_X2RIGHT)
    pf = PlotFrame(pp)
    pf.setSize(664,700) # for slides
    #pf.setSize(569,874) # for slides superzoom
    #pf.setSize(700,418) # for slides
    #pf.setSize(658,628)
    pf.setVisible(True)
    pf.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE)
    pf.setFontSizeForSlide(9.5/11.0,15.0/19.0,16.0/9.0) #big
    #pf.setFontSizeForSlide(0.4,9.5/11.0,16.0/9.0) # superzoom
    #pf.setFontSizeForPrint(8.0,156.3)
    #pf.paintToPng(720.0,2.17,"warps"+str(c)+".png")

def goErrorsIJ():
  #pairs = [(0,4),(4,9),(9,14),(14,17),(17,20)]
  #pairs = [(4,0),(4,9),(4,14),(4,17),(4,20)] # log 4 is nearest to centroid
  #pairs = [(9,15),(15,28),(9,28)]
  #pairs = [(9,14)]
  #pairs = [(0,9)]
  #pairs = [(14,17)]
  pairs = [(0,4)]
  #pairs = [(0,14)]
  sz,f = resample(logs,curve)
  #wlw.setPowError(2.0)
  #wlw.setPowError(1.0)
  #wlw.setPowError(0.5)
  wlw.setPowError(0.25)
  wlw.setMaxShift(250)
  #scale = 40
  #scale =6
  #scale = 2
  #scale = 1.25
  scale = 1.0
  for pair in pairs:
    ia,ib = pair[0],pair[1]
    e = wlw.computeErrorsIJ(f[ia],f[ib])
    ni,nj = len(e[0]),len(e)
    maxerr = max(e)
    #maxerr = -100000000
    norme = zerofloat(nj,ni)
    """
    for i in range(892,1142):
      for j in range(892,1142):
        if e[i][j] > maxerr:
          maxerr = e[i][j]
    for i in range(892,1142):
      for j in range(892,1142):
    """
    for i in range(ni):
      for j in range(nj):
        norme[i][j] = scale*e[i][j]/maxerr
    """
    ee = zerofloat(nj,ni)
    ee[1005][1005] = norme[1005][1005]
    #ee[1005][1006] = norme[1005][1006]
    #ee[1005][1007] = norme[1005][1007]
    #ee[1005][1008] = norme[1005][1008]
    #ee[1005][1009] = norme[1005][1009]
    """
    si = Sampling(ni,1,0)
    sj = Sampling(nj,1,0)
    title = "("+str(ia)+","+str(ib)+")"
    sp = SimplePlot()
    #sp.setSize(700,700)
    sp.setSize(910,700)
    #sp.setTitle(title)
    sp.setVLabel("Log 2 depth index j") ##### make sure to change log #
    sp.setHLabel("Log 1 depth index i") ##### make sure to change log #
    sp.setHFormat("%5f")
    #sp.setHLimits(1335,1350) # for print
    #sp.setVLimits(20,32) # for print
    #sp.setHLimits(1200,2100) #########
    #sp.setVLimits(-250,-75) #########
    #sp.setVLimits(-105,-95) #stencil
    #sp.setHLimits(2725,2740) #stencil
    sp.setVLimits(892,1142)
    sp.setHLimits(892,1142)
    #sp.setVLimits(1005,1025)
    #sp.setHLimits(1005,1025)
    sp.addColorBar("e[i,j]")
    #pv = sp.addPixels(sj,si,transpose(e))
    pv = sp.addPixels(sj,si,transpose(norme))
    #pv = sp.addPixels(sj,si,transpose(ee))
    #pv.setColorModel(cjet)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    #pv.setPercentiles(0,80)
    pv.setClips(0,1)
    sp.setFontSizeForSlide(0.6,9.5/11.0,16.0/9.0)
    #sp.setFontSizeForPrint(8.0,222.0)
    #sp.paintToPng(720.0,3.08,"ijaezoom.png")
    #sp.paintToPng(720.0,3.08,"ijaepath200.png")

    d = wlw.accumulateErrorsIJ(norme)
    sp = SimplePlot()
    sp.setSize(910,700)
    #sp.setTitle(title)
    sp.setVLabel("Log 2 depth index j") ##### make sure to change log #
    sp.setHLabel("Log 1 depth index i") ##### make sure to change log #
    sp.setHFormat("%5f")
    sp.setVLimits(892,1142)
    sp.setHLimits(892,1142)
    sp.addColorBar("Accumulated e[i,j]")
    pv = sp.addPixels(si,sj,transpose(d))
    #pv.setColorModel(cjet)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    #pv.setPercentiles(0,85)
    pv.setClips(10,110)
    sp.setFontSizeForSlide(0.6,9.5/11.0,16.0/9.0)
    sp.paintToPng(720.0,3.08,"aaezoom.png")
    #iw,jw = wlw.findWarping(d)
    #iw = wlw.toFloat(iw)
    #jw = wlw.toFloat(jw)
    #pv = sp.addPoints(iw,jw)
    #pv.setLineColor(Color.WHITE)

def goErrors():
  #pairs = [(0,4),(4,9),(9,14),(14,17),(17,20)]
  #pairs = [(4,0),(4,9),(4,14),(4,17),(4,20)] # log 4 is nearest to centroid
  #pairs = [(9,15),(15,28),(9,28)]
  #pairs = [(9,14)]
  #pairs = [(0,9)]
  #pairs = [(14,17)]
  pairs = [(0,4)]
  sz,f = resample(logs,curve)
  wlw.setPowError(0.25)
  wlw.setMaxShift(250)
  for pair in pairs:
    ia,ib = pair[0],pair[1]
    e = wlw.computeErrors(f[ia],f[ib])
    wlw.interpolateOddErrors(e)
    nl,nk = len(e[0]),len(e)
    lmax = (nl-1)/2
    lmin = -lmax
    maxerr = max(e)
    #maxerr = -100000000
    norme = zerofloat(nl,nk)
    """
    for k in range(1900,2100):
      for l in range(-120,60):
        if e[k][l] > maxerr:
          maxerr = e[k][l]
          print maxerr
    for k in range(1900,2100):
      for l in range(-120,60):
    """
    for k in range(nk):
      for l in range(nl):
        norme[k][l] = e[k][l]/maxerr
    sl = Sampling(nl,1,lmin)
    sk = Sampling(nk,1,0)
    title = "("+str(ia)+","+str(ib)+")"
    sp = SimplePlot()
    sp.setSize(780,555) #zoom
    #sp.setSize(910,700) #big
    #sp.setSize(750,518)
    #sp.setTitle(title)
    sp.setHLabel("Depth index k")
    sp.setVLabel("Lag index l")
    sp.setHFormat("%5f")
    #sp.setHLimits(1335,1350) # for print
    #sp.setVLimits(20,32) # for print
    #sp.setHLimits(1200,2100) #########
    #sp.setVLimits(-250,-75) #########
    #sp.setVLimits(-105,-95) #stencil
    #sp.setHLimits(2725,2740) #stencil
    #sp.setVLimits(-250,250)
    #sp.setHLimits(1784,2284)
    sp.setVLimits(-120,60)
    sp.setHLimits(1900,2100)
    sp.addColorBar("e[k,l]")
    pv = sp.addPixels(sk,sl,transpose(norme))
    #pv = sp.addPixels(sk,sl,transpose(e))
    #pv.setColorModel(cjet)
    pv.setColorModel(ColorMap.GRAY)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setPercentiles(0,100)
    sp.setFontSizeForSlide(0.6,9.5/11.0,16.0/9.0)
    #sp.setFontSizeForPrint(8.0,222.0)
    #sp.paintToPng(720.0,3.08,"stencil.png")
    sp1 = SimplePlot()
    sp1.setSize(780,555)
    #sp1.setSize(910,700)
    #sp1.setSize(750,500)
    #sp1.setTitle(title)
    sp1.setHLabel("Depth index k")
    sp1.setVLabel("Lag index l")
    sp1.setHFormat("%5f")
    #sp1.setHLimits(1200,2100) ##########
    #sp1.setVLimits(-250,-75) ##########
    sp1.setVLimits(-120,60)
    sp1.setHLimits(1900,2100)
    #sp1.setVLimits(-250,250)
    #sp1.setHLimits(1784,2284)
    sp1.addColorBar("e[k,l]")
    pv = sp1.addPixels(sk,sl,transpose(norme))
    #pv = sp1.addPixels(sk,sl,transpose(e))
    #pv.setColorModel(cjet)
    pv.setColorModel(ColorMap.GRAY)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setPercentiles(0,100)
    d = wlw.accumulateErrors(e)
    wlw.interpolateOddErrors(d)
    kw,lw = wlw.findWarping(d)
    kw = wlw.toFloat(kw)
    lw = wlw.toFloat(lw)
    pv1 = sp1.addPoints(kw,lw)
    pv1.setLineColor(Color.RED)
    pv1.setLineStyle(PointsView.Line.DASH)
    pv1.setLineWidth(4.0)
    sp1.setFontSizeForSlide(0.6,9.5/11.0,16.0/9.0)
    sp1.paintToPng(720.0,3.08,"klaepath.png")
    sp = SimplePlot()
    sp.setSize(750,500)
    sp.setTitle(title)
    sp.setHLabel("Depth index k")
    sp.setVLabel("Lag index l")
    sp.setHFormat("%5f")
    pv = sp.addPixels(sk,sl,transpose(d))
    pv.setColorModel(cjet)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setPercentiles(0,90)
    kw,lw = wlw.findWarping(d)
    kw = wlw.toFloat(kw)
    lw = wlw.toFloat(lw)
    pv = sp.addPoints(kw,lw)
    pv.setLineColor(Color.WHITE)

    #pp = PlotPanel(2,1,PlotPanel.Orientation.X1DOWN_X2RIGHT,PlotPanel.AxesPlacement.LEFT_BOTTOM)
    #pp = PlotPanel(1,2,PlotPanel.Orientation.X1DOWN_X2RIGHT)
    pp = PlotPanel(2,1)
    pp.setSize(750,500)
    pp.setVLabel(0,"Lag index l")
    pp.setVLabel(1,"Lag index l")
    pp.setHLabel("Depth index k")
    pp.setHFormat("%5f")
    pp.addColorBar("e[k,l]")
    pp.setHLimits(2550,2850) 
    pp.setVLimits(0,-130,-70) 
    pp.setVLimits(1,-130,-70) 
    pv = pp.addPixels(0,0,sk,sl,transpose(e))
    #pv.setColorModel(cjet)
    pv.setColorModel(ColorMap.GRAY)
    #pv.setClips(0.0,0.5)
    #pv.setClips(0,1)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setPercentiles(0,100)
    pv = pp.addPixels(1,0,sk,sl,transpose(e))
    #pv.setColorModel(cjet)
    pv.setColorModel(ColorMap.GRAY)
    #pv.setClips(0.2,0.7)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv.setPercentiles(0,100)
    pv1 = pp.addPoints(1,0,kw,lw)
    pv1.setLineColor(Color.WHITE)
    pf = PlotFrame(pp)
    pf.setVisible(True)
    pf.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE)
    #pf.setFontSizeForPrint(8.0,222.0)
    #pf.paintToPng(720.0,3.08,"ae.png")

def goResample():
  nlog = len(logs)
  sz,f = resample(logs,curve)
  f = wlw.replaceNulls(f,-0.01)
  c = 0
  for j in range(len(f[0])):
    if f[0][j] != -0.01:
      c +=1
  for j in range(len(f[4])):
    if f[4][j] != -0.01:
      c +=1
  for j in range(len(f[9])):
    if f[9][j] != -0.01:
      c +=1
  for j in range(len(f[14])):
    if f[14][j] != -0.01:
      c +=1
  for j in range(len(f[17])):
    if f[17][j] != -0.01:
      c +=1
  for j in range(len(f[20])):
    if f[20][j] != -0.01:
      c +=1
  print c
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setVLabel("Depth (km)")
  sp.setHLabel("Log index")
  sp.addColorBar() #("Velocity (km/s)")
  pv = sp.addPixels(sz,Sampling(nlog),f)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ajet)
  if curve=="v":
    pv.setClips(2.0,6.0)
  elif curve=="d":
    pv.setClips(2.0,2.8)
  elif curve=="p":
    pv.setClips(0.0,0.5)
  elif curve=="g":
    pv.setClips(50,250)

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
  #sp.setLimits(0.0,0.0,9.0,4.0)
  #pv = sp.addPoints(x,y)
  #pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  #pv.setMarkSize(6)
  if curve=="g":
    x = [x[ 6],x[16],x[18],x[31],x[32],x[33],x[36],
        x[38],x[48],x[55],x[86]] # deepest 11 gamma logs
    y = [y[ 6],y[16],y[18],y[31],y[32],y[33],y[36],
        y[38],y[48],y[55],y[86]] # deepest 11 gamma logs
    pv = sp.addPoints(x,y)
    pv.setLineStyle(PointsView.Line.NONE)
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    pv.setMarkSize(10)
    pv.setMarkColor(Color.RED)
  if curve=="p":
    x = [x[16],x[18],x[19],x[28],x[33],
          x[34],x[35],x[37],x[38],x[39],
          x[45],x[50],x[68]] # deepest 13 porosity logs
    y = [y[16],y[18],y[19],y[28],y[33],
          y[34],y[35],y[37],y[38],y[39],
          y[45],y[50],y[68]] # deepest 13 porosity logs
    pv = sp.addPoints(x,y)
    pv.setLineStyle(PointsView.Line.SOLID)
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    pv.setMarkSize(10)
    pv.setMarkColor(Color.RED)
  if curve=="v":
    """
    x = [x[0],x[4],x[9],x[14],x[17],x[20]]
    y = [y[0],y[4],y[9],y[14],y[17],y[20]]
    pv = sp.addPoints(x,y)
    pv.setLineStyle(PointsView.Line.NONE)
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    pv.setMarkSize(10)
    pv.setMarkColor(Color.RED)
    """
    x1 = [x[0],x[4]]
    y1 = [y[0],y[4]]
    pv = sp.addPoints(x1,y1)
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    pv.setMarkColor(Color.RED)
    pv.setLineColor(Color.RED)
    pv.setLineWidth(2.0)
    pv.setMarkSize(10)
    x1 = [x[0],x[9]]
    y1 = [y[0],y[9]]
    pv = sp.addPoints(x1,y1)
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    pv.setMarkColor(Color.RED)
    pv.setLineColor(Color.RED)
    pv.setLineWidth(2.0)
    pv.setMarkSize(10)
    x1 = [x[0],x[14]]
    y1 = [y[0],y[14]]
    pv = sp.addPoints(x1,y1)
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    pv.setMarkColor(Color.RED)
    pv.setLineColor(Color.RED)
    pv.setLineWidth(2.0)
    pv.setMarkSize(10)
    x1 = [x[0],x[17]]
    y1 = [y[0],y[17]]
    pv = sp.addPoints(x1,y1)
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    pv.setMarkColor(Color.RED)
    pv.setLineColor(Color.RED)
    pv.setLineWidth(2.0)
    pv.setMarkSize(10)
    x1 = [x[0],x[20]]
    y1 = [y[0],y[20]]
    pv = sp.addPoints(x1,y1)
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    pv.setMarkColor(Color.RED)
    pv.setLineColor(Color.RED)
    pv.setLineWidth(2.0)
    pv.setMarkSize(10)
    x1 = [x[4],x[9]]
    y1 = [y[4],y[9]]
    pv = sp.addPoints(x1,y1)
    pv.setLineColor(Color.RED)
    pv.setLineWidth(2.0)
    x1 = [x[4],x[14]]
    y1 = [y[4],y[14]]
    pv = sp.addPoints(x1,y1)
    pv.setLineColor(Color.RED)
    pv.setLineWidth(2.0)
    x1 = [x[4],x[17]]
    y1 = [y[4],y[17]]
    pv = sp.addPoints(x1,y1)
    pv.setLineColor(Color.RED)
    pv.setLineWidth(2.0)
    x1 = [x[4],x[20]]
    y1 = [y[4],y[20]]
    pv = sp.addPoints(x1,y1)
    pv.setLineColor(Color.RED)
    pv.setLineWidth(2.0)
    x1 = [x[9],x[14]]
    y1 = [y[9],y[14]]
    pv = sp.addPoints(x1,y1)
    pv.setLineColor(Color.RED)
    pv.setLineWidth(2.0)
    x1 = [x[9],x[17]]
    y1 = [y[9],y[17]]
    pv = sp.addPoints(x1,y1)
    pv.setLineColor(Color.RED)
    pv.setLineWidth(2.0)
    x1 = [x[9],x[20]]
    y1 = [y[9],y[20]]
    pv = sp.addPoints(x1,y1)
    pv.setLineColor(Color.RED)
    pv.setLineWidth(2.0)
    x1 = [x[14],x[17]]
    y1 = [y[14],y[17]]
    pv = sp.addPoints(x1,y1)
    pv.setLineColor(Color.RED)
    pv.setLineWidth(2.0)
    x1 = [x[14],x[20]]
    y1 = [y[14],y[20]]
    pv = sp.addPoints(x1,y1)
    pv.setLineColor(Color.RED)
    pv.setLineWidth(2.0)
    x1 = [x[17],x[20]]
    y1 = [y[17],y[20]]
    pv = sp.addPoints(x1,y1)
    pv.setLineColor(Color.RED)
    pv.setLineWidth(2.0)

  sp.setFontSizeForSlide(9.5/11.0,1.0,16.0/9.0)
  #sp.setFontSizeForPrint(8.0,222.0)
  sp.paintToPng(720.0,3.08,"sortedmesh.png")

def goMesh():
  mesh = TriMesh()
  for i,log in enumerate(logs):
    node = TriMesh.Node(log.x2[0],log.x3[0])
    node.index = i
    mesh.addNode(node)
  sp = SimplePlot()
  sp.setSize(700,380)
  sp.setHLabel("Crossline (km)")
  sp.setVLabel("Inline (km)")
  sp.setLimits(0.0,0.0,9.0,4.0)
  tmv = TriMeshView(mesh)
  tmv.setLineColor(Color.BLACK)
  tmv.setMarkColor(Color.BLACK)
  tmv.setTriBoundsVisible(True)
  sp.add(tmv)
  #sp.setFontSizeForPrint(8.0,222.0)
  #sp.paintToPng(720.0,3.08,"mesh.png")

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
