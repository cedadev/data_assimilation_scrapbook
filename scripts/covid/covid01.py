import rk
import numpy
import datetime, collections
import matplotlib.pyplot as plt
from utils import get_key_dates, Cases
import model

## python2.7

class BasicException(Exception):
  def __init__(self,msg="none"):
    self.msg = msg
    if msg != "none":
      print (msg) 

class DebugFlowControl(BasicException):
  pass

from scipy.signal import savgol_filter

## y0: numbers infected and not yet showing symtoms
## y1: numbers showing symptoms and infectious
## y2: numbers that are critical

class Addplots(object):
  def __init__(self,dy,y0=0.,filter=None,filter_data = None, centre="max",scale="max"):
    self.dy = dy
    self.y0 = y0
    self.filter = filter
    self.centre = centre
    self.scale = scale
    self.filter_data = filter_data
    self.filter_residues = dict()
    self._extracted = dict()

  def extractCaseData( self, country, c):
    ##x = [t[0] for t in c.data]
    ##y = [t[c.ix[country]] for t in c.data]
    x = c.dd["date"]
    y = c.dd[country]
    ##if country != "England":
      ##x = x[1:]
      ##y = [(y[i+1] - y[i]) for i in range(len(y)-1)]
    if country == 'China':
      yy = sum(y[:43])
      y[:43] = [1.25*z for z in y[:43]]
      y[43] = y[43] - .25*yy
    elif country == 'Spain':
      if abs(y[-1]) > 50000:
        y[-2] = c.data_mask_value
        y[-1] = c.data_mask_value
    
    y = numpy.ma.masked_equal([float(yyy) for yyy in y],float(c.data_mask_value))
    if country == 'Spain':
      print( "SPAIN, MASKED: ",y )

    ##y = [float(yyy) for yyy in y]
    self._extracted[country] = y

    if self.filter != None:
      if self.filter == "savgol":
        yhat = savgol_filter(y, self.filter_data["ww"], self.filter_data["np"]) # window size 51, polynomial order 3
      elif self.filter == "relax_linear":
        yhat = model.relax_linear_fit(y,self.filter_data["weight"],masked = country in ["England","Spain"])
      else:
        raise BasicException( "Bag filter name" )
      self.filter_residues[country] = [y[i] - yhat[i] for i in range(len(y))]
    else:
      yhat = y

    imx = numpy.argmax( yhat )

    if not hasattr( yhat, "mask") or type( yhat.mask ) == numpy.bool_:
      i99 = len(yhat)
    else:
      i99 = max( [i for i in range( imx, len(yhat) ) if not yhat.mask[i] ] )

    imxg = numpy.argmax( [yhat[i+7] - yhat[i] for i in range( len(yhat) - 7 )] ) - 3
    print ( "Date of peak: %s" % x[imx] )
    ioff = 0
    isc = -1
    if self.centre == "max":
      ioff = imx
    elif self.centre == "maxg":
      ioff = imxg
    if self.scale == "max":
      isc = imx
    elif self.scale == "maxg":
      isc = imxg

    xsc = [i-ioff for i in range(len(x))]
    ysc = [z/yhat[isc] + self.y0 for z in yhat]
    print ("YHAT:  ",yhat[imx-5:imx+5] )
    print ("YSC:  ",ysc[imx-5:imx+5] )
    if country == "England":
      print ("ENGLAND: %s, %s, %s" % (len(yhat),i99,type( yhat.mask )))
      print (ysc[imx:i99+1])
    return imx, xsc[:i99+1], ysc[:i99+1]

  def extractRatioData( self, country, c):
    x = [t[0] for t in c.data]
    xt, y = c.ratio(country)  
    x = x[xt[0]:xt[1]+1]

    if self.filter == "savgol":
      yhat = savgol_filter(y, self.filter_data["ww"], self.filter_data["np"]) # window size 51, polynomial order 3
    else:
      yhat = y
    ym = numpy.median( yhat )
    imx = numpy.argmin( yhat - ym )
    print ( "Ratio: Date of median: %s" % x[imx] )
    ##imx = int( len(ym )/2)
    xsc = [i-imx for i in range(len(x))]
    ysc = [z + self.y0 for z in yhat]
    return xt[0], imx, xsc, ysc

  def addplot_by_country(self, country, c, ax,color,filter=True, mode="cases",fit=None, annotate=None, arrow_right=False):

    try:
      if mode == "cases":
        ioff = 0
        imx, xsc, ysc = self.extractCaseData( country, c)
      else:
        ioff, imx, xsc, ysc = self.extractRatioData( country, c)
    except:
      print( "ERROR: no data for %s" % country )
      raise
      return

    ax.plot(xsc,ysc,color=color)
    ax.plot([-50.,75.],[self.y0,]*2,'--',color=color)
    ax.plot( xsc[-1], ysc[-1], "o", color=color )

    if arrow_right:
      ##plt.arrow(xsc[-3], ysc[-3], xsc[-1] - xsc[-3], ysc[-1] - ysc[-3], shape='full', lw=0, length_includes_head=True, head_width=.05)
      ax.annotate("", xy=(xsc[-1], ysc[-1]), xytext=(xsc[-3], ysc[-3]),
            arrowprops=dict(arrowstyle="->"))

    if c.key_dates != None and country in c.key_dates:
      for event in c.key_dates[country]:
        idx, event_type, description = event
        id = c.dd['date'].index( idx )
        mrkr = {'lockdown':'v', 'distancing':'o', 'initial':'s', 'event':"P"}.get( event_type, '-' )
        xd = xsc[id-ioff]
        yd = ysc[id-ioff]
        print( "Marker ",xd,yd,country )
        ax.plot( xd, yd, mrkr, color=color )
 
    bbox_props = dict(fc="cyan", ec="b", lw=1, boxstyle="round")
    if annotate != None:
      ia = annotate['xy_index']
      if ia < 0 or ia > len(xsc)-1:
        ia = len(xsc)/2

      if "annotate_y" in annotate:
        ytarg = annotate["annotate_y"]
        ii = [i for i in range(len(xsc)-1) if (ysc[i]-ytarg)*(ysc[i+1]-ytarg) < 0]
        if len(ii) == 0:
          pass
        else:
          ia = ii[-1]

      xy = (xsc[ia],ysc[ia])
      xytext = annotate['xytext']
      t = ax.annotate(country, ha="left", va="center", rotation=0,
            size=8, xy=xy, xytext=xytext, arrowprops=dict(facecolor=color),
            bbox=bbox_props)
    else:
      t = ax.text(-90., self.y0, country, ha="left", va="center", rotation=0,
            size=8,
            bbox=bbox_props)
    self.y0 += self.dy
    if imx > 6: 
      if (len(xsc)-imx) > 6:
         f3 = numpy.polyfit(xsc[imx-6:imx+6],ysc[imx-6:imx+6],3)
      else:
         f3 = numpy.polyfit(xsc[imx-6:],ysc[imx-6:],2)
      print (country, f3)
    return (xsc,ysc, imx)

mode = "fitfew"

if mode == "many":
  countries = ["New Zealand","China","South Korea","Thailand","South Africa","Greece","France","Germany","United States","Italy","Denmark","Sweden","Austria","Iceland","United Kingdom"]
  colors = ["purple","red","green","orange","pink","black","cyan","purple","red","green","orange","pink","black","cyan", "purple","red","green","orange"]
  title = '%s .. Savgoy-Gorsky[%s,%s]' % (country,ww,np)
  fntag = ""
elif mode == "fitfew":
  countries = ["China","New Zealand","Austria","Germany","France","Italy","Spain","Denmark","Greece","England"]
  ##countries = ["China","New Zealand","Austria","Germany","France","Italy","Spain","Greece"]
  ##countries = ["China","England"]
  wlrf = .10
  fntag = "_wlrf_b%3.3i" % (wlrf*100)
  title = 'Cases Reported,  linear relaxation filter[%s]' % (wlrf)
  colors = ["purple","red","green","pink","orange","blue","brown","magenta","cyan","black"]

a00 = 0.04
a01 = 0.04
a = 0.2
r = .25
a00 = a*(1-r)
a01 = a*r
a11 = a*(1-r)
a12 = a*r
g0 = 24.0
g1 = 0.00
scfac = 1.0
d = model.Demo01('test',g0,g1,a00,a01,a11,a12,scfac, 40,70, model_rank=4,fac4=1.0)
d.run()
d1 = model.Demo01('test',g0,2.1,a00,a01,a11,a12,scfac, 40,70, model_rank=4,fac4=1.0)
d1.run()
d9 = model.Demo01('test',g0,1.9,a00,a01,a11,a12,scfac, 40,70, model_rank=4,fac4=1.0)
d9.run()

xd = [t[0] for t in d.xx]
yd = [t[1][2] for t in d.xx]

imx = numpy.argmax( yd )
print ( "Index of peak: %s [lockdown + %s]" % (imx,imx-40 ) )
xdsc = [i-imx for i in xd]
ydsc = [z/yd[imx] for z in yd]
xld = xdsc[40]
yld = ydsc[40]
fd3 = numpy.polyfit(xdsc[imx-6:imx+6],ydsc[imx-6:imx+6],3)
print ('model', fd3)

xd9 = [t[0] for t in d9.xx]
yd9 = [t[1][2] for t in d9.xx]
imx = numpy.argmax( yd9 )
xdsc9 = [i-imx for i in xd9]
ydsc9 = [z/yd9[imx] for z in yd9]

c = Cases(targets=countries,add_england=True)
dths = Cases(targets=countries,file="total_deaths.csv")
c.setOther( dths )
c.setKeyDates( get_key_dates() )
c.analysis01()

withFilter=True
ww = 11
np = 1
withFilter &= np != 0

f, ax = plt.subplots(1)

ic =0

##adder = Addplots(0.1,0.1,filter="savgolxx",filter_data = {"np":1, "ww":3} )

adder = Addplots(0.0,0.0,filter="relax_linear",filter_data = {"weight":wlrf}, centre="max" )

##plt.ylim(-0.1,2.0)
annotate_dict = None
for country in countries:
  if mode == "fitfew":
    ix = 30 + ic*8
    annotate_dict = {"xytext":(80.,0.3 + ic*0.08), "xy_index":(ic+1)*20, "annotate_y":0.1+ic*0.07 }
  try:
    adder.addplot_by_country( country, c, ax,color=colors[ic],filter=withFilter, annotate=annotate_dict )
  except DebugFlowControl as dbg:
    print ("Caught exception for %s" % country)
  ic += 1

if mode == "many":
  ax.plot(xdsc,ydsc,color="blue")
  ax.plot(xld,yld, "v", color="blue" )
  ax.plot(xdsc9,[y+adder.y0 for y in ydsc9],color="blue")
else:
  if "Germany" in countries:
    ll = len( adder.filter_residues["Germany"] )
    total = [ sum( [v[i] for k, v in adder.filter_residues.items() if k != "Total"] ) for i in range(ll)]
    adder.filter_residues["Total"] = total

plt.savefig("covid%s.png" % fntag, bbox_inches='tight')

#title = 'test, gamma = 0.8 -> 0.1'
plt.ylabel(title)
plt.show()
