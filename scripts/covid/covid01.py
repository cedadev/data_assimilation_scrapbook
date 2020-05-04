import rk
import numpy
import datetime, collections
import matplotlib.pyplot as plt


## python2.7

from scipy.signal import savgol_filter

## y0: numbers infected and not yet showing symtoms
## y1: numbers showing symptoms and infectious
## y2: numbers that are critical

def mint(s):
  if s == "":
    return 0
  else:
    return int(s)

class DateConv(object):
  def __init__(self,ref="2020-01-01"):
      self.oref = self._ord(ref)
      self.ref = ref

  def _ord(self,ds):
      year, month, day = [int(x) for x in ds.split( "-" )]
      date = datetime.date( year, month, day )
      return date.toordinal()

  def convert(self,ds):
    return self._ord(ds) - self.oref

dconv = DateConv()

def get_key_dates():
    """
     -> "A collections.defaultdict(set) instance containing key dates":
       Record key dates.
    """
  ##def __init__(self):
    file = "key_dates.csv"
    ii=open(file)
    results = collections.defaultdict(set)
    for l in ii.readlines():
      parts = [x.strip('"') for x in l.strip().split("\t")]
      country, sdate, key = parts[:3]
      if len(parts) == 4:
         description= parts[3]
      else:
         description = ""
      jdate = dconv.convert( sdate )
      results[country].add( (jdate,key,description) )
    return results

class UKCases(object):
  def __init__(self,file="total_cases.csv",targets=["Reading",]):
#Area name,Area code,Area type,Specimen date,Daily lab-confirmed cases,Previously reported daily cases,Change in daily cases,Cumulative lab-confirmed cases,Previously reported cumulative cases,Change in cumulative cases
    self.other = None

class Cases(object):
  def __init__(self,file="total_cases.csv",targets=["Germany",]):
    self.other = None
    self.key_dates = None
    date_ref = datetime.date( 2020, 1, 1 )
    ii=open(file)
    self.headings = [x.strip('"') for x in ii.readline().strip().split(",")]
    self.data = []
    self.dd = {}
    for l in ii.readlines():
      this = [x.strip('"') for x in l.strip().split(",")]
      this[1:] = [mint(x) for x in this[1:]]
      ##year, month, day = [int(x) for x in this[0].split( "-" )]
      ##date = datetime.date( year, month, day )
      this[0] = dconv.convert( this[0] )
      self.data.append( this )
    ##self.ix = {x:self.headings.index(x) for x in targets}
    self.ix = {x:self.headings.index(x) for x in self.headings}

    self.dd['date'] = [t[0] for t in self.data]
    for h in self.headings[1:]:
      y = [t[self.ix[h]] for t in self.data]
      ix = len(y) - 1
      for i in range( len(y) -1 ):
        y[ix] = y[ix] - y[ix-1]
        if y[ix] < 0:
          print ( "%s:  negative value: %s" % (h,y[ix]) )
        ix += -1
      self.dd[h] = y

  def setKeyDates(self,kd):
    self.key_dates = kd

  def setOther(self,other):
    self.other = other

  def ratio(self,country,sc=10,other=None):
    """ratio of two series ... e.g. cases/deaths"""
    if other != None:
      self.other = other
    assert self.other != None
    ts00 = self.dd[country]
    if filter:
      ts0 = savgol_filter(ts00, ww, np) # window size 51, polynomial order 3
    else:
      ts0 = ts00
    imx = numpy.argmax( ts0 )
    ix0 = imx
    ix9 = imx
    for i in range(imx):
      if ts0[imx-i] < sc:
        break
      ix0 = imx-i
    for i in range(len(ts0)-imx):
      if ts0[imx+i] < sc:
        break
      ix9 = imx+i
    ratio = [self.other.dd[country][x]/float(ts0[x]) for x in range(ix0,ix9+1)]
    print( '%s: %s -- %s' % (country,ix0,ix9) )
    print(  " ".join( ["%4.2f" % x for x in ratio] ) )
    return (ix0,ix9), ratio
      
    
  def analysis01(self):
    
    for ir in range(1,len(self.headings)):
      this = [t[ir] for t in self.data]
      imx = numpy.argmax( this )
      if imx > 10:
        f3 = numpy.polyfit(list( range(-10,1)),this[imx-10:imx+1],2)
        if this[imx] > 1000 and abs(f3[1]) < 20*abs(f3[0]):
          print( self.headings[ir], this[imx], f3 )

class Demo01(object):
  def __init__(self,name,g0,g1,a00, a01,a11,a12, scfac,tc,t9,model_rank=3,fac4=1.):
    self.g0=g0
    self.name=name
    self.tc = tc
    self.t9 = t9
    self.model_rank = model_rank
    if model_rank == 3:
      self.z0 = scfac*numpy.array( [[-a00-a01, g0, 0.],[a01,-.1,0.],[0.,.05,-.2]], dtype='float')
      self.z1 = scfac*numpy.array( [[-a00-a01, g1, 0.],[a01,-.1,0.],[0.,.05,-.2]], dtype='float') 
    else:
      self.z0 = scfac*numpy.array( [[(-a00-a01)*fac4, 0., g0, 0.], [a01*fac4, -a00-a01,0.,0.], [0.,a01,-(a11+a12),0.],[0.,0.,a12,-.2]], dtype='float')
      self.z1 = scfac*numpy.array( [[(-a00-a01)*fac4, 0., g1, 0.], [a01*fac4, -a00-a01,0.,0.], [0.,a01,-(a11+a12),0.],[0.,0.,a12,-.2]], dtype='float') 



    self.w0, self.v0 = numpy.linalg.eig( self.z0 )
    self.w1, self.v1 = numpy.linalg.eig( self.z1 )
    self.f0 = lambda x,y: self.z0.dot(y)
    self.f1 = lambda x,y: self.z1.dot(y)
    self.zz0 = self.z0.transpose().dot( self.z0 )
    self.zz1 = self.z1.transpose().dot( self.z1 )
    self.xx = []

  def fit(self,ioff,obs,w):
    """Fitting (d^2/dt^2 - z^t z) res = w(res-obs)"""
    res = obs[:]
    m0 = 0
    for x in obs:
      m0 += x**2
    m0 = numpy.sqrt( m0/len( obs ) )
    for k in range(50):
      f = [0,]*len(obs)
      s = 0.
      for i in range( 1,len(obs)-1 ):
        s += (res[i]-obs[i])**2 + (res[i-1] - 2.*res[i] + res[i+1] - res[i])**2
        f[i] = res[i-1] - 2.*res[i] + res[i+1] - res[i] - w*(res[i]-obs[i])
      s = numpy.sqrt( s/len( obs ) )
      print (k,s, s/m0, f[25:28] )
      if s < m0*0.001:
        return res
      for i in range( 1,len(obs)-1 ):
        res[i] = res[i] + f[i]/6.
    return res
      
    
  def run(self):
    x = 0.
    h = 0.1
    if self.model_rank == 3:
      y = numpy.array( [1.,0.,0.], dtype='float')
    else:
      y = numpy.array( [1.,0.,0.,0.], dtype='float')

    self.xx.append( (x,y) )
    for k0 in range(self.tc):
      for k1 in range(10):
        x,y = rk.rk4d(self.f0,x,y,h)
      self.xx.append( (x,y) )
      ##print ('%5.2f:: %s' % (x,'%5i, %5i, %5i' % tuple( y ) ) )


    for k0 in range(self.t9):
      for k1 in range(10):
        x,y = rk.rk4d(self.f1,x,y,h)
      self.xx.append( (x,y) )
      ##print ('%5.2f:: %s' % (x,'%5i, %5i, %5i' % tuple( y ) ) )
    self.x = x
    self.y = y
  

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
d = Demo01('test',g0,g1,a00,a01,a11,a12,scfac, 40,70, model_rank=4,fac4=1.0)
d.run()
d1 = Demo01('test',g0,2.1,a00,a01,a11,a12,scfac, 40,70, model_rank=4,fac4=1.0)
d1.run()
d9 = Demo01('test',g0,1.9,a00,a01,a11,a12,scfac, 40,70, model_rank=4,fac4=1.0)
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


countries = ["New Zealand","China","South Korea","Thailand","South Africa","Greece","France","Germany","United States","Italy","Denmark","Sweden","Austria","Iceland","United Kingdom"]
colors = ["purple","red","green","orange","pink","black","cyan","purple","red","green","orange","pink","black","cyan", "purple","red","green","orange"]
c = Cases(targets=countries)
dths = Cases(targets=countries,file="total_deaths.csv")
c.setOther( dths )
c.setKeyDates( get_key_dates() )
c.analysis01()

withFilter=True
ww = 11
np = 1
withFilter &= np != 0

class Addplots(object):
  def __init__(self,dy,y0=0.):
    self.dy = dy
    self.y0 = y0

  def extractCaseData( self, country, c):
    x = [t[0] for t in c.data]
    y = [t[c.ix[country]] for t in c.data]
    x = x[1:]
    y = [(y[i+1] - y[i]) for i in range(len(y)-1)]
    if country == 'China':
      yy = sum(y[:43])
      y[:43] = [1.25*z for z in y[:43]]
      y[43] = y[43] - .25*yy

  #y1 = [numpy.log(z) for z in y]
    if filter:
      yhat = savgol_filter(y, ww, np) # window size 51, polynomial order 3
    else:
      yhat = y
    imx = numpy.argmax( yhat )
    print ( "Date of peak: %s" % x[imx] )
    xsc = [i-imx for i in range(len(x))]
    ysc = [z/yhat[imx] + self.y0 for z in yhat]
    return imx, xsc, ysc

  def extractRatioData( self, country, c):
    x = [t[0] for t in c.data]
    xt, y = c.ratio(country)  
    x = x[xt[0]:xt[1]+1]

    if filter:
      yhat = savgol_filter(y, ww, np) # window size 51, polynomial order 3
    else:
      yhat = y
    ym = numpy.median( yhat )
    imx = numpy.argmin( yhat - ym )
    print ( "Ratio: Date of median: %s" % x[imx] )
    ##imx = int( len(ym )/2)
    xsc = [i-imx for i in range(len(x))]
    ysc = [z + self.y0 for z in yhat]
    return xt[0], imx, xsc, ysc

  def addplot(self, country, c, ax,color,filter=True, mode="cases"):

    try:
      if mode == "cases":
        ioff = 0
        imx, xsc, ysc = self.extractCaseData( country, c)
      else:
        ioff, imx, xsc, ysc = self.extractRatioData( country, c)
    except:
      print( "ERROR: no data for %s" % country )
      return

    ax.plot(xsc,ysc,color=color)
    ax.plot([-50.,75.],[self.y0,]*2,'--',color=color)

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

f, ax = plt.subplots(1)


ax.plot(xdsc,ydsc,color="blue")
ax.plot(xld,yld, "v", color="blue" )
ic =0

adder = Addplots(0.1,0.1)
##plt.ylim(-0.1,2.0)
for country in countries:
  adder.addplot( country, c, ax,color=colors[ic],filter=withFilter )
  ic += 1

ax.plot(xdsc9,[y+adder.y0 for y in ydsc9],color="blue")

plt.savefig("covid.png", bbox_inches='tight')

title = '%s .. Savgoy-Gorsky[%s,%s]' % (country,ww,np)
#title = 'test, gamma = 0.8 -> 0.1'
plt.ylabel(title)
plt.show()
