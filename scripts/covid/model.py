
import numpy, random, scipy
import datetime, collections
import matplotlib.pyplot as plt
import rk


## python2.7

from scipy.signal import savgol_filter

## y0: numbers infected and not yet showing symtoms
## y1: numbers showing symptoms and infectious
## y2: numbers that are critical

def relax_linear_fit(obs,w,relative_w = None, end_weighting={"lfrac":30, "amp":2}, masked=False):
    """Fitting (d^2/dt^2 - z^t z) res = w(res-obs)"""

    if type(obs) == type( [] ):
      res = obs[:]
    else:
      res = obs.copy()
    imx = numpy.argmax( obs )

    if masked:
      mres = res.copy()
      if obs.mask == False:
        i99 = max( [i for i in range( imx, len(obs) )  ] )
        i00 = 0
      else:
        i99 = max( [i for i in range( imx, len(obs) ) if not obs.mask[i] ] )
        if any( [obs.mask[i] for i in range( imx )] ):
          i00 = max( [i for i in range( imx ) if obs.mask[i] ] ) + 1
          mres.mask[:i00] = True
        else:
          i00 = 0
        mres.mask[i99+1:] = True
      res = res[i00:i99+1]
      obs = obs[i00:i99+1]
      imx = numpy.argmax( obs )

    m0 = 0
    for x in obs:
      m0 += x**2
    m0 = numpy.sqrt( m0/len( obs ) )

    if relative_w == None:
      ww = [w,]*len(obs)
    else:
      ww = [x*w for x in relative_w]

    if end_weighting != None:
      nn = int( (len(obs)*end_weighting["lfrac"])/100 )
      sc = [0.]*nn
      for i in range(nn):
        sc[i] = 1. + i*float( end_weighting["amp"] - 1 )/(nn-1)

      for i in range(nn):
        ww[i] = ww[i]*sc[nn-1-i]
        ww[len(obs)-1-i] = ww[len(obs)-1-i]*sc[nn-1-i]


    use_numpy = True

    if use_numpy:
      ab = numpy.zeros( (3,len(obs)), "float" )
      b = numpy.zeros( (len(obs)), "float" )
      ab[0,1:] = - 1.
      ab[2,:-1] = - 1.
      for i in range(len(obs)):
        ab[1,i] = 2. + ww[i]
        b[i] = obs[i]*ww[i]
      ab[1,0] = 1. + ww[i]
      ab[1,len(obs)-1] = 1. + ww[i]
      res = scipy.linalg.solve_banded( (1,1), ab, b )

    else:

      flipper = 0
      for k in range(4000):
        f = [0,]*len(obs)
        s = 0.
        for i0 in range( 1,len(obs)-1 ):
          i = i0*(1-flipper) + (len(obs) -1 - i0)*flipper
          try:
            s += (res[i]-obs[i])**2 + (res[i-1] - 2.*res[i] + res[i+1])**2
  
            f[i] = res[i-1] - 2.*res[i] + res[i+1] - ww[i]*(res[i]-obs[i])
            res[i] = res[i] + 1.0*f[i]/(2.+ww[i])
          except:
            print (flipper, i, i0 )
            raise

        s = numpy.sqrt( s/len( obs ) )
        ##if (k/100)*100 == k:
          ##print (k,s, s/m0, [res[i] - obs[i] for i in range(imx-5,imx+5)] )
        if s < m0*0.001:
          if masked:
            mres[i00:i99+1] = res[:]
            return mres
          else:
            return res

      res[0] = max( [0., 2*res[1] - res[2]] )
      res[-1] = max( [0., 2*res[-2] - res[-3]] )
      flipper = 1 - flipper

    if masked:
      mres[i00:i99+1] = res[:]
      return mres
    else:
      return res

def relax_linear_matrix_fit(obs,ix,zz0,w):
    res = numpy.zeros( (len(obs),4), dtype="float" )
    res[:,ix] = obs[:]
    m0 = 0
    for x in obs:
      m0 += x**2
    m0 = numpy.sqrt( m0/len( obs ) )
    for k in range(200):
      f = [0,]*len(obs)
      s = 0.
      for i in range( 1,len(obs)-1 ):
        s += (res[i]-obs[i])**2 + (res[i-1,ix] - 2.*res[i,ix] + res[i+1,ix])**2
        f[i] = (res[i-1,ix] - 2.*res[i,ix] + res[i+1,ix]) - w*(res[i,ix]-obs[i])
      s = numpy.sqrt( s/len( obs ) )
      print (k,s, s/m0, f[25:28] )
      if s < m0*0.001:
        return res
      for i in range( 1,len(obs)-1 ):
        res[i] = res[i] + 0.9*f[i]/(2.+w)
    return res

class Demo01(object):
  def __init__(self,name,g0,g1,a00, a01,a11,a12, scfac,tc,t9,model_rank=3,fac4=1.):
    self.g0=g0
    self.g1=g1
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

    for k0 in range(self.t9):
      for k1 in range(10):
        x,y = rk.rk4d(self.f1,x,y,h)
      self.xx.append( (x,y) )
      ##print ('%5.2f:: %s' % (x,'%5i, %5i, %5i' % tuple( y ) ) )
    self.x = x
    self.y = y

  def review(self,offset=True,scale=True):
    self.xd = [t[0] for t in self.xx]
    self.yd = [t[1][2] for t in self.xx]

    self.imx = numpy.argmax( self.yd )
    print ( "Index of peak: %s [lockdown + %s]" % (self.imx,self.imx-self.tc ) )
    if offset:
      self.xdsc = [i-self.imx for i in self.xd]
    else:
      self.xdsc = self.xd[:]
    if scale:
      self.ydsc = [z/self.yd[self.imx] for z in self.yd]
    else:
      self.ydsc = self.yd[:]

    self.xld = self.xdsc[self.tc]
    self.yld = self.ydsc[self.tc]

    self.fd3 = numpy.polyfit(self.xdsc[self.imx-6:self.imx+6],self.ydsc[self.imx-6:self.imx+6],3)
    print ('model', self.fd3)

def add_plot(ax, this, col, ixl, bbox_props):
  ax.plot(this.xdsc,this.ydsc,color=col)
  if this.xld != None:
    ax.plot(this.xld,this.yld, "v", color=col )
  if this.g1 != None:
    t = ax.text(this.xdsc[ixl],this.ydsc[ixl], "%4.2f" % this.g1, ha="left", va="center", rotation=0,
        size=8,
        bbox=bbox_props)

def mrun(a,r,g0,gg, ax, color, ixl, bbox_props):
  a00 = a*(1-r)
  a01 = a*r
  a11 = a*(1-r)
  a12 = a*r
  gc = a / ( r**2. )
  g0 = 24.0
  g1 = gg*gc
  scfac = 1.0
  plot_offset = False
  plot_scale = True

  d = model.Demo01('test',g0,g1,a00,a01,a11,a12,scfac, 40,70, model_rank=4,fac4=1.0)
  d.run()
  d.review(offset=plot_offset,scale=plot_scale)
  add_plot(ax, d, color, ixl, bbox_props)


def plot_model():
  a = 0.2
  r = .25
  ixl = 90
##a00 = a*(1-r)
##a01 = a*r
##a11 = a*(1-r)
##a12 = a*r
  g0 = 24.0

##
  f, ax = plt.subplots(1)
  bbox_props = dict(fc="cyan", ec="b", lw=1, boxstyle="round")

  mrun( a, r, g0, 0., ax, "blue", ixl, bbox_props)
  mrun( a, r, g0, 0.9, ax, "blue", ixl, bbox_props)
  mrun( a, r, g0, 1.0, ax, "blue", ixl, bbox_props)
  mrun( a, r*1.5, g0, .0, ax, "green", ixl, bbox_props)
  
  plt.savefig("modelled_pdmc01.png", bbox_inches='tight')

  plt.show()

class Wrap01(object):
  def __init__(self,x,y,g1=None,xld=None):
    self.xdsc = x
    self.ydsc = y
    self.g1 = g1
    self.xld = xld

def plot_relax_linear_fit(case="ex01"):
  colors = ["purple","red","green","orange","pink","black","cyan"]
  if case == "ex01":
    obs = [0.,]*41
    obs[20] = 1.
    relw = [0.,]*41
  ##for i in [10,90,100,110,190]:
    for i in [15,20,25]:
      relw[i] = 1.
    xcoords = range(41)
    ip = 0
    weights = [100.,10., 1., .1, .01]
  elif case == "ex02":
    obs = [0.,]*201
    for i in range( 1,200):
      x = (i-100)*0.01
      obs[i] = x**2 + (random.random() - 0.5)*0.5
    relw = None
    xcoords = range(201)
    ip=2
    weights = [ 1., .1, .01]
    
  else:
    raise

  f, ax = plt.subplots(1)
  bbox_props = dict(fc="cyan", ec="b", lw=1, boxstyle="round")
  if case == "ex02":
    for i in range( 1,200):
      ax.plot(i,obs[i], "o", color="grey" )

  for ww in weights:
    res = relax_linear_fit( obs, ww, relative_w=relw)
    add_plot(ax, Wrap01(xcoords, res), colors[ip], 100, bbox_props)

    ip += 1

  if case == "ex01":
    ax.plot(15,0., "o", color="black" )
    ax.plot(20,1., "o", color="black" )
    ax.plot(25,0., "o", color="black" )

  plt.savefig("relax_linear_fit_%s.png" % case, bbox_inches='tight')
  print( obs )
  plt.show()
  
if __name__ == "__main__":
  plot_relax_linear_fit(case="ex02")
