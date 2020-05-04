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



a = 0.2
r = .25
ixl = 90
##
def add_plot(ax, this, col, ixl, bbox_props):
  ax.plot(this.xdsc,this.ydsc,color=col)
  ax.plot(this.xld,this.yld, "v", color=col )
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

  d = Demo01('test',g0,g1,a00,a01,a11,a12,scfac, 40,70, model_rank=4,fac4=1.0)
  d.run()
  d.review(offset=plot_offset,scale=plot_scale)
  add_plot(ax, d, color, ixl, bbox_props)

##
f, ax = plt.subplots(1)

bbox_props = dict(fc="cyan", ec="b", lw=1, boxstyle="round")
mrun( a, r, g0, 0., ax, "blue", ixl, bbox_props)
mrun( a, r, g0, 0.9, ax, "blue", ixl, bbox_props)
mrun( a, r, g0, 1.0, ax, "blue", ixl, bbox_props)
mrun( a, r*1.5, g0, .0, ax, "green", ixl, bbox_props)

##d1 = Demo01('test',g0,gc,a00,a01,a11,a12,scfac, 40,70, model_rank=4,fac4=1.0)
##d1.run()
##d1.review(offset=plot_offset,scale=plot_scale)

##d9 = Demo01('test',g0,.9*gc,a00,a01,a11,a12,scfac, 40,70, model_rank=4,fac4=1.0)
##d9.run()
##d9.review(offset=plot_offset,scale=plot_scale)


##add_plot(ax, d1, "green", ixl, bbox_props)
##add_plot(ax, d9, "orange", ixl, bbox_props)

plt.savefig("modelled_pdmc01.png", bbox_inches='tight')

if plot_scale:
  title = 'New Symptomatic Cases (scaled)' 
else:
  title = 'New Symptomatic Cases' 
plt.ylabel(title)
plt.show()
