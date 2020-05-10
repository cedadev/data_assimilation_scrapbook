import rk
import numpy
import xlrd
import datetime, collections, os
""" Provides: 
  * get_key_dates
  * Cases 
  * EnglandCases (under construction)
  * mint (for internal use)
  * DateConv (for internal use)
"""

## python2.7

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

class Workbook(object):
  def __init__(self,file):
    assert os.path.isfile(file), 'File %s not found' % file
    self.book = xlrd.open_workbook( file )
    self.sns = self.book.sheet_names()

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

class EnglandCases(object):
  def __init__(self,file="coronavirus-cases_latest.csv",targets=["Reading",],end_truncate=4,mask_value=-99999):
#Area name,Area code,Area type,Specimen date,Daily lab-confirmed cases,Previously reported daily cases,Change in daily cases,Cumulative lab-confirmed cases,Previously reported cumulative cases,Change in cumulative cases
    self.other = None
    self.data_mask_value = mask_value
    cc = collections.defaultdict(dict)
    sa = set()
    for l in open(file).readlines():
      aname,id,atype,date,cases,casesp = [x.strip() for x in l.split(",")[:6]]
      sa.add(atype)
      if atype in ["Region","Nation"]:
        cc[aname][date] = (cases,casesp)

    s1 = sorted( list( cc['England'].keys() ) )
    
    nn = 0
    for k in cc.keys():
      sx = sorted( list( cc[k].keys() ) )
      if s1 != sx:
        print ("Mismatch in date set for %s  vs. England: len %s -- %s" % (k,len(s1),len(sx)) )
        nn += 1

    print ("%s data records for regions %s" % (len(s1),cc.keys()) )
    print (sa )
    self.dd = dict()
    self.by_date = dict()
    self.dd["date"] = s1
    self.dd["England"] = [int(cc["England"][x][0]) for x in s1]
    for k, v in cc["England"].items():
       if k in s1:
         self.by_date[dconv.convert(k)] = int(v[0])

    if end_truncate > 0:
      ks = sorted( list( self.by_date.keys() ) )
      for i in range(1,end_truncate+1):
        self.by_date[ks[-i]] = self.data_mask_value

    wb = Workbook( "data/COVID-19-total-announced-deaths-5-May-2020.xlsx" )
    assert "COVID19 total deaths by region" in wb.sns
    self.sh = wb.book.sheet_by_name( "COVID19 total deaths by region" )
##
## row 15 has dates in xldate format. row 16, col 3 onwards is the daily deaths.
##
##xlrd.xldate_as_datetime( e.sh.row(15)[6].value, 0 )
        

class Cases(object):
  def __init__(self,file="total_cases.csv",targets=["Germany",],add_england=False):
    self.data_mask_value = -99999
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

      this[0] = dconv.convert( this[0] )
      self.data.append( this )
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

    if add_england:
      self.addEnglandData()

  def addEnglandData(self):
    engl = EnglandCases(mask_value=self.data_mask_value)
    this = []
    for d in self.dd["date"]:
      this.append( engl.by_date.get( d, self.data_mask_value ) )
    self.dd["England"] = this

  def setKeyDates(self,kd):
    self.key_dates = kd

  def setOther(self,other):
    self.other = other

  def getDD7(self):
    self.dd7 = dict()
    w = 1./7.
    for k in self.dd.keys():
      this = [0.]*(len(self.dd[k]) - 6) 
      for i in range(len(this)):
        this[i] = w*sum( self.dd[k][i:i+7] )
      self.dd7[k] = this
    self.dd7["date"] = self.dd["date"][3:-3]

  def key_points(self):
    self.kp = dict()
    for k, data in self.dd7.items():
      imx = numpy.argmax( data )
      grad = [data[i+1] - data[i] for i in range(len(data)-1)]
      inf0 = numpy.argmax( grad )
      inf1 = numpy.argmin( grad )
      self.kp[k] = [self.dd7["date"][i] for i in [inf0,imx,inf1]] + [len(data) -1 - inf1 < 3,]

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
      this = numpy.ma.masked_equal([t[ir] for t in self.data],self.data_mask_value)
      imx = numpy.argmax( this )
      if imx > 10:
        f3 = numpy.polyfit(list( range(-10,1)),this[imx-10:imx+1],2)
        if this[imx] > 1000 and abs(f3[1]) < 20*abs(f3[0]):
          print( self.headings[ir], this[imx], f3 )



if __name__ == "__main__":
  ec = EnglandCases()
