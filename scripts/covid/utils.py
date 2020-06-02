import rk
import numpy
import xlrd
import datetime, collections, os
import pandas
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

english_regions = ["South West", "South East", "London", "East of England", "West Midlands",
                "East Midlands", "Yorkshire and The Humber", "North West", "North East"] 

def mint(s):
  if s == "":
    return 0
  elif s.find('.') != -1:
    return float(s)
  else:
    return int(s)

class DateConv(object):
  def __init__(self,ref="2020-01-01"):
      self.oref = self._ord(ref)
      self.ref = ref

  def _ord(self,ds):
      year, month, day = [int(x) for x in ds.split( "-" )]
      date = datetime.date( year, month, day )
      self.date = date
      return date.toordinal()

  def convert(self,ds):
    return self.date.timetuple().tm_yday
    ##return self._ord(ds) - self.oref

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

      y,m,d = [int(xx) for xx in sdate.split('-')]
      jdate = datetime.datetime( y,m,d).timetuple().tm_yday
      ##jdate = dconv.convert( sdate )
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
    self.by_date = collections.defaultdict( dict )
    self.dd["date"] = s1

    for region in ["England",] :
##+ english_regions:
      s1 = sorted( list( cc[region].keys() ) )
      self.dd[region] = [mint(cc[region][x][0]) for x in s1]
      for k, v in cc[region].items():
        if k in s1:
          self.by_date[region][dconv.convert(k)] = mint(v[0])

    if end_truncate > 0:
      ks = sorted( list( self.by_date.keys() ) )
      if len(ks) > end_truncate + 1:
        for i in range(1,end_truncate+1):
          self.by_date[region][ks[-i]] = self.data_mask_value
      else:
        print ("WARNING: request to trincate > data length: %s > %s" % (end_truncate, len(ks)) )

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

  def addEnglandData(self, with_regions=True):
    engl = EnglandCases(mask_value=self.data_mask_value)
    regions = ["England",]
    if with_regions and False:
      regions += english_regions
    
    for region in regions:
      this = []
      for d in self.dd["date"]:
        this.append( engl.by_date[region].get( d, self.data_mask_value ) )
      self.dd[region] = this

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

class ACAPS(object):
  def __init__(self):
    """Database of measures taken by Governments ... just over 16k lines ... on June 1st."""

##
## reading this with pandas is painfully slow
##
##    self.pd = pandas.read_excel( "acaps-_covid19_government_measures_dataset.xlsx", sheet_name='Database' )
    wb = xlrd.open_workbook( "acaps-_covid19_government_measures_dataset.xlsx" )
    s = wb.sheet_by_name( 'Database' )
    self.countries = set()
    kk = 0
    cmap = {  'Czech republic':'Czech Republic', 
              'Moldova Republic Of':'Moldova Republic of',
              'kenya':'Kenya' }

    cc = collections.defaultdict(set)

    self.data = collections.defaultdict( lambda : collections.defaultdict(set) )
    self.regions_by_country =  collections.defaultdict(set)

    ccm = {'regions':3, 'log_type':6, 'category':7, 'measure':8 }


    for i in range(1,s.nrows):
        
        kk += 1
        this = s.row(i)[1].value
        if this != '':
          rec = s.row(i)
          region = s.row(i)[ccm['regions']].value
          for k,v in ccm.items():
            cc[k].add( s.row(i)[v].value )

          if type( rec[12].value ) == type( 0. ):
            dt = xlrd.xldate.xldate_as_datetime(rec[12].value, wb.datemode)
          else:
            print ("Unexpected date type ... %s:" % type( rec[12].value ) )
            print (rec[12].value)
            print (rec)
            continue

          jday = dt.timetuple().tm_yday

          data_key = this
          if region != '':
            data_key += ' (%s)' % region
            self.regions_by_country[this].add( data_key )
          self.data[data_key][jday].add( tuple( [rec[ccm[x]].value for x in ['log_type','category','measure']]) )

          if s.row(i)[7].value in [42,'42']:
             print ('42:  %s, %s' % (i,this) )

          if (kk//200)*200 == kk:
            print( kk,  this  )
          self.countries.add( cmap.get(this,this)  )

    for k in sorted( list( cc.keys() ) ):
      print ( "%s :: %s" % (k,cc[k]) )

        

class JH_Covid(object):
  def __init__(self):
    """Global confirmed cases.
       China, Australia, Canada only provided at the state level .. need to aggregate
       Not in ACAPS:
       'Andorra', 
         'Brunei' --> 'Brunei Darussalam'
         'Congo (Brazzaville)' --> 'Congo'
          'Congo (Kinshasa)' --> 'Congo DR'
          "Cote d'Ivoire" --> "Côte d'Ivoire"
          'Diamond Princess'
          'Czechia' --> 'Czech Republic', 'Czech republic'
          'Holy See'
          'Korea, South' --> 'Korea Republic of'
          'Moldova' --> 'Moldova Republic Of'¸ 'Moldova Republic of'
          'Monaco' 
          'North Macedonia' --> 'North Macedonia Republic Of'
          'Russia' --> 'Russian Federation'
          'Taiwan*' --> 
          'US'  --> United States of America
          'Vietnam' --> 'Viet Nam'
          'Laos'
          'West Bank and Gaza'
          'Kosovo'
          'Burma' --> 'Myanmar'
          'MS Zaandam'
          'Western Sahara'

         Only in ACAPS:
          ['', 'Australia', 'Brunei Darussalam', 'Canada', 'China', 'China, Hong Kong Special Administrative Region', 'Congo', 'Congo DR', 'Czech Republic', 'Czech republic', "Côte d'Ivoire", 'Kiribati', 'Korea DPR', 'Korea Republic of', 'Lao PDR', 'Marshall Islands', 'Micronesia', 'Moldova Republic Of', 'Moldova Republic of', 'Myanmar', 'Nauru', 'North Macedonia Republic Of', 'Palau', 'Palestine', 'Russian Federation', 'Samoa', 'Solomon Islands', 'Tonga', 'Turkmenistan', 'Tuvalu', 'United States of America', 'Vanuatu', 'Viet Nam', 'kenya']
    """
    self.agset = collections.defaultdict( set )
    map_to_acaps = { 'Brunei':'Brunei Darussalam',
          'Congo (Brazzaville)':'Congo',
          'Congo (Kinshasa)':'Congo DR',
          "Cote d'Ivoire":"Côte d'Ivoire",
          'Czechia':'Czech Republic',
          'Korea, South':'Korea Republic of',
          'Moldova':'Moldova Republic of',
          'North Macedonia':'North Macedonia Republic Of',
          'Russia':'Russian Federation',
          'US':'United States of America',
          'Vietnam':'Viet Nam',
          'Burma':'Myanmar' }

    ## starts 22nd Jan.
    ### with ref="2020-01-01", starting at "21"
    ##
    jh = pandas.read_csv( "time_series_covid19_confirmed_global.csv" )
    self.jh = jh
    self.dates = jh.columns[4:]

    mdy = [[int(xx) for xx in x.split('/')] for x in self.dates]
    self.jday = [datetime.datetime( y,m,d).timetuple().tm_yday for m,d,y in mdy]

    self.cases = jh.values[:,4:]
    self.countries = [ ]
    self.regions = [ ]
    self.regions_by_country = collections.defaultdict( set )
    self.lines = [ ]
    self.dd = dict()
    self.dd['date'] = self.jday[1:]
    self.data_mask_value = -1.e20

    for i in jh.index:
      c0 = jh.values[i,1]
      c0 = map_to_acaps.get( c0, c0 )
      this = jh.values[i,4:].tolist()
      if type( jh.values[i,0] ) == type(''):
        r0 = '%s (%s)' % (c0,jh.values[i,0])
        self.regions.append( r0 )
        self.regions_by_country[c0].add( r0 )
        self.lines.append( r0 )
        self.agset[c0].add( r0 )
        self.dd[r0] = [this[j+1] - this[j] for j in range(len(self.dates)-1)]
      else:
        self.dd[c0] = [this[j+1] - this[j] for j in range(len(self.dates)-1)]
        self.countries.append( c0 )
        self.lines.append( c0 )
     
    


if __name__ == "__main__":
  ec = EnglandCases()
