
import utils


def run_key_dates():

##
## import case data
  c = utils.Cases()

# 7-day smoothing

  c.getDD7()

# identify max and inflectia

  c.key_points()

# import lockdown dates
  c.setKeyDates(utils.get_key_dates())

  for k in ["Germany","France","Spain","Italy","Austria","China","New Zealand"]:
    this = list( c.key_dates[k] )
    ll = [x[1] for x in this]
    if "lockdown" in ll:
      ldate = this[ ll.index( "lockdown" ) ][0]
    else:
      ldate = None
    print (k,c.kp[k],ldate)



if __name__ == "__main__":
   run_key_dates()
   
