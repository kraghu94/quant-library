#import heston
#import bsm
import pandas as pd
#import redis1
import datetime
import QuantLib as ql
import black_scholes
import scipy.optimize
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
S0=2811
risk_free_rate=0.022
dividend_rate=0.02
a=pd.read_excel('all spx tickers.xlsx')
#with open('vols.txt', 'w') as f:
#	for i in a.Tickers:
#		f.write(redis1.dPrice(i,'IVOL')+'\n')
with open('vols.txt','r') as f:
	vols=f.readlines()


print(vols)


import numpy as np
data = []
strikes1=[]
expiration_dates1=[]
call_or_put=[]
calculation_date="05/28/19"
#l = []
for i in range(len(a)):
    Ticker = a.Tickers[i]
    expiration_dates = a.Tickers[i][7:-12]
    call_or_put.append(a.Tickers[i][16])
    expiration_dates1.append(expiration_dates)
    strikes = (float(a.Tickers[i][17:-6])) if a.Tickers[i][17:-6] is not None else None
    strikes1.append(strikes)
    #print(vols[i])
    data.append(float(vols[i].strip()))
    #d = {'Ticker': Ticker,'Expiration Dates':expiration_dates,'Strike':strikes, 'I. Vols': vols[i]}
    #l.append(d1)

print(expiration_dates1)
strikes1 = [strikes1[i]for i in range(len(data))if data[i]is not None and call_or_put[i]=='C']
expiration_dates1 = [ql.Date(int(expiration_dates1[i][3:5]),int(expiration_dates1[i][0:2]),int('20'+expiration_dates1[i][6:8]))for i in range(len(expiration_dates1))if data[i]is not None and call_or_put[i]=='C']
data = [data[i] for i in range(len(data))if data[i]is not None and call_or_put[i]=='C']
implied_vols=np.array(data).reshape(len(set(strikes1)),len(set(expiration_dates1))).tolist()
data=np.array(data).reshape(len(set(expiration_dates1)),len(set(strikes1))).tolist()
expiration_dates=expiration_dates1
expiration_dates1=list(sorted(set(expiration_dates1)))
strikes=strikes1
strikes1=list(sorted(set(strikes1)))
#print(data, strikes1,expiration_dates1)

day_count = ql.Actual365Fixed()
calendar = ql.UnitedStates()
#print(data)
calculation_date1 = ql.Date(28, 5, 2019)
#black_var_surface = ql.BlackVarianceSurface(
 # calculation_date1, calendar, 
  # sorted(set(expiration_dates1)), sorted((set(strikes1))), 
   # data, day_count)
#print(abs((datetime.datetime.strptime('01/17/20','%m/%d/%y')-datetime.datetime.strptime('05/21/19','%m/%d/%y')).days)/365)
#print(black_var_surface.blackVol(3, 2125))


flat_ts = ql.YieldTermStructureHandle(
    ql.FlatForward(calculation_date1, risk_free_rate, day_count))
dividend_ts = ql.YieldTermStructureHandle(
    ql.FlatForward(calculation_date1, dividend_rate, day_count))		
#v0=np.arange(start=0,stop=1,step=0.01)
#kappa=np.arange(start=0,stop=1,step=0.01)
#theta=np.arange(start=0,stop=1,step=0.01)
#rho=np.arange(start=0,stop=1,step=0.01)
#sigma=np.arange(start=0,stop=1,step=0.01)
# Set Initial Conditions
v0 = 0.38; kappa = 0.20; theta = 0.6; rho = -0.3; sigma = 1;
#

process = ql.HestonProcess(flat_ts, dividend_ts, 
                           ql.QuoteHandle(ql.SimpleQuote(S0)), 
                           v0, kappa, theta, sigma, rho)
model = ql.HestonModel(process)
engine = ql.AnalyticHestonEngine(model) 

heston_helpers = []
#print(len(expiration_dates1))
#black_var_surface.setInterpolation("bicubic")
#one_year_idx = 4 # 12th row in data is for 1 year expiry
#date = expiration_dates1[one_year_idx]
#print(len(strikes1))
heston_helpers = []

for i, s in enumerate(strikes1):
    for j,t in enumerate(expiration_dates1):
        t = (expiration_dates1[j] - calculation_date1)
        p = ql.Period(t, ql.Days)
        sigma = data[j][i]
        #sigma = black_var_surface.blackVol(t/365.25, s)
        helper = ql.HestonModelHelper(p, calendar, S0, s, 
                                  ql.QuoteHandle(ql.SimpleQuote(sigma/100)),
                                  flat_ts, 
                                  dividend_ts)
        helper.setPricingEngine(engine)
        heston_helpers.append(helper)
lm = ql.LevenbergMarquardt(1e-8, 1e-8, 1e-8)
model.calibrate(heston_helpers, lm, 
                 ql.EndCriteria(500, 50, 1.0e-8,1.0e-8, 1.0e-8))
theta, kappa, sigma, rho, v0 = model.params()
print(theta,kappa,sigma,rho,v0)
#print(black_scholes.black_scholes_heston(ql.Date(17, 1, 2020),2811,2025,0.24978200,0.02,ql.Option.Call,0.022,ql.Date(21, 5, 2019),[v0,kappa,theta,sigma,rho]))
stochastic_vol_invert=[]

for i, s in enumerate(strikes1):
    for j,t in enumerate(expiration_dates1):
    	stochastic_vol_invert.append(scipy.optimize.root(black_scholes.black_scholes_heston,0.1,args=(t,S0,s,0.02,ql.Option.Call,0.022,calculation_date1,[v0,kappa,theta,sigma,rho]))['x'].item())	
stochastic_vol_invert=np.array(stochastic_vol_invert).reshape(len(set(strikes1)),len(set(expiration_dates1))).tolist()
print(stochastic_vol_invert)
black_var_surface = ql.BlackVarianceSurface(calculation_date1, calendar, expiration_dates1, strikes1, stochastic_vol_invert, day_count)
#black_var_surface.setInterpolation("linear")
plot_years = np.arange(0, 1.5, 0.1)
plot_strikes = np.arange(2400, 2975, 1)
fig = plt.figure()
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(plot_strikes, plot_years)
Z = np.array([black_var_surface.blackVol(float(y),float(x)) 
              for xr, yr in zip(X, Y) 
                  for x, y in zip(xr,yr) ]
            ).reshape(len(X), len(X[0]))
surf = ax.plot_surface(X,Y,Z, rstride=1, cstride=1, cmap=cm.seismic, 
               linewidth=0.1)

fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()

#print(theta, kappa, sigma, rho, v0)
avg = 0.0

print  (
    "Strikes", "Market Value", 
    "Model Value", "Relative Error (%)")
print ("="*70)
for i, opt in enumerate(heston_helpers):
    err = (opt.modelValue()/opt.marketValue() - 1.0)
    print (strikes[i], opt.marketValue(), 
        opt.modelValue(), 
        100.0*(opt.modelValue()/opt.marketValue() - 1.0))
    avg += abs(err)
avg = avg*100.0/len(heston_helpers)
print ("-"*70)
print ("Average Abs Error (%%) : %5.3f" % (avg))
