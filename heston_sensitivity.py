import QuantLib as ql # version 1.5
import matplotlib.pyplot as plt

def black_scholes_heston(volatility,maturity_date, spot_price,strike_price,dividend_rate,option_type, risk_free_rate,calculation_date,hestonparams):
	day_count = ql.Actual365Fixed()
	calendar = ql.UnitedStates()
	ql.Settings.instance().evaluationDate = calculation_date
	payoff = ql.PlainVanillaPayoff(option_type, strike_price)
	exercise = ql.EuropeanExercise(maturity_date)
	european_option = ql.VanillaOption(payoff, exercise)
	spot_handle = ql.QuoteHandle(
	ql.SimpleQuote(spot_price))
	flat_ts = ql.YieldTermStructureHandle(ql.FlatForward(calculation_date, risk_free_rate, day_count))
	dividend_yield = ql.YieldTermStructureHandle(
	ql.FlatForward(calculation_date, dividend_rate, day_count))
	volatility=volatility[0]
	flat_vol_ts = ql.BlackVolTermStructureHandle(
	ql.BlackConstantVol(calculation_date, calendar, volatility, day_count))
	bsm_process = ql.BlackScholesMertonProcess(spot_handle, 
										   dividend_yield, 
										   flat_ts, 
										   flat_vol_ts)
	european_option.setPricingEngine(ql.AnalyticEuropeanEngine(bsm_process))
	bs_price = european_option.NPV()
	spot_handle = ql.QuoteHandle(ql.SimpleQuote(spot_price))
	flat_ts = ql.YieldTermStructureHandle(
	ql.FlatForward(calculation_date, risk_free_rate, day_count))
	dividend_yield = ql.YieldTermStructureHandle(
	ql.FlatForward(calculation_date, dividend_rate, day_count))
	heston_process = ql.HestonProcess(flat_ts,
								  dividend_yield,
								  spot_handle,
								  hestonparams[0],
								  hestonparams[1],
								  hestonparams[2],
								  hestonparams[3],
								  hestonparams[4])
	engine = ql.AnalyticHestonEngine(ql.HestonModel(heston_process),0.01, 1000)
	european_option.setPricingEngine(engine)
	h_price = european_option.NPV()
	#return h_price,bs_price
	return bs_price-h_price
#print(black_scholes.black_scholes_heston(0.24683,0.1,t,2811,s,0.02,ql.Option.Call,0.022,calculation_date1,[v0,kappa,theta,sigma,rho]))['x'])
def black_scholes_sensitivity(volatility,maturity_date, spot_price,strike_price,dividend_rate,option_type, risk_free_rate,calculation_date,hestonparams):
	day_count = ql.Actual365Fixed()
	calendar = ql.UnitedStates()
	ql.Settings.instance().evaluationDate = calculation_date
	payoff = ql.PlainVanillaPayoff(option_type, strike_price)
	exercise = ql.EuropeanExercise(maturity_date)
	european_option = ql.VanillaOption(payoff, exercise)
	spot_handle = ql.QuoteHandle(
	ql.SimpleQuote(spot_price))
	flat_ts = ql.YieldTermStructureHandle(ql.FlatForward(calculation_date, risk_free_rate, day_count))
	dividend_yield = ql.YieldTermStructureHandle(
	ql.FlatForward(calculation_date, dividend_rate, day_count))
	volatility=volatility[0]
	flat_vol_ts = ql.BlackVolTermStructureHandle(
	ql.BlackConstantVol(calculation_date, calendar, volatility, day_count))
	bsm_process = ql.BlackScholesMertonProcess(spot_handle, 
										   dividend_yield, 
										   flat_ts, 
										   flat_vol_ts)
	european_option.setPricingEngine(ql.AnalyticEuropeanEngine(bsm_process))
	bs_price = european_option.NPV()
	spot_handle = ql.QuoteHandle(ql.SimpleQuote(spot_price))
	flat_ts = ql.YieldTermStructureHandle(
	ql.FlatForward(calculation_date, risk_free_rate, day_count))
	dividend_yield = ql.YieldTermStructureHandle(
	ql.FlatForward(calculation_date, dividend_rate, day_count))
	heston_process = ql.HestonProcess(flat_ts,
								  dividend_yield,
								  spot_handle,
								  hestonparams[0],
								  hestonparams[1],
								  hestonparams[2],
								  hestonparams[3],
								  hestonparams[4])
	engine = ql.AnalyticHestonEngine(ql.HestonModel(heston_process),0.01, 1000)
	european_option.setPricingEngine(engine)
	h_price = european_option.NPV()
	#return h_price,bs_price
	return bs_price
def heston_sensitivity(volatility,maturity_date, spot_price,strike_price,dividend_rate,option_type, risk_free_rate,calculation_date,hestonparams):
	day_count = ql.Actual365Fixed()
	calendar = ql.UnitedStates()
	ql.Settings.instance().evaluationDate = calculation_date
	payoff = ql.PlainVanillaPayoff(option_type, strike_price)
	exercise = ql.EuropeanExercise(maturity_date)
	european_option = ql.VanillaOption(payoff, exercise)
	spot_handle = ql.QuoteHandle(
	ql.SimpleQuote(spot_price))
	flat_ts = ql.YieldTermStructureHandle(ql.FlatForward(calculation_date, risk_free_rate, day_count))
	dividend_yield = ql.YieldTermStructureHandle(
	ql.FlatForward(calculation_date, dividend_rate, day_count))
	volatility=volatility[0]
	flat_vol_ts = ql.BlackVolTermStructureHandle(
	ql.BlackConstantVol(calculation_date, calendar, volatility, day_count))
	bsm_process = ql.BlackScholesMertonProcess(spot_handle, 
										   dividend_yield, 
										   flat_ts, 
										   flat_vol_ts)
	european_option.setPricingEngine(ql.AnalyticEuropeanEngine(bsm_process))
	bs_price = european_option.NPV()
	spot_handle = ql.QuoteHandle(ql.SimpleQuote(spot_price))
	flat_ts = ql.YieldTermStructureHandle(
	ql.FlatForward(calculation_date, risk_free_rate, day_count))
	dividend_yield = ql.YieldTermStructureHandle(
	ql.FlatForward(calculation_date, dividend_rate, day_count))
	heston_process = ql.HestonProcess(flat_ts,
								  dividend_yield,
								  spot_handle,
								  hestonparams[0],
								  hestonparams[1],
								  hestonparams[2],
								  hestonparams[3],
								  hestonparams[4])
	engine = ql.AnalyticHestonEngine(ql.HestonModel(heston_process),0.01, 1000)
	european_option.setPricingEngine(engine)
	h_price = european_option.NPV()
	#return h_price,bs_price
	return h_price
