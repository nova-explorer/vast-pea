# -*- coding: utf-8 -*-
"""
@author: Olivier Couture
         olivier.couture@list.lu

Description: curve_fit wrappers for different functions and data

TODO:
    - Units could probably be handled better
"""

from scipy.optimize import curve_fit, Bounds
from scipy.special import gamma, polygamma
import numpy as np
import pint

from constants import GAS_CONSTANT

UNIT_REG = pint.UnitRegistry()
PRINT_TEMPLATE_ERROR = "\t{name} : {value:.{dp}f} ± {error:.{dp}f} {unit}"
PRINT_TEMPLATE_ERROR_E = "\t{name} : {value:.{dp}e} ± {error:.{dp}e} {unit}"
PRINT_TEMPLATE_VALUE = "\t{name} : {value:.{dp}f} {unit}"
PRINT_TEMPLATE_VALUE_E = "\t{name} : {value:.{dp}e} {unit}"
PRINT_TEMPLATE_RSQ = "\tR^2 : {:.2f}"

def equation_KWW(x, tau, beta) -> np.ndarray:
    return np.exp(-(x/tau)**beta)

def equation_linear(x, a, b) -> np.ndarray:
    return a * x + b

def equation_exponential_convergence(x, xconv, alpha, gamma):
    return xconv - alpha / x**gamma

def equation_wlffv(T, ffvg, cte):
    return ffvg + cte * T

def fit(x, y, equation, p0=[1,1], sigma=None, absolute_sigma=False, bounds=False, maxfev=1e7):

    if bounds:
        popt, pcov = curve_fit(equation, x, y, p0, sigma=sigma, absolute_sigma=absolute_sigma, bounds=bounds, maxfev=int(maxfev))
    else:
        popt, pcov = curve_fit(equation, x, y, p0, sigma=sigma, absolute_sigma=absolute_sigma, maxfev=int(maxfev))

    fitted = equation(x, *popt)
    perr = np.sqrt(np.diag(pcov))
    r_sq = 1 - ( np.sum((y - fitted)**2) / np.sum((y - np.mean(y))**2) )

    return popt, perr, r_sq, fitted

def find_x_at_y_percent(x:np.ndarray, y:np.ndarray, value:float) -> float:
    idx = ( np.abs(y - y[0]*value) ).argmin()
    return float(x[idx])

def find_decimal_point(x) -> int:
    return int(np.ceil(-np.log10(abs(x))))

def choose_decimal_point(value, error=0) -> int:
    decimal_point = 0
    dp_error = 0
    dp_value = 0

    if error!=0 and value!=0:
        dp_error = find_decimal_point(error)
        dp_value = find_decimal_point(value)
    elif error==0 and value!=0:
        dp_error = 0
        dp_value = find_decimal_point(value)
    elif error!=0 and value == 0:
        dp_error = find_decimal_point(error)
        dp_value = 0

    if dp_error > dp_value:
        decimal_point = dp_error + 1
    else:
        decimal_point = dp_value + 1

    if decimal_point < 0:
        decimal_point = 0

    return decimal_point

def print_value(name, value, error=None, units=None) -> None:

    if find_decimal_point(value) > 6:
        if error:
            str_template = PRINT_TEMPLATE_ERROR_E
        else:
            str_template = PRINT_TEMPLATE_VALUE_E
            error = 0
    else:
        if error:
            str_template = PRINT_TEMPLATE_ERROR
        else:
            str_template = PRINT_TEMPLATE_VALUE
            error = 0

    try:
        print(str_template.format(
            name=name,
            value=value,
            error=error,
            dp=choose_decimal_point(value, error),
            unit=units
        ))
    except OverflowError:
        error = 0
        print(str_template.format(
            name=name,
            value=value,
            error=error,
            dp=choose_decimal_point(value, error),
            unit=units
        ))

class fit_KWW():
    def __init__(self, x:np.ndarray, y:np.ndarray, beta_ini:float, initial_units='fs', final_units='ps', verbose=True):

        if verbose:
            print("Fitting KWW single stretched exponential decay :")

        self.x = x
        self.y = y

        tau_ini = find_x_at_y_percent(self.x, self.y, 1/np.exp(1))
        popt, perr, self.r_sq, self.fitted = fit(self.x, self.y, equation_KWW, [tau_ini, beta_ini])

        self.tau = popt[0]
        self.tau_err = perr[0]
        self.beta = popt[1]
        self.beta_err = perr[1]

        self.tau_ave = self.tau / self.beta * gamma(1/self.beta)
        deriv_tau = 1/self.beta * gamma(1/self.beta)
        deriv_beta = -self.tau * gamma(1/self.beta) * (self.beta + polygamma(0,1/self.beta)) / self.beta**3
        self.tau_ave_err = np.sqrt(deriv_tau**2 * self.tau_err**2 + deriv_beta**2 * self.beta_err**2)

        tau_units = UNIT_REG[initial_units].to(final_units) # pyright: ignore[reportAttributeAccessIssue]
        self.tau *= tau_units.magnitude
        self.tau_err *= tau_units.magnitude
        self.tau_ave *= tau_units.magnitude
        self.tau_ave_err *=tau_units.magnitude

        if verbose:
            print_value(
                name = "Tau",
                value = self.tau,
                error = self.tau_err,
                units = tau_units.units
            )
            print_value(
                name = "Beta",
                value = self.beta,
                error = self.beta_err,
                units = ''
            )
            print(PRINT_TEMPLATE_RSQ.format(self.r_sq))
            print_value(
                name = "Tau_average",
                value = self.tau_ave,
                error = self.tau_ave_err,
                units = tau_units.units
        )

class fit_arrhenius:
    def __init__(self,
        temperature:np.ndarray, tau:np.ndarray,
        initial_units={'temp':'K', 'time':'ps'},
        final_units={'time':'s', 'energy':'kJ/mol'},
        verbose=True) -> None:

        if verbose:
            print("Fitting Arrhenius :")

        temperature = UNIT_REG.Quantity(temperature, initial_units['temp']).to('K')
        tau_units = UNIT_REG[initial_units['time']].to(final_units['time'])
        tau *= tau_units.magnitude

        self.x = 1/temperature.magnitude
        self.y = np.log(tau)
        y_units = 'ln({})'.format(tau_units.units)

        popt, perr, r_sq, fitted = fit(self.x, self.y, equation_linear)

        self.r_sq = r_sq
        self.fitted = fitted

        self.a = popt[0]
        self.a_err = perr[0]

        self.b = popt[1]
        self.b_err = perr[1]

        self.activation_energy = self.a * 8.314
        self.activation_energy_err = self.a_err * 8.314
        act_en_units = UNIT_REG['K*J/mol/K'].to(final_units['energy']) # pyright: ignore[reportAttributeAccessIssue]
        self.activation_energy *= act_en_units.magnitude
        self.activation_energy_err *= act_en_units.magnitude

        self.prefactor = np.exp(self.b)
        self.prefactor_err = abs(self.prefactor*self.b_err)

        if verbose:
            print_value(
                name = "a",
                value = self.a,
                error = self.a_err,
                units = 'K / ' + y_units
            )
            print_value(
                name = "b",
                value = self.b,
                error = self.b_err,
                units = y_units
            )
            print(PRINT_TEMPLATE_RSQ.format(r_sq))
            print_value(
                name = "Activation Energy",
                value = self.activation_energy,
                error = self.activation_energy_err,
                units = act_en_units.units
            )
            print_value(
                name = "Prefactor",
                value = self.prefactor,
                error = self.prefactor_err,
                units = tau_units.units
            )

class fit_arrhenius_yerr:
    def __init__(self,
        temperature:np.ndarray, tau:np.ndarray,
        tau_err:np.ndarray,
        initial_units={'temp':'K', 'time':'ps'},
        final_units={'time':'s', 'energy':'kJ/mol'},
        verbose=True) -> None:

        if verbose:
            print("Fitting Arrhenius :")

        temperature = UNIT_REG.Quantity(temperature, initial_units['temp']).to('K')
        tau_units = UNIT_REG[initial_units['time']].to(final_units['time'])
        tau *= tau_units.magnitude
        tau_err *= tau_units.magnitude

        self.x = 1/temperature.magnitude # pyright: ignore[reportAttributeAccessIssue]
        self.y = np.log(tau)
        y_units = 'ln({})'.format(tau_units.units)
        self.y_err = abs(tau_err/tau)

        popt, perr, r_sq, fitted = fit(self.x, self.y, equation_linear, sigma=self.y_err, absolute_sigma=True)

        self.r_sq = r_sq
        self.fitted = fitted

        self.a = popt[0]
        self.a_err = perr[0]

        self.b = popt[1]
        self.b_err = perr[1]

        self.activation_energy = self.a * 8.314
        self.activation_energy_err = self.a_err * 8.314
        act_en_units = UNIT_REG['K*J/mol/K'].to(final_units['energy']) # pyright: ignore[reportAttributeAccessIssue]
        self.activation_energy *= act_en_units.magnitude
        self.activation_energy_err *= act_en_units.magnitude

        self.prefactor = np.exp(self.b)
        self.prefactor_err = abs(self.prefactor*self.b_err)

        if verbose:
            print_value(
                name = "a",
                value = self.a,
                error = self.a_err,
                units = 'K / ' + y_units
            )
            print_value(
                name = "b",
                value = self.b,
                error = self.b_err,
                units = y_units
            )
            print(PRINT_TEMPLATE_RSQ.format(r_sq))
            print_value(
                name = "Activation Energy",
                value = self.activation_energy,
                error = self.activation_energy_err,
                units = act_en_units.units
            )
            print_value(
                name = "Prefactor",
                value = self.prefactor,
                error = self.prefactor_err,
                units = tau_err.units
            )

class fit_exponential_convergence():
    def __init__(self, x:np.ndarray, y:np.ndarray, alpha_ini=400, gamma_ini=0.2, verbose=True):

        if verbose:
            print("Fitting exponential convergence :")

        if int(x[0]) == 0: # skipping early values that would cause division by 0
            self.x = x[1:]
            self.y = y[1:]
        else:
            self.x = x
            self.y = y

        xconv_ini = max(y)

        popt, perr, self.r_sq, self.fitted = fit(self.x, self.y, equation_exponential_convergence, [xconv_ini, alpha_ini, gamma_ini])

        self.xconv = popt[0]
        self.xconv_err = perr[0]
        self.alpha = popt[1]
        self.alpha_err = perr[1]
        self.gamma = popt[2]
        self.gamma_err = perr[2]

        if verbose:
            print_value(
                name = "x_conv",
                value = self.xconv,
                error = self.xconv_err,
                units = ''
            )
            print_value(
                name = "Alpha",
                value = self.alpha,
                error = self.alpha_err,
                units = ''
            )
            print_value(
                name = "Gamma",
                value = self.gamma,
                error = self.gamma_err,
                units = ''
            )
            print(PRINT_TEMPLATE_RSQ.format(self.r_sq))

    def predict_x(self, y):
        return ( -self.alpha / (y-self.xconv) )**(1/self.gamma)

class fit_tg_window_linear_intersection():
    def __init__(self, x:np.ndarray, y:np.ndarray, tg:float, window=50, verbose=True):
        if verbose:
            print("Fitting extremum linear intersection:")

        sort_index = x.argsort()
        self.x = x[sort_index]
        self.y = y[sort_index]

        window_start = ( np.abs(self.x - (tg-window)) ).argmin()
        window_end = ( np.abs(self.x - (tg+window)) ).argmin()

        self.x_start = self.x[:window_start]
        self.y_start = self.y[:window_start]
        self.x_end = self.x[window_end:]
        self.y_end = self.y[window_end:]
        self.x_mid = self.x[window_start:window_end]
        self.y_mid = self.y[window_start:window_end]

        popt, perr, r_sq, fitted = fit(self.x_start, self.y_start, equation_linear)
        self.rsq_start = r_sq
        self.fitted_start = fitted
        self.a_start = popt[0]
        self.a_err_start = perr[0]
        self.b_start = popt[1]
        self.b_err_start = perr[1]

        popt, perr, r_sq, fitted = fit(self.x_end, self.y_end, equation_linear)
        self.rsq_end = r_sq
        self.fitted_end = fitted
        self.a_end = popt[0]
        self.a_err_end = perr[0]
        self.b_end = popt[1]
        self.b_err_end = perr[1]

        self.intersection = (self.b_end - self.b_start) / (self.a_start - self.a_end)
        self.intersection_err = np.sqrt( self.b_err_start**2 + self.b_err_end**2 + (self.intersection**2)*self.a_err_start**2 + (self.intersection**2)*self.a_err_end**2 ) / (self.a_start-self.a_end)

        if verbose:
            print_value(
                name = "intersection",
                value = self.intersection,
                error = self.intersection_err,
                units = 'K'
            )
            print_value(
                name = "intersection",
                value = self.intersection-273.15,
                error = self.intersection_err,
                units = 'C'
            )
            print("\tStart linear fit")
            print_value(
                name = "a_start",
                value = self.a_start,
                error = self.a_err_start,
                units = 'cm3/gK'
            )
            print_value(
                name = "b_start",
                value = self.b_start,
                error = self.b_err_start,
                units = 'cm3/g'
            )
            print(PRINT_TEMPLATE_RSQ.format(self.rsq_start))

            print("\tEnd linear fit")
            print_value(
                name = "a_end",
                value = self.a_end,
                error = self.a_err_end,
                units = 'cm3/gK'
            )
            print_value(
                name = "b_start",
                value = self.b_end,
                error = self.b_err_end,
                units = 'cm3/g'
            )
            print(PRINT_TEMPLATE_RSQ.format(self.rsq_end))

    def predict(self, region:str):
        if region=="start":
            y = self.x*self.a_start + self.b_start
        if region=="end":
            y = self.x*self.a_end + self.b_end

        return y

class fit_wlffv:
    def __init__(self,
        temperature:np.ndarray, ffv:np.ndarray,
        ffvg_ini=2.5, cte_ini=1, tg=450,
        initial_units={'temp':'K'},
        final_units={'temp':'K'},
        verbose=True) -> None:

        if verbose:
            print("Fitting WLF :")

        temperature = UNIT_REG.Quantity(temperature, initial_units['temp']).to('K')

        self.x = np.array(temperature.magnitude) - tg
        self.y = np.array(ffv) / 100
        ffvg_ini = ffvg_ini / 100
        cte_ini = cte_ini * 3 * 100

        popt, perr, r_sq, fitted = fit(self.x, self.y, equation_wlffv, [ffvg_ini, cte_ini])

        self.r_sq = r_sq
        self.fitted = fitted

        self.ffvg = popt[0] * 100
        self.ffvg_err = perr[0] * 100

        self.cte = popt[1] / 3
        self.cte_err = perr[1] / 3

        # self.tg = popt[2]
        # self.tg_err = perr[2]

        if verbose:
            print_value(
                name = "FFVg",
                value = self.ffvg,
                error = self.ffvg_err,
                units = '%'
            )
            print_value(
                name = "CTE",
                value = self.cte * 1e6,
                error = self.cte_err * 1e6,
                units = '1/K * 1e-6'
            )

            print(PRINT_TEMPLATE_RSQ.format(r_sq))

class fit_diffusion:
    def __init__(self,
        x:np.ndarray, y:np.ndarray,
        initial_units={'dist':'angstrom', 'time':'ns'},
        final_units  ={'dist':'cm',       'time':'s'},
        verbose=True) -> None:

        if verbose:
            print("Fitting Diffusion Coefficient :")

        self.x = x
        self.y = y

        popt, perr, r_sq, fitted = fit(self.x, self.y, equation_linear)

        self.r_sq = r_sq
        self.fitted = fitted

        self.a = popt[0]
        self.a_err = perr[0]

        self.b = popt[1]
        self.b_err = perr[1]

        dist_units = UNIT_REG[initial_units['dist']].to(final_units['dist'])
        time_units = UNIT_REG[initial_units['time']].to(final_units['time'])

        self.diff_coef = self.a / 6 * dist_units.magnitude**2 / time_units.magnitude
        self.diff_err = self.a_err / 6 * dist_units.magnitude**2 / time_units.magnitude

        if verbose:
            print_value(
                name = "a",
                value = self.a,
                error = self.a_err,
                units = '{}^2 / {}'.format(initial_units['dist'], initial_units['time'])
            )
            print_value(
                name = "b",
                value = self.b,
                error = self.b_err,
                units = initial_units['dist']
            )
            print(PRINT_TEMPLATE_RSQ.format(r_sq))
            print_value(
                name = "Diffusion coefficient",
                value = self.diff_coef,
                error = self.diff_err,
                units = '{}^2 / {}'.format(final_units['dist'], final_units['time'])
            )