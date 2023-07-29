from scipy.io import loadmat
from scipy.integrate import solve_ivp
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import basinhopping
import time

# Load the MATLAB files
data_SEIRV_fit = loadmat('data_SEIRV_fit.mat')

# Convert MATLAB arrays to numpy arrays
CRNA2 = data_SEIRV_fit['cRNA2'].flatten()
F2 = data_SEIRV_fit['F2'].flatten()

# Calculate V
V = CRNA2 * F2

# Split 1 week after first day of vaccine (12/11/2020)
split = 78
V = V[:split]

# Define the time span
tspan = np.arange(1, len(V) + 1)

# Define constants from the MATLAB code
tau0 = 189.6
Q0 = 2.5
T0 = 20
A = 3.624836409841919
B = 0.020222716119084
C = 4.466530666828714
D = 16.229757918464635

def getDecay(t):
    # Compute temperature-adjusted decay rate of viral RNA
    T = A * np.sin(B * t - C) + D
    tau = tau0 * Q0 ** (-(T - T0) / 10)
    k = np.log(2) / tau
    return k

def SEIRV(t, y, param):
    # Parameters to be fit
    lambda_, alpha, beta, _ = param

    # Other parameters
    traveltime = 18  # hours
    k = getDecay(t)
    eta = 1 - np.exp(-k * traveltime)
    sigma = 1 / 3
    gamma = 1 / 8

    # Variables
    S, E, I, R, V, E_cumulative = y

    # Compute derivatives
    dS = -lambda_ * S * I
    dE = lambda_ * S * I - sigma * E
    dI = sigma * E - gamma * I
    dR = gamma * I
    dV = alpha * beta * (1 - eta) * I
    dE_cumulative = lambda_ * S * I

    return [dS, dE, dI, dR, dV, dE_cumulative]



# Parameters to be fitted
# lambda: transmission rate per day per person
# alpha: fecal load in gram
# beta: viral shedding in stool in viral RNA copies per gram
# E0: initial exposed population

# Constants
# sigma: rate of movement from E to I (inverse of duration of incubation period)
# gamma: rate of movement from I to R (inverse of duration of infectiousness)
# N0: total population served by DITP
# eta: fraction of viruses that survive the travel time to the wastewater treatment plant
# V0: initial virus concentration in wastewater

sigma = 1 / 3
gamma = 1 / 8
N0 = 2300000
eta = 1 - np.exp(-getDecay(1) * 18)  # traveltime = 18 hours

def seirv_model(y, t, lambda_, alpha, beta, E0):
    S, E, I, R, V, E_cumulative = y
    dS = -lambda_ * S * I
    dE = lambda_ * S * I - sigma * E
    dI = sigma * E - gamma * I
    dR = gamma * I
    dV = alpha * beta * (1 - eta) * I
    dE_cumulative = lambda_ * S * I
    return [dS, dE, dI, dR, dV, dE_cumulative]

def simulate_seirv_model(param, tspan, V0):
    lambda_, alpha, beta, E0 = param
    I0 = V0 / (alpha * beta * (1 - eta))
    R0 = 0
    S0 = N0 - (E0 + I0 + R0)
    ICs = [S0, E0, I0, R0, V0, E0]
    sol = odeint(seirv_model, ICs, tspan, args=(lambda_, alpha, beta, E0))
    return sol

def sse_obj_func(param, tspan, data):
    sol = simulate_seirv_model(param, tspan, data[0])
    cumVirus = sol[:, 4]
    dailyVirus = np.diff(cumVirus)
    temp = np.log10(data[1:]) - np.log10(np.abs(dailyVirus))
    adiff = temp[~np.isnan(temp)]
    return np.sum(adiff ** 2)

# Define the bounds of the parameters
param_bounds = [(0, 1E-4), (51, 796), (4.48526e7, 4.48526e7), (10, 5000)]

# The objective function is the same (sse_obj_func), but we need to flatten it to 1D for basinhopping
flat_sse_obj_func = lambda x: sse_obj_func(x, tspan=tspan, data=V)

# Initialize a random starting point
#np.random.seed(0)
#print(time.time())
init_guess = np.random.uniform(low=[b[0] for b in param_bounds], high=[b[1] for b in param_bounds]) #solution from paper: [9.06e-08, 360, 4.48526e7, 1182]

# Define a function to check if a point is within bounds
def in_bounds(x, bounds):
    return np.all([low <= xi <= high for xi, (low, high) in zip(x, bounds)])

# Define a function for generating a new random point within bounds
def random_within_bounds(x, bounds):
    return np.array([np.random.uniform(low, high) for (low, high) in bounds])

# Define the function for generating a new trial step
def new_trial_step(x):
    return random_within_bounds(x, param_bounds)

# Define the function for accepting or rejecting the new trial step
def accept_trial(f_new, x_new, f_old, x_old):
    return in_bounds(x_new, param_bounds)



# Run basinhopping optimization
res = basinhopping(sse_obj_func, init_guess, niter=25, T=1.0, stepsize=0.5,
                   minimizer_kwargs={'method': 'L-BFGS-B', 'bounds': param_bounds, 'args': (tspan, V)},
                   take_step=new_trial_step, accept_test=accept_trial, seed=0)

# Print the results
res.x, res.fun



# PLOTTING
import matplotlib.pyplot as plt


# Load newRepCases2 from MATLAB file
newRepCases2 = data_SEIRV_fit['newRepCases2'].flatten()

# Define the time span for prediction
tspan_pred = np.arange(1, len(V) + 1 + 30)  # extend for 30 days

# Simulate the SEIRV model with the optimized parameters
sol_opt = simulate_seirv_model(res.x, tspan_pred, V[0])

# Plot the estimated daily incidence and the observed daily incidence
plt.figure(figsize=(12, 6))

# Plot the observed daily incidence
plt.plot(tspan_pred, np.log10(newRepCases2[:split+30]), '.', markersize=10, color='black', label='Observed daily incidence')

# Plot the estimated daily incidence using the optimized parameters
plt.plot(tspan_pred[:-1], np.log10(np.diff(sol_opt[:, 5])), '-', color='green', label='Estimated daily incidence (optimized)')

plt.xlabel('Days')
plt.ylabel('Log10 daily incidence')
plt.legend()
plt.grid(True)
plt.title('Estimated and observed daily incidence')
plt.show()
