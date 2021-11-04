import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#!/usr/bin/env python3
import korali
import numpy as np
from scipy.stats import t

def lStandardNormalGaussian(s):
  v = s["Parameters"][0]
  s["logP(x)"] = -0.5 * v * v
  s["grad(logP(x))"] = [-v]
  s["H(logP(x))"] = [[-1.0]]

def standardNormalGaussianDensity(x):
    mu = 0.0
    sigma = 1.0
    return 1.0 / (np.sqrt(2.0 * np.pi) * sigma) * np.exp(-0.5 * ((x-mu)/sigma)**2)

dim = 1

numSamples = 5000
numExperiments = 4
samples = np.zeros(shape=(numExperiments, numSamples))
for experiment in range(numExperiments):
    # Set Random Seed
    e = korali.Experiment()
    e["Random Seed"] = 0xC0FFEE + experiment

    #######################################################################################
    ############################# File/Console output settings ############################
    #######################################################################################
    e["File Output"]["Frequency"] = 1000
    e["Console Output"]["Frequency"] = 1000
    # e["Console Output"]["Verbosity"] = "Detailed"

    #######################################################################################
    ################################ Variable Initilization ###############################
    #######################################################################################
    e["Variables"][0]["Name"] = "X"
    e["Variables"][0]["Initial Mean"] = np.random.rand()
    e["Variables"][0]["Initial Standard Deviation"] = 1.0

    #######################################################################################
    ################################## Sampler version ####################################
    #######################################################################################
    e["Problem"]["Type"] = "Sampling"
    e["Solver"]["Type"] = "Sampler/HMC"
    e["Solver"]["Version"] = 'Static'
    e["Solver"]["Use Diagonal Metric"] = False #for Version != Static || Version != Riemannian
    e["Solver"]["Inverse Regularization Parameter"] = 1.0 #for Version = Riemannian || Version == Riemannian Const
    e["Problem"]["Probability Function"] = lStandardNormalGaussian

    #######################################################################################
    ############################# Number of samples settings ##############################
    #######################################################################################
    e["Solver"]["Burn In"] = 200
    e["Solver"]["Termination Criteria"]["Max Samples"] = numSamples
    e["Solver"]["Initial Fast Adaption Interval"] = 75
    e["Solver"]["Final Fast Adaption Interval"] = 50
    e["Solver"]["Initial Slow Adaption Interval"] = 25

    #######################################################################################
    ############################# Integration length settings #############################
    #######################################################################################
    # General Step Size settings
    e["Solver"]["Step Size"] = 0.05
    e["Solver"]["Step Size Jitter"] = 0.0
    e["Solver"]["Num Integration Steps"] = 3
    # Adaptive Step Size settings
    e["Solver"]["Use Adaptive Step Size"] = False
    e["Solver"]["Target Acceptance Rate"] = 0.7 #for Use Adaptive Step Size = True
    e["Solver"]["Target Integration Time"] = 1.0 #for Use Adaptive Step Size = True
    e["Solver"]["Adaptive Step Size Speed Constant"] = 0.05 #for Use Adaptive Step Size = True
    e["Solver"]["Adaptive Step Size Stabilization Constant"] = 10.0 #for Use Adaptive Step Size = True
    e["Solver"]["Adaptive Step Size Schedule Constant"] = 0.75 #for Use Adaptive Step Size = True
    # NUTS settings
    e["Solver"]["Use NUTS"] = False
    e["Solver"]["Max Depth"] = 10 #for Use NUTS = True
    # Implicit leapfrog settings
    e["Solver"]["Max Integration Steps"] = 100 #for Version = Riemannian
    
    #######################################################################################
    ################################### Run Experiment ####################################
    #######################################################################################
    k = korali.Engine()
    k.run(e)
    experimentSamples = e["Results"]["Sample Database"]
    experimentSamples = np.reshape(experimentSamples, (-1, dim)).flatten()
    samples[experiment, :] = experimentSamples


#######################################################################################
###################################### Plotting #######################################
#######################################################################################

# Plot a single histogram #############################################################
x_range = 5.0
plt.grid(linestyle=":")
plt.xlim(-x_range, x_range)
plt.ylim(0.0, 0.5)
plt.hist(samples[0, :], bins=50, density=True)

numPoints = 1000
xVals = np.linspace(-5.0, 5.0, numPoints)
stdNormalDistrDensity = np.array([standardNormalGaussianDensity(x) for x in xVals])

plt.xlabel("Position $q$")
plt.ylabel("Rel. Occurrence")
plt.plot(xVals, stdNormalDistrDensity, 'r-')
plt.legend(["PDF", "Samples"])

plt.savefig("samples.png", dpi=200)

plt.clf()

# Plot Convergence ####################################################################
cumSumSamples = np.cumsum(samples, axis=1)
print("cumSumSamples = ", cumSumSamples)
anaExpectation = 0.0
numSamplesRange = np.arange(1, numSamples+1)
empiricalExpectation = np.zeros(shape=(numExperiments, numSamples))
for experiment in range(numExperiments):
    empiricalExpectation[experiment, :] = cumSumSamples[experiment, :] / numSamplesRange
print("empiricalExpectation = ", empiricalExpectation)

errorExpectation = np.abs(anaExpectation - empiricalExpectation)
averageErrorExpectation = np.mean(errorExpectation, axis=0)
print("errorExpectation = ", errorExpectation)
errorIdeal = np.array([1.0/np.sqrt(n) for n in numSamplesRange])

# Least Squares fit
A = np.zeros(shape=(numExperiments*numSamples, 2))
y = np.zeros(shape=(numExperiments*numSamples, 1))
for i in range(numExperiments):
    A[i*numSamples, 0] = 1.0
    for j in range(numSamples):
        A[i*numSamples + j, 1] = np.log(j+1)
        y[i*numSamples + j] = np.log(errorExpectation[i, j])

M = np.matmul(np.transpose(A), A)
ATy = np.matmul(np.transpose(A), y)
print("A = ", A)
print("y = ", y)
print("M = ", M)
print("ATy =", ATy)
x = np.linalg.solve(M, ATy)
print("x = ", x)

LSQFit = np.array([(np.exp(x[0]) * n**x[1]) for n in numSamplesRange])

# Calculate confidence interval
alpha = 0.05
confidenceValue = 1.0 - alpha / 2.0
studentTCoeff = t.ppf(confidenceValue, df=numExperiments-1)
ci = studentTCoeff * np.std(errorExpectation, axis=0) / np.sqrt(numExperiments)

# # for experiment in range(numExperiments):
# #     plt.loglog(numSamplesRange, errorExpectation[experiment, :])
# plt.loglo
g(numSamplesRange, averageErrorExpectation)

fig, ax = plt.subplots()
ax.loglog(numSamplesRange, errorIdeal, 'g--')
ax.loglog(numSamplesRange, LSQFit, 'r--')
ax.loglog(numSamplesRange, averageErrorExpectation)
ax.fill_between(numSamplesRange, (averageErrorExpectation-ci), (averageErrorExpectation+ci), color='b', alpha=.1)
plt.grid(linestyle=":")
plt.xlabel("Number of Samples $N$")
plt.ylabel("Error $|\\mathbb{E}[Q]-\hat{E}[Q]|$")
plt.legend(["i.i.d. conv. $\propto \\frac{1}{\sqrt{N}}$", "LSQ Fit", "$\hat{E}[Error]$", str((100 * (1.0 - alpha))) + "% CI"])
plt.savefig("gaussian_convergence.png", dpi=200)
