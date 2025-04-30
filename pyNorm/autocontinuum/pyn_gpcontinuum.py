import numpy as np
from pyNorm.io import read_inorm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
import matplotlib.pyplot as plt

# Fit a continuum to the data stored in spec['vel'] and spec['flux'] with errors spec['eflux'] using a Gaussian process. The continuum should be smooth over 10s of pixels.  Mask the data around velocities <200 km/s and around other absorption lines. Iteratively identify additional absorption lines and mask them. Return the continuum and the 68% confidence intervals around it. Plot the original spectrum and the continuum fit.

# Read data from IDL save file
# input_filename = 'Data/CIV1548.2i_o.save'
# input_filename = 'Data/sk-67D211_OVI1031.9_fusei_o.save'
input_filename = 'Data/sk-67D211_CIV1548.2_e140mi_o.save'
spec = read_inorm(input_filename)

# Mask really bad pixels and pixels with large velocities
max_vel = 1000.
bdmask = ((np.abs(np.roll(spec['flux'],-1)-spec['flux']) > 5.*np.median(spec['flux'])) | (np.abs((np.roll(spec['eflux'],-1)-spec['eflux'])) > 5.*np.median(spec['eflux']))) | (np.abs(spec['vel']) > max_vel)

# Prepare the data
flx_factor = 10.**-np.floor(np.log10(np.median(spec['flux'][~bdmask])))
X = spec['vel'][~bdmask].reshape(-1, 1)
y = spec['flux'][~bdmask]*flx_factor
sigma = spec['eflux'][~bdmask]*flx_factor

# Mask the data
vel_exclusion = 400.
vel_exclusion = [-100,400]
if len(vel_exclusion) > 1:
    vel_mask = (X > vel_exclusion[0]) & (X < vel_exclusion[1])
else:
    vel_mask = (np.abs(X) < vel_exclusion) 
X_masked = X[~vel_mask].reshape(-1,1)
y_masked = y[~vel_mask[:, 0]]
sigma_masked = sigma[~vel_mask[:, 0]]

# Define the kernel for the Gaussian process
len_scale = 2.  # Length scale of the kernel
kernel = RBF(length_scale=len_scale)

# Create the Gaussian process regressor using uncertainties for weighting.
gp = GaussianProcessRegressor(kernel=kernel, 
                              alpha=sigma_masked ** 2,
                              n_restarts_optimizer=9)

# Fit the Gaussian process to the masked data
gp.fit(X_masked, y_masked)

# Predict the continuum values
continuum, std = gp.predict(X, return_std=True)

# Iteratively mask potential absorption lines
for i in range(10):  # Modify the number of iterations as needed
    # # Find potential absorption lines
    # mask = ((y - continuum) < -1.25 * sigma) | \
    #     ((np.abs(X) < vel_exclusion)[:,0] )

    contiguous_absorption_indices = []
    previous_index = -2

    # Find potential absorption lines
    for i, value in enumerate(y - continuum):
        if value < -0.5 * sigma[i] or \
            vel_mask[i] == True:
            # Check if the indices are contiguous
            if previous_index == i - 1:  
                contiguous_absorption_indices.append(i)
            previous_index = i

    # Create a mask for the contiguous absorption lines
    mask = np.zeros_like(y, dtype=bool)
    mask[contiguous_absorption_indices] = True

    # Mask the potential absorption lines
    X_masked = X[~mask,:].reshape(-1, 1)
    y_masked = y[~mask]
    sigma_masked = sigma[~mask]

    # Fit the Gaussian process to the updated masked data 
    gp = GaussianProcessRegressor(kernel=kernel, 
                                  alpha=sigma_masked ** 2,
                                   n_restarts_optimizer=9)
    gp.fit(X_masked, y_masked)

    # Predict the updated continuum values
    # continuum, std = gp.predict(X, return_std=True)
    continuum, cov = gp.predict(X, return_cov=True)

# Plot the original spectrum and the continuum fit
plt.clf()
plt.plot(X, y, drawstyle='steps-mid')
plt.plot(X, sigma, drawstyle='steps-mid',
         color='gray', alpha=0.5)
plt.plot(X[mask,:], 
         y[mask], 
         'o', markersize=2, color='red',
         label='Masked features')
plt.plot(X, continuum, color='green')

# Plot the 68% confidence intervals around the continuum
plt.fill_between(X.flatten(), 
                 continuum - 1.96 * std, continuum + 1.96 * std, 
                 color='gray', alpha=0.5)

# Limits
plt.ylim(np.array([-0.1, 1.5])*np.max(continuum))
# Add labels and title to the plot
plt.xlabel('Velocity (km/s)')
plt.ylabel('Flux')
plt.title('Continuum Fit')

# Add a legend
plt.legend()

# Plot zero level
plt.axhline(0, color='black', lw=0.5,ls='--')



# Show the plot
plt.show()

# Return the continuum and the 68% confidence intervals
# continuum, continuum - 1.96 * std, continuum + 1.96 * std