import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.interpolate import make_interp_spline  # Import the necessary function

# Generate dynamic profile using the mean values from the RNA-seq
num_samples = 36  # Total number of samples (3 replicates x 2 genotypes x 6 time points)
time_points = ['E12.5', 'P28', 'P56', 'P84', 'P119', 'P150'] * 6
genotypes = ['WT', 'WT', 'WT', 'G93A', 'G93A', 'G93A'] * 6
replicates = [1, 2, 3] * 12
expression_values = np.random.rand(num_samples)  # Expression values 

# Create a DataFrame with the sample data
data = {
    'Time': time_points,
    'Genotype': genotypes,
    'Replicate': replicates,
    'Expression': expression_values
}
df = pd.DataFrame(data)

# Mean expression for each genotype at each time point
mean_expression = df.groupby(['Time', 'Genotype'])['Expression'].mean().unstack()

# Check for and handle invalid values (e.g., replace with mean)
mean_expression.fillna(mean_expression.mean(), inplace=True)

# Convert time points to numerical values
time_numeric = np.arange(len(mean_expression.index))

# Spline interpolation
spl_wt = make_interp_spline(time_numeric, mean_expression['WT'], k=3)  # k=3 for cubic spline
spl_g93a = make_interp_spline(time_numeric, mean_expression['G93A'], k=3)
time_smooth = np.linspace(time_numeric.min(), time_numeric.max(), 300)  
expression_smooth_wt = spl_wt(time_smooth)
expression_smooth_g93a = spl_g93a(time_smooth)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(time_smooth, expression_smooth_wt, color='green', label='WT')
plt.plot(time_smooth, expression_smooth_g93a, color='red', label='G93A')
plt.xticks(time_numeric, mean_expression.index)  # Set x-ticks back to original time points
plt.xlabel('Time')
plt.ylabel('Mean Expression')
plt.title('Mean Expression over Time for WT and G93A')
plt.legend()
plt.grid(True)
plt.show()

#####################

import matplotlib.pyplot as plt
import numpy as np

# Time points
time_points = np.array([12.5, 28, 56, 84, 119, 150])

# --- Late Gradual Reduction Response ---
y_late_gradual = np.array([1, 1, 0.8, 0.6, 0.4, 0.2])  # Adjust these parameters for the curve

plt.figure(figsize=(5, 3))
plt.plot(time_points, y_late_gradual, color='black', linewidth=3)
plt.title('Late Gradual Reduction Response')
plt.xticks(time_points)
plt.yticks([])  # Hide y-axis ticks
plt.ylim(0, 1.2)  # Set y-axis limits
plt.show()

# --- Early Rapid Reduction Response ---
y_early_rapid = np.array([1, 1, 0.2, 0.2, 0.2, 0.2])  

plt.figure(figsize=(5, 3))
plt.plot(time_points, y_early_rapid, color='black', linewidth=3)
plt.title('Early Rapid Reduction Response')
plt.xticks(time_points)
plt.yticks([])
plt.ylim(0, 1.2)
plt.show()

# --- Long Oscillatory Response ---
y_long_oscillatory = np.array([1, 1, 0.2, 0.2, 0.2, 1])  

plt.figure(figsize=(5, 3))
plt.plot(time_points, y_long_oscillatory, color='black', linewidth=3)
plt.title('Long Oscillatory Response')
plt.xticks(time_points)
plt.yticks([])
plt.ylim(0, 1.2)
plt.show()

# --- Late Gradual Response (^) ---
y_late_gradual_up = np.array([0.2, 0.3, 0.4, 0.6, 0.8, 1]) 

plt.figure(figsize=(5, 3))
plt.plot(time_points, y_late_gradual_up, color='black', linewidth=3)
plt.title('Late Gradual Response (^)')
plt.xticks(time_points)
plt.yticks([])
plt.ylim(0, 1.2)
plt.show()

# --- Rapid Long Response (v) ---
y_rapid_long_down = np.array([0.2, 1, 1, 0.8, 0.6, 0.4]) 

plt.figure(figsize=(5, 3))
plt.plot(time_points, y_rapid_long_down, color='black', linewidth=3)
plt.title('Rapid Long Response (v)')
plt.xticks(time_points)
plt.yticks([])
plt.ylim(0, 1.2)
plt.show()

# --- Rapid Long Response (^) ---
y_rapid_long_up = np.array([0.2, 0.8, 0.8, 0.8, 1, 1])  

plt.figure(figsize=(5, 3))
plt.plot(time_points, y_rapid_long_up, color='black', linewidth=3)
plt.title('Rapid Long Response (^)')
plt.xticks(time_points)
plt.yticks([])
plt.ylim(0, 1.2)
plt.show()