import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.interpolate import make_interp_spline  # Import the necessary function

# Generate sample RNA-seq data with random values
num_samples = 36  # Total number of samples (3 replicates x 2 genotypes x 6 time points)
time_points = ['E12.5', 'P28', 'P56', 'P84', 'P119', 'P150'] * 6
genotypes = ['WT', 'WT', 'WT', 'G93A', 'G93A', 'G93A'] * 6
replicates = [1, 2, 3] * 12
expression_values = np.random.rand(num_samples)  # Random expression values between 0 and 1

# Create a DataFrame with the sample data
data = {
    'Time': time_points,
    'Genotype': genotypes,
    'Replicate': replicates,
    'Expression': expression_values
}
df = pd.DataFrame(data)

# Calculate mean expression for each genotype at each time point
mean_expression = df.groupby(['Time', 'Genotype'])['Expression'].mean().unstack()

# Check for and handle invalid values (e.g., replace with mean)
mean_expression.fillna(mean_expression.mean(), inplace=True)

# Convert time points to numerical values for smoothing
time_numeric = np.arange(len(mean_expression.index))

# Smooth the lines using spline interpolation
spl_wt = make_interp_spline(time_numeric, mean_expression['WT'], k=3)  # k=3 for cubic spline
spl_g93a = make_interp_spline(time_numeric, mean_expression['G93A'], k=3)
time_smooth = np.linspace(time_numeric.min(), time_numeric.max(), 300)  # More points for smoother curve
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