import numpy as np

Given_SNR = 61.96 
# Signal parameters
amplitude = 1  # Amplitude of the sinusoidal signal
frequency = 10  # Frequency of the sinusoidal signal in Hz
sampling_frequency = 1000  # Sampling frequency in Hz
duration = 1  # Duration of the signal in seconds
num_samples = int(sampling_frequency * duration)
time = np.arange(0, duration, 1/sampling_frequency)

# Generating the full-scale sinusoidal input signal
input_signal = amplitude * np.sin(2 * np.pi * frequency * time)

# Calculate Signal Power (Ps) by averaging the square of the signal
signal_power = np.mean(input_signal**2)

# Initialize the simulated resolution as None
resolution_simulated = None

# Iterate through different quantization bit values
for n in range(1, 20):
    # Calculate the number of quantization levels
    quantization_levels = 2**n
    
    # Calculate the quantization step size
    quantization_step = 2 * amplitude / quantization_levels
    
    # Quantize the input signal
    quantized_signal = np.round(input_signal / quantization_step) * quantization_step
    
    # Calculate the quantization error signal
    quantization_error_signal = input_signal - quantized_signal
    
    # Calculate Signal Power (Pe) of quantization error signal
    signal_power_quantization_error = np.mean(quantization_error_signal**2)
    
    # Calculate Signal Power in decibels (dB)
    SNR_simulated = 10 * np.log10(signal_power / signal_power_quantization_error)
    
    # Check if the simulated SNR is close to the given SNR with a tolerance of 1 dB
    if np.isclose(Given_SNR, SNR_simulated, atol=1):
        resolution_simulated = n
        break

print(f"Simulated Resolution: {resolution_simulated} bits")

