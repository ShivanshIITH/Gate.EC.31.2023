import numpy as np
#include <stdio.h>
#include <math.h>

// Constants
#define PI 3.14159265358979323846

// Signal parameters
double Given_SNR = 61.96;
double amplitude = 1.0;            // Amplitude of the sinusoidal signal
double frequency = 10.0;           // Frequency of the sinusoidal signal in Hz
double sampling_frequency = 1000.0; // Sampling frequency in Hz
double duration = 1.0;             // Duration of the signal in seconds
int num_samples;

int main() {
    double time_step;
    int resolution_simulated = -1;

    // Calculate the number of samples
    num_samples = (int)(sampling_frequency * duration);
    time_step = 1.0 / sampling_frequency;

    // Generate the full-scale sinusoidal input signal and calculate Signal Power (Ps)
    double signal_power = 0.0;
    for (int i = 0; i < num_samples; i++) {
        double time = i * time_step;
        double input_signal = amplitude * sin(2 * PI * frequency * time);
        signal_power += input_signal * input_signal;
    }

    // Iterate through different quantization bit values
    for (int n = 1; n < 20; n++) {
        // Calculate the number of quantization levels
        int quantization_levels = 1 << n;

        // Calculate the quantization step size
        double quantization_step = 2 * amplitude / quantization_levels;

        // Quantize the input signal and calculate Signal Power (Pe) of quantization error signal
        double signal_power_quantization_error = 0.0;
        for (int i = 0; i < num_samples; i++) {
            double time = i * time_step;
            double input_signal = amplitude * sin(2 * PI * frequency * time);
            double quantized_signal = round(input_signal / quantization_step) * quantization_step;
            double quantization_error_signal = input_signal - quantized_signal;
            signal_power_quantization_error += quantization_error_signal * quantization_error_signal;
        }

        // Calculate Signal Power in decibels (dB)
        double SNR_simulated = 10 * log10(signal_power / signal_power_quantization_error);

        // Check if the simulated SNR is close to the given SNR with a tolerance of 1 dB
        if (fabs(Given_SNR - SNR_simulated) < 1.0) {
            resolution_simulated = n;
            break;
        }
    }

    printf("Simulated Resolution: %d bits\n", resolution_simulated);

    return 0;
}

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

