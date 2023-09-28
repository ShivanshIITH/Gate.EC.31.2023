# Define the SNR in dB
snr_db = 61.96

# Constants
quantization_noise_db = 1.76
noise_bandwidth_factor_db = 6.02

# Calculate the resolution in bits
resolution_bits = (snr_db - quantization_noise_db) / noise_bandwidth_factor_db

print(f"The resolution of the ADC is  {resolution_bits:.2f} bits")

