import numpy as np
data = np.fromfile("results_outAddlie.bin", dtype=np.float32)

# Reshape into (rows, cols*2)
data = data.reshape(rows, cols * 2)
print("shape: {data.shape}")