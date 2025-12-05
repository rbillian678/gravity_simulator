CUDA N-Body Simulation (Shared Memory Version)

This project implements a simple N-body gravitational simulation using CUDA.
Each body interacts with every other body, and the computation is parallelized across GPU threads.

How it Works
	•	Each thread updates one body.
	•	Bodies are processed in tiles of 256 and stored in shared memory.
	•	This reduces repeated global memory reads and improves performance.
	•	After forces are accumulated, each body updates its acceleration, velocity, and position.

The simulation loads initial body data from a CSV file:
name, x, y, z, vx, vy, vz, ax, ay, az, mass

⸻

Build Instructions

Compile using NVCC:

nvcc -std=c++17 -arch=sm_75 -O2 1000bodies_opt.cu -o 1000bodies_opt


Run Instructions

Place 1000bodies.csv in the same folder, then run:

./1000bodies_opt

The program prints the positions of all bodies at daily intervals.

