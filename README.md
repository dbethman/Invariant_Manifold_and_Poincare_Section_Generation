# Invariant Manifold and Poincaré Section Generation (WORK IN PROGRESS)

A Python-based Astrodynamics toolset for generating, propagating, and visualizing invariant manifolds and Poincaré sections for periodic orbits in the Circular Restricted Three-Body Problem (CR3BP). 

This project utilizes the [ASSET](https://github.com/astrodynamics/asset) (Astrodynamics Software and Science Enabling Toolkit) library for optimal control, numerical integration, and state transition matrix (STM) calculations.

##  Features

* **Periodic Orbit Targeting:** Uses multiple shooting and optimal control collocation (LGL3) to converge orbits in the Earth-Moon system.
* **Manifold Generation:** Extracts stable and unstable eigenvectors from the Monodromy matrix to perturb the initial state and propagate invariant manifolds forward and backward in time.
* **Lunar Impact Filtering:** Automatically detects and filters out manifold branches that impact the lunar surface (radius < 1700 km).
* **3D System Visualization:** Plots the converged orbit and manifolds alongside dimensionally-accurate Earth and Moon spheres, enforcing a true 1:1:1 aspect ratio.
* **Poincaré Sections:** Computes continuous-trajectory plane crossings (e.g., at the Earth's X-coordinate) using linear interpolation for exact crossing states.
* **3D Velocity Mapping:** Projects Poincaré crossings into 3D space $(Y, Z, V)$ to map speed characteristics at the crossing plane, complete with Geosynchronous Orbit (GEO) footprint references.

## 📦 Dependencies

This project requires Python 3.x and the following packages:
* `numpy`
* `matplotlib`
* `asset_asrl` (ASSET)

## 🛠️ Usage

### 1. Generate the Manifold Data
Run the main solver to converge the orbit, compute the manifolds, and save the trajectory data to a Pickle file.

```bash
python Manifold_Functions.py
