import numpy as np
import matplotlib.pyplot as plt

# === Rotor and Flight Parameters ===
V = 20                  # Forward velocity (m/s)
Omega = 30              # Rotor angular speed (rad/s)
R = 5                   # Rotor radius (m)
n_radial = 50           # Radial steps

# Blade pitch settings (in radians)
theta_0 = np.radians(6)      # Collective
theta_1c = np.radians(2)     # Longitudinal cyclic
theta_1s = np.radians(1)     # Lateral cyclic

# === Blade azimuth and radius grid ===
psi = np.linspace(0, 2*np.pi, 360)  # Azimuth angle
r = np.linspace(0.05, R, n_radial)  # Avoid r = 0 to prevent singularities
r_grid, psi_grid = np.meshgrid(r, psi)

# === Compute Pitch Distribution θ(ψ) ===
theta_grid = theta_0 + theta_1c * np.cos(psi_grid) + theta_1s * np.sin(psi_grid)

# === Compute inflow angle φ(ψ, r) ===
phi_grid = np.arctan2(V * np.sin(psi_grid), Omega * r_grid)

# === Compute Angle of Attack α = θ - φ ===
alpha_grid = theta_grid - phi_grid
alpha_deg = np.rad2deg(alpha_grid)

# === Transform to nondimensional x-y (rotor disk coordinates) ===
x = (r_grid / R) * np.cos(psi_grid)
y = (r_grid / R) * np.sin(psi_grid)

# === Plot ===
fig, axs = plt.subplots(1, 2, figsize=(14, 6))

# -- Blade flapping angle --
beta = np.rad2deg(theta_0 + theta_1c * np.cos(psi) + theta_1s * np.sin(psi))
axs[0].plot(np.degrees(psi), beta, label="Blade Flapping Angle", color="tab:orange")
axs[0].set_title("Blade Flapping Angle over Revolution", fontsize=14)
axs[0].set_xlabel("Azimuth Angle ψ (deg)")
axs[0].set_ylabel("Flapping Angle β (deg)")
axs[0].grid(True)
axs[0].legend()

# -- AoA contour plot --
levels = np.arange(0, 18, 2)
contour = axs[1].contour(x, y, alpha_deg, levels=levels, cmap='jet')
axs[1].clabel(contour, inline=True, fontsize=10)
axs[1].set_title("Curves of Constant Angle of Attack", fontsize=14, style='italic')
axs[1].set_xlabel("Blade position (nondimensional)")
axs[1].set_ylabel("Blade position (nondimensional)")
axs[1].axvline(0, color='k', linestyle='--')
axs[1].text(-0.6, 0.8, "Retreating side", fontsize=10)
axs[1].text(0.6, 0.8, "Advancing side", fontsize=10)
axs[1].set_aspect('equal')
axs[1].grid(False)

plt.tight_layout()
plt.show()
