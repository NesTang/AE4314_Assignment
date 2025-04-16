import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import aircraft_data as acd

# === Rotor and Flight Parameters ===
V = 20                  # Forward flight speed (m/s)
Omega = 30              # Rotor angular speed (rad/s)
R = 5                   # Rotor radius (m)
n_radial = 50           # Radial steps
n_azimuth = 360         # Azimuth resolution

# === Blade pitch input (in radians) ===
theta_0 = np.radians(6)      # Collective pitch
theta_1c = np.radians(2)     # Longitudinal cyclic
theta_1s = np.radians(1)     # Lateral cyclic

# === Body rotational rates ===
q = np.radians(20)  # Body pitch rate (rad/s)
p = np.radians(10)  # Body roll rate (rad/s)

# === Effective cyclic pitch due to body motion ===
theta_1c_eff = theta_1c + (q * R) / Omega
theta_1s_eff = theta_1s - (p * R) / Omega

# === Blade azimuth and radial grid ===
psi = np.linspace(0, 2*np.pi, n_azimuth)  # Azimuth angle
r = np.linspace(0.05, R, n_radial)        # Radial positions (avoid 0)

# Create meshgrid for calculation
r_grid, psi_grid = np.meshgrid(r, psi)

# === Compute local pitch distribution θ(ψ) with effective cyclic ===
theta_grid = theta_0 + theta_1c_eff * np.cos(psi_grid) + theta_1s_eff * np.sin(psi_grid)

# === Compute inflow angle φ(ψ, r) ===
phi_grid = np.arctan2(V * np.sin(psi_grid), Omega * r_grid)

# === Compute angle of attack α = θ - φ ===
alpha_grid = theta_grid - phi_grid
alpha_deg = np.rad2deg(alpha_grid)

# === Transform to nondimensional x-y plane (rotor disk) ===
x = (r_grid / R) * np.cos(psi_grid)
y = (r_grid / R) * np.sin(psi_grid)

# === Compute flapping angle beta(ψ) ===
# (for now, approximated as pitch input)
beta = np.rad2deg(theta_0 + theta_1c_eff * np.cos(psi) + theta_1s_eff * np.sin(psi))

# === --- Part 3: Fourier Series Coefficients --- ===
beta_rad = np.radians(beta)  # Convert to radians

# Integrate using Simpson’s rule for better accuracy
beta_1c = (2 / np.pi) * simpson(beta_rad * np.cos(psi), x=psi)
beta_1s = (2 / np.pi) * simpson(beta_rad * np.sin(psi), x=psi)

# Disc tilt angles in radians (approx.)
tilt_longitudinal = beta_1c
tilt_lateral = beta_1s

print("=== Disc Tilt Angles ===")
print(f"Longitudinal Tilt (β₁c): {np.degrees(tilt_longitudinal):.3f} deg")
print(f"Lateral Tilt (β₁s): {np.degrees(tilt_lateral):.3f} deg")

# === Plotting ===
fig, axs = plt.subplots(1, 2, figsize=(14, 6))

# -- Blade flapping angle vs azimuth --
beta = np.rad2deg(theta_0 + theta_1c_eff * np.cos(psi) + theta_1s_eff * np.sin(psi))
axs[0].plot(np.degrees(psi), beta, color="tab:orange", label="Flapping Angle")
axs[0].set_title("Blade Flapping Angle over One Revolution", fontsize=14)
axs[0].set_xlabel("Azimuth Angle ψ (degrees)")
axs[0].set_ylabel("Flapping Angle β (degrees)")
axs[0].grid(True)
axs[0].legend()

# -- Angle of Attack contour plot --
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
