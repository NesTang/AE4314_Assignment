import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
import aircraft_data as acd

# === Rotor and Flight Parameters ===
V = 20                  # Forward flight speed (m/s)
V_tip = 700          # Tip speed (ft/s) given in NASA report
V_tip = V_tip * 0.3048  # Convert ft/s to m/s
R = acd.HelicopterLynx().rotor_d / 2     # Rotor radius (m)
Omega = V_tip / R              # Rotor angular speed (rad/s)
n_radial = 50
n_azimuth = 360

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

# === Azimuth and radial grids ===
psi = np.linspace(0, 2*np.pi, n_azimuth)
r = np.linspace(0.05, R, n_radial)
r_grid, psi_grid = np.meshgrid(r, psi)

# === Compute pitch angle at each location on the disc ===
theta_grid = theta_0 + theta_1c_eff * np.cos(psi_grid) + theta_1s_eff * np.sin(psi_grid)

# === Compute inflow angle φ ===
phi_grid = np.arctan2(V * np.sin(psi_grid), Omega * r_grid)

# === Compute angle of attack α = θ - φ ===
alpha_grid = theta_grid - phi_grid
alpha_deg = np.rad2deg(alpha_grid)

# === Compute blade flapping angle beta (as function of psi only) ===
beta_psi = theta_0 + theta_1c_eff * np.cos(psi) + theta_1s_eff * np.sin(psi)
beta_deg_psi = np.rad2deg(beta_psi)

# === Fourier series (Part 3) ===
beta_rad = beta_psi  # already in radians
beta_1c = (2 / np.pi) * simpson(beta_rad * np.cos(psi), x=psi)
beta_1s = (2 / np.pi) * simpson(beta_rad * np.sin(psi), x=psi)

# Disc tilt angles
tilt_long = np.degrees(beta_1c)
tilt_lat = np.degrees(beta_1s)

print("=== Disc Tilt Angles ===")
print(f"Longitudinal Tilt (β₁c): {tilt_long:.3f} deg")
print(f"Lateral Tilt (β₁s): {tilt_lat:.3f} deg")

# === Rotor disk X-Y for plotting ===
x = (r_grid / R) * np.cos(psi_grid - np.pi/2)
y = (r_grid / R) * np.sin(psi_grid - np.pi/2)

# === Flapping angle 2D grid for contour plotting ===
psi2D, r2D = np.meshgrid(psi, np.linspace(0.05, 1, n_radial))
beta_grid = theta_0 + theta_1c_eff * np.cos(psi2D) + theta_1s_eff * np.sin(psi2D)
beta_deg_grid = np.rad2deg(beta_grid)
x_beta = r2D * np.cos(psi2D - np.pi/2)
y_beta = r2D * np.sin(psi2D - np.pi/2)

# === Plotting both as disk contour plots ===
fig, axs = plt.subplots(1, 2, figsize=(14, 6))

# -- Part 1: Flapping Angle Contour --
cont1 = axs[0].contour(x_beta, y_beta, beta_deg_grid, levels=12, cmap='viridis')
axs[0].clabel(cont1, inline=True, fontsize=9)
axs[0].set_title("Blade Flapping Angle β (deg)", fontsize=14)
axs[0].set_xlabel("Blade position (nondimensional)")
axs[0].set_ylabel("Blade position (nondimensional)")
axs[0].text(-0.6, 0.8, "Retreating side", fontsize=10)
axs[0].text(0.6, 0.8, "Advancing side", fontsize=10)
axs[0].set_aspect('equal')
axs[0].grid(False)

# -- Part 2: Angle of Attack Contour --
cont2 = axs[1].contour(x, y, alpha_deg, levels=np.arange(0, 18, 2), cmap='jet')
axs[1].clabel(cont2, inline=True, fontsize=9)
axs[1].set_title("Angle of Attack α (deg)", fontsize=14)
axs[1].set_xlabel("Blade position (nondimensional)")
axs[1].set_ylabel("Blade position (nondimensional)")
axs[1].text(-0.6, 0.8, "Retreating side", fontsize=10)
axs[1].text(0.6, 0.8, "Advancing side", fontsize=10)
axs[1].set_aspect('equal')
axs[1].grid(False)

plt.tight_layout()
plt.show()
