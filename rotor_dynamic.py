import numpy as np
import matplotlib.pyplot as plt
import aircraft_data as ad

# Constants
V = ad.HelicopterLynx().V_max    # Forward flight velocity (m/s)
Omega = 30          # Rotor angular velocity (rad/s)
R = ad.HelicopterLynx().rotor_d / 2               # Rotor radius (m)

# Blade pitch settings (convert to radians)
theta_0 = np.deg2rad(6)   # Collective pitch
theta_1c = np.deg2rad(2)  # Longitudinal cyclic
theta_1s = np.deg2rad(1)  # Lateral cyclic

# Blade azimuth positions (0 to 2*pi)
psi = np.linspace(0, 2*np.pi, 360)

# 1) Blade flapping angle beta(psi)
# Approximation: beta1c ≈ (1/2) * theta_1c, beta1s ≈ (1/2) * theta_1s
beta0 = 0
beta1c = 0.5 * theta_1c
beta1s = 0.5 * theta_1s
beta = beta0 + beta1c * np.cos(psi) + beta1s * np.sin(psi)

# 2) Angle of attack variation α(ψ)
# Grid resolution
n_radial = 50
n_azimuthal = 360
r = np.linspace(0.05, R, n_radial)  # Avoid r=0 to prevent division by zero
psi = np.linspace(0, 2 * np.pi, n_azimuthal)
r_grid, psi_grid = np.meshgrid(r, psi)

# Blade pitch distribution θ(ψ)
theta_grid = theta_0 + theta_1c * np.cos(psi_grid) + theta_1s * np.sin(psi_grid)

# Inflow angle φ ≈ V*sin(ψ)/(Ω*r)
phi_grid = np.arctan2(V * np.sin(psi_grid), Omega * r_grid)

# Angle of attack α = θ - φ
alpha_grid = theta_grid - phi_grid
alpha_deg = np.rad2deg(alpha_grid)

# Convert to rotor disc X-Y coordinates
x = (r_grid / R) * np.cos(psi_grid)
y = (r_grid / R) * np.sin(psi_grid)

# 3) Fourier series components of alpha to get disc tilt
a0 = np.mean(alpha_grid)
a1c = 2 * np.mean(alpha_grid * np.cos(psi))
a1s = 2 * np.mean(alpha_grid * np.sin(psi))

# Convert to degrees for interpretation
disc_tilt_longitudinal_deg = np.rad2deg(a1c)
disc_tilt_lateral_deg = np.rad2deg(a1s)

# Output tilt angles
print(f"Disc Tilt (Longitudinal): {disc_tilt_longitudinal_deg:.2f} degrees")
print(f"Disc Tilt (Lateral): {disc_tilt_lateral_deg:.2f} degrees")

# Blade Flapping
plt.figure(figsize=(7, 7))
plt.plot(np.rad2deg(psi), np.rad2deg(beta), color='orange')
plt.title('Blade Flapping Angle vs Azimuth')
plt.xlabel('Azimuth Angle (degrees)')
plt.ylabel('Flapping Angle β (degrees)')
plt.grid(True)
plt.show()

# Plotting the contours of constant angle of attack
plt.figure(figsize=(7, 7))
contour = plt.contour(x, y, alpha_deg, levels=np.arange(0, 18, 2), cmap='jet')
plt.clabel(contour, inline=True, fontsize=10)
plt.title("Curves of constant angle of attack", fontsize=14, style='italic')
plt.xlabel("Blade position (nondimensional)")
plt.ylabel("Blade position (nondimensional)")
plt.axvline(0, color='k', linestyle='--')
plt.text(-0.6, 0, "Retreating side", fontsize=10, ha='center')
plt.text(0.6, 0, "Advancing side", fontsize=10, ha='center')
plt.axis("equal")
plt.grid(True)
plt.show()
