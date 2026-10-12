#!/usr/bin/env python3
"""Reproduce elementary checks for manuscript.tex.

Requires Python 3 and NumPy. These are consistency checks of the stated
classical models, not simulations of general relativity or physical tests
of the continuum conjecture. No files are changed by this script.
"""

import math

import numpy as np


def report(label, value, tolerance=None):
    print(f"{label}: {value:.12g}")
    if tolerance is not None:
        assert value < tolerance, (label, value, tolerance)


print("Analytic-model consistency checks; deterministic seed 20261012")
rng = np.random.default_rng(20261012)

# Retarded geometry: solve for the emission distance, then evaluate the
# direction and Jacobian independently of the simplified D expression.
beta = rng.uniform(0.0, 0.995, 2000)
s = rng.normal(size=(2000, 3))
xi = s[:, 0]
b2 = s[:, 1] ** 2 + s[:, 2] ** 2
D = np.sqrt(xi**2 + (1 - beta**2) * b2)
R = (beta * xi + D) / (1 - beta**2)
retarded_vector = s.copy()
retarded_vector[:, 0] += beta * R
n = retarded_vector / np.linalg.norm(retarded_vector, axis=1)[:, None]
kappa = 1 - beta * n[:, 0]
report("retarded distance relative residual",
       np.max(np.abs(np.linalg.norm(retarded_vector, axis=1) / R - 1)), 2e-12)
report("kappa R = D relative residual",
       np.max(np.abs(kappa * R / D - 1)), 2e-12)

# Numerical derivative of momentum (five-point stencil), in c=m=1 units.
v = np.linspace(-0.95, 0.95, 101)
step = 2e-5
p = lambda w: w / np.sqrt(1 - w**2)
dp = (p(v - 2 * step) - 8 * p(v - step)
      + 8 * p(v + step) - p(v + 2 * step)) / (12 * step)
report("dp/dv vs gamma^3 maximum relative residual",
       np.max(np.abs(dp * (1 - v**2) ** 1.5 - 1)), 2e-9)

# Check the acoustic metric density, including the conformal prefactor.
metric_error = 0.0
for _ in range(40):
    density = rng.uniform(0.5, 3.0)
    cs = rng.uniform(0.5, 3.0)
    vel = rng.normal(size=3) * cs * 0.3
    metric = np.eye(4)
    metric[0, 0] = -(cs**2 - vel @ vel)
    metric[0, 1:] = metric[1:, 0] = -vel
    metric *= density / cs
    density_from_metric = np.sqrt(-np.linalg.det(metric)) * np.linalg.inv(metric)
    expected = np.empty((4, 4))
    expected[0, 0] = -1
    expected[0, 1:] = expected[1:, 0] = -vel
    expected[1:, 1:] = cs**2 * np.eye(3) - np.outer(vel, vel)
    expected *= density / cs**2
    metric_error = max(metric_error, np.max(np.abs(density_from_metric - expected)))
report("acoustic sqrt(-g) g^-1 maximum absolute residual", metric_error, 2e-13)

# Independently integrate the TT angular projector to check the 8*pi/5
# coefficient linking far-zone flux to total quadrupole power.
random_matrix = rng.normal(size=(3, 3))
B = (random_matrix + random_matrix.T) / 2
B -= np.eye(3) * np.trace(B) / 3
cos_theta, weights = np.polynomial.legendre.leggauss(32)
azimuth = np.arange(64) * 2 * math.pi / 64
angular_integral = 0.0
for zz, weight in zip(cos_theta, weights):
    for angle in azimuth:
        normal = np.array([math.sqrt(1 - zz**2) * math.cos(angle),
                           math.sqrt(1 - zz**2) * math.sin(angle), zz])
        transverse = np.eye(3) - np.outer(normal, normal)
        projected = transverse @ B @ transverse
        projected -= transverse * np.trace(transverse @ B) / 2
        angular_integral += weight * np.sum(projected**2) * 2 * math.pi / 64
expected_integral = 8 * math.pi / 5 * np.sum(B**2)
report("TT angular integral relative residual",
       abs(angular_integral / expected_integral - 1), 2e-12)

# Evaluate the moving kink PDE by finite differences of the field itself.
K, ell, speed = 1.7, 0.8, 2.3
def kink(t, x, velocity):
    gamma = 1 / math.sqrt(1 - (velocity / speed) ** 2)
    return 4 * np.arctan(np.exp(gamma * (x - velocity * t) / ell))

def second_derivative(fun, point, h):
    return (-fun(point + 2 * h) + 16 * fun(point + h) - 30 * fun(point)
            + 16 * fun(point - h) - fun(point - 2 * h)) / (12 * h**2)

pde_error = 0.0
energy_error = 0.0
momentum_error = 0.0
z = np.linspace(-25, 25, 60001)
integrate = getattr(np, "trapezoid", None)
if integrate is None:
    integrate = np.trapz
for fraction in (0.0, 0.2, 0.8, 0.95, 0.99):
    velocity = fraction * speed
    gamma = 1 / math.sqrt(1 - fraction**2)
    t = 0.13
    xx = velocity * t + np.linspace(-4, 4, 71) * ell / gamma
    # Steps resolve the Lorentz-contracted profile equally at each speed.
    dx = ell / gamma * 2e-3
    dt = dx / speed
    phi_tt = second_derivative(lambda tt: kink(tt, xx, velocity), t, dt)
    phi_xx = second_derivative(lambda xx_: kink(t, xx_, velocity), xx, dx)
    residual = phi_tt / speed**2 - phi_xx + np.sin(kink(t, xx, velocity)) / ell**2
    pde_error = max(pde_error, np.max(np.abs(residual)) * ell**2 / gamma**2)

    # Integrate the unsimplified energy and momentum densities.
    phi = 4 * np.arctan(np.exp(z))
    phi_x = 2 * gamma / ell / np.cosh(z)
    phi_t = -velocity * phi_x
    energy_density = K * (phi_t**2 / speed**2 + phi_x**2) / 2
    energy_density += K / ell**2 * (1 - np.cos(phi))
    momentum_density = -K / speed**2 * phi_t * phi_x
    x_grid = z * ell / gamma
    E_numeric = integrate(energy_density, x_grid)
    P_numeric = integrate(momentum_density, x_grid)
    E0 = 8 * K / ell
    energy_error = max(energy_error, abs(E_numeric / (gamma * E0) - 1))
    if velocity:
        momentum_error = max(momentum_error,
                             abs(P_numeric / (gamma * E0 * velocity / speed**2) - 1))
    else:
        assert P_numeric == 0
report("dimensionless kink PDE finite-difference residual", pde_error, 3e-8)
report("kink integrated energy maximum relative residual", energy_error, 2e-11)
report("kink integrated momentum maximum relative residual", momentum_error, 2e-11)

# Differentiate the periodic quadrupole numerically using a spectral method,
# rather than inserting its closed-form third derivative.
G = 6.67430e-11
c = 299792458.0
mu, length, amplitude, omega = 1.0, 0.1, 0.001, 2 * math.pi
N = 256
period = 2 * math.pi / omega
times = np.arange(N) * period / N
separation = length + amplitude * np.cos(omega * times)
Q = mu * separation[:, None] ** 2 * np.array([2 / 3, -1 / 3, -1 / 3])[None, :]
frequencies = 2 * math.pi * np.fft.fftfreq(N, d=period / N)
Q3 = np.fft.ifft((1j * frequencies[:, None])**3 * np.fft.fft(Q, axis=0), axis=0).real
power_numeric = G / (5 * c**5) * np.mean(np.sum(Q3**2, axis=1))
power_formula = 4 * G * mu**2 * omega**6 / (15 * c**5)
power_formula *= length**2 * amplitude**2 + 4 * amplitude**4
report("quadrupole power numerical vs analytic relative residual",
       abs(power_numeric / power_formula - 1), 2e-10)
report("oscillating-pair power (W)", power_formula)
report("one-kilogram, 0.1-m-radius compactness", G / (0.1 * c**2))
report("uniform-sphere binding mass fraction", -3 * G / (5 * 0.1 * c**2))

print("All checks passed. These results do not test the vacuum-continuum conjecture.")
