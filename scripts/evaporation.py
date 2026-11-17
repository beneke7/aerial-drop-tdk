"""Vectorized Denner (2026) isolated-water-drop heat and mass transfer.

This module implements the heat/mass-transfer model in Denner (2026),
Eqs. (12)--(33), in SI units.  It returns state derivatives; callers own the
time integrator and must stop a trajectory when its radius reaches zero or a
separate breakup/impact event occurs.  The model uses Denner's constant liquid
water density and heat capacity and assumes a spatially uniform drop
temperature, a spherical drop, uniform ambient conditions, and no radiation.

The functions accept scalar or broadcast-compatible NumPy inputs. Ambient air
properties are evaluated at the ambient temperature and relative humidity,
while water surface tension and latent heat use the instantaneous drop
temperature.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]

R_UNIVERSAL_J_MOL_K = 8.3145
M_WATER_KG_MOL = 18.0153e-3
RHO_WATER_KG_M3 = 997.0
CP_WATER_J_KG_K = 4181.0
T_CRITICAL_WATER_K = 647.096
P_AMBIENT_PA = 101325.0


@dataclass(frozen=True)
class AirProperties:
    """Denner Eq. (27) air properties at the supplied ambient state."""

    density_kg_m3: FloatArray
    dynamic_viscosity_pa_s: FloatArray
    thermal_conductivity_w_m_k: FloatArray
    specific_heat_j_kg_k: FloatArray


@dataclass(frozen=True)
class TransferRates:
    """Instantaneous drop transfer rates and diagnostic dimensionless groups.

    All fields are NumPy arrays, including zero-dimensional arrays for scalar
    inputs. Positive ``radius_rate_m_s`` and ``mass_rate_kg_s`` mean
    condensation; negative values mean evaporation. Positive
    ``convective_heat_out_w`` means heat leaves the drop.
    """

    radius_rate_m_s: FloatArray
    temperature_rate_k_s: FloatArray
    mass_rate_kg_s: FloatArray
    convective_heat_out_w: FloatArray
    latent_mass_term_w: FloatArray
    sherwood_number: FloatArray
    nusselt_number: FloatArray
    reynolds_number: FloatArray
    schmidt_number: FloatArray
    prandtl_number: FloatArray
    wet_bulb_temperature_k: FloatArray


def _array(name: str, value: ArrayLike) -> FloatArray:
    result = np.asarray(value, dtype=np.float64)
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    return result


def _broadcast_state(
    radius_m: ArrayLike,
    drop_temperature_k: ArrayLike,
    relative_speed_m_s: ArrayLike,
    ambient_temperature_k: ArrayLike,
    relative_humidity: ArrayLike,
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray, FloatArray]:
    names = (
        "radius_m",
        "drop_temperature_k",
        "relative_speed_m_s",
        "ambient_temperature_k",
        "relative_humidity",
    )
    values = np.broadcast_arrays(
        *(
            _array(name, value)
            for name, value in zip(
                names,
                (
                    radius_m,
                    drop_temperature_k,
                    relative_speed_m_s,
                    ambient_temperature_k,
                    relative_humidity,
                ),
                strict=True,
            )
        )
    )
    radius, drop_t, speed, ambient_t, rh = values
    if np.any(radius <= 0.0):
        raise ValueError("radius_m must be positive; terminate empty drops outside this RHS")
    if np.any(speed < 0.0):
        raise ValueError("relative_speed_m_s must be non-negative")
    if np.any((rh < 0.0) | (rh > 1.0)):
        raise ValueError("relative_humidity must be a fraction in [0, 1]")
    if np.any((ambient_t < 273.15) | (ambient_t > 323.15)):
        raise ValueError(
            "ambient_temperature_k must be within the liquid-water model range 0 to 50 degC"
        )
    if np.any((drop_t < 273.15) | (drop_t > 373.15)):
        raise ValueError(
            "drop_temperature_k must be within the model's 0 to 100 degC property range"
        )
    # Stull (2011), which Denner uses for Eq. (26), reports validity at
    # 5--99% RH. Exactly saturated conditions are an endpoint-only physical
    # case; the wet-bulb Nu correction itself is undefined there.
    if np.any((rh < 0.05) | ((rh > 0.99) & (rh < 1.0))):
        raise ValueError(
            "Denner's Stull wet-bulb closure is supported here for RH 0.05--0.99; RH=1 is endpoint-only"
        )
    return tuple(values)  # type: ignore[return-value]


def air_properties(ambient_temperature_k: ArrayLike, relative_humidity: ArrayLike) -> AirProperties:
    """Return air density, viscosity, conductivity and heat capacity (Eq. 27).

    The polynomial coefficients are Denner (2026), Table 1. The independent
    variables are dry-bulb temperature in degrees Celsius and RH as a fraction.
    The fit is based on Tsilingiris (2018), whose comparative review spans
    0--100 degC; Denner uses it for ambient air in the firefighting cases.
    """
    temperature = _array("ambient_temperature_k", ambient_temperature_k)
    rh = _array("relative_humidity", relative_humidity)
    temperature, rh = np.broadcast_arrays(temperature, rh)
    if np.any((temperature < 273.15) | (temperature > 323.15)):
        raise ValueError(
            "ambient_temperature_k must be within Denner's moist-air fit range 0 to 50 degC"
        )
    if np.any((rh < 0.0) | (rh > 1.0)):
        raise ValueError("relative_humidity must be a fraction in [0, 1]")
    theta = temperature - 273.15

    def polynomial(a1: float, b1: float, b2: float, c1: float, c2: float) -> FloatArray:
        return a1 + theta * (b1 + b2 * rh) + theta**2 * (c1 + c2 * rh)

    return AirProperties(
        density_kg_m3=polynomial(1.2932, -4.6688e-3, -8.6965e-4, 1.6165e-5, 2.2437e-5),
        dynamic_viscosity_pa_s=polynomial(
            1.7192e-5, 4.9467e-8, -2.2425e-9, -4.2412e-11, -3.2389e-10
        ),
        thermal_conductivity_w_m_k=polynomial(
            2.4031e-2, 7.5904e-5, -3.1199e-6, -4.8653e-8, -1.3015e-7
        ),
        specific_heat_j_kg_k=polynomial(1006.4, 3.9654e-2, 2.0101, 4.15208e-4, -1.673e-1),
    )


def wet_bulb_temperature_k(
    ambient_temperature_k: ArrayLike, relative_humidity: ArrayLike
) -> FloatArray:
    """Stull (2011) wet-bulb correlation as used by Denner Eq. (26).

    Stull's fit is stated valid for -20 to 50 degC and RH 5--99%, except for
    combined low-humidity/cold cases. For exactly saturated air this returns
    dry-bulb temperature as a physical endpoint, not as Stull's fit value.
    """
    temperature = _array("ambient_temperature_k", ambient_temperature_k)
    rh = _array("relative_humidity", relative_humidity)
    temperature, rh = np.broadcast_arrays(temperature, rh)
    if np.any((temperature < 273.15) | (temperature > 323.15)):
        raise ValueError(
            "ambient_temperature_k must be within the liquid-water model range 0 to 50 degC"
        )
    if np.any((rh < 0.05) | ((rh > 0.99) & (rh < 1.0)) | (rh > 1.0)):
        raise ValueError(
            "Stull wet-bulb correlation requires RH 0.05--0.99, with RH=1 supported only as endpoint"
        )

    theta_c = temperature - 273.15
    rh_percent = 100.0 * rh
    calculated_c = (
        theta_c * np.arctan(0.151977 * np.sqrt(rh_percent + 8.313659))
        + np.arctan(theta_c + rh_percent)
        - np.arctan(rh_percent - 1.676331)
        + 0.00391838 * rh_percent**1.5 * np.arctan(0.023101 * rh_percent)
        - 4.686035
    )
    # Stull's empirical fit is outside its stated range at RH=1 and has a
    # small residual there. Assign the physical saturated endpoint explicitly.
    calculated_c = np.where(rh == 1.0, theta_c, calculated_c)
    return calculated_c + 273.15


def saturation_vapor_pressure_pa(temperature_k: ArrayLike) -> FloatArray:
    """Huang (2018) saturation pressure over liquid water, Denner Eq. (30).

    The supplied correlation is used for liquid water from 0--100 degC, within
    the range relevant to this model. It returns pascals.
    """
    temperature = _array("temperature_k", temperature_k)
    if np.any((temperature < 273.15) | (temperature > 373.15)):
        raise ValueError(
            "temperature_k must be within the liquid-water property range 273.15--373.15 K"
        )
    return np.exp(34.494 - 4924.99 / (temperature - 36.05)) * (temperature - 168.15) ** -1.57


def water_surface_tension_n_m(drop_temperature_k: ArrayLike) -> FloatArray:
    """Water-air surface tension from Denner Eq. (28), based on Vargaftik et al."""
    temperature = _array("drop_temperature_k", drop_temperature_k)
    if np.any((temperature < 273.15) | (temperature > 373.15)):
        raise ValueError("drop_temperature_k must be within 273.15--373.15 K")
    reduced = (T_CRITICAL_WATER_K - temperature) / T_CRITICAL_WATER_K
    return 0.2358 * reduced**1.256 * (1.0 - 0.625 * reduced)


def latent_heat_vaporization_j_kg(drop_temperature_k: ArrayLike) -> FloatArray:
    """Water latent heat from Denner Eq. (29), based on Xiang (1997)."""
    temperature = _array("drop_temperature_k", drop_temperature_k)
    if np.any((temperature < 273.15) | (temperature > 373.15)):
        raise ValueError("drop_temperature_k must be within 273.15--373.15 K")
    reduced_critical = (T_CRITICAL_WATER_K - temperature) / T_CRITICAL_WATER_K
    reduced_drop = (T_CRITICAL_WATER_K - temperature) / temperature
    ratio = (
        6.8477
        * (reduced_drop**0.325 + 1.7383 * reduced_critical**0.835)
        / (1.0 + 1.3913 * reduced_critical)
    )
    return ratio * R_UNIVERSAL_J_MOL_K * T_CRITICAL_WATER_K / M_WATER_KG_MOL


def water_vapor_diffusivity_m2_s(ambient_temperature_k: ArrayLike) -> FloatArray:
    """Water-vapor diffusivity in air from Denner Eq. (31), at 101325 Pa."""
    temperature = _array("ambient_temperature_k", ambient_temperature_k)
    if np.any((temperature < 273.15) | (temperature > 323.15)):
        raise ValueError("ambient_temperature_k must be within Denner's 0 to 50 degC model range")
    return 2.178e-5 * (temperature / 273.15) ** 1.81


def _sherwood_or_nusselt(reynolds: FloatArray, diffusivity_number: FloatArray) -> FloatArray:
    """Beard/Pruppacher-Rasmussen piecewise sphere correlation (Eq. 21/23)."""
    transfer_parameter = np.sqrt(reynolds) * diffusivity_number ** (1.0 / 3.0)
    return np.where(
        transfer_parameter < 1.4,
        2.0 + 0.316 * reynolds * diffusivity_number ** (2.0 / 3.0),
        1.56 + 0.616 * transfer_parameter,
    )


def transfer_rates(
    radius_m: ArrayLike,
    drop_temperature_k: ArrayLike,
    relative_speed_m_s: ArrayLike,
    ambient_temperature_k: ArrayLike,
    relative_humidity: ArrayLike,
    *,
    apply_wet_bulb_correction: bool = True,
) -> TransferRates:
    """Evaluate Denner heat and mass transfer rates for isolated drops.

    Parameters use SI units except ``relative_humidity``, a fraction in
    [0, 1]. Inputs are broadcast together. Radius must be positive; terminate
    a drop at zero radius in the caller. Ambient temperatures are limited to
    0--50 degC, computed wet-bulb temperature must be at least 0 degC, and
    drop temperatures are limited to 0--100 degC. The corrected closure is
    supported at RH 5--99%. At exactly RH=100%, ``apply_wet_bulb_correction``
    must be false: Eq. (25) is 0/0 and Stull's fit is only validated through
    99%. This explicitly evaluates only the base Nu0 correlation.

    Returned derivatives are ``dr/dt`` [m/s] and ``dT/dt`` [K/s]. ``mass_rate``
    follows ``rho_w 4 pi r^2 dr/dt`` [kg/s]. The two signed energy rates are
    the outward convective heat flux and the latent mass term used in Denner
    Eq. (17), both [W].
    """
    radius, drop_t, speed, ambient_t, rh = _broadcast_state(
        radius_m,
        drop_temperature_k,
        relative_speed_m_s,
        ambient_temperature_k,
        relative_humidity,
    )
    air = air_properties(ambient_t, rh)
    twb = wet_bulb_temperature_k(ambient_t, rh)
    if np.any(twb < 273.15):
        raise ValueError("wet-bulb temperature below 0 degC is outside this liquid-water model")
    if apply_wet_bulb_correction and np.any(rh == 1.0):
        raise ValueError(
            "Eq. (25) is singular at RH=1; use RH<=0.99 or disable the correction explicitly"
        )
    psat_air = saturation_vapor_pressure_pa(ambient_t)
    psat_surface = saturation_vapor_pressure_pa(drop_t)
    diffusivity = water_vapor_diffusivity_m2_s(ambient_t)
    reynolds = 2.0 * radius * air.density_kg_m3 * speed / air.dynamic_viscosity_pa_s
    schmidt = air.dynamic_viscosity_pa_s / (air.density_kg_m3 * diffusivity)
    prandtl = air.dynamic_viscosity_pa_s * air.specific_heat_j_kg_k / air.thermal_conductivity_w_m_k
    sherwood = _sherwood_or_nusselt(reynolds, schmidt)
    nusselt_0 = _sherwood_or_nusselt(reynolds, prandtl)

    # Eqs. (16) and (25): evaluate the mass-flux radius rate at Td = Twb for
    # the wet-bulb consistency correction to Nu. Eq. (25) is not evaluated at
    # RH=1; the corrected path was rejected above because its limit is not
    # defined by Stull's fit outside 99% RH.
    twb_radius_rate = (
        sherwood
        * diffusivity
        * M_WATER_KG_MOL
        / (2.0 * radius * RHO_WATER_KG_M3 * R_UNIVERSAL_J_MOL_K)
        * (rh * psat_air / ambient_t - saturation_vapor_pressure_pa(twb) / twb)
    )
    denominator = air.thermal_conductivity_w_m_k * (twb - ambient_t)
    delta_nu = np.zeros_like(nusselt_0)
    if apply_wet_bulb_correction:
        delta_nu = (
            2.0
            * radius
            * RHO_WATER_KG_M3
            * latent_heat_vaporization_j_kg(twb)
            * twb_radius_rate
            / denominator
            - nusselt_0
        )
    nusselt = nusselt_0 + np.maximum(0.0, delta_nu)

    radius_rate = (
        sherwood
        * diffusivity
        * M_WATER_KG_MOL
        / (2.0 * radius * RHO_WATER_KG_M3 * R_UNIVERSAL_J_MOL_K)
        * (rh * psat_air / ambient_t - psat_surface / drop_t)
    )
    latent = latent_heat_vaporization_j_kg(drop_t)
    temperature_rate = (
        3.0
        / (radius * CP_WATER_J_KG_K)
        * (
            nusselt
            * air.thermal_conductivity_w_m_k
            * (ambient_t - drop_t)
            / (2.0 * radius * RHO_WATER_KG_M3)
            + latent * radius_rate
        )
    )
    mass_rate = 4.0 * np.pi * radius**2 * RHO_WATER_KG_M3 * radius_rate
    heat_transfer_coefficient = nusselt * air.thermal_conductivity_w_m_k / (2.0 * radius)
    convective_heat_out = heat_transfer_coefficient * 4.0 * np.pi * radius**2 * (drop_t - ambient_t)
    latent_mass_term = latent * mass_rate
    return TransferRates(
        radius_rate_m_s=radius_rate,
        temperature_rate_k_s=temperature_rate,
        mass_rate_kg_s=mass_rate,
        convective_heat_out_w=convective_heat_out,
        latent_mass_term_w=latent_mass_term,
        sherwood_number=sherwood,
        nusselt_number=nusselt,
        reynolds_number=reynolds,
        schmidt_number=schmidt,
        prandtl_number=prandtl,
        wet_bulb_temperature_k=twb,
    )


def characteristic_evaporation_time_s(
    initial_radius_m: ArrayLike,
    ambient_temperature_k: ArrayLike,
    relative_humidity: ArrayLike,
    sherwood_number: ArrayLike = 2.0,
) -> FloatArray:
    """Stationary-drop d-squared evaporation time, Denner Eq. (33).

    Infinite time is returned when the Eq. (33) driving term is non-positive
    (saturated or condensation conditions). This is a reference scale, not a
    finite-flight lifetime prediction.
    """
    radius, ambient_t, rh, sherwood = np.broadcast_arrays(
        _array("initial_radius_m", initial_radius_m),
        _array("ambient_temperature_k", ambient_temperature_k),
        _array("relative_humidity", relative_humidity),
        _array("sherwood_number", sherwood_number),
    )
    if np.any(radius <= 0.0) or np.any(sherwood <= 0.0):
        raise ValueError("initial_radius_m and sherwood_number must be positive")
    twb = wet_bulb_temperature_k(ambient_t, rh)
    if np.any(twb < 273.15):
        raise ValueError("wet-bulb temperature below 0 degC is outside this liquid-water model")
    drive = (
        saturation_vapor_pressure_pa(twb) / twb
        - rh * saturation_vapor_pressure_pa(ambient_t) / ambient_t
    )
    denominator = sherwood * water_vapor_diffusivity_m2_s(ambient_t) * M_WATER_KG_MOL * drive
    numerator = radius**2 * RHO_WATER_KG_M3 * R_UNIVERSAL_J_MOL_K
    return np.divide(
        numerator,
        denominator,
        out=np.full_like(numerator, np.inf),
        where=denominator > 0.0,
    )


def d2_law_surface_area_fraction(
    time_s: ArrayLike,
    initial_radius_m: ArrayLike,
    ambient_temperature_k: ArrayLike,
    relative_humidity: ArrayLike,
    sherwood_number: ArrayLike = 2.0,
) -> FloatArray:
    """Eq. (32) surface-area fraction for the wet-bulb terminal limit.

    Values are zero after the analytical lifetime. The equation is a
    quasi-steady reference (drop already at wet-bulb temperature), not the
    transient heat-capacity solution of Eqs. (16)/(20).
    """
    time, radius, ambient_t, rh, sherwood = np.broadcast_arrays(
        _array("time_s", time_s),
        _array("initial_radius_m", initial_radius_m),
        _array("ambient_temperature_k", ambient_temperature_k),
        _array("relative_humidity", relative_humidity),
        _array("sherwood_number", sherwood_number),
    )
    if np.any(time < 0.0) or np.any(radius <= 0.0) or np.any(sherwood <= 0.0):
        raise ValueError("time_s must be non-negative; radius and Sherwood number must be positive")
    twb = wet_bulb_temperature_k(ambient_t, rh)
    if np.any(twb < 273.15):
        raise ValueError("wet-bulb temperature below 0 degC is outside this liquid-water model")
    drive = (
        saturation_vapor_pressure_pa(twb) / twb
        - rh * saturation_vapor_pressure_pa(ambient_t) / ambient_t
    )
    slope = (
        sherwood
        * water_vapor_diffusivity_m2_s(ambient_t)
        * M_WATER_KG_MOL
        * drive
        / (radius**2 * RHO_WATER_KG_M3 * R_UNIVERSAL_J_MOL_K)
    )
    return np.maximum(0.0, 1.0 - time * slope)
