import csv
import math
import os
import sys

import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import mission_params as mp
import guidance

G0 = mp.G0
R = mp.KERBIN_RADIUS
MU = mp.KERBIN_MU
OMEGA = 2 * math.pi / mp.KERBIN_ROTATION_PERIOD
ROTATION_IN_PLANE = math.cos(math.radians(mp.TARGET_INCLINATION)) / math.cos(math.radians(mp.KSC_LATITUDE))

_atm_table = None


def current_g_h(h, mu=MU, r=R):
    return mu / ((r + h) ** 2)


def first_cosmic_velocity(h=0.0, mu=MU, r=R):
    return math.sqrt(mu / (r + h))


def tsiolkovsky(isp, m0, m1):
    return isp * G0 * math.log(m0 / m1)


def drag_force(cda, ro, v):
    return cda * ro * v ** 2 / 2


def isp_at_pressure(isp_sl, isp_vac, p):
    return isp_vac + (isp_sl - isp_vac) * p / mp.KERBIN_P0


def use_ksp_atmosphere(csv_path):
    global _atm_table
    rows = read_ksp_flight(csv_path)
    pairs = sorted((r["alt"], r["pressure"], r["density"]) for r in rows
                   if r["alt"] < mp.KERBIN_ATM_HEIGHT and r["pressure"] > 0)
    if len(pairs) < 10:
        return False
    h, p, ro = (np.array(c) for c in zip(*pairs))
    _atm_table = (h, p, ro)
    return True


def find_p_atm(h):
    if h >= mp.KERBIN_ATM_HEIGHT:
        return 0.0
    if _atm_table is not None:
        return float(np.interp(h, _atm_table[0], _atm_table[1]))
    return mp.KERBIN_P0 * math.exp(-max(h, 0.0) / mp.KERBIN_SCALE_HEIGHT)


def find_ro(h):
    if h >= mp.KERBIN_ATM_HEIGHT:
        return 0.0
    if _atm_table is not None:
        return float(np.interp(h, _atm_table[0], _atm_table[2]))
    return mp.KERBIN_RHO0 * find_p_atm(h) / mp.KERBIN_P0


def orbit_from_state(r, vr, vt, mu=MU):
    energy = (vr ** 2 + vt ** 2) / 2 - mu / r
    h = r * vt
    e = math.sqrt(max(0.0, 1 + 2 * energy * h ** 2 / mu ** 2))
    if energy >= 0:
        return math.inf, h ** 2 / mu / (1 + e) - R, math.inf, e
    a = -mu / (2 * energy)
    return a * (1 + e) - R, a * (1 - e) - R, a, e


def orbital_period(a, mu=MU):
    return 2 * math.pi * math.sqrt(a ** 3 / mu)


def delta_v_budget(vostok):
    b, c, u = vostok.boosters, vostok.core, vostok.upper
    m0 = vostok.total_mass
    burned_1 = (b.mass_flow + c.mass_flow) * b.burn_time
    m1 = m0 - burned_1
    isp_1 = (b.thrust_vac + c.thrust_vac) / ((b.mass_flow + c.mass_flow) * G0)
    m2 = m1 - b.dry_mass
    core_left = c.fuel - c.mass_flow * b.burn_time
    m3 = m2 - core_left - vostok.fairing_mass
    m4 = m3 - c.dry_mass
    m5 = m4 - u.fuel
    return [
        ("Боковые блоки + блок А", tsiolkovsky(isp_1, m0, m1)),
        ("Блок А", tsiolkovsky(c.isp_vac, m2, m3)),
        ("Блок Е", tsiolkovsky(u.isp_vac, m4, m5)),
    ]


def derivatives(state, engines, direction, cda):
    r, phi, vr, vt, m = state
    h = r - R
    p = find_p_atm(h)
    ro = find_ro(h)

    thrust = 0.0
    mass_flow = 0.0
    for flow, isp_sl, isp_vac in engines:
        thrust += flow * G0 * isp_at_pressure(isp_sl, isp_vac, p)
        mass_flow += flow

    vt_air = vt - OMEGA * r * ROTATION_IN_PLANE
    v_air = math.hypot(vr, vt_air)
    drag = drag_force(cda, ro, v_air)
    dr, dt = direction
    ar = thrust / m * dr - mu_over_r2(r) + vt ** 2 / r
    at = thrust / m * dt - vr * vt / r
    if v_air > 0:
        ar -= drag / m * vr / v_air
        at -= drag / m * vt_air / v_air
    gx = thrust / m * dr - (drag / m * vr / v_air if v_air > 0 else 0)
    gy = thrust / m * dt - (drag / m * vt_air / v_air if v_air > 0 else 0)
    extra = {"thrust": thrust, "drag": drag, "g_load": math.hypot(gx, gy) / G0, "v_air": v_air,
             "q": ro * v_air ** 2 / 2}
    return (vr, vt / r, ar, at, -mass_flow), extra


def mu_over_r2(r):
    return MU / r ** 2


def rk4_step(state, dt, engines, direction, cda):
    k1, extra = derivatives(state, engines, direction, cda)
    s2 = [x + dt / 2 * k for x, k in zip(state, k1)]
    k2, _ = derivatives(s2, engines, direction, cda)
    s3 = [x + dt / 2 * k for x, k in zip(state, k2)]
    k3, _ = derivatives(s3, engines, direction, cda)
    s4 = [x + dt * k for x, k in zip(state, k3)]
    k4, _ = derivatives(s4, engines, direction, cda)
    new = [x + dt / 6 * (a + 2 * b + 2 * c + d) for x, a, b, c, d in zip(state, k1, k2, k3, k4)]
    return new, extra


def pitch_direction(pitch_deg):
    p = math.radians(pitch_deg)
    return math.sin(p), math.cos(p)


class Log:
    def __init__(self):
        self.rows = {k: [] for k in ("t", "alt", "speed", "speed_surface", "g_load", "mass",
                                     "pitch", "q", "apoapsis", "periapsis")}

    def add(self, t, state, extra, pitch):
        r, phi, vr, vt, m = state
        apo, peri, _, _ = orbit_from_state(r, vr, vt)
        for k, v in (("t", t), ("alt", r - R), ("speed", math.hypot(vr, vt)),
                     ("speed_surface", extra["v_air"]), ("g_load", extra["g_load"]),
                     ("mass", m), ("pitch", pitch), ("q", extra["q"]),
                     ("apoapsis", apo), ("periapsis", peri)):
            self.rows[k].append(v)

    def arrays(self):
        out = {k: np.array(v) for k, v in self.rows.items()}
        out["accel"] = np.gradient(out["speed"], out["t"]) if len(out["t"]) > 1 else out["speed"] * 0
        return out


def simulate_ascent(vostok, dt=0.02, log_every=0.5, t_max=1200.0):
    t = 0.0
    state = [R, 0.0, 0.0, OMEGA * R * ROTATION_IN_PLANE, vostok.total_mass]
    events = []
    log = Log()
    next_log = 0.0
    cutoff = False
    t_cutoff = None

    while t < t_max:
        b, c, u = vostok.boosters, vostok.core, vostok.upper
        r, phi, vr, vt, m = state

        if vostok.has_boosters and b.fuel <= 0:
            vostok.has_boosters = False
            events.append((t, "Отделение боковых блоков"))
        if vostok.has_fairing and t >= mp.T_FAIRING:
            vostok.has_fairing = False
            events.append((t, "Сброс обтекателя"))
        if vostok.has_core and (t >= mp.T_CORE_SEP or (not vostok.has_boosters and c.fuel <= 0)):
            vostok.has_core = False
            events.append((t, "Отделение блока А, запуск блока Е"))
        state[4] = vostok.total_mass

        engines = []
        active = []
        if vostok.has_boosters and b.fuel > 0:
            engines.append((b.mass_flow, b.isp_sl, b.isp_vac))
            active.append(b)
        if vostok.has_core and c.fuel > 0:
            engines.append((c.mass_flow, c.isp_sl, c.isp_vac))
            active.append(c)
        if not vostok.has_core and not cutoff and u.fuel > 0:
            engines.append((u.mass_flow, u.isp_sl, u.isp_vac))
            active.append(u)

        if vostok.has_core:
            pitch = guidance.ascent_pitch(t)
        else:
            thrust = u.thrust_vac if engines else 0.0
            pitch = guidance.upper_stage_pitch(r - R, vr, vt, r, thrust / m)

        new_state, extra = rk4_step(state, dt, engines, pitch_direction(pitch), vostok.cda)
        for block in active:
            block.fuel = max(0.0, block.fuel - block.mass_flow * dt)

        if t >= next_log:
            log.add(t, state, extra, pitch)
            next_log += log_every

        state = new_state
        t += dt

        if not vostok.has_core and not cutoff:
            apo, peri, _, _ = orbit_from_state(state[0], state[2], state[3])
            if apo >= mp.TARGET_APOAPSIS or u.fuel <= 0:
                cutoff = True
                t_cutoff = t
                events.append((t, "Выключение блока Е" if u.fuel > 0 else "Блок Е: топливо кончилось"))
        if cutoff and t >= t_cutoff + mp.SPACECRAFT_SEP_DELAY:
            vostok.has_upper = False
            state[4] = vostok.total_mass
            events.append((t, "Отделение корабля"))
            break

    log.add(t, state, {"g_load": 0.0, "v_air": 0.0, "q": 0.0}, 0.0)
    return log.arrays(), events, state, t


def simulate_descent(state, t0, spacecraft_mass=mp.KSP_SPACECRAFT_MASS,
                     capsule_mass=mp.KSP_DESCENT_MODULE_MASS, log_every=1.0):
    r, phi, vr, vt, _ = state
    state = [r, phi, vr, vt, spacecraft_mass]
    _, _, a, _ = orbit_from_state(r, vr, vt)
    period = orbital_period(a)
    t_retro = t0 + mp.RETRO_ORBIT_FRACTION * period

    retro_thrust, retro_isp_sl, retro_isp_vac = mp.KSP_RETRO_ENGINE
    retro_flow = retro_thrust / (retro_isp_vac * G0)
    events = [(t0, "Корабль на орбите, период %.1f мин" % (period / 60))]
    log = Log()
    t = t0
    next_log = t0
    phase = "orbit"
    t_burn_end = None
    phi_retro = None
    cda = 0.0
    max_g = 0.0
    info = {"period": period}

    while True:
        r, phi, vr, vt, m = state
        h = r - R
        engines = []
        direction = (0.0, 0.0)

        if phase == "orbit" and t >= t_retro:
            phase = "burn"
            phi_retro = phi
            info["m_before_burn"] = m
            events.append((t, "Включение ТДУ"))
        if phase == "burn":
            _, peri, _, _ = orbit_from_state(r, vr, vt)
            if peri <= mp.REENTRY_PERIAPSIS:
                phase = "after_burn"
                t_burn_end = t
                info["retro_burn_time"] = t - t_retro
                info["retro_dv"] = tsiolkovsky(retro_isp_vac, info["m_before_burn"], m)
                events.append((t, "Выключение ТДУ: %.0f с, Δv = %.0f м/с"
                               % (info["retro_burn_time"], info["retro_dv"])))
            else:
                engines = [(retro_flow, retro_isp_sl, retro_isp_vac)]
                v = math.hypot(vr, vt)
                direction = (-vr / v, -vt / v)
        if phase == "after_burn" and t >= t_burn_end + mp.MODULE_SEP_DELAY:
            phase = "descent"
            state[4] = capsule_mass
            cda = mp.KSP_CDA_CAPSULE
            events.append((t, "Разделение СА и приборного отсека"))
        if phase == "descent" and h <= mp.EJECTION_ALTITUDE and "ejection" not in info:
            info["ejection"] = (t, math.hypot(vr, vt - OMEGA * r * ROTATION_IN_PLANE))
            info["max_g_reentry"] = max_g
            events.append((t, "Высота 7 км: катапультирование космонавта, парашют СА взведён"))
        if phase == "descent" and h <= mp.CHUTE_ALTITUDE and "chute" not in info:
            info["chute"] = t
            events.append((t, "Раскрытие основного парашюта"))
        if "chute" in info:
            x = min(1.0, (t - info["chute"]) / mp.CHUTE_OPENING_TIME)
            cda = mp.KSP_CDA_CAPSULE + mp.KSP_CDA_CHUTE_SEMI + x * (mp.KSP_CDA_CHUTE - mp.KSP_CDA_CHUTE_SEMI)
        elif "ejection" in info:
            x = min(1.0, (t - info["ejection"][0]) / mp.CHUTE_OPENING_TIME)
            cda = mp.KSP_CDA_CAPSULE + x * mp.KSP_CDA_CHUTE_SEMI

        in_air = h < mp.KERBIN_ATM_HEIGHT + 1000
        dt = 0.02 if h < mp.EJECTION_ALTITUDE else (0.05 if (in_air or phase == "burn") else 1.0)
        new_state, extra = rk4_step(state, dt, engines, direction, cda)
        if phase == "descent":
            max_g = max(max_g, extra["g_load"])
        if t >= next_log:
            log.add(t, state, extra, 0.0)
            next_log += log_every if not in_air else 0.25
        state = new_state
        t += dt

        if state[0] - R <= 0:
            v_land = math.hypot(state[2], state[3] - OMEGA * state[0] * ROTATION_IN_PLANE)
            events.append((t, "Приземление, скорость %.1f м/с" % v_land))
            info["landing_time"] = t
            info["landing_speed"] = v_land
            info["range_from_retro_km"] = (state[1] - phi_retro) * R / 1000
            break
        if t > t0 + 3 * period:
            events.append((t, "Не вошёл в атмосферу за 3 витка"))
            break

    info.setdefault("max_g_reentry", max_g)
    info["max_g"] = max_g
    return log.arrays(), events, info


def read_ksp_flight(csv_path):
    with open(csv_path, encoding="utf-8") as f:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]
