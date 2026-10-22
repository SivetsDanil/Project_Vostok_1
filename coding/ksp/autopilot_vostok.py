import math
import os
import sys
import time

import krpc
import matplotlib.pyplot as plt

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import mission_params as mp
import guidance
import time_warp
from telemetry import Telemetry

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "mat_model", "data")


def smooth(value, target, smooth_factor=0.1):
    return value + (target - value) * smooth_factor


def event(tel, text):
    print(f"[T+{tel.t:7.1f}] {text}")
    tel.add_event(text)


def fuel_in_next_decouple_stage(vessel):
    res = vessel.resources_in_decouple_stage(stage=vessel.control.current_stage - 1, cumulative=False)
    return res.amount("LiquidFuel")


def launch_and_ascent(sc, vessel, tel):
    body = vessel.orbit.body
    flight_inertial = vessel.flight(body.non_rotating_reference_frame)
    ap = vessel.auto_pilot
    base_heading = mp.launch_heading()
    print(f"[ascent] курс старта {base_heading:.1f}° (наклонение {mp.TARGET_INCLINATION}°)")

    ap.reference_frame = vessel.surface_reference_frame
    ap.target_pitch_and_heading(90, base_heading)
    ap.engage()

    vessel.control.throttle = 1.0
    vessel.control.activate_next_stage()
    time.sleep(1.5)
    vessel.control.activate_next_stage()
    tel.t0 = sc.ut
    event(tel, "Подъём! Поехали!")

    boosters_attached = True
    fairing_attached = True
    core_attached = True
    pitch = 90.0
    t_cutoff = None

    while True:
        t = tel.t
        alt = flight_inertial.mean_altitude
        r = alt + body.equatorial_radius

        if boosters_attached and fuel_in_next_decouple_stage(vessel) < 0.1:
            vessel.control.activate_next_stage()
            boosters_attached = False
            event(tel, "Отделение боковых блоков")
        elif not boosters_attached and fairing_attached and t >= mp.T_FAIRING:
            vessel.control.activate_next_stage()
            fairing_attached = False
            event(tel, "Сброс обтекателя")
        elif not fairing_attached and core_attached and \
                (t >= mp.T_CORE_SEP or fuel_in_next_decouple_stage(vessel) < 0.1):
            vessel.control.activate_next_stage()
            core_attached = False
            event(tel, "Отделение блока А, запуск блока Е")

        if core_attached:
            pitch = guidance.ascent_pitch(t)
        else:
            a_thrust = vessel.thrust / vessel.mass
            target = guidance.upper_stage_pitch(alt, flight_inertial.vertical_speed,
                                                flight_inertial.horizontal_speed, r, a_thrust)
            pitch = smooth(pitch, target, 0.1)

        heading = base_heading
        if vessel.orbit.speed > 1000:
            error = math.degrees(vessel.orbit.inclination) - mp.TARGET_INCLINATION
            heading += max(-5.0, min(5.0, error))
        ap.target_pitch_and_heading(pitch, heading)

        if not core_attached and t > mp.T_CORE_SEP + 5:
            no_fuel = vessel.thrust == 0
            if vessel.orbit.apoapsis_altitude >= mp.TARGET_APOAPSIS or no_fuel:
                vessel.control.throttle = 0.0
                t_cutoff = t
                event(tel, "Выключение блока Е" + (": кончилось топливо" if no_fuel else ""))
                break

        time.sleep(0.05)

    while tel.t < t_cutoff + mp.SPACECRAFT_SEP_DELAY:
        time.sleep(0.1)
    vessel.control.activate_next_stage()
    event(tel, "Отделение корабля — выход на орбиту")
    orbit = vessel.orbit
    print(f"[orbit] апоцентр {orbit.apoapsis_altitude / 1000:.1f} км, "
          f"перицентр {orbit.periapsis_altitude / 1000:.1f} км, "
          f"наклонение {math.degrees(orbit.inclination):.2f}°, период {orbit.period / 60:.1f} мин")
    return t_cutoff


def coast_and_retro(sc, vessel, tel):
    ap = vessel.auto_pilot
    ap.reference_frame = vessel.orbital_reference_frame
    ap.target_direction = (0, -1, 0)
    ap.engage()

    retro_ut = sc.ut + mp.RETRO_ORBIT_FRACTION * vessel.orbit.period
    print(f"[orbit] торможение через {(retro_ut - sc.ut) / 60:.1f} мин (доля витка как у Гагарина)")
    tel.period = 2.0
    time_warp.safe_warp_to(sc, vessel, retro_ut - 30)
    ap.reference_frame = vessel.orbital_reference_frame
    ap.target_direction = (0, -1, 0)
    ap.engage()
    tel.period = 0.1
    while sc.ut < retro_ut:
        time.sleep(0.05)

    m0 = vessel.mass
    t_start = tel.t
    vessel.control.throttle = 1.0
    event(tel, "Включение ТДУ")
    time.sleep(0.5)
    while vessel.orbit.periapsis_altitude > mp.REENTRY_PERIAPSIS and vessel.thrust > 0:
        time.sleep(0.05)
    vessel.control.throttle = 0.0
    burn_time = tel.t - t_start
    dv = mp.KSP_RETRO_ENGINE[2] * mp.G0 * math.log(m0 / vessel.mass)
    event(tel, f"Выключение ТДУ: {burn_time:.1f} с, dv = {dv:.0f} м/с")

    time.sleep(mp.MODULE_SEP_DELAY)
    vessel.control.activate_next_stage()
    event(tel, "Разделение СА и приборного отсека")
    return burn_time, dv


def descent(sc, vessel, tel):
    ap = vessel.auto_pilot
    body = vessel.orbit.body
    flight_srf = vessel.flight(body.reference_frame)

    orbit = vessel.orbit
    try:
        ta = -orbit.true_anomaly_at_radius(body.equatorial_radius + mp.KERBIN_ATM_HEIGHT + 2000)
        entry_ut = orbit.ut_at_true_anomaly(ta)
        if entry_ut < sc.ut:
            entry_ut += orbit.period
    except Exception:
        entry_ut = sc.ut + orbit.time_to_periapsis - 120
    tel.period = 2.0
    time_warp.safe_warp_to(sc, vessel, entry_ut)
    tel.period = 0.1

    ap.reference_frame = vessel.surface_velocity_reference_frame
    ap.target_direction = (0, -1, 0)
    ap.engage()
    event(tel, "Вход в атмосферу")

    max_g = 0.0
    ejected = False
    holding = True
    while True:
        alt = flight_srf.mean_altitude
        max_g = max(max_g, vessel.flight().g_force)
        if holding and alt < 30000:
            ap.disengage()
            holding = False
        if not ejected and alt <= mp.EJECTION_ALTITUDE:
            ejected = True
            vessel.control.activate_next_stage()
            event(tel, "Высота 7 км: катапультирование (выйдите в EVA и раскройте парашют вручную)")
        if vessel.situation in (sc.VesselSituation.landed, sc.VesselSituation.splashed):
            event(tel, f"Приземление, скорость {flight_srf.speed:.1f} м/с")
            break
        time.sleep(0.1)
    return max_g


def print_summary(tel, t_cutoff, burn_time, dv, max_g_reentry):
    print("\n================ ИТОГИ: KSP vs «Восток-1» ================")
    for t, name in tel.events:
        ref = [t_real for t_real, real_name in mp.REAL_EVENTS if name.startswith(real_name)]
        print(f"T+{t:7.1f} с  {name}" + (f"   (реально T+{ref[0]:.0f} с)" if ref else ""))
    max_g_ascent = max((r["g_load"] for r in tel.rows if 0 <= r["t"] <= t_cutoff), default=0)
    print(f"Перегрузка на выведении: {max_g_ascent:.1f} g (реально ~{mp.REAL_ASCENT_MAX_G:.0f} g)")
    print(f"ТДУ: {burn_time:.0f} с, {dv:.0f} м/с (реально {mp.REAL_RETRO_BURN_TIME:.0f} с, ~{mp.REAL_RETRO_DV:.0f} м/с)")
    print(f"Перегрузка на спуске: {max_g_reentry:.1f} g (реально {mp.REAL_REENTRY_MAX_G[0]:.0f}-"
          f"{mp.REAL_REENTRY_MAX_G[1]:.0f} g)")
    print(f"Полёт: {tel.t / 60:.0f} мин (реально {mp.REAL_FLIGHT_TIME / 60:.0f} мин)")


def plot(tel):
    t = [r["t"] for r in tel.rows]
    for key, label, title in (("alt", "Высота, м", "График высоты"),
                              ("speed", "Скорость, м/с", "График скорости"),
                              ("g_load", "Перегрузка, g", "График перегрузки")):
        plt.plot(t, [r[key] for r in tel.rows])
        plt.xlabel("Время, с")
        plt.ylabel(label)
        plt.title(title)
        plt.grid(True)
        plt.show()


if __name__ == "__main__":
    conn = krpc.connect(name="Vostok-1 autopilot")
    sc = conn.space_center
    vessel = sc.active_vessel

    tel = Telemetry()
    tel.start()
    try:
        t_cutoff = launch_and_ascent(sc, vessel, tel)
        burn_time, dv = coast_and_retro(sc, vessel, tel)
        max_g_reentry = descent(sc, vessel, tel)
    finally:
        tel.stop()
        os.makedirs(DATA_DIR, exist_ok=True)
        tel.save(os.path.join(DATA_DIR, "ksp_flight.csv"), os.path.join(DATA_DIR, "ksp_events.csv"))

    print_summary(tel, t_cutoff, burn_time, dv, max_g_reentry)
    plot(tel)
