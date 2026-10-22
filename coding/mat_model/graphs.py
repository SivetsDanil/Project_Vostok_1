import os
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import mission_params as mp
from objects import Block, Vostok
import calc

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
KSP_FLIGHT = os.path.join(DATA_DIR, "ksp_flight.csv")


def build_vostok():
    return Vostok(
        boosters=Block("Б, В, Г, Д", *mp.KSP_BOOSTER),
        core=Block("А", *mp.KSP_CORE),
        upper=Block("Е", *mp.KSP_UPPER),
        fairing_mass=mp.KSP_FAIRING_MASS,
        spacecraft_mass=mp.KSP_SPACECRAFT_MASS,
        cda_packet=mp.KSP_CDA_PACKET,
        cda_core=mp.KSP_CDA_CORE,
    )


def load_ksp():
    if not os.path.exists(KSP_FLIGHT):
        print("Данных KSP нет — строим только мат. модель (запустите ksp/autopilot_vostok.py)")
        return None
    rows = calc.read_ksp_flight(KSP_FLIGHT)
    ksp = {k: np.array([r[k] for r in rows if r["t"] >= 0]) for k in rows[0]}
    if len(ksp["t"]) > 1:
        ksp["accel"] = np.gradient(ksp["speed"], ksp["t"])
    if calc.use_ksp_atmosphere(KSP_FLIGHT):
        print("Атмосфера Кербина для модели взята из телеметрии KSP")
    return ksp


def mark_real_events(t_max):
    for t, name in mp.REAL_EVENTS:
        if t <= t_max:
            plt.axvline(t, color="gray", linestyle=":", linewidth=1)
            plt.text(t, plt.ylim()[1], " " + name, rotation=90, va="top", fontsize=7, color="gray")


def compare(model, ksp, key, ylabel, title, t_max, scale=1.0):
    plt.figure()
    mask = model["t"] <= t_max
    plt.plot(model["t"][mask], model[key][mask] * scale, "b-", linewidth=1, label="Мат. модель")
    if ksp is not None and key in ksp:
        mask = ksp["t"] <= t_max
        plt.plot(ksp["t"][mask], ksp[key][mask] * scale, "r--", linewidth=1, label="KSP")
    plt.xlabel("Время, с")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    mark_real_events(t_max)
    plt.legend()
    plt.show()


if __name__ == "__main__":
    ksp = load_ksp()
    vostok = build_vostok()

    print("Стартовая масса: %.1f т" % (vostok.total_mass / 1000))
    print("Время работы: боковые блоки %.0f с, блок А %.0f с, блок Е %.0f с"
          % (vostok.boosters.burn_time, vostok.core.burn_time, vostok.upper.burn_time))
    for name, dv in calc.delta_v_budget(vostok):
        print("Запас Δv, %s: %.0f м/с" % (name, dv))

    ascent, events, state, t_sep = calc.simulate_ascent(vostok)
    apo, peri, a, e = calc.orbit_from_state(state[0], state[2], state[3])
    descent, descent_events, info = calc.simulate_descent(state, t_sep)

    print("\n--- Циклограмма: модель / реальный «Восток-1» ---")
    real = dict((name, t) for t, name in mp.REAL_EVENTS)
    for t, name in events + descent_events:
        ref = [t_real for real_name, t_real in real.items() if name.startswith(real_name)]
        print("T+%7.1f с  %s" % (t, name) + ("   (реально T+%.0f с)" % ref[0] if ref else ""))

    print("\n--- Орбита: модель / реальная ---")
    print("Перицентр %.0f км / %.0f км" % (peri / 1000, mp.REAL_PERIGEE / 1000))
    print("Апоцентр  %.0f км / %.0f км" % (apo / 1000, mp.REAL_APOGEE / 1000))
    print("Период    %.1f мин / %.1f мин (Кербин меньше Земли)" % (info["period"] / 60, mp.REAL_PERIOD / 60))
    print("Перегрузка на выведении %.1f g / ~%.0f g" % (ascent["g_load"].max(), mp.REAL_ASCENT_MAX_G))
    print("ТДУ: %.0f с, %.0f м/с / %.0f с, ~%.0f м/с" % (info.get("retro_burn_time", 0), info.get("retro_dv", 0),
                                                        mp.REAL_RETRO_BURN_TIME, mp.REAL_RETRO_DV))
    print("Перегрузка на спуске %.1f g / %.0f-%.0f g" % (info["max_g_reentry"], *mp.REAL_REENTRY_MAX_G))
    print("Перегрузка при раскрытии парашюта %.1f g" % info["max_g"])
    if "landing_time" in info:
        print("Полёт %.0f мин / %.0f мин" % (info["landing_time"] / 60, mp.REAL_FLIGHT_TIME / 60))

    t_max = t_sep + 20
    compare(ascent, ksp, "speed", "Скорость, м/с", "Зависимость скорости от времени", t_max)
    compare(ascent, ksp, "alt", "Высота, км", "Зависимость высоты от времени", t_max, scale=1e-3)
    compare(ascent, ksp, "accel", "Ускорение, м/с²", "Зависимость ускорения от времени", t_max)
    compare(ascent, ksp, "g_load", "Перегрузка, g", "Перегрузка на выведении", t_max)

    plt.figure()
    plt.plot(descent["t"], descent["alt"] / 1000, "b-", linewidth=1, label="Мат. модель")
    if ksp is not None:
        mask = ksp["t"] >= t_sep
        plt.plot(ksp["t"][mask], ksp["alt"][mask] / 1000, "r--", linewidth=1, label="KSP")
    plt.xlabel("Время, с")
    plt.ylabel("Высота, км")
    plt.title("Виток, торможение и спуск")
    plt.grid(True)
    plt.legend()
    plt.show()

    if ksp is not None:
        plt.figure()
        plt.plot(ksp["lon"], ksp["lat"], "r.", markersize=1)
        plt.xlabel("Долгота, °")
        plt.ylabel("Широта, °")
        plt.title("Трасса полёта над Кербином (KSP)")
        plt.grid(True)
        plt.show()
