import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import mission_params as mp
from objects import Block, Vostok
import calc

FUEL_MARGIN = 0.06


def booster_parts(limiter):
    engine, count, _, _, structure = mp.KSP_BOOSTER_PARTS
    thrust, _, isp_vac, _ = mp.ENGINES[engine]
    propellant = thrust * limiter / (isp_vac * mp.G0) * mp.T_BOOSTER_BURN
    return engine, count, limiter, propellant, structure


def evaluate(x, dt=0.05):
    lb, lu, fuel, p1 = x
    upper = list(mp.KSP_UPPER_PARTS)
    upper[2], upper[3] = lu, fuel
    mp.PITCH_AT_BOOSTER_SEP = p1
    vostok = Vostok(Block("Б, В, Г, Д", *mp.ksp_block(booster_parts(lb))),
                    Block("А", *mp.KSP_CORE), Block("Е", *mp.ksp_block(tuple(upper))),
                    mp.KSP_FAIRING_MASS, mp.KSP_SPACECRAFT_MASS, mp.KSP_CDA_PACKET, mp.KSP_CDA_CORE)
    b, c = vostok.boosters, vostok.core
    twr = (b.thrust_vac * b.isp_sl / b.isp_vac + c.thrust_vac * c.isp_sl / c.isp_vac) / (vostok.total_mass * mp.G0)
    m0 = vostok.total_mass
    _, events, state, _ = calc.simulate_ascent(vostok, dt=dt)
    apo, peri, _, _ = calc.orbit_from_state(state[0], state[2], state[3])
    t_cut = next((t for t, name in events if name.startswith("Выключение")), None)
    return dict(twr=twr, m0=m0, apo=apo, peri=peri, t_cut=t_cut, fuel_left=vostok.upper.fuel / fuel)


def cost(x):
    if x[0] <= 0.05 or x[1] <= 0.01 or x[2] <= 100:
        return 1e9
    r = evaluate(x)
    if r["t_cut"] is None or r["t_cut"] <= mp.T_CORE_SEP + 1:
        return 1e4 + abs(min(r["apo"], 5e6) - mp.TARGET_APOAPSIS) / 1000
    j = ((r["t_cut"] - mp.REAL_T_UPPER_CUTOFF) / 5.0) ** 2
    j += ((r["peri"] - mp.TARGET_PERIAPSIS) / 3000.0) ** 2
    j += ((r["fuel_left"] - FUEL_MARGIN) / 0.02) ** 2
    j += ((r["twr"] - 1.39) / 0.1) ** 2
    return j


def nelder_mead(f, x0, steps, iterations=120):
    n = len(x0)
    simplex = [list(x0)] + [[x0[j] + (steps[i] if j == i else 0) for j in range(n)] for i in range(n)]
    values = [f(x) for x in simplex]
    for it in range(iterations):
        order = sorted(range(n + 1), key=lambda i: values[i])
        simplex = [simplex[i] for i in order]
        values = [values[i] for i in order]
        centroid = [sum(x[j] for x in simplex[:-1]) / n for j in range(n)]
        worst = simplex[-1]
        reflected = [c + (c - w) for c, w in zip(centroid, worst)]
        fr = f(reflected)
        if fr < values[0]:
            expanded = [c + 2 * (c - w) for c, w in zip(centroid, worst)]
            fe = f(expanded)
            simplex[-1], values[-1] = (expanded, fe) if fe < fr else (reflected, fr)
        elif fr < values[-2]:
            simplex[-1], values[-1] = reflected, fr
        else:
            contracted = [c + 0.5 * (w - c) for c, w in zip(centroid, worst)]
            fc = f(contracted)
            if fc < values[-1]:
                simplex[-1], values[-1] = contracted, fc
            else:
                best = simplex[0]
                simplex = [best] + [[b + 0.5 * (x - b) for b, x in zip(best, xx)] for xx in simplex[1:]]
                values = [values[0]] + [f(x) for x in simplex[1:]]
        if it % 10 == 0:
            print("  %3d  cost %.3f  x = %s" % (it, values[0], ["%.3f" % v for v in simplex[0]]), flush=True)
    i = min(range(n + 1), key=lambda k: values[k])
    return simplex[i]


if __name__ == "__main__":
    x0 = [mp.KSP_BOOSTER_PARTS[2], mp.KSP_UPPER_PARTS[2], mp.KSP_UPPER_PARTS[3], mp.PITCH_AT_BOOSTER_SEP]
    x = nelder_mead(cost, x0, steps=[0.03, 0.04, 300.0, 5.0])
    r = evaluate(x, dt=0.02)
    print("\nТяговооружённость на старте %.2f (реально 1,39), стартовая масса %.1f т" % (r["twr"], r["m0"] / 1000))
    print("Выключение блока Е T+%.0f с, орбита %.0f x %.0f км, запас топлива %.0f %%"
          % (r["t_cut"], r["peri"] / 1000, r["apo"] / 1000, r["fuel_left"] * 100))
    lb, lu, fuel, p1 = x
    print("\nВписать в mission_params.py:")
    print("KSP_BOOSTER_PARTS = (%r, %d, %.3f, %.0f, %.0f)" % booster_parts(lb))
    print("KSP_UPPER_PARTS = (%r, %d, %.3f, %.0f, %.0f)"
          % (mp.KSP_UPPER_PARTS[0], 1, lu, fuel, mp.KSP_UPPER_PARTS[4]))
    print("PITCH_AT_BOOSTER_SEP = %.1f" % p1)
