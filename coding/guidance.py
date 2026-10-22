import math

import mission_params as mp


def smoothstep(x):
    x = max(0.0, min(1.0, x))
    return 3 * x * x - 2 * x * x * x


def ascent_pitch(t):
    if t < mp.T_VERTICAL:
        return 90.0
    if t < mp.T_BOOSTER_BURN:
        s = smoothstep((t - mp.T_VERTICAL) / (mp.T_BOOSTER_BURN - mp.T_VERTICAL))
        return 90.0 + (mp.PITCH_AT_BOOSTER_SEP - 90.0) * s
    x = min(1.0, (t - mp.T_BOOSTER_BURN) / (mp.T_CORE_SEP - mp.T_BOOSTER_BURN))
    return mp.PITCH_AT_BOOSTER_SEP + (mp.PITCH_AT_CORE_SEP - mp.PITCH_AT_BOOSTER_SEP) * x


def upper_stage_pitch(alt, v_vertical, v_horizontal, r, a_thrust,
                      target_alt=mp.TARGET_PERIAPSIS, mu=mp.KERBIN_MU):
    vr_min, vr_max = mp.UPPER_VR_LIMITS
    vr_cmd = max(vr_min, min(vr_max, (target_alt - alt) / mp.UPPER_TAU))
    a_vertical = (vr_cmd - v_vertical) / 10.0 + mu / r ** 2 - v_horizontal ** 2 / r
    if a_thrust <= 0:
        return 0.0
    s = max(-0.5, min(0.95, a_vertical / a_thrust))
    return math.degrees(math.asin(s))
