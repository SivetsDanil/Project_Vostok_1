import time


def recover_physics(sc, vessel=None, wait=0.6):
    try:
        if sc.warp_rate > 0:
            sc.warp_to(sc.ut)
    except:
        pass
    time.sleep(wait)
    if vessel is not None:
        try:
            vessel.control.throttle = 0.0
            vessel.auto_pilot.disengage()
            time.sleep(0.1)
            vessel.auto_pilot.engage()
        except:
            pass
    time.sleep(0.1)


def safe_warp_to(sc, vessel, target_ut):
    now = sc.ut
    rem = target_ut - now
    print(f"[warp] now={now:.2f}, target={target_ut:.2f}, rem={rem:.2f}")
    if rem <= 0:
        print("[warp] Целевой момент уже прошёл.")
        return
    if rem > 600:
        t = target_ut - 600
        print(f"[warp] coarse warp → UT {t:.2f}")
        try:
            sc.warp_to(t)
        except Exception as e:
            print("[warp] coarse error:", e)
        recover_physics(sc)

    now = sc.ut
    rem = target_ut - now
    if rem > 120:
        t = target_ut - 120
        print(f"[warp] medium warp → UT {t:.2f}")
        try:
            sc.warp_to(t)
        except Exception as e:
            print("[warp] medium error:", e)
        recover_physics(sc)

    now = sc.ut
    rem = target_ut - now
    if rem > 8:
        t = target_ut - 8
        print(f"[warp] fine warp → UT {t:.2f}")
        try:
            sc.warp_to(t)
        except Exception as e:
            print("[warp] fine error:", e)
        recover_physics(sc, wait=0.8)

    print("[warp] финальная синхронизация в realtime...")
    timeout = 12
    t0 = time.time()

    while sc.ut < target_ut:
        if time.time() - t0 > timeout:
            print("[warp] WARNING: таймаут ожидания UT!")
            break
        time.sleep(0.05)

    recover_physics(sc, vessel, wait=0.4)
    print(f"[warp] Достигнут UT: {sc.ut:.2f} (цель: {target_ut:.2f})")
