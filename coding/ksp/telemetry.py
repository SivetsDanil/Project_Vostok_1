import csv
import math
import threading
import time

import krpc

FIELDS = ["t", "alt", "speed", "speed_surface", "vertical_speed", "g_load", "mass", "thrust",
          "pitch", "heading", "pressure", "density", "drag", "lat", "lon",
          "apoapsis", "periapsis", "inclination", "stage"]


class Telemetry(threading.Thread):
    def __init__(self, period=0.1):
        super().__init__(daemon=True)
        self.period = period
        self.t0 = None
        self.rows = []
        self.events = []
        self._stop_flag = False

        self.conn = krpc.connect(name="Vostok telemetry")
        sc = self.conn.space_center
        vessel = sc.active_vessel
        body = vessel.orbit.body
        flight = vessel.flight()
        flight_srf = vessel.flight(body.reference_frame)
        orbit = vessel.orbit
        add = self.conn.add_stream

        self.ut = add(getattr, sc, "ut")
        self.streams = {
            "alt": add(getattr, flight, "mean_altitude"),
            "speed": add(getattr, orbit, "speed"),
            "speed_surface": add(getattr, flight_srf, "speed"),
            "vertical_speed": add(getattr, flight_srf, "vertical_speed"),
            "g_load": add(getattr, flight, "g_force"),
            "mass": add(getattr, vessel, "mass"),
            "thrust": add(getattr, vessel, "thrust"),
            "pitch": add(getattr, flight, "pitch"),
            "heading": add(getattr, flight, "heading"),
            "pressure": add(getattr, flight, "static_pressure"),
            "density": add(getattr, flight, "atmosphere_density"),
            "drag": add(getattr, flight, "drag"),
            "lat": add(getattr, flight, "latitude"),
            "lon": add(getattr, flight, "longitude"),
            "apoapsis": add(getattr, orbit, "apoapsis_altitude"),
            "periapsis": add(getattr, orbit, "periapsis_altitude"),
            "inclination": add(getattr, orbit, "inclination"),
            "stage": add(getattr, vessel.control, "current_stage"),
        }

    @property
    def t(self):
        return self.ut() - self.t0 if self.t0 is not None else 0.0

    def run(self):
        while not self._stop_flag:
            if self.t0 is not None:
                try:
                    row = {"t": self.t}
                    for name, stream in self.streams.items():
                        row[name] = stream()
                    row["drag"] = math.sqrt(sum(x * x for x in row["drag"]))
                    row["inclination"] = math.degrees(row["inclination"])
                    self.rows.append(row)
                except Exception as e:
                    print("[telemetry]", e)
            time.sleep(self.period)

    def add_event(self, text):
        self.events.append((self.t, text))

    def stop(self):
        self._stop_flag = True
        self.join(timeout=2)

    def save(self, flight_csv, events_csv):
        with open(flight_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(self.rows)
        with open(events_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["t", "event"])
            writer.writerows(self.events)
        print(f"[telemetry] сохранено {len(self.rows)} строк -> {flight_csv}")
