G0 = 9.80665


class Block:
    def __init__(self, name, count, full_mass, dry_mass, thrust_vac, isp_sl, isp_vac):
        self.name = name
        self.count = count
        self.full_mass = full_mass * count
        self.dry_mass = dry_mass * count
        self.thrust_vac = thrust_vac * count
        self.isp_sl = isp_sl
        self.isp_vac = isp_vac
        self.fuel = self.full_mass - self.dry_mass

    @property
    def mass_flow(self):
        return self.thrust_vac / (self.isp_vac * G0)

    @property
    def burn_time(self):
        return self.fuel / self.mass_flow

    @property
    def mass(self):
        return self.dry_mass + self.fuel


class Vostok:
    def __init__(self, boosters, core, upper, fairing_mass, spacecraft_mass,
                 cda_packet, cda_core):
        self.boosters = boosters
        self.core = core
        self.upper = upper
        self.fairing_mass = fairing_mass
        self.spacecraft_mass = spacecraft_mass
        self.cda_packet = cda_packet
        self.cda_core = cda_core

        self.has_boosters = True
        self.has_fairing = True
        self.has_core = True
        self.has_upper = True

    @property
    def total_mass(self):
        m = self.spacecraft_mass
        if self.has_fairing:
            m += self.fairing_mass
        for block, attached in ((self.boosters, self.has_boosters),
                                (self.core, self.has_core),
                                (self.upper, self.has_upper)):
            if attached:
                m += block.mass
        return m

    @property
    def cda(self):
        return self.cda_packet if self.has_boosters else self.cda_core
