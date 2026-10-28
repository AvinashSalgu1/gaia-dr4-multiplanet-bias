"""Physical constants (SI), matching gaiamock's values where it defines them."""

G_SI = 6.6743e-11                 # m^3 kg^-1 s^-2
MSUN_KG = 1.98840987069805e30     # kg
AU_M = 1.4959787e11               # m
DAY_S = 86400.0                   # s
MJUP_MSUN = 9.547919e-4           # M_Jup / M_sun (IAU 2015 nominal GM ratio)

# Gaia data-release reference epochs (JD) and observation windows,
# as used by gaiamock.rescale_times_astrometry / get_gost_one_position.
DR_REF_JD = {"dr3": 2457389.0, "dr4": 2457936.875, "dr5": 2458818.5}
N_CCD_PER_TRANSIT = 8
