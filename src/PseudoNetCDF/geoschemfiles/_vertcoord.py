import numpy as np

# GEOS-5 Ap [hPa] for 72 levels (73 edges)
geos5_native_a = np.array([0.000000e+00, 4.804826e-02, 6.593752e+00, 1.313480e+01, 1.961311e+01, 2.609201e+01, 3.257081e+01, 3.898201e+01, 4.533901e+01, 5.169611e+01, 5.805321e+01, 6.436264e+01, 7.062198e+01, 7.883422e+01, 8.909992e+01, 9.936521e+01, 1.091817e+02, 1.189586e+02, 1.286959e+02, 1.429100e+02, 1.562600e+02, 1.696090e+02, 1.816190e+02, 1.930970e+02, 2.032590e+02, 2.121500e+02, 2.187760e+02, 2.238980e+02, 2.243630e+02, 2.168650e+02, 2.011920e+02, 1.769300e+02, 1.503930e+02, 1.278370e+02, 1.086630e+02, 9.236572e+01, 7.851231e+01, 6.660341e+01, 5.638791e+01, 4.764391e+01, 4.017541e+01, 3.381001e+01, 2.836781e+01, 2.373041e+01, 1.979160e+01, 1.645710e+01, 1.364340e+01, 1.127690e+01, 9.292942e+00, 7.619842e+00, 6.216801e+00, 5.046801e+00, 4.076571e+00, 3.276431e+00, 2.620211e+00, 2.084970e+00, 1.650790e+00, 1.300510e+00, 1.019440e+00, 7.951341e-01, 6.167791e-01, 4.758061e-01, 3.650411e-01, 2.785261e-01, 2.113490e-01, 1.594950e-01, 1.197030e-01, 8.934502e-02, 6.600001e-02, 4.758501e-02, 3.270000e-02, 2.000000e-02, 1.000000e-02])

# Bp [unitless] for 72 levels (73 edges)
geos5_native_b = np.array([1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01, 9.203870e-01, 8.989080e-01, 8.774290e-01, 8.560180e-01, 8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01, 7.493782e-01, 7.211660e-01, 6.858999e-01, 6.506349e-01, 6.158184e-01, 5.810415e-01, 5.463042e-01, 4.945902e-01, 4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01, 2.467411e-01, 2.003501e-01, 1.562241e-01, 1.136021e-01, 6.372006e-02, 2.801004e-02, 6.960025e-03, 8.175413e-09, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00])

# The reduced vertical coordinate skips edges
geos5_reduced_idx = np.array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 38, 40, 42, 44, 48, 52, 56, 60, 64, 68, 72])
geos5_reduced_a = geos5_native_a[geos5_reduced_idx]
geos5_reduced_b = geos5_native_b[geos5_reduced_idx]

# GEOS-4 AP [hPa] for 55 levels (56 edges)
geos4_native_a = np.array([0.000000e0, 0.000000e0, 12.704939e0, 35.465965e0, 66.098427e0, 101.671654e0, 138.744400e0, 173.403183e0, 198.737839e0, 215.417526e0, 223.884689e0, 224.362869e0, 216.864929e0, 201.192093e0, 176.929993e0, 150.393005e0, 127.837006e0, 108.663429e0, 92.365662e0, 78.512299e0, 66.603378e0, 56.387939e0, 47.643932e0, 40.175419e0, 33.809956e0, 28.367815e0, 23.730362e0, 19.791553e0, 16.457071e0, 13.643393e0, 11.276889e0, 9.292943e0, 7.619839e0, 6.216800e0, 5.046805e0, 4.076567e0, 3.276433e0, 2.620212e0, 2.084972e0, 1.650792e0, 1.300508e0, 1.019442e0, 0.795134e0, 0.616779e0, 0.475806e0, 0.365041e0, 0.278526e0, 0.211349e0, 0.159495e0, 0.119703e0, 0.089345e0, 0.066000e0, 0.047585e0, 0.032700e0, 0.020000e0, 0.010000e0])

# GEOS-4 BP [unitless] for 55 levels (56 edges)
geos4_native_b = np.array([1.000000e0, 0.985110e0, 0.943290e0, 0.867830e0, 0.764920e0, 0.642710e0, 0.510460e0, 0.378440e0, 0.270330e0, 0.183300e0, 0.115030e0, 0.063720e0, 0.028010e0, 0.006960e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0, 0.000000e0])

# The reduced vertical coordinate skips edges
geos4_reduced_idx = np.array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 23, 25, 27, 31, 35, 39, 43, 47, 51, 55])

geos4_reduced_a = geos4_native_a[geos4_reduced_idx]
geos4_reduced_b = geos4_native_b[geos4_reduced_idx]

geos_hyai = {'GEOS-5-NATIVE': geos5_native_a,
             'GEOS-5-REDUCED': geos5_reduced_a,
             'GEOS-4-NATIVE': geos4_native_a,
             'GEOS-4-REDUCED': geos4_reduced_a,
             'MERRA-NATIVE': geos5_native_a,
             'MERRA-REDUCED': geos5_reduced_a
            }
geos_hyam = dict([(k, np.convolve(v, [.5, .5], 'valid')) for k, v in geos_hyai.iteritems()])

geos_hybi = {'GEOS-5-NATIVE': geos5_native_b,
             'GEOS-5-REDUCED': geos5_reduced_b,
             'GEOS-4-NATIVE': geos4_native_b,
             'GEOS-4-REDUCED': geos4_reduced_b,
             'MERRA-NATIVE': geos5_native_b,
             'MERRA-REDUCED': geos5_reduced_b
            }
h0 = 7.6 # km
geos_hybm = dict([(k, np.convolve(v, [.5, .5], 'valid')) for k, v in geos_hybi.iteritems()])

geos_etai = {}
geos_etai_pressure = {}
geos_etai_height = {}
geos_etam = {}
geos_etam_pressure = {}
geos_etam_height = {}
geos_eta_slp = 1013.25
for k, ai in geos_hyai.iteritems():
    bi = geos_hybi[k]
    am = geos_hyam[k]
    bm = geos_hybm[k]
    pressurei = geos_etai_pressure[k] = ai + bi * geos_eta_slp
    pressurem = geos_etam_pressure[k] = am + bm * geos_eta_slp
    ptop = pressurei[-1]
    geos_etai[k] = (pressurei - ptop) / (geos_eta_slp - ptop)
    geos_etam[k] = (pressurem - ptop) / (geos_eta_slp - ptop)
    
    geos_etai_height[k] = -np.log(pressurei / geos_eta_slp) * h0
    geos_etam_height[k] = -np.log(pressurem / geos_eta_slp) * h0

def std_alt():
    T_0 = 233.
    T_Delta = 75.
    g = 9.80665
    R = 287.05
    H = 10000.
    z = np.log((np.exp(np.log(pressurem / geos_eta_slp) / (-g/R * H/T_0)) * (T_0 + T_Delta) - T_Delta) / T_0) * H
    return z
    #p = geos_eta_slp * np.exp(-1 * g/R * H/T_0 * np.log( (np.exp(z/H)*T_0 + T_Delta)/(T_0 + T_Delta)))
