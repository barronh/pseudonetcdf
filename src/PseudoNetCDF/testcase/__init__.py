__all__ = ['camxfiles_paths', 'net_balance_paths', 'geoschemfiles_paths', 'icarttfiles_paths', 'all_paths', 'self_described_paths']

from os.path import join, realpath, abspath

camxfiles_paths = dict(wind = 'camxfiles/wind/test.wind',
                 landuse = 'camxfiles/landuse/test.landuse',
                 temperature = 'camxfiles/temperature/test.temperature',
                 cloud_rain = 'camxfiles/cloud_rain/test.cloud_rain',
                 vertical_diffusivity = 'camxfiles/vertical_diffusivity/test.vertical_diffusivity',
                 humidity = 'camxfiles/humidity/test.humidity',
                 height_pressure = 'camxfiles/height_pressure/test.height_pressure',
                 uamiv = 'camxfiles/uamiv/test.uamiv',
                 point_source = 'camxfiles/point_source/test.point_source',
                 lateral_boundary = 'camxfiles/lateral_boundary/test.lateral_boundary')

for key, val in camxfiles_paths.items():
    camxfiles_paths[key] = abspath(join(*__path__ + val.split('/')))

net_balance_paths = dict(
                 mrg_file = 'net_balance/test.mrg_file',
                 net_file = 'net_balance/test.net_file',
                 )

for key, val in net_balance_paths.items():
    net_balance_paths[key] = join(*__path__ + val.split('/'))

geoschemfiles_paths = dict(bpch = 'geoschemfiles/test.bpch')

for key, val in geoschemfiles_paths.items():
    geoschemfiles_paths[key] = join(*__path__ + val.split('/'))

icarttfiles_paths = dict(ffi1001 = join(*__path__ +  ['icarttfiles', 'test.ffi1001']))

all_paths = dict()
all_paths.update(camxfiles_paths, **geoschemfiles_paths)
all_paths.update(icarttfiles_paths, **net_balance_paths)

self_described_paths = dict()
for k in 'uamiv point_source lateral_boundary humidity vertical_diffusivity'.split():
    self_described_paths[k] = camxfiles_paths[k]

self_described_paths['bpch'] = geoschemfiles_paths['bpch']
self_described_paths['bpch2'] = geoschemfiles_paths['bpch']
self_described_paths['ffi1001'] = icarttfiles_paths['ffi1001']
