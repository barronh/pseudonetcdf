__all__ = ['camxfiles_paths', 'net_balance_paths', 'geoschemfiles_paths']

from os.path import join

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

for key, val in camxfiles_paths.iteritems():
    camxfiles_paths[key] = join(*__path__ + val.split('/'))

net_balance_paths = dict(
                 mrg_file = 'net_balance/test.mrg_file',
                 net_file = 'net_balance/test.net_file',
                 )

for key, val in net_balance_paths.iteritems():
    net_balance_paths[key] = join(*__path__ + val.split('/'))

geoschemfiles_paths = dict(bpch = 'geoschemfiles/test.bpch')

for key, val in geoschemfiles_paths.iteritems():
    geoschemfiles_paths[key] = join(*__path__ + val.split('/'))
