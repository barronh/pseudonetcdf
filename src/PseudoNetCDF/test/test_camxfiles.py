__all__ = ['test_cloud_rain_Memmap',
           'test_height_pressure_Memmap',
           'test_height_pressure_Read',
           'test_humidity_Memmap',
           'test_humidity_Read',
           'test_ipr_Memmap',
           'test_ipr_Read',
           'test_irr_Read',
           'test_irr_Memmap',
           'test_landuse_Memmap',
           'test_lateral_boundary_Memmap',
           'test_one3d_Memmap',
           'test_one3d_Read',
           'test_point_source_Memmap',
           'test_point_source_Read',
           'test_temperature_Memmap',
           'test_temperature_Read',
           'test_uamiv_Memmap',
           'test_uamiv_Read',
           'test_uamiv_Write',
           'test_vertical_diffusivity_Memmap',
           'test_vertical_diffusivity_Read',
           'test_wind_Memmap',
           'test_wind_Read']

import PseudoNetCDF.camxfiles as cx
test_cloud_rain_Memmap = cx.cloud_rain.Memmap.TestMemmap
# test_finst_Memmap = cx.finst.Memmap.TestMemmap
test_height_pressure_Memmap = cx.height_pressure.Memmap.TestMemmap
test_height_pressure_Read = cx.height_pressure.Read.TestRead
test_humidity_Memmap = cx.humidity.Memmap.TestMemmap
test_humidity_Read = cx.humidity.Read.TestRead
test_ipr_Memmap = cx.ipr.Memmap.TestMemmap
test_ipr_Read = cx.ipr.Read.TestRead
test_irr_Read = cx.irr.Read.TestRead
test_irr_Memmap = cx.irr.Memmap.TestMemmap
test_landuse_Memmap = cx.landuse.Memmap.TestMemmap
test_lateral_boundary_Memmap = cx.lateral_boundary.Memmap.TestMemmap
test_one3d_Memmap = cx.one3d.Memmap.TestMemmap
test_one3d_Read = cx.one3d.Read.TestRead
test_point_source_Memmap = cx.point_source.Memmap.TestMemmap
test_point_source_Read = cx.point_source.Read.TestRead
test_temperature_Memmap = cx.temperature.Memmap.TestMemmap
test_temperature_Read = cx.temperature.Read.TestRead
test_uamiv_Memmap = cx.uamiv.Memmap.TestMemmap
test_uamiv_Read = cx.uamiv.Read.TestuamivRead
test_uamiv_Write = cx.uamiv.Write.TestMemmaps
test_vertical_diffusivity_Memmap = cx.vertical_diffusivity.Memmap.TestMemmap
test_vertical_diffusivity_Read = cx.vertical_diffusivity.Read.TestRead
test_wind_Memmap = cx.wind.Memmap.TestMemmap
test_wind_Read = cx.wind.Read.TestRead
