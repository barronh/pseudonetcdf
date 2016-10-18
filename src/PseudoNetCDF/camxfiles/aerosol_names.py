__all__ = ['aerosol_names']
__doc__ = """
mech_aerosol_options keys:
    camx6_3_chemparam_2_cf
	camx6_3_chemparam_6_cmu
	camx6_3_chemparam_2_cf_vbs
	camx6_3_chemparam_6_none
	camx6_3_chemparam_2_none
	camx6_3_chemparam_4_none
	camx6_3_chemparam_6_cf_hg
	camx6_3_chemparam_3_none
	camx6_3_chemparam_6_cf
	camx6_3_chemparam_6_cf_vbs
	camx6_3_chemparam_4_cf
	camx6_3_chemparam_4_none_srfmod
	camx6_3_chemparam_5_cf
	camx6_3_chemparam_5_none
	camx6_3_chemparam_3_cf
	camx6_3_chemparam_6_inert
"""
mech_aerosol_options = {
    "camx6_3_chemparam_2_cf": ('PNO3', 'PSO4', 'PNH4', 'POA', 'SOA1', 'SOA2', 'SOA3', 'SOA4', 'SOA5', 'SOA6', 'SOA7', 'SOAH', 'SOPA', 'SOPB', 'PEC', 'FPRM', 'FCRS', 'CPRM', 'CCRS', 'NA', 'PCL', 'PH2O'),
    "camx6_3_chemparam_2_cf_vbs": ('PNO3', 'PSO4', 'PNH4', 'PAS0', 'PAS1', 'PAS2', 'PAS3', 'PAS4', 'PBS0', 'PBS1', 'PBS2', 'PBS3', 'PBS4', 'PAP0', 'PAP1', 'PAP2', 'PAP3', 'PAP4', 'PCP0', 'PCP1', 'PCP2', 'PCP3', 'PCP4', 'PFP0', 'PFP1', 'PFP2', 'PFP3', 'PFP4', 'PEC', 'FPRM', 'FCRS', 'CPRM', 'CCRS', 'NA', 'PCL', 'PH2O'),
    "camx6_3_chemparam_2_none": (),
    "camx6_3_chemparam_3_cf": ('PNO3', 'PSO4', 'PNH4', 'POA', 'SOA1', 'SOA2', 'SOA3', 'SOA4', 'SOA5', 'SOA6', 'SOA7', 'SOAH', 'SOPA', 'SOPB', 'PEC', 'FPRM', 'FCRS', 'CPRM', 'CCRS', 'NA', 'PCL', 'PH2O'),
    "camx6_3_chemparam_3_none": (),
    "camx6_3_chemparam_4_cf": ('PNO3', 'PSO4', 'PNH4', 'POA', 'SOA1', 'SOA2', 'SOA3', 'SOA4', 'SOA5', 'SOA6', 'SOA7', 'SOAH', 'SOPA', 'SOPB', 'PEC', 'FPRM', 'FCRS', 'CPRM', 'CCRS', 'NA', 'PCL', 'PH2O'),
    "camx6_3_chemparam_4_none": (),
    "camx6_3_chemparam_4_none_srfmod": (),
    "camx6_3_chemparam_5_cf": ('PNO3', 'PSO4', 'PNH4', 'POA', 'SOA1', 'SOA2', 'SOA3', 'SOA4', 'SOA5', 'SOA6', 'SOA7', 'SOAH', 'SOPA', 'SOPB', 'PEC', 'FPRM', 'FCRS', 'CPRM', 'CCRS', 'NA', 'PCL', 'PH2O'),
    "camx6_3_chemparam_5_none": (),
    "camx6_3_chemparam_6_cf": ('PNO3', 'PSO4', 'PNH4', 'POA', 'SOA1', 'SOA2', 'SOA3', 'SOA4', 'SOA5', 'SOA6', 'SOA7', 'SOAH', 'SOPA', 'SOPB', 'PEC', 'FPRM', 'FCRS', 'CPRM', 'CCRS', 'NA', 'PCL', 'PH2O'),
    "camx6_3_chemparam_6_cf_hg": ('PNO3', 'PSO4', 'PNH4', 'POA', 'SOA1', 'SOA2', 'SOA3', 'SOA4', 'SOA5', 'SOA6', 'SOA7', 'SOAH', 'SOPA', 'SOPB', 'PEC', 'FPRM', 'FCRS', 'CPRM', 'CCRS', 'NA', 'PCL', 'PH2O', 'HGP', 'HGIIP', 'HGIIPC'),
    "camx6_3_chemparam_6_cf_vbs": ('PNO3', 'PSO4', 'PNH4', 'PAS0', 'PAS1', 'PAS2', 'PAS3', 'PAS4', 'PBS0', 'PBS1', 'PBS2', 'PBS3', 'PBS4', 'PAP0', 'PAP1', 'PAP2', 'PAP3', 'PAP4', 'PCP0', 'PCP1', 'PCP2', 'PCP3', 'PCP4', 'PFP0', 'PFP1', 'PFP2', 'PFP3', 'PFP4', 'PEC', 'FPRM', 'FCRS', 'CPRM', 'CCRS', 'NA', 'PCL', 'PH2O'),
    "camx6_3_chemparam_6_cmu": ('PNO3', 'PSO4', 'PNH4', 'POA', 'SOA1', 'SOA2', 'SOA3', 'SOA4', 'SOA5', 'SOA6', 'SOA7', 'SOAH', 'SOPA', 'SOPB', 'PEC', 'CRST', 'NA', 'PCL', 'PH2O'),
    "camx6_3_chemparam_6_inert": ('FCRS', 'CCRS'),
    "camx6_3_chemparam_6_none": (),
    }
aerosol_names = ('PEC', 'PNH4', 'PCP0', 'SOA1', 'PBS0', 'SOA5', 'FCRS', 'PH2O', 'SOPB', 'PBS3', 'PFP4', 'PCL', 'SOA3', 'HGIIPC', 'SOPA', 'PAS0', 'PFP1', 'PAP4', 'CPRM', 'SOA7', 'SOA6', 'NA', 'PAP1', 'PBS1', 'HGIIP', 'PCP4', 'CRST', 'SOAH', 'FPRM', 'SOA2', 'PAS4', 'PAS2', 'PAS3', 'PAS1', 'POA', 'PFP0', 'SOA4', 'HGP', 'PAP2', 'PCP3', 'CCRS', 'PCP1', 'PFP3', 'PBS2', 'PSO4', 'PNO3', 'PFP2', 'PCP2', 'PBS4', 'PAP3', 'PAP0')
