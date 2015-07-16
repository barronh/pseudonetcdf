#!/usr/bin/env python
import sys
import os
import re
import warnings

from collections import defaultdict
import textwrap
import json
from datetime import datetime, timedelta


import numpy as np
from numpy import exp, log, pi, zeros

from netCDF4 import Dataset
netcdf = Dataset
from PseudoNetCDF import getvarpnc, slice_dim, extract
from PseudoNetCDF.conventions.ioapi import add_cf_from_ioapi
from PseudoNetCDF.pncgen import pncgen
from PseudoNetCDF.geoschemfiles import bpch

def geticoncdl():
    return """netcdf DEFAULT_ICON {
dimensions:
	TSTEP = UNLIMITED ; // (1 currently)
	DATE-TIME = 2 ;
	VAR = 106 ;
	LAY = %(NLAYS)d ;
	ROW = %(NROWS)d ;
	COL = %(NCOLS)d ;
variables:
	int TFLAG(TSTEP, VAR, DATE-TIME) ;
		TFLAG:units = "<YYYYDDD,HHMMSS>" ;
		TFLAG:long_name = "TFLAG           " ;
		TFLAG:var_desc = "Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS                                " ;
	float NO2(TSTEP, LAY, ROW, COL) ;
		NO2:long_name = "NO2             " ;
		NO2:units = "ppmV            " ;
		NO2:coordinates = "lon lat" ;
		NO2:var_desc = "Variable NO2                                                                    " ;
	float NO(TSTEP, LAY, ROW, COL) ;
		NO:long_name = "NO              " ;
		NO:units = "ppmV            " ;
		NO:coordinates = "lon lat" ;
		NO:var_desc = "Variable NO                                                                     " ;
	float O3(TSTEP, LAY, ROW, COL) ;
		O3:long_name = "O3              " ;
		O3:units = "ppmV            " ;
		O3:coordinates = "lon lat" ;
		O3:var_desc = "Variable O3                                                                     " ;
	float NO3(TSTEP, LAY, ROW, COL) ;
		NO3:long_name = "NO3             " ;
		NO3:units = "ppmV            " ;
		NO3:coordinates = "lon lat" ;
		NO3:var_desc = "Variable NO3                                                                    " ;
	float OH(TSTEP, LAY, ROW, COL) ;
		OH:long_name = "OH              " ;
		OH:units = "ppmV            " ;
		OH:coordinates = "lon lat" ;
		OH:var_desc = "Variable OH                                                                     " ;
	float HO2(TSTEP, LAY, ROW, COL) ;
		HO2:long_name = "HO2             " ;
		HO2:units = "ppmV            " ;
		HO2:coordinates = "lon lat" ;
		HO2:var_desc = "Variable HO2                                                                    " ;
	float N2O5(TSTEP, LAY, ROW, COL) ;
		N2O5:long_name = "N2O5            " ;
		N2O5:units = "ppmV            " ;
		N2O5:coordinates = "lon lat" ;
		N2O5:var_desc = "Variable N2O5                                                                   " ;
	float HNO3(TSTEP, LAY, ROW, COL) ;
		HNO3:long_name = "HNO3            " ;
		HNO3:units = "ppmV            " ;
		HNO3:coordinates = "lon lat" ;
		HNO3:var_desc = "Variable HNO3                                                                   " ;
	float HONO(TSTEP, LAY, ROW, COL) ;
		HONO:long_name = "HONO            " ;
		HONO:units = "ppmV            " ;
		HONO:coordinates = "lon lat" ;
		HONO:var_desc = "Variable HONO                                                                   " ;
	float PNA(TSTEP, LAY, ROW, COL) ;
		PNA:long_name = "PNA             " ;
		PNA:units = "ppmV            " ;
		PNA:coordinates = "lon lat" ;
		PNA:var_desc = "Variable PNA                                                                    " ;
	float H2O2(TSTEP, LAY, ROW, COL) ;
		H2O2:long_name = "H2O2            " ;
		H2O2:units = "ppmV            " ;
		H2O2:coordinates = "lon lat" ;
		H2O2:var_desc = "Variable H2O2                                                                   " ;
	float NTR(TSTEP, LAY, ROW, COL) ;
		NTR:long_name = "NTR             " ;
		NTR:units = "ppmV            " ;
		NTR:coordinates = "lon lat" ;
		NTR:var_desc = "Variable NTR                                                                    " ;
	float ROOH(TSTEP, LAY, ROW, COL) ;
		ROOH:long_name = "ROOH            " ;
		ROOH:units = "ppmV            " ;
		ROOH:coordinates = "lon lat" ;
		ROOH:var_desc = "Variable ROOH                                                                   " ;
	float FORM(TSTEP, LAY, ROW, COL) ;
		FORM:long_name = "FORM            " ;
		FORM:units = "ppmV            " ;
		FORM:coordinates = "lon lat" ;
		FORM:var_desc = "Variable FORM                                                                   " ;
	float ALD2(TSTEP, LAY, ROW, COL) ;
		ALD2:long_name = "ALD2            " ;
		ALD2:units = "ppmV            " ;
		ALD2:coordinates = "lon lat" ;
		ALD2:var_desc = "Variable ALD2                                                                   " ;
	float PAR(TSTEP, LAY, ROW, COL) ;
		PAR:long_name = "PAR             " ;
		PAR:units = "ppmV            " ;
		PAR:coordinates = "lon lat" ;
		PAR:var_desc = "Variable PAR                                                                    " ;
	float CO(TSTEP, LAY, ROW, COL) ;
		CO:long_name = "CO              " ;
		CO:units = "ppmV            " ;
		CO:coordinates = "lon lat" ;
		CO:var_desc = "Variable CO                                                                     " ;
	float MEPX(TSTEP, LAY, ROW, COL) ;
		MEPX:long_name = "MEPX            " ;
		MEPX:units = "ppmV            " ;
		MEPX:coordinates = "lon lat" ;
		MEPX:var_desc = "Variable MEPX                                                                   " ;
	float FACD(TSTEP, LAY, ROW, COL) ;
		FACD:long_name = "FACD            " ;
		FACD:units = "ppmV            " ;
		FACD:coordinates = "lon lat" ;
		FACD:var_desc = "Variable FACD                                                                   " ;
	float C2O3(TSTEP, LAY, ROW, COL) ;
		C2O3:long_name = "C2O3            " ;
		C2O3:units = "ppmV            " ;
		C2O3:coordinates = "lon lat" ;
		C2O3:var_desc = "Variable C2O3                                                                   " ;
	float PAN(TSTEP, LAY, ROW, COL) ;
		PAN:long_name = "PAN             " ;
		PAN:units = "ppmV            " ;
		PAN:coordinates = "lon lat" ;
		PAN:var_desc = "Variable PAN                                                                    " ;
	float PACD(TSTEP, LAY, ROW, COL) ;
		PACD:long_name = "PACD            " ;
		PACD:units = "ppmV            " ;
		PACD:coordinates = "lon lat" ;
		PACD:var_desc = "Variable PACD                                                                   " ;
	float AACD(TSTEP, LAY, ROW, COL) ;
		AACD:long_name = "AACD            " ;
		AACD:units = "ppmV            " ;
		AACD:coordinates = "lon lat" ;
		AACD:var_desc = "Variable AACD                                                                   " ;
	float PANX(TSTEP, LAY, ROW, COL) ;
		PANX:long_name = "PANX            " ;
		PANX:units = "ppmV            " ;
		PANX:coordinates = "lon lat" ;
		PANX:var_desc = "Variable PANX                                                                   " ;
	float OLE(TSTEP, LAY, ROW, COL) ;
		OLE:long_name = "OLE             " ;
		OLE:units = "ppmV            " ;
		OLE:coordinates = "lon lat" ;
		OLE:var_desc = "Variable OLE                                                                    " ;
	float ETH(TSTEP, LAY, ROW, COL) ;
		ETH:long_name = "ETH             " ;
		ETH:units = "ppmV            " ;
		ETH:coordinates = "lon lat" ;
		ETH:var_desc = "Variable ETH                                                                    " ;
	float IOLE(TSTEP, LAY, ROW, COL) ;
		IOLE:long_name = "IOLE            " ;
		IOLE:units = "ppmV            " ;
		IOLE:coordinates = "lon lat" ;
		IOLE:var_desc = "Variable IOLE                                                                   " ;
	float TOL(TSTEP, LAY, ROW, COL) ;
		TOL:long_name = "TOL             " ;
		TOL:units = "ppmV            " ;
		TOL:coordinates = "lon lat" ;
		TOL:var_desc = "Variable TOL                                                                    " ;
	float CRES(TSTEP, LAY, ROW, COL) ;
		CRES:long_name = "CRES            " ;
		CRES:units = "ppmV            " ;
		CRES:coordinates = "lon lat" ;
		CRES:var_desc = "Variable CRES                                                                   " ;
	float OPEN(TSTEP, LAY, ROW, COL) ;
		OPEN:long_name = "OPEN            " ;
		OPEN:units = "ppmV            " ;
		OPEN:coordinates = "lon lat" ;
		OPEN:var_desc = "Variable OPEN                                                                   " ;
	float MGLY(TSTEP, LAY, ROW, COL) ;
		MGLY:long_name = "MGLY            " ;
		MGLY:units = "ppmV            " ;
		MGLY:coordinates = "lon lat" ;
		MGLY:var_desc = "Variable MGLY                                                                   " ;
	float XYL(TSTEP, LAY, ROW, COL) ;
		XYL:long_name = "XYL             " ;
		XYL:units = "ppmV            " ;
		XYL:coordinates = "lon lat" ;
		XYL:var_desc = "Variable XYL                                                                    " ;
	float ISOP(TSTEP, LAY, ROW, COL) ;
		ISOP:long_name = "ISOP            " ;
		ISOP:units = "ppmV            " ;
		ISOP:coordinates = "lon lat" ;
		ISOP:var_desc = "Variable ISOP                                                                   " ;
	float SO2(TSTEP, LAY, ROW, COL) ;
		SO2:long_name = "SO2             " ;
		SO2:units = "ppmV            " ;
		SO2:coordinates = "lon lat" ;
		SO2:var_desc = "Variable SO2                                                                    " ;
	float SULF(TSTEP, LAY, ROW, COL) ;
		SULF:long_name = "SULF            " ;
		SULF:units = "ppmV            " ;
		SULF:coordinates = "lon lat" ;
		SULF:var_desc = "Variable SULF                                                                   " ;
	float ETHA(TSTEP, LAY, ROW, COL) ;
		ETHA:long_name = "ETHA            " ;
		ETHA:units = "ppmV            " ;
		ETHA:coordinates = "lon lat" ;
		ETHA:var_desc = "Variable ETHA                                                                   " ;
	float BENZENE(TSTEP, LAY, ROW, COL) ;
		BENZENE:long_name = "BENZENE         " ;
		BENZENE:units = "ppmV            " ;
		BENZENE:coordinates = "lon lat" ;
		BENZENE:var_desc = "Variable BENZENE                                                                " ;
	float ASO4J(TSTEP, LAY, ROW, COL) ;
		ASO4J:long_name = "ASO4J           " ;
		ASO4J:units = "micrograms/m**3 " ;
		ASO4J:coordinates = "lon lat" ;
		ASO4J:var_desc = "Variable ASO4J                                                                  " ;
	float ASO4I(TSTEP, LAY, ROW, COL) ;
		ASO4I:long_name = "ASO4I           " ;
		ASO4I:units = "micrograms/m**3 " ;
		ASO4I:coordinates = "lon lat" ;
		ASO4I:var_desc = "Variable ASO4I                                                                  " ;
	float AALKJ(TSTEP, LAY, ROW, COL) ;
		AALKJ:long_name = "AALKJ           " ;
		AALKJ:units = "micrograms/m**3 " ;
		AALKJ:coordinates = "lon lat" ;
		AALKJ:var_desc = "Variable AALKJ                                                                  " ;
	float AXYL1J(TSTEP, LAY, ROW, COL) ;
		AXYL1J:long_name = "AXYL1J          " ;
		AXYL1J:units = "micrograms/m**3 " ;
		AXYL1J:coordinates = "lon lat" ;
		AXYL1J:var_desc = "Variable AXYL1J                                                                 " ;
	float AXYL2J(TSTEP, LAY, ROW, COL) ;
		AXYL2J:long_name = "AXYL2J          " ;
		AXYL2J:units = "micrograms/m**3 " ;
		AXYL2J:coordinates = "lon lat" ;
		AXYL2J:var_desc = "Variable AXYL2J                                                                 " ;
	float AXYL3J(TSTEP, LAY, ROW, COL) ;
		AXYL3J:long_name = "AXYL3J          " ;
		AXYL3J:units = "micrograms/m**3 " ;
		AXYL3J:coordinates = "lon lat" ;
		AXYL3J:var_desc = "Variable AXYL3J                                                                 " ;
	float ATOL1J(TSTEP, LAY, ROW, COL) ;
		ATOL1J:long_name = "ATOL1J          " ;
		ATOL1J:units = "micrograms/m**3 " ;
		ATOL1J:coordinates = "lon lat" ;
		ATOL1J:var_desc = "Variable ATOL1J                                                                 " ;
	float ATOL2J(TSTEP, LAY, ROW, COL) ;
		ATOL2J:long_name = "ATOL2J          " ;
		ATOL2J:units = "micrograms/m**3 " ;
		ATOL2J:coordinates = "lon lat" ;
		ATOL2J:var_desc = "Variable ATOL2J                                                                 " ;
	float ATOL3J(TSTEP, LAY, ROW, COL) ;
		ATOL3J:long_name = "ATOL3J          " ;
		ATOL3J:units = "micrograms/m**3 " ;
		ATOL3J:coordinates = "lon lat" ;
		ATOL3J:var_desc = "Variable ATOL3J                                                                 " ;
	float ABNZ1J(TSTEP, LAY, ROW, COL) ;
		ABNZ1J:long_name = "ABNZ1J          " ;
		ABNZ1J:units = "micrograms/m**3 " ;
		ABNZ1J:coordinates = "lon lat" ;
		ABNZ1J:var_desc = "Variable ABNZ1J                                                                 " ;
	float ABNZ2J(TSTEP, LAY, ROW, COL) ;
		ABNZ2J:long_name = "ABNZ2J          " ;
		ABNZ2J:units = "micrograms/m**3 " ;
		ABNZ2J:coordinates = "lon lat" ;
		ABNZ2J:var_desc = "Variable ABNZ2J                                                                 " ;
	float ABNZ3J(TSTEP, LAY, ROW, COL) ;
		ABNZ3J:long_name = "ABNZ3J          " ;
		ABNZ3J:units = "micrograms/m**3 " ;
		ABNZ3J:coordinates = "lon lat" ;
		ABNZ3J:var_desc = "Variable ABNZ3J                                                                 " ;
	float ATRP1J(TSTEP, LAY, ROW, COL) ;
		ATRP1J:long_name = "ATRP1J          " ;
		ATRP1J:units = "micrograms/m**3 " ;
		ATRP1J:coordinates = "lon lat" ;
		ATRP1J:var_desc = "Variable ATRP1J                                                                 " ;
	float ATRP2J(TSTEP, LAY, ROW, COL) ;
		ATRP2J:long_name = "ATRP2J          " ;
		ATRP2J:units = "micrograms/m**3 " ;
		ATRP2J:coordinates = "lon lat" ;
		ATRP2J:var_desc = "Variable ATRP2J                                                                 " ;
	float AISO1J(TSTEP, LAY, ROW, COL) ;
		AISO1J:long_name = "AISO1J          " ;
		AISO1J:units = "micrograms/m**3 " ;
		AISO1J:coordinates = "lon lat" ;
		AISO1J:var_desc = "Variable AISO1J                                                                 " ;
	float AISO2J(TSTEP, LAY, ROW, COL) ;
		AISO2J:long_name = "AISO2J          " ;
		AISO2J:units = "micrograms/m**3 " ;
		AISO2J:coordinates = "lon lat" ;
		AISO2J:var_desc = "Variable AISO2J                                                                 " ;
	float ASQTJ(TSTEP, LAY, ROW, COL) ;
		ASQTJ:long_name = "ASQTJ           " ;
		ASQTJ:units = "micrograms/m**3 " ;
		ASQTJ:coordinates = "lon lat" ;
		ASQTJ:var_desc = "Variable ASQTJ                                                                  " ;
	float ACORS(TSTEP, LAY, ROW, COL) ;
		ACORS:long_name = "ACORS           " ;
		ACORS:units = "micrograms/m**3 " ;
		ACORS:coordinates = "lon lat" ;
		ACORS:var_desc = "Variable ACORS                                                                  " ;
	float ASOIL(TSTEP, LAY, ROW, COL) ;
		ASOIL:long_name = "ASOIL           " ;
		ASOIL:units = "micrograms/m**3 " ;
		ASOIL:coordinates = "lon lat" ;
		ASOIL:var_desc = "Variable ASOIL                                                                  " ;
	float NUMATKN(TSTEP, LAY, ROW, COL) ;
		NUMATKN:long_name = "NUMATKN         " ;
		NUMATKN:units = "#/m**3          " ;
		NUMATKN:coordinates = "lon lat" ;
		NUMATKN:var_desc = "Variable NUMATKN                                                                " ;
	float NUMACC(TSTEP, LAY, ROW, COL) ;
		NUMACC:long_name = "NUMACC          " ;
		NUMACC:units = "#/m**3          " ;
		NUMACC:coordinates = "lon lat" ;
		NUMACC:var_desc = "Variable NUMACC                                                                 " ;
	float NUMCOR(TSTEP, LAY, ROW, COL) ;
		NUMCOR:long_name = "NUMCOR          " ;
		NUMCOR:units = "#/m**3          " ;
		NUMCOR:coordinates = "lon lat" ;
		NUMCOR:var_desc = "Variable NUMCOR                                                                 " ;
	float SRFATKN(TSTEP, LAY, ROW, COL) ;
		SRFATKN:long_name = "SRFATKN         " ;
		SRFATKN:units = "m**2/m**3       " ;
		SRFATKN:coordinates = "lon lat" ;
		SRFATKN:var_desc = "Variable SRFATKN                                                                " ;
	float SRFACC(TSTEP, LAY, ROW, COL) ;
		SRFACC:long_name = "SRFACC          " ;
		SRFACC:units = "m**2/m**3       " ;
		SRFACC:coordinates = "lon lat" ;
		SRFACC:var_desc = "Variable SRFACC                                                                 " ;
	float SRFCOR(TSTEP, LAY, ROW, COL) ;
		SRFCOR:long_name = "SRFCOR          " ;
		SRFCOR:units = "m**2/m**3       " ;
		SRFCOR:coordinates = "lon lat" ;
		SRFCOR:var_desc = "Variable SRFCOR                                                                 " ;
	float AISO3J(TSTEP, LAY, ROW, COL) ;
		AISO3J:long_name = "AISO3J          " ;
		AISO3J:units = "micrograms/m**3 " ;
		AISO3J:coordinates = "lon lat" ;
		AISO3J:var_desc = "Variable AISO3J                                                                 " ;
	float AOLGAJ(TSTEP, LAY, ROW, COL) ;
		AOLGAJ:long_name = "AOLGAJ          " ;
		AOLGAJ:units = "micrograms/m**3 " ;
		AOLGAJ:coordinates = "lon lat" ;
		AOLGAJ:var_desc = "Variable AOLGAJ                                                                 " ;
	float AOLGBJ(TSTEP, LAY, ROW, COL) ;
		AOLGBJ:long_name = "AOLGBJ          " ;
		AOLGBJ:units = "micrograms/m**3 " ;
		AOLGBJ:coordinates = "lon lat" ;
		AOLGBJ:var_desc = "Variable AOLGBJ                                                                 " ;
	float NH3(TSTEP, LAY, ROW, COL) ;
		NH3:long_name = "NH3             " ;
		NH3:units = "ppmV            " ;
		NH3:coordinates = "lon lat" ;
		NH3:var_desc = "Variable NH3                                                                    " ;
	float SV_ALK(TSTEP, LAY, ROW, COL) ;
		SV_ALK:long_name = "SV_ALK          " ;
		SV_ALK:units = "ppmV            " ;
		SV_ALK:coordinates = "lon lat" ;
		SV_ALK:var_desc = "Variable SV_ALK                                                                 " ;
	float SV_XYL1(TSTEP, LAY, ROW, COL) ;
		SV_XYL1:long_name = "SV_XYL1         " ;
		SV_XYL1:units = "ppmV            " ;
		SV_XYL1:coordinates = "lon lat" ;
		SV_XYL1:var_desc = "Variable SV_XYL1                                                                " ;
	float SV_XYL2(TSTEP, LAY, ROW, COL) ;
		SV_XYL2:long_name = "SV_XYL2         " ;
		SV_XYL2:units = "ppmV            " ;
		SV_XYL2:coordinates = "lon lat" ;
		SV_XYL2:var_desc = "Variable SV_XYL2                                                                " ;
	float SV_TOL1(TSTEP, LAY, ROW, COL) ;
		SV_TOL1:long_name = "SV_TOL1         " ;
		SV_TOL1:units = "ppmV            " ;
		SV_TOL1:coordinates = "lon lat" ;
		SV_TOL1:var_desc = "Variable SV_TOL1                                                                " ;
	float SV_TOL2(TSTEP, LAY, ROW, COL) ;
		SV_TOL2:long_name = "SV_TOL2         " ;
		SV_TOL2:units = "ppmV            " ;
		SV_TOL2:coordinates = "lon lat" ;
		SV_TOL2:var_desc = "Variable SV_TOL2                                                                " ;
	float SV_BNZ1(TSTEP, LAY, ROW, COL) ;
		SV_BNZ1:long_name = "SV_BNZ1         " ;
		SV_BNZ1:units = "ppmV            " ;
		SV_BNZ1:coordinates = "lon lat" ;
		SV_BNZ1:var_desc = "Variable SV_BNZ1                                                                " ;
	float SV_BNZ2(TSTEP, LAY, ROW, COL) ;
		SV_BNZ2:long_name = "SV_BNZ2         " ;
		SV_BNZ2:units = "ppmV            " ;
		SV_BNZ2:coordinates = "lon lat" ;
		SV_BNZ2:var_desc = "Variable SV_BNZ2                                                                " ;
	float SV_TRP1(TSTEP, LAY, ROW, COL) ;
		SV_TRP1:long_name = "SV_TRP1         " ;
		SV_TRP1:units = "ppmV            " ;
		SV_TRP1:coordinates = "lon lat" ;
		SV_TRP1:var_desc = "Variable SV_TRP1                                                                " ;
	float SV_TRP2(TSTEP, LAY, ROW, COL) ;
		SV_TRP2:long_name = "SV_TRP2         " ;
		SV_TRP2:units = "ppmV            " ;
		SV_TRP2:coordinates = "lon lat" ;
		SV_TRP2:var_desc = "Variable SV_TRP2                                                                " ;
	float SV_ISO1(TSTEP, LAY, ROW, COL) ;
		SV_ISO1:long_name = "SV_ISO1         " ;
		SV_ISO1:units = "ppmV            " ;
		SV_ISO1:coordinates = "lon lat" ;
		SV_ISO1:var_desc = "Variable SV_ISO1                                                                " ;
	float SV_ISO2(TSTEP, LAY, ROW, COL) ;
		SV_ISO2:long_name = "SV_ISO2         " ;
		SV_ISO2:units = "ppmV            " ;
		SV_ISO2:coordinates = "lon lat" ;
		SV_ISO2:var_desc = "Variable SV_ISO2                                                                " ;
	float SV_SQT(TSTEP, LAY, ROW, COL) ;
		SV_SQT:long_name = "SV_SQT          " ;
		SV_SQT:units = "ppmV            " ;
		SV_SQT:coordinates = "lon lat" ;
		SV_SQT:var_desc = "Variable SV_SQT                                                                 " ;
	float AMGJ(TSTEP, LAY, ROW, COL) ;
		AMGJ:long_name = "AMGJ            " ;
		AMGJ:units = "micrograms/m**3 " ;
		AMGJ:coordinates = "lon lat" ;
		AMGJ:var_desc = "AMGJ            " ;
	float ACAJ(TSTEP, LAY, ROW, COL) ;
		ACAJ:long_name = "ACAJ            " ;
		ACAJ:units = "micrograms/m**3 " ;
		ACAJ:coordinates = "lon lat" ;
		ACAJ:var_desc = "ACAJ            " ;
	float ANH4J(TSTEP, LAY, ROW, COL) ;
		ANH4J:long_name = "ANH4J           " ;
		ANH4J:units = "micrograms/m**3 " ;
		ANH4J:coordinates = "lon lat" ;
		ANH4J:var_desc = "ANH4J           " ;
	float AFEJ(TSTEP, LAY, ROW, COL) ;
		AFEJ:long_name = "AFEJ            " ;
		AFEJ:units = "micrograms/m**3 " ;
		AFEJ:coordinates = "lon lat" ;
		AFEJ:var_desc = "AFEJ            " ;
	float ACLK(TSTEP, LAY, ROW, COL) ;
		ACLK:long_name = "ACLK            " ;
		ACLK:units = "micrograms/m**3 " ;
		ACLK:coordinates = "lon lat" ;
		ACLK:var_desc = "ACLK            " ;
	float APNCOMI(TSTEP, LAY, ROW, COL) ;
		APNCOMI:long_name = "APNCOMI         " ;
		APNCOMI:units = "micrograms/m**3 " ;
		APNCOMI:coordinates = "lon lat" ;
		APNCOMI:var_desc = "APNCOMI         " ;
	float AMNJ(TSTEP, LAY, ROW, COL) ;
		AMNJ:long_name = "AMNJ            " ;
		AMNJ:units = "micrograms/m**3 " ;
		AMNJ:coordinates = "lon lat" ;
		AMNJ:var_desc = "AMNJ            " ;
	float APNCOMJ(TSTEP, LAY, ROW, COL) ;
		APNCOMJ:long_name = "APNCOMJ         " ;
		APNCOMJ:units = "micrograms/m**3 " ;
		APNCOMJ:coordinates = "lon lat" ;
		APNCOMJ:var_desc = "APNCOMJ         " ;
	float AIRDEN(TSTEP, LAY, ROW, COL) ;
		AIRDEN:long_name = "AIRDEN          " ;
		AIRDEN:units = "molec/cm3       " ;
		AIRDEN:coordinates = "lon lat" ;
		AIRDEN:var_desc = "AIRDEN          " ;
	float ANO3I(TSTEP, LAY, ROW, COL) ;
		ANO3I:long_name = "ANO3I           " ;
		ANO3I:units = "micrograms/m**3 " ;
		ANO3I:coordinates = "lon lat" ;
		ANO3I:var_desc = "ANO3I           " ;
	float ANO3J(TSTEP, LAY, ROW, COL) ;
		ANO3J:long_name = "ANO3J           " ;
		ANO3J:units = "micrograms/m**3 " ;
		ANO3J:coordinates = "lon lat" ;
		ANO3J:var_desc = "ANO3J           " ;
	float ANO3K(TSTEP, LAY, ROW, COL) ;
		ANO3K:long_name = "ANO3K           " ;
		ANO3K:units = "micrograms/m**3 " ;
		ANO3K:coordinates = "lon lat" ;
		ANO3K:var_desc = "ANO3K           " ;
	float ASIJ(TSTEP, LAY, ROW, COL) ;
		ASIJ:long_name = "ASIJ            " ;
		ASIJ:units = "micrograms/m**3 " ;
		ASIJ:coordinates = "lon lat" ;
		ASIJ:var_desc = "ASIJ            " ;
	float AECJ(TSTEP, LAY, ROW, COL) ;
		AECJ:long_name = "AECJ            " ;
		AECJ:units = "micrograms/m**3 " ;
		AECJ:coordinates = "lon lat" ;
		AECJ:var_desc = "AECJ            " ;
	float AOTHRJ(TSTEP, LAY, ROW, COL) ;
		AOTHRJ:long_name = "AOTHRJ          " ;
		AOTHRJ:units = "micrograms/m**3 " ;
		AOTHRJ:coordinates = "lon lat" ;
		AOTHRJ:var_desc = "AOTHRJ          " ;
	float ANH4I(TSTEP, LAY, ROW, COL) ;
		ANH4I:long_name = "ANH4I           " ;
		ANH4I:units = "micrograms/m**3 " ;
		ANH4I:coordinates = "lon lat" ;
		ANH4I:var_desc = "ANH4I           " ;
	float ATIJ(TSTEP, LAY, ROW, COL) ;
		ATIJ:long_name = "ATIJ            " ;
		ATIJ:units = "micrograms/m**3 " ;
		ATIJ:coordinates = "lon lat" ;
		ATIJ:var_desc = "ATIJ            " ;
	float AKJ(TSTEP, LAY, ROW, COL) ;
		AKJ:long_name = "AKJ             " ;
		AKJ:units = "micrograms/m**3 " ;
		AKJ:coordinates = "lon lat" ;
		AKJ:var_desc = "AKJ             " ;
	float ASEACAT(TSTEP, LAY, ROW, COL) ;
		ASEACAT:long_name = "ASEACAT         " ;
		ASEACAT:units = "micrograms/m**3 " ;
		ASEACAT:coordinates = "lon lat" ;
		ASEACAT:var_desc = "ASEACAT         " ;
	float ASO4K(TSTEP, LAY, ROW, COL) ;
		ASO4K:long_name = "ASO4K           " ;
		ASO4K:units = "micrograms/m**3 " ;
		ASO4K:coordinates = "lon lat" ;
		ASO4K:var_desc = "ASO4K           " ;
	float ISPD(TSTEP, LAY, ROW, COL) ;
		ISPD:long_name = "ISPD            " ;
		ISPD:units = "ppmV            " ;
		ISPD:coordinates = "lon lat" ;
		ISPD:var_desc = "ISPD            " ;
	float ALDX(TSTEP, LAY, ROW, COL) ;
		ALDX:long_name = "ALDX            " ;
		ALDX:units = "ppmV            " ;
		ALDX:coordinates = "lon lat" ;
		ALDX:var_desc = "ALDX            " ;
	float AALJ(TSTEP, LAY, ROW, COL) ;
		AALJ:long_name = "AALJ            " ;
		AALJ:units = "micrograms/m**3 " ;
		AALJ:coordinates = "lon lat" ;
		AALJ:var_desc = "AALJ            " ;
	float ANAJ(TSTEP, LAY, ROW, COL) ;
		ANAJ:long_name = "ANAJ            " ;
		ANAJ:units = "micrograms/m**3 " ;
		ANAJ:coordinates = "lon lat" ;
		ANAJ:var_desc = "ANAJ            " ;
	float APOCJ(TSTEP, LAY, ROW, COL) ;
		APOCJ:long_name = "APOCJ           " ;
		APOCJ:units = "micrograms/m**3 " ;
		APOCJ:coordinates = "lon lat" ;
		APOCJ:var_desc = "APOCJ           " ;
	float APOCI(TSTEP, LAY, ROW, COL) ;
		APOCI:long_name = "APOCI           " ;
		APOCI:units = "micrograms/m**3 " ;
		APOCI:coordinates = "lon lat" ;
		APOCI:var_desc = "APOCI           " ;
	float AECI(TSTEP, LAY, ROW, COL) ;
		AECI:long_name = "AECI            " ;
		AECI:units = "micrograms/m**3 " ;
		AECI:coordinates = "lon lat" ;
		AECI:var_desc = "AECI            " ;
	float ACLJ(TSTEP, LAY, ROW, COL) ;
		ACLJ:long_name = "ACLJ            " ;
		ACLJ:units = "micrograms/m**3 " ;
		ACLJ:coordinates = "lon lat" ;
		ACLJ:var_desc = "ACLJ            " ;

// global attributes:
		:IOAPI_VERSION = "$Id: @(#) ioapi library version 3.1 $                                           " ;
		:EXEC_ID = "BCON_V5g_Darwin13_x86_64gfortran                                                " ;
		:FTYPE = 2 ;
		:CDATE = 2014031 ;
		:CTIME = 230334 ;
		:WDATE = 2014031 ;
		:WTIME = 230334 ;
		:SDATE = %(SDATE)d ;
		:STIME = %(STIME)d ;
		:TSTEP = %(TSTEP)d ;
		:NTHIK = 1 ;
		:NCOLS = %(NCOLS)d ;
		:NROWS = %(NROWS)d ;
		:NLAYS = %(NLAYS)d ;
		:NVARS = 106 ;
		:GDTYP = %(GDTYP)d ;
		:P_ALP = %(P_ALP)f ;
		:P_BET = %(P_BET)f ;
		:P_GAM = %(P_GAM)f ;
		:XCENT = %(XCENT)f ;
		:YCENT = %(YCENT)f ;
		:XORIG = %(XORIG)f ;
		:YORIG = %(YORIG)f ;
		:XCELL = %(XCELL)f ;
		:YCELL = %(YCELL)f ;
		:VGTYP = 7 ;
		:VGTOP = %(VGTOP)f ;
		:VGLVLS = %(VGLVLS)s ;
		:GDNAM = "%(GDNAM)16s" ;
		:UPNAM = "%(UPNAM)16s" ;
		:VAR-LIST = "NO2             NO              O3              NO3             OH              HO2             N2O5            HNO3            HONO            PNA             H2O2            NTR             ROOH            FORM            ALD2            PAR             CO              MEPX            FACD            C2O3            PAN             PACD            AACD            PANX            OLE             ETH             IOLE            TOL             CRES            OPEN            MGLY            XYL             ISOP            SO2             SULF            ETHA            BENZENE         ASO4J           ASO4I           AALKJ           AXYL1J          AXYL2J          AXYL3J          ATOL1J          ATOL2J          ATOL3J          ABNZ1J          ABNZ2J          ABNZ3J          ATRP1J          ATRP2J          AISO1J          AISO2J          ASQTJ           ACORS           ASOIL           NUMATKN         NUMACC          NUMCOR          SRFATKN         SRFACC          SRFCOR          AISO3J          AOLGAJ          AOLGBJ          NH3             SV_ALK          SV_XYL1         SV_XYL2         SV_TOL1         SV_TOL2         SV_BNZ1         SV_BNZ2         SV_TRP1         SV_TRP2         SV_ISO1         SV_ISO2         SV_SQT          AMGJ            ACAJ            ANH4J           AFEJ            ACLK            APNCOMI         AMNJ            APNCOMJ         AIRDEN          ANO3I           ANO3J           ANO3K           ASIJ            AECJ            AOTHRJ          ANH4I           ATIJ            AKJ             ASEACAT         ASO4K           ISPD            ALDX            AALJ            ANAJ            APOCJ           APOCI           AECI            ACLJ            " ;
		:FILEDESC = "BCON output file BNDY_CONC_1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    " ;
		:HISTORY = "" ;
}
"""
def getbconcdl():
    return """netcdf BCON_V5g_CMAQ-BENCHMARK_profile {
dimensions:
	TSTEP = UNLIMITED ;
	DATE-TIME = 2 ;
	VAR = 78 ;
	LAY = %(NLAYS)d ;
	PERIM = %(PERIM)d ;
variables:
    int TFLAG(TSTEP, VAR, DATE-TIME) ;
            TFLAG:units = "<YYYYDDD,HHMMSS>" ;
            TFLAG:long_name = "TFLAG           " ;
            TFLAG:var_desc = "Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS                                " ;
	float NO2(TSTEP, LAY, PERIM) ;
		NO2:long_name = "NO2             " ;
		NO2:units = "ppmV            " ;
		NO2:var_desc = "Variable NO2                                                                    " ;
	float NO(TSTEP, LAY, PERIM) ;
		NO:long_name = "NO              " ;
		NO:units = "ppmV            " ;
		NO:var_desc = "Variable NO                                                                     " ;
	float O3(TSTEP, LAY, PERIM) ;
		O3:long_name = "O3              " ;
		O3:units = "ppmV            " ;
		O3:var_desc = "Variable O3                                                                     " ;
	float NO3(TSTEP, LAY, PERIM) ;
		NO3:long_name = "NO3             " ;
		NO3:units = "ppmV            " ;
		NO3:var_desc = "Variable NO3                                                                    " ;
	float OH(TSTEP, LAY, PERIM) ;
		OH:long_name = "OH              " ;
		OH:units = "ppmV            " ;
		OH:var_desc = "Variable OH                                                                     " ;
	float HO2(TSTEP, LAY, PERIM) ;
		HO2:long_name = "HO2             " ;
		HO2:units = "ppmV            " ;
		HO2:var_desc = "Variable HO2                                                                    " ;
	float N2O5(TSTEP, LAY, PERIM) ;
		N2O5:long_name = "N2O5            " ;
		N2O5:units = "ppmV            " ;
		N2O5:var_desc = "Variable N2O5                                                                   " ;
	float HNO3(TSTEP, LAY, PERIM) ;
		HNO3:long_name = "HNO3            " ;
		HNO3:units = "ppmV            " ;
		HNO3:var_desc = "Variable HNO3                                                                   " ;
	float HONO(TSTEP, LAY, PERIM) ;
		HONO:long_name = "HONO            " ;
		HONO:units = "ppmV            " ;
		HONO:var_desc = "Variable HONO                                                                   " ;
	float PNA(TSTEP, LAY, PERIM) ;
		PNA:long_name = "PNA             " ;
		PNA:units = "ppmV            " ;
		PNA:var_desc = "Variable PNA                                                                    " ;
	float H2O2(TSTEP, LAY, PERIM) ;
		H2O2:long_name = "H2O2            " ;
		H2O2:units = "ppmV            " ;
		H2O2:var_desc = "Variable H2O2                                                                   " ;
	float NTR(TSTEP, LAY, PERIM) ;
		NTR:long_name = "NTR             " ;
		NTR:units = "ppmV            " ;
		NTR:var_desc = "Variable NTR                                                                    " ;
	float ROOH(TSTEP, LAY, PERIM) ;
		ROOH:long_name = "ROOH            " ;
		ROOH:units = "ppmV            " ;
		ROOH:var_desc = "Variable ROOH                                                                   " ;
	float FORM(TSTEP, LAY, PERIM) ;
		FORM:long_name = "FORM            " ;
		FORM:units = "ppmV            " ;
		FORM:var_desc = "Variable FORM                                                                   " ;
	float ALD2(TSTEP, LAY, PERIM) ;
		ALD2:long_name = "ALD2            " ;
		ALD2:units = "ppmV            " ;
		ALD2:var_desc = "Variable ALD2                                                                   " ;
	float PAR(TSTEP, LAY, PERIM) ;
		PAR:long_name = "PAR             " ;
		PAR:units = "ppmV            " ;
		PAR:var_desc = "Variable PAR                                                                    " ;
	float CO(TSTEP, LAY, PERIM) ;
		CO:long_name = "CO              " ;
		CO:units = "ppmV            " ;
		CO:var_desc = "Variable CO                                                                     " ;
	float MEPX(TSTEP, LAY, PERIM) ;
		MEPX:long_name = "MEPX            " ;
		MEPX:units = "ppmV            " ;
		MEPX:var_desc = "Variable MEPX                                                                   " ;
	float FACD(TSTEP, LAY, PERIM) ;
		FACD:long_name = "FACD            " ;
		FACD:units = "ppmV            " ;
		FACD:var_desc = "Variable FACD                                                                   " ;
	float C2O3(TSTEP, LAY, PERIM) ;
		C2O3:long_name = "C2O3            " ;
		C2O3:units = "ppmV            " ;
		C2O3:var_desc = "Variable C2O3                                                                   " ;
	float PAN(TSTEP, LAY, PERIM) ;
		PAN:long_name = "PAN             " ;
		PAN:units = "ppmV            " ;
		PAN:var_desc = "Variable PAN                                                                    " ;
	float PACD(TSTEP, LAY, PERIM) ;
		PACD:long_name = "PACD            " ;
		PACD:units = "ppmV            " ;
		PACD:var_desc = "Variable PACD                                                                   " ;
	float AACD(TSTEP, LAY, PERIM) ;
		AACD:long_name = "AACD            " ;
		AACD:units = "ppmV            " ;
		AACD:var_desc = "Variable AACD                                                                   " ;
	float PANX(TSTEP, LAY, PERIM) ;
		PANX:long_name = "PANX            " ;
		PANX:units = "ppmV            " ;
		PANX:var_desc = "Variable PANX                                                                   " ;
	float OLE(TSTEP, LAY, PERIM) ;
		OLE:long_name = "OLE             " ;
		OLE:units = "ppmV            " ;
		OLE:var_desc = "Variable OLE                                                                    " ;
	float ETH(TSTEP, LAY, PERIM) ;
		ETH:long_name = "ETH             " ;
		ETH:units = "ppmV            " ;
		ETH:var_desc = "Variable ETH                                                                    " ;
	float IOLE(TSTEP, LAY, PERIM) ;
		IOLE:long_name = "IOLE            " ;
		IOLE:units = "ppmV            " ;
		IOLE:var_desc = "Variable IOLE                                                                   " ;
	float TOL(TSTEP, LAY, PERIM) ;
		TOL:long_name = "TOL             " ;
		TOL:units = "ppmV            " ;
		TOL:var_desc = "Variable TOL                                                                    " ;
	float CRES(TSTEP, LAY, PERIM) ;
		CRES:long_name = "CRES            " ;
		CRES:units = "ppmV            " ;
		CRES:var_desc = "Variable CRES                                                                   " ;
	float OPEN(TSTEP, LAY, PERIM) ;
		OPEN:long_name = "OPEN            " ;
		OPEN:units = "ppmV            " ;
		OPEN:var_desc = "Variable OPEN                                                                   " ;
	float MGLY(TSTEP, LAY, PERIM) ;
		MGLY:long_name = "MGLY            " ;
		MGLY:units = "ppmV            " ;
		MGLY:var_desc = "Variable MGLY                                                                   " ;
	float XYL(TSTEP, LAY, PERIM) ;
		XYL:long_name = "XYL             " ;
		XYL:units = "ppmV            " ;
		XYL:var_desc = "Variable XYL                                                                    " ;
	float ISOP(TSTEP, LAY, PERIM) ;
		ISOP:long_name = "ISOP            " ;
		ISOP:units = "ppmV            " ;
		ISOP:var_desc = "Variable ISOP                                                                   " ;
	float SO2(TSTEP, LAY, PERIM) ;
		SO2:long_name = "SO2             " ;
		SO2:units = "ppmV            " ;
		SO2:var_desc = "Variable SO2                                                                    " ;
	float SULF(TSTEP, LAY, PERIM) ;
		SULF:long_name = "SULF            " ;
		SULF:units = "ppmV            " ;
		SULF:var_desc = "Variable SULF                                                                   " ;
	float ETHA(TSTEP, LAY, PERIM) ;
		ETHA:long_name = "ETHA            " ;
		ETHA:units = "ppmV            " ;
		ETHA:var_desc = "Variable ETHA                                                                   " ;
	float BENZENE(TSTEP, LAY, PERIM) ;
		BENZENE:long_name = "BENZENE         " ;
		BENZENE:units = "ppmV            " ;
		BENZENE:var_desc = "Variable BENZENE                                                                " ;
	float ASO4J(TSTEP, LAY, PERIM) ;
		ASO4J:long_name = "ASO4J           " ;
		ASO4J:units = "micrograms/m**3 " ;
		ASO4J:var_desc = "Variable ASO4J                                                                  " ;
	float ASO4I(TSTEP, LAY, PERIM) ;
		ASO4I:long_name = "ASO4I           " ;
		ASO4I:units = "micrograms/m**3 " ;
		ASO4I:var_desc = "Variable ASO4I                                                                  " ;
	float AALKJ(TSTEP, LAY, PERIM) ;
		AALKJ:long_name = "AALKJ           " ;
		AALKJ:units = "micrograms/m**3 " ;
		AALKJ:var_desc = "Variable AALKJ                                                                  " ;
	float AXYL1J(TSTEP, LAY, PERIM) ;
		AXYL1J:long_name = "AXYL1J          " ;
		AXYL1J:units = "micrograms/m**3 " ;
		AXYL1J:var_desc = "Variable AXYL1J                                                                 " ;
	float AXYL2J(TSTEP, LAY, PERIM) ;
		AXYL2J:long_name = "AXYL2J          " ;
		AXYL2J:units = "micrograms/m**3 " ;
		AXYL2J:var_desc = "Variable AXYL2J                                                                 " ;
	float AXYL3J(TSTEP, LAY, PERIM) ;
		AXYL3J:long_name = "AXYL3J          " ;
		AXYL3J:units = "micrograms/m**3 " ;
		AXYL3J:var_desc = "Variable AXYL3J                                                                 " ;
	float ATOL1J(TSTEP, LAY, PERIM) ;
		ATOL1J:long_name = "ATOL1J          " ;
		ATOL1J:units = "micrograms/m**3 " ;
		ATOL1J:var_desc = "Variable ATOL1J                                                                 " ;
	float ATOL2J(TSTEP, LAY, PERIM) ;
		ATOL2J:long_name = "ATOL2J          " ;
		ATOL2J:units = "micrograms/m**3 " ;
		ATOL2J:var_desc = "Variable ATOL2J                                                                 " ;
	float ATOL3J(TSTEP, LAY, PERIM) ;
		ATOL3J:long_name = "ATOL3J          " ;
		ATOL3J:units = "micrograms/m**3 " ;
		ATOL3J:var_desc = "Variable ATOL3J                                                                 " ;
	float ABNZ1J(TSTEP, LAY, PERIM) ;
		ABNZ1J:long_name = "ABNZ1J          " ;
		ABNZ1J:units = "micrograms/m**3 " ;
		ABNZ1J:var_desc = "Variable ABNZ1J                                                                 " ;
	float ABNZ2J(TSTEP, LAY, PERIM) ;
		ABNZ2J:long_name = "ABNZ2J          " ;
		ABNZ2J:units = "micrograms/m**3 " ;
		ABNZ2J:var_desc = "Variable ABNZ2J                                                                 " ;
	float ABNZ3J(TSTEP, LAY, PERIM) ;
		ABNZ3J:long_name = "ABNZ3J          " ;
		ABNZ3J:units = "micrograms/m**3 " ;
		ABNZ3J:var_desc = "Variable ABNZ3J                                                                 " ;
	float ATRP1J(TSTEP, LAY, PERIM) ;
		ATRP1J:long_name = "ATRP1J          " ;
		ATRP1J:units = "micrograms/m**3 " ;
		ATRP1J:var_desc = "Variable ATRP1J                                                                 " ;
	float ATRP2J(TSTEP, LAY, PERIM) ;
		ATRP2J:long_name = "ATRP2J          " ;
		ATRP2J:units = "micrograms/m**3 " ;
		ATRP2J:var_desc = "Variable ATRP2J                                                                 " ;
	float AISO1J(TSTEP, LAY, PERIM) ;
		AISO1J:long_name = "AISO1J          " ;
		AISO1J:units = "micrograms/m**3 " ;
		AISO1J:var_desc = "Variable AISO1J                                                                 " ;
	float AISO2J(TSTEP, LAY, PERIM) ;
		AISO2J:long_name = "AISO2J          " ;
		AISO2J:units = "micrograms/m**3 " ;
		AISO2J:var_desc = "Variable AISO2J                                                                 " ;
	float ASQTJ(TSTEP, LAY, PERIM) ;
		ASQTJ:long_name = "ASQTJ           " ;
		ASQTJ:units = "micrograms/m**3 " ;
		ASQTJ:var_desc = "Variable ASQTJ                                                                  " ;
	float ACORS(TSTEP, LAY, PERIM) ;
		ACORS:long_name = "ACORS           " ;
		ACORS:units = "micrograms/m**3 " ;
		ACORS:var_desc = "Variable ACORS                                                                  " ;
	float ASOIL(TSTEP, LAY, PERIM) ;
		ASOIL:long_name = "ASOIL           " ;
		ASOIL:units = "micrograms/m**3 " ;
		ASOIL:var_desc = "Variable ASOIL                                                                  " ;
	float NUMATKN(TSTEP, LAY, PERIM) ;
		NUMATKN:long_name = "NUMATKN         " ;
		NUMATKN:units = "#/m**3          " ;
		NUMATKN:var_desc = "Variable NUMATKN                                                                " ;
	float NUMACC(TSTEP, LAY, PERIM) ;
		NUMACC:long_name = "NUMACC          " ;
		NUMACC:units = "#/m**3          " ;
		NUMACC:var_desc = "Variable NUMACC                                                                 " ;
	float NUMCOR(TSTEP, LAY, PERIM) ;
		NUMCOR:long_name = "NUMCOR          " ;
		NUMCOR:units = "#/m**3          " ;
		NUMCOR:var_desc = "Variable NUMCOR                                                                 " ;
	float SRFATKN(TSTEP, LAY, PERIM) ;
		SRFATKN:long_name = "SRFATKN         " ;
		SRFATKN:units = "m**2/m**3       " ;
		SRFATKN:var_desc = "Variable SRFATKN                                                                " ;
	float SRFACC(TSTEP, LAY, PERIM) ;
		SRFACC:long_name = "SRFACC          " ;
		SRFACC:units = "m**2/m**3       " ;
		SRFACC:var_desc = "Variable SRFACC                                                                 " ;
	float SRFCOR(TSTEP, LAY, PERIM) ;
		SRFCOR:long_name = "SRFCOR          " ;
		SRFCOR:units = "m**2/m**3       " ;
		SRFCOR:var_desc = "Variable SRFCOR                                                                 " ;
	float AISO3J(TSTEP, LAY, PERIM) ;
		AISO3J:long_name = "AISO3J          " ;
		AISO3J:units = "micrograms/m**3 " ;
		AISO3J:var_desc = "Variable AISO3J                                                                 " ;
	float AOLGAJ(TSTEP, LAY, PERIM) ;
		AOLGAJ:long_name = "AOLGAJ          " ;
		AOLGAJ:units = "micrograms/m**3 " ;
		AOLGAJ:var_desc = "Variable AOLGAJ                                                                 " ;
	float AOLGBJ(TSTEP, LAY, PERIM) ;
		AOLGBJ:long_name = "AOLGBJ          " ;
		AOLGBJ:units = "micrograms/m**3 " ;
		AOLGBJ:var_desc = "Variable AOLGBJ                                                                 " ;
	float NH3(TSTEP, LAY, PERIM) ;
		NH3:long_name = "NH3             " ;
		NH3:units = "ppmV            " ;
		NH3:var_desc = "Variable NH3                                                                    " ;
	float SV_ALK(TSTEP, LAY, PERIM) ;
		SV_ALK:long_name = "SV_ALK          " ;
		SV_ALK:units = "ppmV            " ;
		SV_ALK:var_desc = "Variable SV_ALK                                                                 " ;
	float SV_XYL1(TSTEP, LAY, PERIM) ;
		SV_XYL1:long_name = "SV_XYL1         " ;
		SV_XYL1:units = "ppmV            " ;
		SV_XYL1:var_desc = "Variable SV_XYL1                                                                " ;
	float SV_XYL2(TSTEP, LAY, PERIM) ;
		SV_XYL2:long_name = "SV_XYL2         " ;
		SV_XYL2:units = "ppmV            " ;
		SV_XYL2:var_desc = "Variable SV_XYL2                                                                " ;
	float SV_TOL1(TSTEP, LAY, PERIM) ;
		SV_TOL1:long_name = "SV_TOL1         " ;
		SV_TOL1:units = "ppmV            " ;
		SV_TOL1:var_desc = "Variable SV_TOL1                                                                " ;
	float SV_TOL2(TSTEP, LAY, PERIM) ;
		SV_TOL2:long_name = "SV_TOL2         " ;
		SV_TOL2:units = "ppmV            " ;
		SV_TOL2:var_desc = "Variable SV_TOL2                                                                " ;
	float SV_BNZ1(TSTEP, LAY, PERIM) ;
		SV_BNZ1:long_name = "SV_BNZ1         " ;
		SV_BNZ1:units = "ppmV            " ;
		SV_BNZ1:var_desc = "Variable SV_BNZ1                                                                " ;
	float SV_BNZ2(TSTEP, LAY, PERIM) ;
		SV_BNZ2:long_name = "SV_BNZ2         " ;
		SV_BNZ2:units = "ppmV            " ;
		SV_BNZ2:var_desc = "Variable SV_BNZ2                                                                " ;
	float SV_TRP1(TSTEP, LAY, PERIM) ;
		SV_TRP1:long_name = "SV_TRP1         " ;
		SV_TRP1:units = "ppmV            " ;
		SV_TRP1:var_desc = "Variable SV_TRP1                                                                " ;
	float SV_TRP2(TSTEP, LAY, PERIM) ;
		SV_TRP2:long_name = "SV_TRP2         " ;
		SV_TRP2:units = "ppmV            " ;
		SV_TRP2:var_desc = "Variable SV_TRP2                                                                " ;
	float SV_ISO1(TSTEP, LAY, PERIM) ;
		SV_ISO1:long_name = "SV_ISO1         " ;
		SV_ISO1:units = "ppmV            " ;
		SV_ISO1:var_desc = "Variable SV_ISO1                                                                " ;
	float SV_ISO2(TSTEP, LAY, PERIM) ;
		SV_ISO2:long_name = "SV_ISO2         " ;
		SV_ISO2:units = "ppmV            " ;
		SV_ISO2:var_desc = "Variable SV_ISO2                                                                " ;
	float SV_SQT(TSTEP, LAY, PERIM) ;
		SV_SQT:long_name = "SV_SQT          " ;
		SV_SQT:units = "ppmV            " ;
		SV_SQT:var_desc = "Variable SV_SQT                                                                 " ;

// global attributes:
		:IOAPI_VERSION = "$Id: @(#) ioapi library version 3.1 $                                           " ;
		:EXEC_ID = "BCON_V5g_Darwin13_x86_64gfortran                                                " ;
		:FTYPE = 2 ;
		:CDATE = 2014031 ;
		:CTIME = 230334 ;
		:WDATE = 2014031 ;
		:WTIME = 230334 ;
		:SDATE = %(SDATE)d ;
		:STIME = %(STIME)d ;
		:TSTEP = %(TSTEP)d ;
		:NTHIK = 1 ;
		:NCOLS = %(NCOLS)d ;
		:NROWS = %(NROWS)d ;
		:NLAYS = %(NLAYS)d ;
		:NVARS = 78 ;
		:GDTYP = %(GDTYP)d ;
		:P_ALP = %(P_ALP)f ;
		:P_BET = %(P_BET)f ;
		:P_GAM = %(P_GAM)f ;
		:XCENT = %(XCENT)f ;
		:YCENT = %(YCENT)f ;
		:XORIG = %(XORIG)f ;
		:YORIG = %(YORIG)f ;
		:XCELL = %(XCELL)f ;
		:YCELL = %(YCELL)f ;
		:VGTYP = 7 ;
		:VGTOP = %(VGTOP)f ;
		:VGLVLS = %(VGLVLS)s ;
		:GDNAM = "%(GDNAM)16s" ;
		:UPNAM = "%(UPNAM)16s" ;
		:VAR-LIST = "NO2             NO              O3              NO3             OH              HO2             N2O5            HNO3            HONO            PNA             H2O2            NTR             ROOH            FORM            ALD2            PAR             CO              MEPX            FACD            C2O3            PAN             PACD            AACD            PANX            OLE             ETH             IOLE            TOL             CRES            OPEN            MGLY            XYL             ISOP            SO2             SULF            ETHA            BENZENE         ASO4J           ASO4I           AALKJ           AXYL1J          AXYL2J          AXYL3J          ATOL1J          ATOL2J          ATOL3J          ABNZ1J          ABNZ2J          ABNZ3J          ATRP1J          ATRP2J          AISO1J          AISO2J          ASQTJ           ACORS           ASOIL           NUMATKN         NUMACC          NUMCOR          SRFATKN         SRFACC          SRFCOR          AISO3J          AOLGAJ          AOLGBJ          NH3             SV_ALK          SV_XYL1         SV_XYL2         SV_TOL1         SV_TOL2         SV_BNZ1         SV_BNZ2         SV_TRP1         SV_TRP2         SV_ISO1         SV_ISO2         SV_SQT          " ;
		:FILEDESC = "BCON output file BNDY_CONC_1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    " ;
		:HISTORY = "" ;
}
"""

messages = ""
def formatwarning(message, category, filename, lineno, line = 0):
    global messages
    strout = "\n***********\n%s:%s: %s:\n\t%s\n***********\n" % (filename, lineno, category.__name__, message)
    messages += strout
    return strout

warnings.formatwarning = formatwarning

warn = warnings.warn

def get_template(option, gcversion):
    # First add the preamble
    out =r"""{
    "comment1": "comments are added like this",
    "comment2": "Each line has the form \"Species\": {\"expression\": \"some expression\", \"unit\": \"some unit\"},",
    "comment3": "Variables with the output unit (micrograms/m**3) will be multiplied by air density (moles/m**3) and molar mass in the script. It should not be done in the expression.",
    "comment4": "Set \"manual_unit\": true to bypass script unit corrections",
    "comment5": "If you do not have PSURF and/or TMPU, you can use values from the standard atmosphere by replacing the 'expression' value with the 'for_stdatm_use' value",
    "comment6": "more comments are added like this",
    "PRESS": {
            "expression": "(hyam[:].reshape(1, -1).T + hybm[:].reshape(1, -1).T * PSURF[:][:, [0]].T).T * 100.",
            "outunit": "Pa"
        },
    "AIRMOLDEN": {
            "expression": "(hyam[:].reshape(1, -1).T + hybm[:].reshape(1, -1).T * PSURF[:][:, [0]].T).T * 100. / 8.3144621 / TMPU[:]",
            "outunit": "moles/m**3",
            "for_stdatm_use": "((hyam[:, None] + hybm[:, None] * 1013.25) * 100 / 8.3144621 / np.maximum(216.6, 288.15 * ((hyam[:] + hybm[:] * 1013.25) / 1013.25)**C)[:NLAYS, None] * np.ones_like(O3).T).T"
            
        }, 
    "AIRMASSDEN": {
            "expression": "0.0289645 * PSURF[:] * 100 / 8.3144621 / TMPU[:]",
            "outunit": "kg/m**3",
            "for_stdatm_use_as_expr": "(0.0289645 * ((hyam[:, None] + hybm[:, None] * 1013.25) * 100 / 8.3144621 / np.array([287.7,286.9,286.0,285.2,284.3,283.4,282.6,281.7,280.7,279.8,278.9,277.9,276.8,275.3,273.6,271.8,270.0,268.2,265.8,262.8,259.7,256.4,252.9,249.2,245.2,241.0,236.4,230.5,223.6,216.8,216.6,216.6,216.6,216.6,216.6,216.6,216.6,217.5,219.8,222.0,225.3,233.6,250.5,269.2,260.0,237.7,214.3])[:NLAYS, None] * np.ones_like(O3).T).T"
        }, """
    # Next add the gas-phase
    if 'cb05' in option:
        out += r"""
    "CMAQSPECIES": {
        "AIRDEN": {
            "expression": "PSURF[:] * 100 / 8.3144621 / TMPU[:] * 6.022e23 / 1e6",
            "outunit": "molec/cm3",
            "manual_unit": true,
            "comment": "Full calculated unit",
            "for_stdatm_use_as_expr": "((hyam[:, None] + hybm[:, None] * 1013.25) * 100 / 8.3144621 / np.array([287.7,286.9,286.0,285.2,284.3,283.4,282.6,281.7,280.7,279.8,278.9,277.9,276.8,275.3,273.6,271.8,270.0,268.2,265.8,262.8,259.7,256.4,252.9,249.2,245.2,241.0,236.4,230.5,223.6,216.8,216.6,216.6,216.6,216.6,216.6,216.6,216.6,217.5,219.8,222.0,225.3,233.6,250.5,269.2,260.0,237.7,214.3])[:NLAYS, None] * np.ones_like(O3).T).T * 6.022e23 / 1e6"
        },
        "O3": {
            "expression": "O3",
            "outunit": "ppmV"
        },
        "NO": {
            "expression": "NO",
            "outunit": "ppmV"
        },
        "NO2": {
            "expression": "NO2",
            "outunit": "ppmV"
        },
        "NO3": {
            "expression": "NO3",
            "outunit": "ppmV"
        },
        "NH3": {
            "expression": "NH3",
            "outunit": "ppmV"
        },
        "N2O5": {
            "expression": "N2O5",
            "outunit": "ppmV"
        },
        "HONO": {
            "expression": "HNO2",
            "outunit": "ppmV"
        },
        "PNA": {
            "expression": "HNO4",
            "outunit": "ppmV"
        },
        "SO2": {
            "expression": "SO2",
            "outunit": "ppmV"
        },
        "CO": {
            "expression": "CO",
            "outunit": "ppmV"
        },
        "FORM": {
            "expression": "CH2O",
            "outunit": "ppmV"
        },
        "ALD2": {
            "expression": "ALD2",
            "outunit": "ppmV"
        },
        "ALDX": {
            "expression": "RCHO",
            "outunit": "ppmV"
        },
        "OLE": {
            "expression": "0.5 * PRPE",
            "outunit": "ppmV",
            "comment1": "Half the carbon is assumed to be propene and half butene",
            "comment2": "propene maps to 1 PAR and 1 OLE."
        },
        "IOLE": {
            "expression": "0.5 * 1./4. * 3. * PRPE",
            "outunit": "ppmV",
            "comment1": "Half the carbon is assumed to be propene and half butene",
            "comment2": "IOLE is a 4 carbon species.",
            "comment2": "PRPE is a 3 carbon species."
        },
        "ETHA": {
            "expression": "C2H6",
            "outunit": "ppmV"
        },
        "PAR": {
            "expression": "1.5 * C3H8 + 4. * ALK4 + 3. * ACET + 4. * MEK + 1. * BENZ + 0.5 * PRPE",
            "outunit": "ppmV",
            "comment1": "Propane maps as 1.5 PAR and 1.5 unreactive",
            "comment2": "ALK4 is a 4C species in GEOS-Chem and maps as 4 PAR",
            "comment3": "Official mapping of acetone for CB05 is 3 PAR, but should likely be 2 or less",
            "comment4": "Methyl ethyl ketone maps a 4 PAR.",
            "comment5": "Benzene maps as 1PAR and 1BENZENE for SOA",
            "comment6": "Half the carbon is assumed to be propene which is 1 PAR + 1 OLE"
        },
        "TOL": {
            "expression": "TOLU",
            "outunit": "ppmV"
        },
        "XYL": {
            "expression": "XYLE",
            "outunit": "ppmV"
        },
        "ISOP": {
            "expression": "ISOP",
            "outunit": "ppmV"
        },
        "ISPD": {
            "expression": "MACR + MVK",
            "outunit": "ppmV"
        },
        "HNO3": {
            "expression": "HNO3",
            "outunit": "ppmV"
        },
        "NTR": {
            "expression": "R4N2",
            "outunit": "ppmV"
        },
        "H2O2": {
            "expression": "H2O2",
            "outunit": "ppmV"
        },
        "MEPX": {
            "expression": "MP",
            "outunit": "ppmV"
        },
        "ROOH": {
            "expression": "RIP",
            "outunit": "ppmV"
        },
        "PAN": {
            "expression": "PAN",
            "outunit": "ppmV"
        },
        "PANX": {
            "expression": "PPN + PMN",
            "outunit": "ppmV"
        },
        "PACD": {
            "expression": "MAP",
            "outunit": "ppmV"
        },"""
        
        if 'ae6' in option:
            if gcversion > 8:
                out += r"""
        "BENZENE": {
            "expression": "BENZ",
            "outunit": "ppmV"
        },
        "ABNZ1J": {
            "expression": "0.12 * (ASOAN + ASOA1 + ASOA2 + ASOA3)",
            "outunit": "micrograms/m**3"
        },
        "ABNZ2J": {
            "expression": "0.04 * (ASOAN + ASOA1 + ASOA2 + ASOA3)",
            "outunit": "micrograms/m**3"
        },
        "ABNZ3J": {
            "expression": "0.32 * (ASOAN + ASOA1 + ASOA2 + ASOA3)",
            "outunit": "micrograms/m**3"
        },
        "AALJ": {
            "expression": "0.05695 * DST1",
            "outunit": "micrograms/m**3"
        },
        "ACAJ": {
            "expression": "0.0118 * SALA + 0.07940 * DST1",
            "outunit": "micrograms/m**3"
        },
        "ACLJ": {
            "expression": "0.00945 * DST1 + 0.5538 * SALA",
            "outunit": "micrograms/m**3"
        },
        "ACLK": {
            "expression": "0.01190 * DST2 + 0.01190 * DST3 + 0.01190 * DST4 + 0.5538 * SALC",
            "outunit": "micrograms/m**3"
        },
        "AECI": {
            "expression": "0.001 * BCPI + 0.001 * BCPO",
            "outunit": "micrograms/m**3"
        },
        "AECJ": {
            "expression": "0.999 * BCPI + 0.999 * BCPO",
            "outunit": "micrograms/m**3"
        },
        "AFEJ": {
            "expression": "0.03355 * DST1",
            "outunit": "micrograms/m**3"
        },
        "AISO1J": {
            "expression": "0.75 * (ISOA1 + ISOA2 + ISOA3)",
            "outunit": "micrograms/m**3"
        },
        "AISO2J": {
            "expression": "0.25 * (ISOA1 + ISOA2 + ISOA3)",
            "outunit": "micrograms/m**3"
        },
        "AKJ": {
            "expression": "0.0114 * SALA + 0.03770 * DST1",
            "outunit": "micrograms/m**3"
        },
        "AMGJ": {
            "expression": "0.0368 * SALA",
            "outunit": "micrograms/m**3"
        },
        "AMNJ": {
            "expression": "0.00115 * DST1",
            "outunit": "micrograms/m**3"
        },
        "ANAJ": {
            "expression": "0.3086 * SALA + 0.03935 * DST1",
            "outunit": "micrograms/m**3"
        },
        "ANH4I": {
            "expression": "0.01 * NH4",
            "outunit": "micrograms/m**3"
        },
        "ANH4J": {
            "expression": "0.00005 * DST1 + 0.99 * NH4",
            "outunit": "micrograms/m**3"
        },
        "ANO3I": {
            "expression": "0.01 * NIT",
            "outunit": "micrograms/m**3"
        },
        "ANO3J": {
            "expression": "0.00020 * DST1 + 0.99 * NIT",
            "outunit": "micrograms/m**3"
        },
        "ANO3K": {
            "expression": "0.0016 * DST2 + 0.0016 * DST3 + 0.0016 * DST4 + NITs",
            "outunit": "micrograms/m**3"
        },
        "AOTHRJ": {
            "expression": "0.50219 * DST1",
            "outunit": "micrograms/m**3"
        },
        "APNCOMI": {
            "expression": "0.4 * 0.001 * OCPI + 0.4 * 0.001 * OCPO",
            "outunit": "micrograms/m**3"
        },
        "APNCOMJ": {
            "expression": "0.4 * 0.999 * OCPI + 0.4 * 0.999 * OCPO",
            "outunit": "micrograms/m**3"
        },
        "APNCOMJ": {
            "expression": "0.0043 * DST1",
            "outunit": "micrograms/m**3"
        },
        "APOCI": {
            "expression": "0.001 * OCPI + 0.001 * OCPO",
            "outunit": "micrograms/m**3"
        },
        "APOCJ": {
            "expression": "0.999 * OCPI + 0.999 * OCPO + 0.01075 * DST1",
            "outunit": "micrograms/m**3"
        },
        "ASEACAT": {
            "expression": "0.3685 * SALC",
            "outunit": "micrograms/m**3"
        },
        "ASIJ": {
            "expression": "0.19435 * DST1",
            "outunit": "micrograms/m**3"
        },
        "ASO4I": {
            "expression": "0.01 * SO4",
            "outunit": "micrograms/m**3"
        },
        "ASO4J": {
            "expression": "0.99 * SO4 + 0.0225 * DST1 + 0.0776 * SALA",
            "outunit": "micrograms/m**3"
        },
        "ASO4K": {
            "expression": "0.0776 * SALC + 0.02655 * DST2 + 0.02655 * DST3 + 0.02655 * DST4 + SO4s",
            "outunit": "micrograms/m**3"
        },
        "ASOIL": {
            "expression": "0.95995 * DST2 + 0.95995 * DST3 + 0.95995 * DST4",
            "outunit": "micrograms/m**3"
        },
        "ASQTJ": {
            "expression": "0.33 * TSOA0 + 0.33 * TSOA1 + 0.33 * TSOA2 + 0.33 * TSOA3",
            "outunit": "micrograms/m**3"
        },
        "ATIJ": {
            "expression": "0.0028 * DST1",
            "outunit": "micrograms/m**3"
        },
        "ATOL1J": {
            "expression": "0.04 * (ASOAN + ASOA1 + ASOA2 + ASOA3)",
            "outunit": "micrograms/m**3"
        },
        "ATOL2J": {
            "expression": "0.04 * (ASOAN + ASOA1 + ASOA2 + ASOA3)",
            "outunit": "micrograms/m**3"
        },
        "ATOL3J": {
            "expression": "0.29 * (ASOAN + ASOA1 + ASOA2 + ASOA3)",
            "outunit": "micrograms/m**3"
        },
        "ATRP1J": {
            "expression": "0.33 * TSOA0 + 0.33 * TSOA1 + 0.33 * TSOA2 + 0.33 * TSOA3",
            "outunit": "micrograms/m**3"
        },
        "ATRP2J": {
            "expression": "0.34 * TSOA0 + 0.34 * TSOA1 + 0.34 * TSOA2 + 0.34 * TSOA3",
            "outunit": "micrograms/m**3"
        },
        "AXYL1J": {
            "expression": "0.03 * (ASOAN + ASOA1 + ASOA2 + ASOA3)",
            "outunit": "micrograms/m**3"
        },
        "AXYL2J": {
            "expression": "0.01 * (ASOAN + ASOA1 + ASOA2 + ASOA3)",
            "outunit": "micrograms/m**3"
        },
        "AXYL3J": {
            "expression": "0.11 * (ASOAN + ASOA1 + ASOA2 + ASOA3)",
            "outunit": "micrograms/m**3"
        },
        "SV_BNZ1": {
            "expression": "0.06 * (ASOG1 + ASOG2 + ASOG3)",
            "outunit": "ppmV"
        },
        "SV_BNZ2": {
            "expression": "0.23 * (ASOG1 + ASOG2 + ASOG3)",
            "outunit": "ppmV"
        },
        "SV_ISO1": {
            "expression": "0.75 * (ISOG1 + ISOG2 + ISOG3)",
            "outunit": "ppmV"
        },
        "SV_ISO2": {
            "expression": "0.25 * (ISOG1 + ISOG2 + ISOG3)",
            "outunit": "ppmV"
        },
        "SV_SQT": {
            "expression": "0.33 * TSOG0 + 0.33 * TSOG1 + 0.33 * TSOG2 + 0.33 * TSOG3",
            "outunit": "ppmV"
        },
        "SV_TOL1": {
            "expression": "0.23 * (ASOG1 + ASOG2 + ASOG3)",
            "outunit": "ppmV"
        },
        "SV_TOL2": {
            "expression": "0.23 * (ASOG1 + ASOG2 + ASOG3)",
            "outunit": "ppmV"
        },
        "SV_TRP1": {
            "expression": "0.33 * TSOG0 + 0.33 * TSOG1 + 0.33 * TSOG2 + 0.33 * TSOG3",
            "outunit": "ppmV"
        },
        "SV_TRP2": {
            "expression": "0.34 * TSOG0 + 0.34 * TSOG1 + 0.34 * TSOG2 + 0.34 * TSOG3",
            "outunit": "ppmV"
        },
        "SV_XYL1": {
            "expression": "0.19 * (ASOG1 + ASOG2 + ASOG3)",
            "outunit": "ppmV"
        },
        "SV_XYL2": {
            "expression": "0.06 * (ASOG1 + ASOG2 + ASOG3)",
            "outunit": "ppmV"
        }"""
            else:
                out += r"""
        "AALJ": {
            "expression": "0.05695 * DST1",
            "outunit": "micrograms/m**3"
        }
        "AALKJ": {
            "expression": "AALKJ",
            "outunit": "micrograms/m**3"
        }
        "ABNZ1J": {
            "expression": "0.12 * SOA5",
            "outunit": "micrograms/m**3"
        }
        "ABNZ2J": {
            "expression": "0.04 * SOA5",
            "outunit": "micrograms/m**3"
        }
        "ABNZ3J": {
            "expression": "0.32 * SOA5",
            "outunit": "micrograms/m**3"
        }
        "ACAJ": {
            "expression": "0.0118 * SALA",
            "outunit": "micrograms/m**3"
        }
        "ACAJ": {
            "expression": "0.07940 * DST1",
            "outunit": "micrograms/m**3"
        }
        "ACLJ": {
            "expression": "0.00945 * DST1",
            "outunit": "micrograms/m**3"
        }
        "ACLJ": {
            "expression": "0.5538 * SALA",
            "outunit": "micrograms/m**3"
        }
        "ACLK": {
            "expression": "0.01190 * DST2",
            "outunit": "micrograms/m**3"
        }
        "ACLK": {
            "expression": "0.01190 * DST3",
            "outunit": "micrograms/m**3"
        }
        "ACLK": {
            "expression": "0.01190 * DST4",
            "outunit": "micrograms/m**3"
        }
        "ACLK": {
            "expression": "0.5538 * SALC",
            "outunit": "micrograms/m**3"
        }
        "ACORS": {
            "expression": "ACORS",
            "outunit": "micrograms/m**3"
        }
        "AECI": {
            "expression": "0.001 * BCPI",
            "outunit": "micrograms/m**3"
        }
        "AECI": {
            "expression": "0.001 * BCPO",
            "outunit": "micrograms/m**3"
        }
        "AECJ": {
            "expression": "0.999 * BCPI",
            "outunit": "micrograms/m**3"
        }
        "AECJ": {
            "expression": "0.999 * BCPO",
            "outunit": "micrograms/m**3"
        }
        "AFEJ": {
            "expression": "0.03355 * DST1",
            "outunit": "micrograms/m**3"
        }
        "AISO1J": {
            "expression": "0.75 * SOA4",
            "outunit": "micrograms/m**3"
        }
        "AISO2J": {
            "expression": "0.25 * SOA4",
            "outunit": "micrograms/m**3"
        }
        "AISO3J": {
            "expression": "AISO3J",
            "outunit": "micrograms/m**3"
        }
        "AKJ": {
            "expression": "0.0114 * SALA",
            "outunit": "micrograms/m**3"
        }
        "AKJ": {
            "expression": "0.03770 * DST1",
            "outunit": "micrograms/m**3"
        }
        "AMGJ": {
            "expression": "0.0368 * SALA",
            "outunit": "micrograms/m**3"
        }
        "AMNJ": {
            "expression": "0.00115 * DST1",
            "outunit": "micrograms/m**3"
        }
        "ANAJ": {
            "expression": "0.3086 * SALA",
            "outunit": "micrograms/m**3"
        }
        "ANAJ": {
            "expression": "0.03935 * DST1",
            "outunit": "micrograms/m**3"
        }
        "ANH4I": {
            "expression": "0.01 * NH4",
            "outunit": "micrograms/m**3"
        }
        "ANH4J": {
            "expression": "0.00005 * DST1",
            "outunit": "micrograms/m**3"
        }
        "ANH4J": {
            "expression": "0.99 * NH4",
            "outunit": "micrograms/m**3"
        }
        "ANO3I": {
            "expression": "0.01 * NIT",
            "outunit": "micrograms/m**3"
        }
        "ANO3J": {
            "expression": "0.00020 * DST1",
            "outunit": "micrograms/m**3"
        }
        "ANO3J": {
            "expression": "0.99 * NIT",
            "outunit": "micrograms/m**3"
        }
        "ANO3K": {
            "expression": "0.0016 * DST2",
            "outunit": "micrograms/m**3"
        }
        "ANO3K": {
            "expression": "0.0016 * DST3",
            "outunit": "micrograms/m**3"
        }
        "ANO3K": {
            "expression": "0.0016 * DST4",
            "outunit": "micrograms/m**3"
        }
        "ANO3K": {
            "expression": "NITs",
            "outunit": "micrograms/m**3"
        }
        "AOLGAJ": {
            "expression": "AOLGAJ",
            "outunit": "micrograms/m**3"
        }
        "AOLGBJ": {
            "expression": "AOLGBJ",
            "outunit": "micrograms/m**3"
        }
        "AOTHRJ": {
            "expression": "0.50219 * DST1",
            "outunit": "micrograms/m**3"
        }
        "APNCOMI": {
            "expression": "0.4 * 0.001 * OCPI",
            "outunit": "micrograms/m**3"
        }
        "APNCOMI": {
            "expression": "0.4 * 0.001 * OCPO",
            "outunit": "micrograms/m**3"
        }
        "APNCOMJ": {
            "expression": "0.4 * 0.999 * OCPI",
            "outunit": "micrograms/m**3"
        }
        "APNCOMJ": {
            "expression": "0.4 * 0.999 * OCPO",
            "outunit": "micrograms/m**3"
        }
        "APNCOMJ": {
            "expression": "0.0043 * DST1",
            "outunit": "micrograms/m**3"
        }
        "APOCI": {
            "expression": "0.001 * OCPI",
            "outunit": "micrograms/m**3"
        }
        "APOCI": {
            "expression": "0.001 * OCPO",
            "outunit": "micrograms/m**3"
        }
        "APOCJ": {
            "expression": "0.999 * OCPI",
            "outunit": "micrograms/m**3"
        }
        "APOCJ": {
            "expression": "0.999 * OCPO",
            "outunit": "micrograms/m**3"
        }
        "APOCJ": {
            "expression": "0.01075 * DST1",
            "outunit": "micrograms/m**3"
        }
        "ASEACAT": {
            "expression": "0.3685 * SALC",
            "outunit": "micrograms/m**3"
        }
        "ASIJ": {
            "expression": "0.19435 * DST1",
            "outunit": "micrograms/m**3"
        }
        "ASO4I": {
            "expression": "0.01 * SO4",
            "outunit": "micrograms/m**3"
        }
        "ASO4J": {
            "expression": "0.99 * SO4",
            "outunit": "micrograms/m**3"
        }
        "ASO4J": {
            "expression": "0.0225 * DST1",
            "outunit": "micrograms/m**3"
        }
        "ASO4J": {
            "expression": "0.0776 * SALA",
            "outunit": "micrograms/m**3"
        }
        "ASO4K": {
            "expression": "0.0776 * SALC",
            "outunit": "micrograms/m**3"
        }
        "ASO4K": {
            "expression": "0.02655 * DST2",
            "outunit": "micrograms/m**3"
        }
        "ASO4K": {
            "expression": "0.02655 * DST3",
            "outunit": "micrograms/m**3"
        }
        "ASO4K": {
            "expression": "0.02655 * DST4",
            "outunit": "micrograms/m**3"
        }
        "ASO4K": {
            "expression": "SO4s",
            "outunit": "micrograms/m**3"
        }
        "ASOIL": {
            "expression": "0.95995 * DST2",
            "outunit": "micrograms/m**3"
        }
        "ASOIL": {
            "expression": "0.95995 * DST3",
            "outunit": "micrograms/m**3"
        }
        "ASOIL": {
            "expression": "0.95995 * DST4",
            "outunit": "micrograms/m**3"
        }
        "ASQTJ": {
            "expression": "SOA3",
            "outunit": "micrograms/m**3"
        }
        "ATIJ": {
            "expression": "0.0028 * DST1",
            "outunit": "micrograms/m**3"
        }
        "ATOL1J": {
            "expression": "0.04 * SOA5",
            "outunit": "micrograms/m**3"
        }
        "ATOL2J": {
            "expression": "0.04 * SOA5",
            "outunit": "micrograms/m**3"
        }
        "ATOL3J": {
            "expression": "0.29 * SOA5",
            "outunit": "micrograms/m**3"
        }
        "ATRP1J": {
            "expression": "0.33 * SOA1",
            "outunit": "micrograms/m**3"
        }
        "ATRP1J": {
            "expression": "0.33 * SOA2",
            "outunit": "micrograms/m**3"
        }
        "ATRP2J": {
            "expression": "0.67 * SOA1",
            "outunit": "micrograms/m**3"
        }
        "ATRP2J": {
            "expression": "0.67 * SOA2",
            "outunit": "micrograms/m**3"
        }
        "AXYL1J": {
            "expression": "0.03 * SOA5",
            "outunit": "micrograms/m**3"
        }
        "AXYL2J": {
            "expression": "0.01 * SOA5",
            "outunit": "micrograms/m**3"
        }
        "AXYL3J": {
            "expression": "0.11 * SOA5",
            "outunit": "micrograms/m**3"
        }
        "NH3": {
            "expression": "NH3",
            "outunit": "micrograms/m**3"
        }
        "NUMACC": {
            "expression": "NUMACC",
            "outunit": "micrograms/m**3"
        }
        "NUMATKN": {
            "expression": "NUMATKN",
            "outunit": "micrograms/m**3"
        }
        "NUMCOR": {
            "expression": "NUMCOR",
            "outunit": "micrograms/m**3"
        }
        "SRFACC": {
            "expression": "SRFACC",
            "outunit": "micrograms/m**3"
        }
        "SRFATKN": {
            "expression": "SRFATKN",
            "outunit": "micrograms/m**3"
        }
        "SRFCOR": {
            "expression": "SRFCOR",
            "outunit": "micrograms/m**3"
        }
        "SULF": {
            "expression": "SULF",
            "outunit": "micrograms/m**3"
        }
        "SV_ALK": {
            "expression": "SV_ALK",
            "outunit": "ppmV"
        }
        "SV_BNZ1": {
            "expression": "0.06 * SOG5",
            "outunit": "ppmV"
        }
        "SV_BNZ2": {
            "expression": "0.23 * SOG5",
            "outunit": "ppmV"
        }
        "SV_ISO1": {
            "expression": "0.75 * SOG4",
            "outunit": "ppmV"
        }
        "SV_ISO2": {
            "expression": "0.25 * SOG4",
            "outunit": "ppmV"
        }
        "SV_SQT": {
            "expression": "SOG3",
            "outunit": "ppmV"
        }
        "SV_TOL1": {
            "expression": "0.23 * SOG5",
            "outunit": "ppmV"
        }
        "SV_TOL2": {
            "expression": "0.23 * SOG5",
            "outunit": "ppmV"
        }
        "SV_TRP1": {
            "expression": "0.33 * SOG1",
            "outunit": "ppmV"
        }
        "SV_TRP1": {
            "expression": "0.33 * SOG2",
            "outunit": "ppmV"
        }
        "SV_TRP2": {
            "expression": "0.67 * SOG1",
            "outunit": "ppmV"
        }
        "SV_TRP2": {
            "expression": "0.67 * SOG2",
            "outunit": "ppmV"
        }
        "SV_XYL1": {
            "expression": "0.19 * SOG5",
            "outunit": "ppmV"
        }
        "SV_XYL2": {
            "expression": "0.06 * SOG5",
            "outunit": "ppmV"
        }"""

    
    out += r"""
    }
}
"""
    
    return out

class speciesstruct(object):
    def __init__(self, name, ind, mode, density, version, found):
        """
        name - text string identifying aerosol
        ind - variable number (still needed?)
        mode - 0: Aitken, 1: Accumulation, 2: Coarse
        density - aerosol density (kg/m^3)
        version - 0: both, 5: AERO5, 6: AERO6
        """
        self.name = name
        self.mode = mode
        self.density = density
        self.version = version
        self.found = found

def repair_ae(f, myioo):
    if myioo is None:
        status = warn = error = warnings.warn
    else:
        warn = myioo.warn
        status = myioo.status
        error = myioo.error
    
    RHOSO4  = 1.8e3 # bulk density of aerosol sulfate
    RHONH4  = 1.8e3 # bulk density of aerosol ammonium
    RHONO3  = 1.8e3 # bulk density of aerosol nitrate
    RHOORG  = 1.3e3 # bulk density for aerosol organics following Carlton et al. 2010
    RHOSOIL = 2.6e3 # bulk density for aerosol soil dust
    RHOSEAS = 2.2e3 # bulk density for marine aerosol
    RHOANTH = 2.2e3 # bulk density for anthropogenic aerosol
    SGINIAT = 1.7   # initial sigma-G for Aitken mode
    SGINIAC = 2.0   # initial sigma-G for accumulation mode
    SGINICO = 2.2   # initial sigma-G for coarse mode
    DGINIAT = 0.01E-6  # geometric mean diameter for Aitken mode [ m ]
    DGINIAC = 0.07E-6  # geometric mean diameter for accum  mode [ m ]
    DGINICO = 1.0E-6   # geometric mean diameter for coarse mode [ m ]
    CONMIN  = 1.0E-30  # minimum concentration [ ug/m**3 ]
    nspecies = 57      # number of aerosol species treated

    #...conversion factors for number and surface area
    NUMFAC = {}
    NUMFAC['ATKN'] = 1.0 / ( ( DGINIAT ** 3.0 ) * exp( ( 9.0 / 2.0 ) * ( ( log( SGINIAT ) ) ** 2.0 ) ) )
    NUMFAC['ACC']  = 1.0 / ( ( DGINIAC ** 3.0 ) * exp( ( 9.0 / 2.0 ) * ( ( log( SGINIAC ) ) ** 2.0 ) ) )
    NUMFAC['COR']  = 1.0 / ( ( DGINICO ** 3.0 ) * exp( ( 9.0 / 2.0 ) * ( ( log( SGINICO ) ) ** 2.0 ) ) )
    SRFFAC = {}
    SRFFAC['ATKN'] = pi / ( DGINIAT * exp( ( 5.0 / 2.0 ) * ( ( log( SGINIAT ) ) ** 2.0 ) ) )
    SRFFAC['ACC']  = pi / ( DGINIAC * exp( ( 5.0 / 2.0 ) * ( ( log( SGINIAC ) ) ** 2.0 ) ) )
    SRFFAC['COR']  = pi / ( DGINICO * exp( ( 5.0 / 2.0 ) * ( ( log( SGINICO ) ) ** 2.0 ) ) )

    bcspcs = [speciesstruct (    'ACLI',   0, 'ATKN', RHOSEAS, 0, False),
             speciesstruct (    'AECI',   0, 'ATKN', RHOANTH, 0, False),
             speciesstruct (    'ANAI',   0, 'ATKN', RHOSEAS, 0, False),
             speciesstruct (   'ANH4I',   0, 'ATKN',  RHONH4, 0, False),
             speciesstruct (   'ANO3I',   0, 'ATKN',  RHONO3, 0, False),
             speciesstruct (   'ASO4I',   0, 'ATKN',  RHOSO4, 0, False),
             speciesstruct (    'A25I',   0, 'ATKN', RHOANTH, 5, False),
             speciesstruct ( 'AORGPAI',   0, 'ATKN',  RHOORG, 5, False),
             speciesstruct (  'AOTHRI',   0, 'ATKN', RHOANTH, 6, False),
             speciesstruct ( 'APNCOMI',   0, 'ATKN',  RHOORG, 6, False),
             speciesstruct (   'APOCI',   0, 'ATKN',  RHOORG, 6, False),
             speciesstruct (   'AALKJ',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'ABNZ1J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'ABNZ2J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'ABNZ3J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (    'ACLJ',   0,  'ACC', RHOSEAS, 0, False),
             speciesstruct (    'AECJ',   0,  'ACC', RHOANTH, 0, False),
             speciesstruct (  'AISO1J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'AISO2J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'AISO3J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (    'ANAJ',   0,  'ACC', RHOSEAS, 0, False),
             speciesstruct (   'ANH4J',   0,  'ACC',  RHONH4, 0, False),
             speciesstruct (   'ANO3J',   0,  'ACC',  RHONO3, 0, False),
             speciesstruct (  'AOLGAJ',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'AOLGBJ',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'AORGCJ',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (   'ASO4J',   0,  'ACC',  RHOSO4, 0, False),
             speciesstruct (   'ASQTJ',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'ATOL1J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'ATOL2J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'ATOL3J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'ATRP1J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'ATRP2J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'AXYL1J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'AXYL2J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (  'AXYL3J',   0,  'ACC',  RHOORG, 0, False),
             speciesstruct (    'A25J',   0,  'ACC', RHOANTH, 5, False),
             speciesstruct ( 'AORGPAJ',   0,  'ACC',  RHOORG, 5, False),
             speciesstruct (    'AALJ',   0,  'ACC', RHOSOIL, 6, False),
             speciesstruct (    'ACAJ',   0,  'ACC', RHOSOIL, 6, False),
             speciesstruct (    'AFEJ',   0,  'ACC', RHOSOIL, 6, False),
             speciesstruct (     'AKJ',   0,  'ACC', RHOSOIL, 6, False),
             speciesstruct (    'AMGJ',   0,  'ACC', RHOSEAS, 6, False),
             speciesstruct (    'AMNJ',   0,  'ACC', RHOSOIL, 6, False),
             speciesstruct (  'AOTHRJ',   0,  'ACC', RHOSOIL, 6, False),
             speciesstruct ( 'APNCOMJ',   0,  'ACC',  RHOORG, 6, False),
             speciesstruct (   'APOCJ',   0,  'ACC',  RHOORG, 6, False),
             speciesstruct (    'ASIJ',   0,  'ACC', RHOSOIL, 6, False),
             speciesstruct (    'ATIJ',   0,  'ACC', RHOSOIL, 6, False),
             speciesstruct (    'ACLK',   0,  'COR', RHOSEAS, 0, False),
             speciesstruct (   'ACORS',   0,  'COR', RHOANTH, 0, False),
             speciesstruct (   'ANH4K',   0,  'COR',  RHONH4, 0, False),
             speciesstruct (   'ANO3K',   0,  'COR',  RHONO3, 0, False),
             speciesstruct (   'ASO4K',   0,  'COR',  RHOSO4, 0, False),
             speciesstruct (   'ASOIL',   0,  'COR', RHOSOIL, 0, False),
             speciesstruct (    'ANAK',   0,  'COR', RHOSEAS, 5, False),
             speciesstruct ( 'ASEACAT',   0,  'COR', RHOSEAS, 6, False)]

    for spc in bcspcs:
        try:
            data_shape = f.variables[spc.name].shape
            break
        except:
            pass
    else:
        warn("There are no aerosol species")
        return
    not_found = defaultdict(lambda: [])
    for spc in bcspcs:
        try:
            spc.found = spc.name in f.variables.keys()
        except:
            spc.found = False
        if not spc.found:
            not_found[spc.version].append(spc.name)
    version_check = []
    for k, v in not_found.iteritems():
        if k != 0:
            version_check.append(k)
        warn('Some variables from %d were not found: %s' % (k, v))
    
    if len(version_check) > 1:
        warn('Some variables from aerosol versions %s were not found' % ' and '.join([str(v) for v in version_check]))
    moment3 = dict([(k, zeros(data_shape, dtype = 'f')) for k in 'ATKN ACC COR'.split()])
    
    for spc in bcspcs:
        try:
            bcval = f.variables[spc.name]
        except KeyError:
            continue
        v = 1.0e-9*6.0/( pi*spc.density ) * bcval[:]
        moment3[spc.mode] += v
    
    for modek, modv in moment3.iteritems():
        numkey = 'NUM' + modek
        srfkey = 'SRF' + modek
        if numkey in f.variables.keys():
            numvar = f.variables[numkey]
        else:
            bcspc = [spc.name for spc in bcspcs if spc.found][0]
            dims = f.variables[bcspc].dimensions
            numvar = f.createVariable(numkey, 'f', dims)
            numvar.units = '#/m**3'.ljust(16);
            numvar.long_name = numkey.ljust(16);
            numvar.var_desc = numkey.ljust(16);

        if srfkey in f.variables.keys():
            srfvar = f.variables[srfkey]
        else:
            bcspc = [spc.name for spc in bcspcs if spc.found][0]
            dims = f.variables[bcspc].dimensions
            srfvar = f.createVariable(srfkey, 'f', dims)
            srfvar.units = 'm**2/m**3'.ljust(16);
            srfvar.long_name = srfkey.ljust(16);
            srfvar.var_desc = srfkey.ljust(16);
        numvar[:] = NUMFAC[modek] * modv
        srfvar[:] = SRFFAC[modek] * modv
    f.sync()

def makeregriddedbdy(metbdy, args, spcs):
    lonlatcoords = '/'.join(['%f,%f' % (lon, lat) for lon, lat in zip(metbdy.variables['longitude'][:].ravel(), metbdy.variables['latitude'][:].ravel())])
    out = []
    for ND49, ND49_REGRID_BDY in zip(args.ND49, args.ND49_REGRID_BDY):
        outf = extract(getvarpnc(ND49, None), lonlatcoords, method = args.extractmethod)
        if args.persistintermediate:
            pncgen(outf, ND49_REGRID_BDY, verbose = False)
            outf = Dataset(ND49_REGRID_BDY)
        out.append(outf)
    return out

def makeregriddedcro(metcro, args, spcs):
    lonlatcoords = '/'.join(['%f,%f' % (lon, lat) for lon, lat in zip(metcro.variables['longitude'][:].ravel(), metcro.variables['latitude'][:].ravel())])
    out = []
    for ND49, ND49_REGRID_CRO in zip(args.ND49, args.ND49_REGRID_CRO):
        outf = extract(slice_dim(getvarpnc(ND49, None), 'time,0'), lonlatcoords, method = args.extractmethod)
        if args.persistintermediate:
            pncgen(outf, ND49_REGRID_CRO, verbose = False)
            outf = Dataset(ND49_REGRID_CRO)
        out.append(outf)
        
    return out

def getdefault(oldcon, vark, noutstep):
    defval = np.ma.filled(oldcon.variables[vark][:], 0)
    if defval.shape[0] == 1:
        defval = defval.repeat(noutstep, 0)
    elif defval.shape[0] > noutstep:
        warn('Default (I,B)CON has time dimension that is greater than output (%d>%d); only first %d are used.' % (defval.shape[0], noutstep, noutstep))
        defval = defval[0:noutstep]
    return defval


def makebcon(args):
    global messages
    messages += ' '.join(sys.argv[:]) + '\n'
    mappings_file = json.load(file(args.mapping, mode = 'r'))
    mappings = mappings_file['CMAQSPECIES']
    spcs = ','.join(reduce(tuple.__add__, [compile(mapping['expression'], 'mapping', 'eval').co_names for mapping in mappings.values() + [mappings_file['AIRMOLDEN']]]))

    if args.METBDY3D is None:
        dobcon = False
    else:
        dobcon = True

    if args.METCRO3D is None:
        doicon = False
    else:
        doicon = True
    
    if dobcon:
        metbdyfiles, metbdyargs = pncparse(has_ofile = False, plot_options = False, interactive = False, args = args.METBDY3D.split(' '), parser = None)
        metbdy = getvarpnc(metbdyfiles[0], ['TA', 'PRES'])
        add_cf_from_ioapi(metbdy)
        regridded_nd49_bdy = makeregriddedbdy(metbdy, args, spcs)
        

        metbdyprops = dict([(propk, getattr(metbdy, propk)) for propk in metbdy.ncattrs()])
        metbdyprops['VGLVLS'] = 'f, '.join(['%f' % vgl for vgl in metbdyprops['VGLVLS']]) + 'f'
        metbdyprops['PERIM'] = len(metbdy.dimensions['PERIM'])
        metbdyprops['GDTYP'] = metbdy.GDTYP

        if args.BCON == 'dummybcon.nc':
            bcontmp = file('bcon.cdl.tmp', 'w')
            bcontmp.write(getbconcdl() % metbdyprops)
            bcontmp.close()
            os.system("ncgen -o dummybcon.nc bcon.cdl.tmp")
            args.BCON = 'dummybcon.nc'
            oldbcon = Dataset(args.BCON, 'r+')
            oldbcon.variables['TFLAG'][0] = 0  
            for k, v in oldbcon.variables.iteritems():
                v[0] = 1e-32      
        else:
            bconfiles, bconargs = pncparse(has_ofile = False, plot_options = False, interactive = False, args = args.BCON.split(' '), parser = None)
            oldbcon = bconfiles[0]
        newbcon = Dataset(args.NEWBCON, mode = 'w', format = 'NETCDF3_CLASSIC')
        for propk in oldbcon.ncattrs():
            propv = getattr(oldbcon, propk)
            setattr(newbcon, propk, propv)

        newbcon.createDimension('TSTEP', None)
        newbcon.createDimension('DATE-TIME', 2)
        for dimk in ('LAY', 'PERIM'):
            newbcon.createDimension(dimk, len(oldbcon.dimensions[dimk]))
        for vark, oldv in oldbcon.variables.iteritems():
            if vark == 'TFLAG': continue
            newv = newbcon.createVariable(vark, oldv.dtype.char, oldv.dimensions)
            for propk in oldv.ncattrs():
                propv = getattr(oldv, propk)
                setattr(newv, propk, propv)
        nbconoutsteps = sum([len(tmpf.dimensions['time']) for tmpf in regridded_nd49_bdy])
        if len(newbcon.dimensions['PERIM']) != metbdyprops['PERIM']:
            raise ValueError('Default (I,B)CON has different perimeter dimension (%d) than the output file (%d).' % (len(newbcon.dimensions['PERIM']), metbdyprops['PERIM']))
        obconsteps = len(oldbcon.dimensions['TSTEP'])
        if obconsteps != 1 and obconsteps < nbconoutsteps:
            raise ValueError('Default BCON has time dimension that is less than output (%d<%d) and not 1. Use time-independent file or file with same times.' % (obconsteps, nbconoutsteps))

    
    if doicon:
        metcrofiles, metcroargs = pncparse(has_ofile = False, plot_options = False, interactive = False, args = args.METCRO3D.split(' '), parser = None)
        metcro = slice_dim(getvarpnc(metcrofiles[0], ['TA', 'PRES']), 'TSTEP,0')
        add_cf_from_ioapi(metcro)
        regridded_nd49_cro = makeregriddedcro(metcro, args, spcs)
        metcroprops = metbdyprops.copy()
        metcroprops['GDTYP'] = metcro.GDTYP
        if args.ICON == 'dummyicon.nc':
            icontmp = file('icon.cdl.tmp', 'w')
            icontmp.write(geticoncdl() % metcroprops)
            icontmp.close()
            os.system("ncgen -o dummyicon.nc icon.cdl.tmp")
            args.ICON = 'dummyicon.nc'
            oldicon = Dataset(args.ICON, 'r+')
            oldicon.variables['TFLAG'][0] = 0        
            for k, v in oldicon.variables.iteritems():
                v[0] = 1e-32
        else:
            iconfiles, iconargs = pncparse(has_ofile = False, plot_options = False, interactive = False, args = args.ICON.split(' '), parser = None)
            oldicon = iconfiles[0]
        oldicon = Dataset(args.ICON, mode = 'r+', format = 'NETCDF3_CLASSIC')
        newicon = Dataset(args.NEWICON, mode = 'w', format = 'NETCDF3_CLASSIC')
        for propk in oldicon.ncattrs():
            propv = getattr(oldicon, propk)
            setattr(newicon, propk, propv)
    
        newicon.createDimension('TSTEP', None)
        for dimk in ('DATE-TIME', 'LAY', 'ROW', 'COL'):
            newicon.createDimension(dimk, len(oldicon.dimensions[dimk]))
        for vark, oldv in oldicon.variables.iteritems():
            if vark == 'TFLAG': continue
            newv = newicon.createVariable(vark, oldv.dtype.char, oldv.dimensions)
            for propk in oldv.ncattrs():
                propv = getattr(oldv, propk)
                setattr(newv, propk, propv)
        if len(newicon.dimensions['ROW']) != metcroprops['NROWS']:
            raise ValueError('Default (I,B)CON has different ROW dimension (%d) than the output file (%d).' % (len(newicon.dimensions['ROW']), metcroprops['NROWS']))
        if len(newicon.dimensions['COL']) != metcroprops['NCOLS']:
            raise ValueError('Default (I,B)CON has different COL dimension (%d) than the output file (%d).' % (len(newicon.dimensions['COL']), metcroprops['NCOLS']))
        

    minval = 1e-32
    addedkeys = []
    infiles = []
    if dobcon:
        infiles += [(newbcon, ('TSTEP', 'LAY', 'PERIM'))]
    if doicon:
        infiles += [(newicon, ('TSTEP', 'LAY', 'ROW', 'COL'))]

    for vark, varo in mappings.iteritems():
        for newcon, dims in infiles:
            if vark not in newcon.variables.keys():
                newv = newcon.createVariable(vark, 'f', dims)
                newv.units = varo['outunit'].ljust(16)
                newv.var_desc = newv.long_name = vark.ljust(16)
                addedkeys.append(vark)

    infiles = []
    if dobcon:
        infiles += [(metbdy, regridded_nd49_bdy, oldbcon, newbcon, nbconoutsteps)]
    if doicon:
        infiles +=  [(metcro, regridded_nd49_cro, oldicon, newicon, 1)]
    general_messages = messages
    messages = ""
    for metfile, regridded_nd49, oldcon, newcon, noutstep in infiles:
        if args.sigmaeta:
            cpress = np.convolve([0.5, 0.5], metfile.VGLVLS, mode = 'valid')
        else:
            cpressv = metfile.variables['PRES']
            assert(cpressv.units.strip() == 'Pa')
        oldkeys = []
        for vark, varo in newcon.variables.iteritems():
            toff = 0
            if vark in ('TFLAG',): continue
            if vark in mappings:
                sys.stdout.write(vark + ', ')
                sys.stdout.flush()
                expr, ounit = eval('expression, outunit', None, mappings[vark])
                # Temporary version for units
                nd49 = regridded_nd49[0]
                coexpr = compile(expr, 'expr', 'eval')
                inunits = [nd49.variables[cn].units.strip() for cn in coexpr.co_names if cn in nd49.variables]
                if len(inunits) == 0:
                    warn("Cannot evaluate any parts of '%s' = '%s'; using default file" % (vark, expr))
                    oldkeys.append(vark)
                    minout = getdefault(oldcon, vark, noutstep)
                else:
                    inC = dict([(cn, nd49.variables[cn].carbon) for cn in coexpr.co_names if cn in nd49.variables])
                    inkgpermole = dict([(cn, nd49.variables[cn].kgpermole) for cn in coexpr.co_names if cn in nd49.variables])
                    def conv(matcho):
                        found = {}
                        found.update(matcho.groupdict())
                        spc = found['spc']
                        if spc not in nd49.variables:
                            return spc
                        found['carbon'] = inC[spc]
                        found['kgpermole'] = inkgpermole[spc]
                        if 'micrograms' in ounit:
                            out =  '%(spc)s[:] / %(carbon)f * %(kgpermole)f' % found
                        elif found['carbon'] == 1.:
                            out =  '%(spc)s[:]' % found
                        else:
                            out =  '%(spc)s[:] / %(carbon)f' % found
                        return out
                    for cn in coexpr.co_names:
                        expr = re.compile(r'(?P<spc>\b' + cn + r'\b)').sub(conv, expr)
            
                    inunit = inunits[0]
                    manual_unit = mappings[vark].get('manual_unit', False)
                    ounit = ounit.strip()
                    if args.verbose:
                        sys.stdout.write('\n%s %s\n' % (expr, ounit))
                        sys.stdout.flush()
                    out = np.zeros((noutstep,) + varo[:].shape[1:], dtype = 'f')
                    for ndi, nd49 in enumerate(regridded_nd49):
                        if not (ounit == varo.units.strip()):
                            sys.stdout.write('\n%s %s %s\n' % (vark, varo.units.strip(), ounit))
                            sys.stdout.flush()
                        try:
                            temp_val = eval(expr, None, nd49.variables)[:]
                        except:
                            warn("Cannot evaluate part of '%s' = '%s'" % (vark, expr))
                            oldkeys.append(vark)
                            out = getdefault(oldcon, vark, noutstep)
                            break
                    
                        if manual_unit:
                            pass
                        else:
                            assert((np.array(inunits) == inunit).all())
                            if inunit == ounit:
                                pass
                            elif ounit == 'ppmV' and inunit in ('ppbv', 'ppbC'):
                                # ppbC has automatically been converted to ppbv
                                temp_val *= 1e-3
                            elif ounit == 'micrograms/m**3' and inunit == 'ppbv':
                                airmolden = eval(mappings_file['AIRMOLDEN']['expression'], None, nd49.variables)
                                temp_val *= airmolden
                            else:
                                raise ValueError('Error: in unit/outunit combo unknown: "%s", "%s"' % (inunit, ounit))
                        toff
                        oldrm = 0
                        for ti, temp_hour in enumerate(temp_val):
                            if args.sigmaeta:
                                rpress = ((nd49.variables['etam_pressure'] - newcon.VGTOP / 100) / (1013.25 - newcon.VGTOP / 100))
                                while rpress.ndim < out[0].ndim:
                                    rpress = np.expand_dims(rpress, -1)
                                rpress = rpress * np.ones((rpress.shape[0],) + out[0, 0].shape, dtype = 'f')
                            else:
                                cpress = cpressv[ti, :]
                                cpress = cpress.reshape(cpress.shape[0], -1)
                                assert((np.diff(cpress, axis = 0).mean(1) < 0).all())
                                rpressv = eval(mappings_file['PRESS']['expression'], None, nd49.variables)
                                if mappings_file['PRESS']['outunit'] == 'Pa':
                                    rpress = rpressv[ti,:]
                                elif mappings_file['PRESS']['outunit'] == 'hPa':
                                    rpress = rpressv[ti,:] * 100
                                else:
                                    warn("Assuming pressure is Pa, but got %s" % mappings_file['PSURF']['outunit'])
                                rpress = rpress.reshape(rpress.shape[0], cpress.shape[1])
                                assert((np.diff(rpress, axis = 0).mean(1) < 0).all())
                        
                            if args.sigmaeta:
                                while cpress.ndim < out[0].ndim:
                                    cpress = np.expand_dims(cpress, -1)
                                cpress = cpress * np.ones_like(out[0])
                        
                            outvals = np.zeros(cpress.shape, dtype = 'f')
                            for pi, (cp, rp, column) in enumerate(zip(cpress.T, rpress.T, temp_hour.T)):
                                outvals[:, pi] = np.interp(cp[::-1], rp[::-1], column[::-1])[::-1]
                            out[toff + ti, :, :] = outvals.reshape(*out.shape[1:])
                        toff = toff + ti + 1
                    if args.verbose:
                        sys.stdout.write('\n%s %.5e %.5e %.5e %.5e %.5e\n' % tuple(['Raw'] + np.percentile(out, [0, 25, 50, 75, 100]).tolist()))
                        sys.stdout.flush()
                    minout = np.maximum(out, minval)
                    if args.verbose:
                        sys.stdout.write('\n%s %.5e %.5e %.5e %.5e %.5e\n' % tuple(['Min'] + np.percentile(minout, [0, 25, 50, 75, 100]).tolist()))
                        sys.stdout.flush()
            else:
                oldkeys.append(vark)
                minout = getdefault(oldcon, vark, noutstep)
            varo[0:minout.shape[0], :, :] = minout[:]

        #varlist = getattr(newcon, 'VAR-LIST')
        #nvars = getattr(newcon, 'NVARS') + len(addedkeys)
        #setattr(newcon, 'VAR-LIST', varlist + ''.join([nk.ljust(16) for nk in addedkeys]))
        allkeys = [nk.ljust(16) for nk in newcon.variables.keys() if nk != 'TFLAG']
        nvars = len(allkeys)
        newvarlist = ''.join(allkeys)
        setattr(newcon, 'VAR-LIST', newvarlist)
        setattr(newcon, 'NVARS', nvars)
        newcon.createDimension('VAR', nvars)
        tflag = newcon.createVariable('TFLAG', 'i', ('TSTEP', 'VAR', 'DATE-TIME'))
        tflag.units = "<YYYYDDD,HHMMSS>" ;
        tflag.long_name = "TFLAG           " ;
        tflag.var_desc = "Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS                                " ;
        if args.timeindependent:
            newcon.SDATE = np.int32(-635)
            newcon.STIME = np.int32(0)
            newcon.TSTEP = np.int32(0)
            tflag[:] = 0.
        else:
            gtime = regridded_nd49[0].variables['time']
            gsdate = (datetime.strptime(gtime.units.split(' since ')[-1].replace(' UTC', ''), '%Y-%m-%d %H:%M:%S')  + timedelta(hours = int(gtime[0]))).strftime('%Y-%m-%d %H:%M:%S')
            if args.sdate is None:
                args.sdate = gsdate
            if args.tstep is None:
                try:
                    gtstep = gtime[1] - gtime[0]
                    args.tstep = gtstep * 10000
                except IndexError:
                    raise ValueError('Unable to determine timestep from file; use --tstep option to assign manually')
            sdate = datetime.strptime(args.sdate, '%Y-%m-%d %H:%M:%S')
            newcon.SDATE = np.int32(sdate.strftime('%Y%j'))
            newcon.STIME = np.int32(sdate.strftime('%H%M%S'))
            newcon.TSTEP = np.int32(args.tstep)
            sdate = datetime.strptime(str(newcon.SDATE) + ' ' + '%06d' % newcon.STIME, '%Y%j %H%M%S')
            dt = datetime.strptime('%06d' % newcon.TSTEP, '%H%M%S') - datetime(1900, 1, 1)
            tnow = sdate
            for ti in range(varo.shape[0]):
                jdate = int(tnow.strftime('%Y%j'))
                itime = int(tnow.strftime('%H%M%S'))
                tflag[ti, :, :] = np.array([[jdate, itime]]).repeat(newcon.NVARS, 0)
                tnow = tnow + dt
        newcon.sync()
        repair_ae(newcon, None)
        if len(addedkeys) > 0:
            warn('CON did not originally have:\n' + ', '.join(addedkeys) + '\n')
        if len(oldkeys) > 0:
            warn('Using old CON for:\n' + ', '.join(oldkeys) + '\n')
        setattr(newcon, 'HISTORY', getattr(newcon, 'HISTORY', '') + general_messages + messages)
        newcon.sync()

if __name__ == '__main__':
    from PseudoNetCDF.pncparse import getparser, pncparse
    
    parser = getparser(has_ofile = False, plot_options = False, interactive = False)
    parser.add_argument('--METBDY3D', default = None, type = str, help='path (or PseudoNetCDF commands) to a MCIP METBDY3D file')
    parser.add_argument('--METCRO3D', default = None, help='path (or PseudoNetCDF commands) to a MCIP METCRO3D file')
    parser.add_argument('--BCON', default = "dummybcon.nc", type = str, help='path (or PseudoNetCDF commands) to an old BCON file (default dummybcon.nc - created on the fly)')
    parser.add_argument('--ICON', default = "dummyicon.nc", help='path (or PseudoNetCDF commands) to an old ICON file')
    parser.add_argument('--timeindependent', default = False, action = 'store_true', help = 'Start date/time is any start date time')
    parser.add_argument('-d', '--sdate', default = None, help = 'Start date in YYYY-MM-DD HH:MM:SS format')
    parser.add_argument('--sigmaeta', action = 'store_true', help = 'Use sigma from CMAQ and calculate sigma for GEOS-Chem from pure eta levels for interpolation')
    parser.add_argument('--verbose', action = 'store_true', default = False, help = 'Add more printing details')
    parser.add_argument('--persist', dest = 'persistintermediate', action = 'store_true', default = False, help = 'Save interpolated files (may have speed benefits)')
    parser.add_argument('--tstep', default = None, help = 'Time increment between ND49 files')
    parser.add_argument('--mapping', default = 'mappings.json', help = 'Path to mappings file (i.e., json formatted dictionary); use --template to get a mapping')
    parser.add_argument('--outbfolder', default = ".", help = 'Path to output folder for BCON.')
    parser.add_argument('--outifolder', default = ".", help = 'Path to output folder for ICON.')
    parser.add_argument('--template', dest = 'template', default = None, choices = ['cb05tucl_ae6_aq', 'saprc07tc_ae6_aq'], type = str, help = 'Print template mappings.json to stdout and exit (still requires ND49 argument)')
    parser.add_argument('--gc-version', dest = 'gcversion', default = 9, choices = [8, 9, 10], type = int, help = 'GEOS-Chem version is used to determine SOA mapping.')
    parser.set_defaults(format = 'bpch,nogroup=True')
    parser.description = """pncgeos2cmaq makes boundary conditions from GEOS-Chem for CMAQ."""
    parser.epilog = """
Requirements:
    ifile - path to either a ND49 or TPCORE
    mapping - path to a table of mappings (see --template for defaults)
    one or both of:
        --METBDY3D - path to METBDY3D MCIP output file
        --METCRO3D - path to METCRO3D MCIP output file

Example:
    $ pncgeos2cmaq.py inputs/ts20120301.bpch --template=cb04tucl_ae6_aq > mappings.json
    $ pncgeos2cmaq.py inputs/ts20120301.bpch --METBDY3D inputs/METBDY3D_20120301 --METCRO3D inputs/METCRO3D_20120301 --mapping=mappings.json  --outifolder=icon/ --outbfolder=bcon/

"""
    
    args = parser.parse_args()
    if not args.template is None:
        print get_template(args.template, args.gcversion)
        exit()
    if not os.path.exists(args.mapping):
        parser.print_help()
        raise ValueError("mappings.json must exist and can be created using the --template option")
    ifiles, args = pncparse(has_ofile = False, parser = parser)
    args.ND49 = ifiles
    args.ND49PATH = args.ipath
    if args.METCRO3D is None and args.METBDY3D is None:
        parser.print_help()
        raise ValueError('METCRO3D or METBDY3D are required (or both)')
    args.ND49_REGRID_BDY = [os.path.join(args.outbfolder, os.path.basename(p).replace('.nc', '') + '.BDY.nc') for p in args.ND49PATH]
    args.ND49_REGRID_CRO = [os.path.join(args.outifolder, os.path.basename(p).replace('.nc', '') + '.CRO.nc') for p in args.ND49PATH]
    args.NEWBCON = os.path.join(args.outbfolder, os.path.basename(args.ND49PATH[0].replace('.nc', '') + '.BCON.nc'))
    args.NEWICON = os.path.join(args.outifolder, os.path.basename(args.ND49PATH[0].replace('.nc', '') + '.ICON.nc'))
    makebcon(args)