# Calculate Mean Emissions of Formaldehyde from a CAMx emissions file #

Step 1: [Download the "other" inputs from the CAMx testcase](http://camx.com/down/testcase.php)

Step 2: untar the outputs

Step 3: open a python terminal and type the follow (the path may be different)

```
path = "emiss.stl.12kmsmall.20020603.a1.bin"
from PseudoNetCDF.camxfiles.Memmaps import uamiv

emiss_file = uamiv(path)
form = emiss_file.variables['FORM']

print "Formaldehyde mean",form.units,form[:].mean()
print "Formaldehyde Hourly mean", form.units
for form_hour in form:
    print form_hour.mean()
```