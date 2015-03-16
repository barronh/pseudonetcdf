# Increase Emissions of Formaldehyde from a CAMx emissions file #

Step 1: [Download the "other" inputs from the CAMx testcase](http://camx.com/down/testcase.php)

Step 2: untar the outputs

Step 3: copy emiss.stl.12kmsmall.20020603.a1.bin to emiss.stl.12kmsmall.20020603.a1.bin.copy
Step 4: open a python terminal and type the follow (the path may be different)

```
path = "emiss.stl.12kmsmall.20020603.a1.bin.copy"
from PseudoNetCDF.camxfiles.Memmaps import uamiv

emiss_file = uamiv(path,'r+')
form = emiss_file.variables['FORM']

print "Formaldehyde mean",form.units,form[:].mean()
form[:] = form[:] * 2. # double formaldehyde
emiss_file.sync()
print "Formaldehyde mean",form.units,form[:].mean()
```


Wow!  That was easy.