# AGNREF
Composite XSPEC model, consisting of a disc geometry as in AGNSED (Kubota &amp; Done 2018) and an additional outflow component giving an additional thermal re-processor component and reflected component. If you use this in your work please reference Hagen & Done (in prep.)


Requirements
--------------
* Working installation of Heasoft (including XSPEC!) This has been tested on Heasoft v.6.29 and v.6.30, with XSPEC v.12.12.0 and v.12.12.1


Installation
--------------
1. Clone this repository (or just download the source code. You'll need the agnref.f and agnref.dat files)
2. cd into the src directore (or if you've put the relevant file somewhere else cd into there).
3. Open xspec
4. Run: `initpackage agnref agnref.dat .` in xspec (this will compile the code)
5. Now load the model with: `lmod agnref .` in xspec

If you don't want to type lmod everytime you open xspec, you can cd into ~/.xspec and edit/create and xspec.rc file. In this file put the line
`lmod agnref /path/to/agnref`
Xspec will now load the model automatically upon opening

  
Model Description
------------------
For a full descripton see section 4.3 in Hagen & Done (in prep.).
However, briefly: The model considers a geometry as outlined in Kubota & Done (2018) for AGNSED (i.e a standard outer disc, a warm Comptonizing region where the disc fails to thermalise, and a hot Comptonizing corona giving the X-ray emission). We then also consider a bi-conical outflow, which gives a reflected component (modeled with rdblur*pexmon (Fabian et al 1989; Nandra et al 2007)), and a re-processed thermal component (modelled with bbody). The relative normalisations of these components are set by the fraction of X-ray power intercepted by the outflow.
