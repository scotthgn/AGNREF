# AGNREF
Composite XSPEC model, consisting of a disc geometry as in AGNSED (Kubota &amp; Done 2018) and an additional outflow component giving an additional thermal re-processor component and reflected component. If you use this in your work please reference Hagen & Done (2022, submitted).
https://ui.adsabs.harvard.edu/abs/2022arXiv221004924H/abstract


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


Model Parameters
-----------------
Par 1 : Mass - Black hole mass in solar masses <br/>
Par 2 : Dist - Co-Moving distance in Mpc <br/>
Par 3 : $\log(\dot{m})$ - Mass accretion rate, where $\dot{m} = \dot{M}/\dot{M}_{\text{edd}}$ is in units of the Eddington mass accretion rate<br/>
Par 4 : astar - Dimensionless black hole spin <br/>
Par 5 : cos_inc - cosine of the inclination angle, measured from z-axis with the disc in the x-y plane <br/>
Par 6 : kTe_hot - Electron temperature of hot Comptonisation region in keV - if negative, then only returns hot Compton component <br/>
Par 7 : kTe_warm - Electron temperature of warm Comtonisation region in keV - if negative, then only returns warm Compton component <br/>
Par 8 : gamma_hot - Photon index for hot Comptonisation region <br/>
Par 9 : gamma_warm - Photon index for warm Comptonisation region - if negative, then only returns the outer disc component <br/>
Par 10 : r_hot - Outer radius of hot Comptonisation region in Rg <br/>
Par 11 : r_warm - Outer radius of warm Comptonisation region in Rg <br/>
Par 12 : log_rout - Outer radius of the outer disc in Rg - If negative, then uses the self-gravity radius, calculated as in Laor and Netzer (1989) <br/>
Par 13 : hmax - Scale height of the corona in Rg <br/>
Par 14 : fcov - Covering fraction of the outflow, as seen by the black hole, in $\Omega/4\pi$ <br/>
Par 15 : kT_wind - Black-body temperature of wind/outflow component in keV <br/>
Par 16 : Awind - Albedo of wind/outflow <br/>
Par 17 : inc_blue - flag, 1=>include blurring of reflected component (using rdblur), 2=>no blurring <br/>
Par 18 : rin_blur - Convineince parameter. Comes from the way rdblur is written (assumes a standard disc in Keplerian orbit). Needed to determine the 'amount' of blurring <br/>
Par 19 : Redshift <br/>
par 20 : Normalisation - should be fixed to 1!! - The code calculates energetics self-consistently <br/>
