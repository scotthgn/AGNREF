c      Fortran version of AGNREF
c      Ref. Hagen and Done (in prep.)
c
c      Combines AGNSED (Kubota & Done 2018) with PEXMON (Nandra er al.
c      2007), and a thermal re-processor component. We assume a geometry
c      where both the reflection and additional thermal component
c      originate from the same region - taken to be some smooth wind
c      covering a fraction of the sky. Hence these components are tied
c      together.
c      We also tie the normalisation of these to AGNSED, by assuming
c      that all reflection/re-processing comes from illumination by the
c      central X-ray corona. Hence we can write
c
c      Lref = A*Lx*fc*0.5
c      Lth = (1-A)*Lx*fc*0.5
c
c      where Lx is the luminosity of the central X-ray corona, and fc is
c      the covering fraction (in units d\Omega/4\pi) of the
c      reflector/re-processor. Note that Lx is calculated
c      self-consistently as in AGNSED
c      The factor 0.5 comes from the fact that we only see the wind
c      emission from one side - however the covering fraction includes
c      both sides
c
c      The model consists of four regions: a standard outer disc, a warm
c      Comptonization region where the disc has not thermalised properly,
c      a spherical inner corona modelled as a hot Comptonization region,
c      and a distant re-processor/reflector
c 
c      The Comptonization is modelled using NTHCOMP (Zdziarski, Johnson &
c      Magdziarz 1996; Zycki, Done & Smith 1999)
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      
       subroutine agnref(ear,ne,param,ifl,photar,photer)
       implicit none

       integer npars
       parameter(npars=19)
      
       integer ne, ifl
       real ear(0:ne), photar(ne), param(npars), photer(ne)
       real oldpar(npars)

       integer Nnew !Defining own energy grid for model calc.
       parameter(Nnew=2000) !Number of bin edges
       real enew(0:Nnew), ph(Nnew) !ph is photon array on new grid
       real dloge, newemin, newemax

       real fstart(ne), fend(ne)
       integer istart(ne), iend(ne)

       logical parchange, echange
       save oldpar

       integer i, n             !Iteration indeces
       double precision zfac

c      ear: Energy array
c      ne: Size of energy array
c      param: Parameter value array
c      ifl: Spectrum number
c      photar: flux array (size ne-1)
c
c      param(1):  Mass, Msol
c      param(2):  Distance, Mpc
c      param(3):  log mdot, Mass accretion rate,  L/Ledd
c      param(4):  astar, BH spin
c      param(5):  cos inc, Inclination
c      param(6):  kTe_hot, hot corona temperature, keV
c      param(7):  kTe_warm, warm compton temperature, keV
c      param(8):  gamma_hot, hot corona photon index
c      param(9):  gamma_warm, warm compton photon index
c      param(10): r_hot, hot corona radius, Rg
c      param(11): r_warm, warm compton outer radius, Rg
c      param(12): log r_out, outer disc radius, Rg
c      param(13): hmax, scale height of corona, Rg
c      param(14): fcov, Covering fraction of re-processor, dOmega/4pi
c      param(15): kT_wind, Temperature of re-processor, keV
c      param(16): Awind, Albedo of reprocessor
c      param(17): blur, 0=off/1=on - whether to convolve pexmon with rdblur
c      param(18): rin_blur, effective radius for rdblur, Rg
c      param(19): redshift
       

c      checking if parameters have changed
       parchange=.false.
       do i=1,npars,1
          if (param(i).ne.oldpar(i)) parchange=.true.
       end do
       
c      Checking if we need to change energy grid
c      The default extends from 1e-4 to 1e3 keV
       zfac = (1.0 + param(19))
       if (ear(0).eq.0.0) then
          newemin=ear(1) - (ear(2)-ear(1))/10.0
       else
          newemin=min(1.0e-4, ear(0)*zfac)
       end if
       newemax=max(1.0e3, ear(ne)*zfac)

       if ((enew(0).ne.newemin).or.(enew(Nnew).ne.newemax)) then
          echange=.true.
       end if
       
       
c      Calculating new energy grid if necessary
       if (echange) then
          dloge=log10(newemax/newemin)/float(Nnew)
          enew(0) = newemin
          do n=1, Nnew, 1
             enew(n)=10**(log10(enew(0))+dloge*float(n))
          end do
       end if

c      Call model if parameters or energy grid have changed
       if (parchange.or.echange) then
          
          call calc_refspec(enew, Nnew, param, ifl, ph)

          !Redshift correct energy bins
          do n=1, Nnew, 1
             enew(n) = enew(n)/zfac
             ph(n) = ph(n)/zfac
          end do

          do i=1, npars, 1
             oldpar(i)=param(i)
          end do
       end if

c      Re-bin onto original xspec grid
       call inibin(Nnew, enew, ne, ear, istart, iend, fstart, fend, 0)
       call erebin(Nnew, ph, ne, istart, iend, fstart, fend, photar)

       return
       end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c      Main subroutine for calculating disc spec, warm Comp spec, and
c      hot Comp spec
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


       subroutine calc_refspec(es, nn, param, ifl, ph)

c      Function to calculate the SED. Calculates the spectral
c      contribution from each region (disc, warm comp, hot comp,
c      and reprocessor/outflow)
c      Temperature found through the Novikov-Thorne equations
c       


       implicit none
       integer nn, ifl
       real es(0:nn), ph(nn), param(*)
       real ph_ann(nn), phw_ann(nn), phh_ann(nn), phh_shape(nn) !disc, warm, hot
       real ph_ref(nn), ph_rep(nn)
       double precision M, mdot
       double precision rout, rw, rh, rsg, risco
       double precision astar, hmax
       double precision kTh, kTw
       double precision gammah, gammaw
       double precision cos_inc, dist
       double precision fcov_wind, kTwind, Awind, rin_blur
       integer inc_blur

       double precision pi, G, c, h, sigma_sb,  mp, kB
       double precision eta, Ledd, Mdot_edd, Rg
       double precision dr_dex, nrd, nrw, nrh, nrseed
       double precision dlog_rd, dlog_rw, dlog_rh, dlog_rseed
       double precision rmid, dr
       double precision kkev, kevhz, keverg
       double precision dflux, en, dnphot, tr, tw
       double precision calc_risco, calc_rsg, calc_fcol
       double precision efficiency, nttemp, repnttemp
       double precision lphseed, fcov, theta0
       double precision lseed_ann, lseed_cann, ldiss_ann
       double precision tseed_h, ysb
       double precision Lx_tot, fx_tot

       real pxpar(12), ntwpar(5), nthpar(5), rdpar(4), bbpar(2)
       real pxpherr(nn), ntwpherr(nn), nthpherr(nn)
       real phref_err(nn), phrd_err(nn), phbb_err(nn)
       double precision dflux_rep

       double precision normw_ann, normh_ann, norm_rep, norm_ref
       double precision ntw_out, nth_out, px_out, rep_out

       logical return_disc, return_warm, return_hot
       logical return_ref, return_rep
       integer i, n

       
c      Constants
       pi = 4.0*atan(1.0)
       G = 6.67d-8 * 1.99d33   !Grav const. cm^-3 s^-1 Msol^-1
       c = 3.0d10              !Speed of light. cm/s
       h = 6.62617d-27         !Plank const. erg s
       sigma_sb = 5.670367d-5  !Stefan-Boltzmann const. erg s^-1 cm^-2 K^-4
       mp = 1.67d-24           !Proton mass. g
       kB = 1.38d-16           !Boltzmann const. erg/K

c      Unit conversion constants
       kkev = 1.16048d7          !K/keV
       kevhz = 2.417965d17      !Hz/keV
       keverg = 1.60218d-9    !erg/keV

       dr_dex = 30              !Radial resolution, Nbins per decade (ish)
      

c      Read parameters
       M = dble(param(1))
       dist = dble(param(2)) !Mpc
       mdot = dble(10**(param(3)))
       astar = dble(param(4))
       cos_inc = dble(param(5))
       kTh = dble(abs(param(6)))
       kTw = dble(abs(param(7)))
       gammah = dble(param(8))
       gammaw = dble(abs(param(9)))
       rw = dble(param(11))
       hmax = dble(param(13))
       fcov_wind = dble(param(14))
       kTwind = dble(abs(param(15)))
       Awind = dble(abs(param(16)))
       inc_blur = int(param(17))
       rin_blur = dble(param(18))
       
       
c      Getting accretion attributes
       risco = calc_risco(astar)
       eta = efficiency(risco)
       rsg = calc_rsg(M, mdot)
       
       Ledd = 1.39d38 * M       !erg/s
       Mdot_edd = Ledd/(eta * c**2) !g/s
       Rg = (G*M)/c**2          !cm
       
       dist = dist*1.0d6        !pc
       dist = dist*3.086d18    !cm
       
c      Filling parameter arrays (for those that dont change!!!)       
c      warm nthcomp params
       ntwpar(1) = gammaw
       ntwpar(2) = kTw
       ntwpar(4) = 0.0 !BB assumed from each annulus
       ntwpar(5) = 0.0 !Redshift - dealt with within this script

c      hot nthomp params
       nthpar(1) = gammah
       nthpar(2) = kTh
       nthpar(4) = 0.0
       nthpar(5) = 0.0

c      pexmon pars
       pxpar(1) = gammah       !photon index
       pxpar(2) = 1.0d5        !fold E, keV
       pxpar(3) = -1           !return only reflected
       pxpar(4) = 0.0          !redshift, delt with elsewhere
       pxpar(5) = 1.0          !Abundances, set to solar
       pxpar(6) = 1.0          !Fe abundance, set to solar
       pxpar(7) = acos(cos_inc) * (180.0/pi) !Inclination, deg
       !pexmon normalisation will be tied to Lhot, later in script

c      rdblur pars
       rdpar(1) = -3.0          !Power law index
       rdpar(2) = rin_blur      !Effective bluring radius
       rdpar(3) = 1.0d6         !Some large radius...
       rdpar(4) = acos(cos_inc) * (180.0/pi)

c      Black-Body pars
       bbpar(1) = kTwind
       bbpar(2) = 1.0
       
c      Checking switching parameters
       return_disc=.true.
       return_warm=.true.
       return_hot=.true.
       return_ref=.true.
       return_rep=.true.
       
       if (param(12).lt.0.0) then !checkig outer radius
          rout = rsg
       else
          rout = dble(10**param(12))
       end if

       if ((param(10).lt.0.0).or.(param(10).lt.risco)) then !Checking inner disc radius
          rh = risco
          return_hot=.false.
          call xwrite('r_hot < risco!!! Setting r_hot = risco', 10)
       else
          rh = dble(param(10))
       end if

       if ((param(11).lt.0.0).or.(param(11).lt.risco)) then !Checking warm region
          rw = risco
          return_warm=.false.
          return_hot=.false.
          call xwrite('r_warm < risco!! Setting rw = risco', 10)
       else
          rw = dble(param(11))
       end if

       
       if (param(9).lt.0.0) then !if gammaw<0 show ONLY disc
          return_warm=.false.
          return_hot=.false.
          return_ref=.false.
          return_rep=.false.
       end if

       if (param(7).lt.0.0) then !if kTw<0 show ONLY warm Comp
          return_disc=.false.
          return_hot=.false.
          return_ref=.false.
          return_rep=.false.
       end if

       if (param(6).lt.0.0) then !if kTh<0 show ONLY hot comp
          return_disc=.false.
          return_warm=.false.
          return_ref=.false.
          return_rep=.false.
       end if

       if (param(15).lt.0.0) then !if kTwind<0 show ONLY thermal component
          return_disc=.false.
          return_warm=.false.
          return_hot=.false.
          return_ref=.false.
       end if

       if (param(16).lt.0.0) then !if Awind<0 show ONLY reflected component
          return_disc=.false.
          return_warm=.false.
          return_hot=.false.
          return_rep=.false.
       end if
         

       if (rw.le.rh) then
          call xwrite('r_warm <= r_hot!!! Re-setting r_w = r_h', 10)
          rw = rh
          return_warm=.false.
       end if

       if (rw.ge.rout) then
          call xwrite('r_warm >= r_out!!! Re-setting r_w = r_out', 10)
          rw = rout
          return_disc=.false.
       end if


       if (hmax.gt.rh) then
          call xwrite('hmax > r_hot!!! Re-setting hmax = r_hot', 10)
          hmax = rh
       end if


c-----------------------------------------------------------------------
c      Section for calculating hot compton region (ie corona)
c-----------------------------------------------------------------------
c      First finding the total seed photon luminosity/flux
c      Integrated flux from entire disc seen by corona
       
c      This section gets calculated even if return_hot false, since
c      we need the xray luminosity for normalisation of the reflected
c      and re-processed components + disc temperature.
c      The switch for including hot component in SED is instead further
c      down the script
       nrseed = (log10(rout) - log10(rh)) * dr_dex
       nrseed = ceiling(nrseed) !round up to nearest integer
       dlog_rseed = (log10(rout) - log10(rh))/nrseed

       lphseed = 0.0
       do i=1, int(nrseed), 1
          rmid = 10**(log10(rh)+float(i-1)*dlog_rseed+dlog_rseed/2.0)
          dr = 10**(log10(rmid) + dlog_rseed/2.0)
          dr = dr - 10**(log10(rmid) - dlog_rseed/2.0)

          if ((rmid+dr/2.0+dr).gt.rout) then !Ensuring bin within region
             dr = rout - 10**(log10(rmid) - dlog_rseed/2.0)
          end if

          if (hmax.le.rmid) then
             theta0 = asin(hmax/rmid)
             fcov = theta0 - 0.5*sin(2.0*theta0) !corona covering fraction seen from rmid
          else
             fcov = 0.0
          end if

          tr = nttemp(rmid, M, mdot, Mdot_edd, astar, risco) !Temp. K

          lseed_ann = sigma_sb * tr**4.0 !erg/s/cm^2
          lseed_ann = lseed_ann * 4*pi*rmid*dr*Rg**2 * (fcov/pi) !erg/s
          lseed_ann = lseed_ann/(4*pi * dist**2) !erg/s/cm^2

          lphseed = lphseed + lseed_ann !total seed photon flux erg/s/cm^2
       end do

c      Now finding seed photon temperature
c      Assumed to be ~inner disc temp
       tseed_h = nttemp(rh, M, mdot, Mdot_edd, astar, risco) !K
       if (rh.lt.rw) then
          ysb = (gammaw*(4.0/9.0))**(-4.5) !Compton y-param for warm region
          tseed_h = (tseed_h/kkev) * exp(ysb) !in keV
       else
          tseed_h = (tseed_h/kkev) 
       end if


          
       !calling nthcomp now, since intrinsic shape does not change!!!
       nthpar(3) = tseed_h
       call donthcomp(es, nn, nthpar, ifl, phh_shape, nthpherr) 
       !total flux output from nthcom
       nth_out = 0.0
       do n=1, nn, 1
          nth_out = nth_out + phh_shape(n)*es(n)*kevhz*h !erg/s/cm^2
       end do

c      Now calculating emission from each coronal annulus
       nrh = (log10(rh) - log10(risco)) * dr_dex
       nrh = ceiling(nrh)    !round up to integer
       dlog_rh = (log10(rh) - log10(risco))/nrh !actual spacing

       fx_tot = 0.0 !total xray flux - erg/s/cm^2
       do i=1, int(nrh), 1
          rmid = 10**(log10(risco)+float(i-1)*dlog_rh+dlog_rh/2.0)
          dr = 10**(log10(rmid) + dlog_rh/2.0)
          dr = dr - 10**(log10(rmid) - dlog_rh/2.0)

          if ((rmid+dr/2.0+dr).gt.rh) then !Ensuring bin within region
             dr = rh - 10**(log10(rmid) - dlog_rh/2.0)
          end if

          tr = nttemp(rmid, M, mdot, Mdot_edd, astar, risco) !K
          ldiss_ann = sigma_sb * tr**4.0
          ldiss_ann = ldiss_ann * 4.0*pi*rmid*dr * Rg**2.0 !erg/s
          ldiss_ann = ldiss_ann/(4.0*pi * dist**2.0) !erg/s/cm^2

          !Assuming equall amounts of emission from each coronal annulus
          lseed_cann = lphseed/(nrh-1.0) 
          normh_ann = ldiss_ann + lseed_cann !erg/s/cm^2

          fx_tot = fx_tot + normh_ann !erg/s/cm^2
          !Now applying normalisation
          do n=1, nn, 1
             if (nth_out.eq.0) then
                phh_ann(n) = 0.0
             else
                phh_ann(n) = phh_shape(n) * (normh_ann/nth_out) !photons/s/cm^2
             end if
          end do

          if (return_hot) then
c            Adding to total output array
             do n=1, nn, 1
                if (i.eq.1) then
                   ph(n) = phh_ann(n) !if no disc or warm comp start from scratch
                else
                   ph(n) = ph(n) + phh_ann(n)
                end if
             end do
          end if
       end do
       Lx_tot = fx_tot * 4.0*pi*dist**2 !Xray lum - erg/s
       

c-----------------------------------------------------------------------
c      Section for calculating the disc region
c-----------------------------------------------------------------------
       if (return_disc) then
          nrd = (log10(rout) - log10(rw)) * dr_dex !nr of bins in region
          nrd = ceiling(nrd) !round up to nearest integer value
          dlog_rd = (log10(rout) - log10(rw))/nrd
          do i=1, int(nrd), 1
             rmid = 10**(log10(rw)+float(i-1)*dlog_rd+dlog_rd/2.0)
             dr = 10**(log10(rmid)+dlog_rd/2.0)
             dr = dr - 10**(log10(rmid)-dlog_rd/2.0)

             if ((rmid+dr/2.0 + dr).gt.rout) then !Makes big bin at edge
                dr = rout - 10**(log10(rmid)-dlog_rd/2.0)
             end if

             tr=repnttemp(rmid,M,mdot,Mdot_edd,astar,risco,hmax,Lx_tot) !Temp. K
             
c            Calculating BB emission over each energy
             do n=1, nn, 1
                en = dble(log10(es(n-1)) + log10(es(n)))
                en = en/2.0 !Evaluate in mid bin
                en = 10**en
                
                if (en.lt.30.0*tr) then
                   dflux = (pi*2.0 * h * (en*kevhz)**3.0)/(c**2.0)
                   dflux = dflux * 1.0/(exp((h*en*kevhz)/(kB*tr)) - 1.0)
                   dflux = dflux * 4.0*pi*rmid*dr*(cos_inc/0.5) *Rg**2.0 !multiplying by area

                   dflux = dflux/(4.0*pi*dist**2.0) !erg/s/cm^2/Hz
                   
                   ph_ann(n) = dflux/(h*kevhz*en) !photons/s/cm^2/Hz
                   ph_ann(n) = ph_ann(n)*kevhz*(es(n)-es(n-1)) !photons/s/cm^2                
                else
                   ph_ann(n) = 0.0d0
                end if         
             end do


c            Adding into main photon array
             do n=1, nn, 1
                if (i.eq.1) then
                   if (return_hot) then
                      ph(n) = ph(n) + ph_ann(n)
                   else
                      ph(n) = ph_ann(n) !if no hot start from scratch
                   end if
                else
                   ph(n)=ph(n)+ph_ann(n)
                end if
             end do 
          end do
       end if


c-----------------------------------------------------------------------
c      Section for calculating the warm Compton region
c-----------------------------------------------------------------------
       if (return_warm) then
          nrw = (log10(rw) - log10(rh)) * dr_dex !nr bins in region
          nrw = ceiling(nrw)    !round up to nearest integer
          dlog_rw = (log10(rw) - log10(rh))/nrw !actual bin size

          do i=1, int(nrw), 1
             rmid = 10**(log10(rh)+float(i-1)*dlog_rw + dlog_rw/2.0)
             dr = 10**(log10(rmid) + dlog_rw/2.0)
             dr = dr - 10**(log10(rmid) - dlog_rw/2.0)

             if ((rmid+dr/2.0+dr).gt.rw) then !Ensuring bin within region
                dr = rw - 10**(log10(rmid) - dlog_rw/2.0)
             end if

             tr=repnttemp(rmid,M,mdot,Mdot_edd,astar,risco,hmax,Lx_tot) !Temp. K
             ntwpar(3) = tr/kkev !setting seed temp to disc annulus temp (keV)

c            Calling nthcomp
             call donthcomp(es, nn, ntwpar, ifl, phw_ann, ntwpherr) !photons/s/cm^s/

             normw_ann = sigma_sb*tr**4.0
             normw_ann = normw_ann*4.0*pi*rmid*dr*(cos_inc/0.5)*Rg**2.0 !erg/s/, emission from annulus
             normw_ann = normw_ann/(4.0*pi*dist**2.0) !erg/s/cm^2/
             
c            Finding total flux output from nthcomp
             ntw_out = 0.0
             do n=1, nn, 1
                ntw_out = ntw_out + phw_ann(n)*es(n)*kevhz*h !erg/s/cm^2
             end do

c            Now applying normalisation
             do n=1, nn, 1
                if (ntw_out.eq.0) then
                   phw_ann(n) = 0.0
                else
                   phw_ann(n) = phw_ann(n) * (normw_ann/ntw_out) !photons/s/cm^2
                end if
             end do
             
c            Adding to total output array
             do n=1, nn, 1
                if(i.eq.1) then
                   if ((return_disc).or.(return_hot)) then 
                      ph(n) = ph(n) + phw_ann(n)
                   else
                      ph(n) = phw_ann(n) !if no disc or hot we start from scratch!
                   end if
                else
                   ph(n) = ph(n) + phw_ann(n)
                end if
             end do     
          end do
       end if


c-----------------------------------------------------------------------
c      Section for calculating reflected component (i.e pexmon)
c-----------------------------------------------------------------------
       if (return_ref) then
          !Calling pexmon - no longer splitting into annuli
          call pexmon(es, nn, pxpar, ifl, ph_ref, phref_err)

          !Applying rdblur if included
          if (inc_blur.eq.1) then
             call rdblur(es, nn, rdpar, ifl, ph_ref, phrd_err)
          end if
          
          !Finding total pexmon output
          px_out = 0.0
          do n=1, nn, 1
             px_out = px_out + ph_ref(n)*es(n)*kevhz*h !erg/s/cm^2
          end do
 
          !Applying normalisation - Lref=Awind*Lx*fcov_wind
          norm_ref = Awind*fx_tot*fcov_wind*0.5 !erg/s/cm^2
          do n=1, nn, 1
             if (px_out.eq.0.0) then
                ph_ref(n) = 0.0
             else
                ph_ref(n)=ph_ref(n) * (norm_ref/px_out)
             end if
          end do
  
          !Adding to total output array
          do n=1, nn, 1
             if ((return_disc).or.(return_warm).or.(return_hot)) then
                ph(n) = ph(n) + ph_ref(n)
             else
                ph(n) = ph_ref(n) !if only reflected component
             end if
          end do
       end if


c-----------------------------------------------------------------------
c      Section for calculating re-processed thermal component
c-----------------------------------------------------------------------
       if (return_rep) then
          !Calculating black-body shape of wind
          call bbody(es, nn, bbpar, ifl, ph_rep, phbb_err)

          !Finding black-body output flux
          rep_out = 0.0
          do n=1, nn, 1
             rep_out = rep_out + ph_rep(n)*es(n)*kevhz*h !erg/s/cm^2
          end do

          !Applying normalisation
          norm_rep = (1.0-Awind) * fcov_wind * fx_tot*0.5 !erg/s/cm^2
          do n=1, nn, 1
             if (rep_out.eq.0.0) then
                ph_rep(n) = 0.0
             else
                ph_rep(n) = ph_rep(n) * (norm_rep/rep_out)
             end if
          end do

          !Adding to total output array
          do n=1, nn, 1
             if ((return_disc).or.(return_warm).or.(return_hot)) then
                ph(n) = ph(n) + ph_rep(n)
             else if (return_ref) then
                ph(n) = ph(n) + ph_rep(n)
             else
                ph(n) = ph_rep(n) !if only reprocessed component
             end if
          end do
       end if
       
       end
             






       
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c      Functions for calculating disc properties
c      e.g risco, eta (efficiency), rsg, disc T, etc.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

       function calc_risco(astar)
c      Calculates the innermost stable circular orbit/marginally stable
c      orbit
       implicit none

       double precision astar, calc_risco
       double precision Z1_1, Z1_2, Z1, Z2
       
       Z1_1 = (1.0 - astar**2.0)**(1.0/3.0)
       Z1_2 = (1.0 + astar)**(1.0/3.0) + (1.0 - astar)**(1.0/3.0)
       Z1 = 1.0 + Z1_1 * Z1_2

       Z2 = sqrt(3.0 * astar**2.0 + Z1**2.0)
    
       if(astar.ge.0.0) then
          calc_risco = 3.0 + Z2 - sqrt((3.0-Z1) * (3.0+Z1+2.0*Z2))
       else
          calc_risco = 3.0 + Z2 + sqrt((3.0-Z1) * (3.0+Z1+2.0*Z2))
       end if
       
       return
       end


       function efficiency(risc)
c      Calculates the accretion efficiency, s.t L_bol=eta*Mdot*c^2
c      Uses GR case where eta = 1 - sqrt(1 - 2/(3*risco))
c      Taken from: The Physics and Evolution of Active Galactic Nuclei,
c      H. Netzer, 2013, p.38
       
       implicit none
       double precision risc, efficiency

       efficiency = 1.0 - sqrt(1.0 - 2.0/(3.0*risc))
       return
       end
       

       function calc_rsg(M, mdot)
c      Calculates the self gravity radius following Laor & Netzer 1989
c
c      Note! Assumes \alpha=0.1

       implicit none
       double precision M, mdot, calc_rsg
       double precision m9, alpha

       alpha = 0.1
       m9 = M/1.0d9
       calc_rsg = 2150*m9**(-2.0/9.0)*mdot**(4.0/9.0)*alpha**(2.0/9.0)

       return
       end


       function NTpars(r, astar, risc)
c      Calculates the Novikov-Thorne parameters at a given radius

       implicit none
       double precision r, astar, risc
       double precision pi, y, yisc, y1, y2, y3
       double precision B, C1, C2, C, NTpars
       double precision C2_1, C2_2, C2_3

       pi = 4.0*atan(1.0)

       y = sqrt(r)
       yisc = sqrt(risc)
       y1 = 2.0*cos((acos(astar) - pi)/3.0)
       y2 = 2.0*cos((acos(astar) + pi)/3.0)
       y3 = -2.0*cos(acos(astar)/3.0)

       B = 1.0 - (3.0/r) + ((2.0*astar)/(r**(3.0/2.0)))

       C1 = 1 - (yisc/y) - ((3.0*astar)/(2.0*y)) * log(y/yisc)

       C2_1 = (3.0*(y1-astar)**2.0)/(y*y1 * (y1-y2) * (y1-y3))
       C2_1 = C2_1 * log((y-y1)/(yisc-y1))
       
       C2_2 = (3.0*(y2-astar)**2.0)/(y*y2 * (y2-y1) * (y2-y3))
       C2_2 = C2_2 * log((y-y2)/(yisc-y2))

       C2_3 = (3.0*(y3-astar)**2.0)/(y*y3 * (y3-y1) * (y3-y2))
       C2_3 = C2_3 * log((y-y3)/(yisc-y3))

       C2 = C2_1 + C2_2 + C2_3
       C = C1 - C2
       NTpars = C/B

       return
       end


       function nttemp(r, M, mdot, Mdot_edd, astar, risc)
c      Function to calculate Novikov-Thorne temperature at radius r

      
       implicit none

       double precision r, M, mdot, Mdot_edd, astar, risc
       double precision pi, G, sigma_sb, c, Rg, Rt
       double precision NTpars
       double precision nttemp

       pi = 4.0*atan(1.0)
       G = 6.67d-8 * 1.99d33 !cm^3 s^-1 Msol^-1
       sigma_sb = 5.670367d-5 !erg cm^-2 s^-1 K^-4
       c = 3.0d10 !cm/s

       Rg = (G*M)/c**2.0 !cm
       Rt = NTpars(r, astar, risc)

       nttemp = (3.0*G*M*mdot*Mdot_edd)/(8*pi*sigma_sb*(r*Rg)**3) !K^4
       nttemp = nttemp*Rt
       nttemp = nttemp**(1.0/4.0) !K

       return
       end

      
       function repnttemp(r, M, mdot, Mdot_edd, astar, risc, hmax, Lx)
c      Disc temperature inclusing re-processing

       implicit none

       double precision r, M, mdot, Mdot_edd, astar, risc, hmax, Lx
       double precision pi, G, sigma_sb, c, Rg, Rt
       double precision repnttemp, nttemp, reptemp
       double precision frep, Rphys, Hphys

       pi = 4.0*atan(1.0)
       G = 6.67d-8 * 1.99d33 !cm^3 s^-1 Msol^-1
       sigma_sb = 5.670367d-5 !erg cm^-2 s^-1 K^-4
       c = 3.0d10               !cm/s

       Rg = (G*M)/c**2.0  !cm
       Rphys = r*Rg !cm
       Hphys = hmax*Rg !cm
       
       frep = Lx/(4.0*pi * (Rphys**2.0 + Hphys**2.0)) !erg/s/cm^2
       frep = frep * Hphys/sqrt(Rphys**2.0 + Hphys**2.0)
       frep = frep * (1.0-0.3)    !Disc albedo assumed as 0.3

       reptemp = frep/sigma_sb

       repnttemp = nttemp(r, M, mdot, Mdot_edd, astar, risc)
       repnttemp = repnttemp**4.0
       repnttemp = repnttemp + reptemp
       repnttemp = repnttemp**(1.0/4.0) !K

       return
       end
       
