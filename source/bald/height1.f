c@height1.tex
c--------1---------2---------3---------4---------5---------6---------7-c
c
        subroutine height1 (
     &   rminor,   rmajor,   kappa95,  delta95,
     &   current,  btor,     dense,    hydmass,
     &   zeff,     q,        t_guess,  t_ped)
 
c Revision History
c ----------------
c      date          description
c
c    28 June 2001    First Version by Thawatchai Onjun
c
c..Inputs:
c
c    rminor    :   Minor radius  [m]
c    rmajor    :   Major radius  [m]
c    kappa95   :   Elongation at 95\ 0.000000lux surface
c    delta95   :   Triangularity at 95\ 0.000000lux surface
c    current   :   Plasma current [MA]
c    btor      :   Magnetic field [T]
c    dense     :   Pedestal density [particle/m^3]
c    hydmass   :   Average hydrogenic mass  [amu]
c    zeff      :   Effective charge
c    q         :   Safety factor
c    t_guess   :   Guessed pedestal temperature [keV]
c
c..Outputs:
c
c    t_ped     :   Predicted pedestal temperature [keV]
c
      IMPLICIT NONE
c
      real
     &   rminor,   rmajor,   kappa95,  delta95,
     &   current,  btor,     dense,    hydmass,
     &   zeff,     q,        t_guess
c
      real
     &   zpi,      znuhat_e, znuhat_i, zep,
     &   x,        c1,       c2,       c3,
     &   c4,       D,        c_b,      alpha_c,
     &   t_ped
c
c..physical constants
c
      zpi     = atan2 ( 0.0, -1.0 )
c
c..calculate collisionality
c	
      call collision(
     &   rminor,   rmajor,   dense,    hydmass,
     &   zeff,     q,        t_guess,  znuhat_e,
     &   znuhat_i)
c
c..calculate inverse aspect ratio
c
      zep  = rminor/rmajor
c
c..calculate trapped particle fraction
c
      x  = sqrt(2*zep)
c
c..calculate the collisionality effect in bootstrap current
c  using Wesson formula
c
      c1 = (4+2.6*x)/(1+1.02*sqrt(znuhat_e)+1.07*znuhat_e)	
     &     /(1+1.07*(zep**(3./2.))*znuhat_e)
c
c..by assuming that T_i = T_e, c2 = c1
c
      c2 = c1
c
      c3 = ((7+6.5*x)/(1+0.57*sqrt(znuhat_e)+0.61*znuhat_e)
     &     /(1+0.61*(zep**(3./2.))*znuhat_e)) - 2.5*c1
c
      c4 = ((((-1.17/(1+0.46*x))+0.35*sqrt(znuhat_i))
     &     /(1+0.7*sqrt(znuhat_i)))
     &     +(2.1*(znuhat_i**2)*(zep**3)))
     &     /(1+(znuhat_e**2)*(zep**3))
     &     /(1+(znuhat_i**2)*(zep**3))
     &     *c2
c
      D  = 2.4+5.4*x+2.6*(x**2)
c
      c_b = (c1+c2+(c3+c4)/2)/(2.0*sqrt(2.0)*zpi*D)
c
c..calculate the maximum normalized pressure gradient from
c  ballooning instability
c
      alpha_c = (0.8*(1+(kappa95**2)*(1+2.0*(delta95**2))))
     &             /(1+ ((0.2*c_b*(1+(kappa95**2)*(1+2.0*(delta95**2))))
     &             /sqrt(zep)))
c
c..calcalate the temperature at the top of pedestal based on
c  ExB suppression of turbulence
c
      t_ped   = (8.689439E+25)*((btor/q)**2)
     &           *((sqrt(hydmass)/rmajor)**(2./3.))
     &           *((alpha_c/dense)**(4./3.))
c
      return
      end
c
