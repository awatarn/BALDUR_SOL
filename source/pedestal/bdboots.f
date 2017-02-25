c
c This subroutine is used to calculate the collisionality effect in 
c bootstrap current.
c
c@bootstrap.f
c--------1---------2---------3---------4---------5---------6---------7-c
c
        subroutine bdboots (
     &   zrminor,   zrmajor,   nuhat_i, nuhat_e,
     &   c_b)
c
c  Input:
c  ------
c zrminor     minor radius [m]
c zrmajor     major radius [m]
c nuhat_e     electron normalized collisionality 
c nuhat_i     ion normalized collisionality
c     
c  Output:
c  -------
c c_b         factor of colliosnality that modifies the bootstrap current  
c
      IMPLICIT NONE
c
      real
     &   zrminor,   zrmajor,   nuhat_e, nuhat_i,
     &   pi,        esp,       x,       c1,
     &   c2,        c3,        c4,      D,
     &   c_b
c
c..physical constants
c
      pi = atan2 ( 0.0, -1.0 )
c
c..calculate inverse aspect ratio
c
      esp = zrminor/zrmajor
c
c..calculate trapped particle fraction
c
      x  = sqrt(2*esp)
c
c..calculate the collisionality effect in bootstrap current
c  using Wesson formula 
c  (J. Wesson, Tokamaks (Clarendon, Oxford, England, 1997))
c
      c1 = (4+2.6*x)/(1+1.02*sqrt(nuhat_e)+1.07*nuhat_e)	
     &     /(1+1.07*(esp**(3./2.))*nuhat_e)
c
c..by assuming that T_i = T_e, c2 = c1
c
      c2 = c1
c
      c3 = ((7+6.5*x)/(1+0.57*sqrt(nuhat_e)+0.61*nuhat_e)
     &     /(1+0.61*(esp**(3./2.))*nuhat_e)) - 2.5*c1
c
      c4 = ((((-1.17/(1+0.46*x))+0.35*sqrt(nuhat_i))
     &     /(1+0.7*sqrt(nuhat_i)))
     &     +(2.1*(nuhat_i**2)*(esp**3)))
     &     /(1+(nuhat_e**2)*(esp**3))
     &     /(1+(nuhat_i**2)*(esp**3))
     &     *c2
c
      D  = 2.4+5.4*x+2.6*(x**2)
c
      c_b = (c1+c2+(c3+c4)/2)/(2.0*sqrt(2.0)*pi*D)
c
      return
      end
c
