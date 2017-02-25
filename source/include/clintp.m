c@clintp for the island package DISLAND2, Bateman, PPPL, 7-nov-91
c
      common /cintrp/
     &   xintrp, hamji, xmai, xmbi, xmci, xmdi, xmei, xmfi
     & , hixixi(kpthta), hithth(kpthta), hizeze(kpthta)
     & , hixith(kpthta), hixize(kpthta), hithze(kpthta)
     & , xinprf, uj0zei, duj0zi, dpnrmi, b0uzei, db0zei
     & , rinv2i, drnv2i, qprofi, dqprfi
c
c  Values interpolated in sbrtn intham and intprf.
c
c       The following metric elements in Hamada coordinates are 
c  interpolated to xi = xintrp in sbrtn intham and
c  stored as a function of theta(jt):
c
c  hixixi(jt) = ( d R / d xi )^2 + ( d Y / d xi )^2 
c
c  hithth(jt) = ( d R / d theta )^2 + ( d Y / d theta )^2
c
c  hizeze(jt) = 1 / R^2
c
c  hixith(jt) = (d R / d xi)(d R / d theta)
c                   + (d Y / d xi)(d Y / d theta)
c
c  hixize(jt) = d Z / d xi
c
c  hithze(jt) = d Z / d theta
