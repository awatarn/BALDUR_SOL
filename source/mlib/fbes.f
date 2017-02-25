        function fbes(kbes,pzx)
        dimension bes(10)
        data nmax/10/
        bes(nmax)=0.
        if((pzx.eq.0.).and.(kbes.eq.1)) go to 50
        if((pzx.eq.0.).and.(kbes.eq.0)) go to 40
        bes(nmax-1)=1.
c       compute trial values
        do 20 n=nmax-2,1,-1
 20     bes(n)=(2*n/pzx)*bes(n+1)-bes(n+2)
c       normalize
        sum=bes(1)
        do 300 ni=3,(nmax-1),2
 300    sum=sum+2*bes(ni)
        fac=1./sum
        if(kbes.eq.0) fbes=bes(1)*fac
        if(kbes.eq.1) fbes=bes(2)*fac
        return
 40     fbes=1.
        return
 50     fbes=0.
        return
        end
