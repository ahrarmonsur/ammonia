      program h2odim
c modification of Toby's code to evaluate the integral over angles using direct MC integration
      implicit double precision(a-h,o-z)

      dimension com1(3),com2(3),Eulan1(3),Eulan2(3),commin(3),Eulmin(3)
      dimension com2plus(3),com2minus(3)

      parameter (nlgrid=20,ncgrid=20,one=1.d0,nrmax=10000)
      parameter (boltz=1.9872065d-3,pi=3.14159265358979323846d0)
c      parameter (massH=1.00794d0,massO=15.9994d0)
c      parameter (massH=1.008d0,massO=16.d0)
      parameter (massH=1.008d0,massO=15.9998d0,dr=0.001d0)

c      parameter (boltz=0.0019872041d0,pi=3.14159265358979323846d0)
      dimension glgrid(nlgrid),wgtgl(nlgrid),gcgrid(ncgrid),wc(ncgrid)
      dimension thetagrid(nlgrid),phigrid(ncgrid)
      dimension rthetaphi(nlgrid*ncgrid*nrmax),rotmat1(3,3)
      dimension rotmatglobal(3,3,nlgrid*ncgrid*ncgrid),rotmat2(3,3)
      dimension sum(nrmax),Emin(nrmax)
      dimension dn1(3),dn2(3),dmu(nrmax)
      dimension sum2(nrmax)
      character*30 argum
      dimension ROwf(3),
     +     R1wf(3),R2wf(3),RMwf(3), rCM(3)

      open(8,file="pmf_nojacobian_kcalpermol")
      open(9,file="pmf_nojacobian_kJpermol")
      open(10,file="pmf_withjacobian_kcalpermol")
      open(11,file="pmf_withjacobian_kJpermol")
      open(12,file="density_nojacobian")
      open(13,file="density_withjacobian")
      open(14,file="Emin")
      open(15,file="dadr")
      open(16,file="mu")


      dr2=2.d0*dr
c TIP4P geometry
      rOH=0.9572d0
      aHOH=104.52d0

c SCP
c      rOH=1.0000d0
c      aHOH=109.47
      
      rCM(1)=0.d0
      rCM(2)=0.d0
      rCM(3)=massH*rOH*cos((aHOH/2.d0)*pi/180.d0)
     +       *2.d0/(2.d0*massH+massO)

      ROwf(1)=0.d0
      ROwf(2)=0.d0
      ROwf(3)=rCM(3)
      RMwf(1)=0.d0
      RMwf(2)=0.d0
      RMwf(3)=ROwf(3)-.15d0

      R1wf(1)=rOH*cos((90.d0-aHOH/2.d0)*pi/180.d0)
      R2wf(1)=-rOH*cos((90.d0-aHOH/2.d0)*pi/180.d0)
      R1wf(2)=0.d0
      R2wf(2)=0.d0
      
      R1wf(3)=-rOH*cos((aHOH/2.d0)*pi/180.d0)+rCM(3)
      R2wf(3)=R1wf(3)

c      print *,ROwf(3),RMwf(3),R1wf(1),R1wf(3)

cE2c ... all coordinates and Euler angles are in space-fixed frame
c     write(6,*)'punch in COM coordinates of water 1 in Angs'
c     read(5,*)(com1(i),i=1,3)
c     write(6,*)'punch in COM coordinates of water 2 in Angs'
c     read(5,*)(com2(i),i=1,3)
c     write(6,*)'punch in Euler angles of water 1 in Radian'
c     read(5,*)(Eulan1(i),i=1,3)
c     write(6,*)'punch in Euler angles of water 2 in Radian'
c     read(5,*)(Eulan2(i),i=1,3)

c ... this is the test of the potential calculation
c     call caleng(com1,com2,Eulan1,Eulan2,E2h2o)

c     stop
c ... calculate the integral below

      call getarg(1,argum)
      read(argum,*)rmin
      call getarg(2,argum)
      read(argum,*)rmax
      call getarg(3,argum)
      read(argum,*)rstep
      call getarg(4,argum)
      read(argum,*) mc_max
      call getarg(5,argum)
      read(argum,*) tempt


      nstep=(rmax-rmin)/real(rstep)
  
c      print *, "rotmat built"
c     write(6,*)(acos(gcgrid(i)),i=1,ncgrid)
c     stop 'temporary stop'

c ... set r between the two H2O to be 3 angs
c     r=3.0
c ... set temperature to be 5K
      beta=1.0d0/(boltz*tempt)
c ... set the zero energy level to be -6 kcal/mol and all potentials will be subtracted by this
c ... level before being put in the exponential
c     E0=-6.0
      E0=-26.0875718462/4.184d0
      E0=0.d0
c      E0=-27.9
c ... set COM of 1st H2O to be at the origin and all three Euler angles to be 0
      do i=1,3
        com1(i)=0.0
        Eulan1(i)=0.0
      enddo
    
      do 10 ir=1,nstep
         r=rmin+(ir-1.d0)*rstep
         rplus=r+dr
         rminus=r-dr
c     ... Emin is to catch the minimum energy
         Emin(ir)=0.0d0
         gradmin=0.d0
c     ... phi and theta are the two polar angles describing the position of 2nd H2O in the first H2O frame
c     ... ph,th,and ch are the three Euler angles describing the orientation of 2nd H2O
c     ... the following 5-fold nested loop is to perform the integration over all the five angles
         dmu(ir)=0.d0
         sum(ir)=0.d0
         sum2(ir)=0.d0
         sumgrad=0.0d0
c     ... for loop over phi, C2v symmetry is used and only integrate over 0 to Pi/2.
c     ... After the loop, the integral is
         do mc_step=1,mc_max

c mc move theta,th=0,Pi phi,ph,ch=0,pi/2
c           theta=acos(2.d0*rand()-1.d0)
c           phi=2.d0*pi*rand() 

           th1=acos(2.d0*rand()-1.d0)
c           th1=pi*rand()
           ph1=2.d0*pi*rand() 
           ch1=2.d0*pi*rand() 

           th2=acos(2.d0*rand()-1.d0)
c           th2=pi*rand()
           ph2=2.d0*pi*rand() 
           ch2=2.d0*pi*rand()            

           com2(1)=0.d0
           com2(2)=0.d0
           com2(3)=r

           com2plus(1)=0.d0
           com2plus(2)=0.d0
           com2plus(3)=rplus

           com2minus(1)=0.d0
           com2minus(2)=0.d0
           com2minus(3)=rminus

           cp=cos(ph1)
           sp=sin(ph1)
           ct=cos(th1)
           st=sin(th1)
           ck=cos(ch1)
           sk=sin(ch1)
   
           dn1(1)=cp*st
           dn1(2)=sp*st
           dn1(3)=ct
               
           rotmat1(1,1)=cp*ct*ck-sp*sk
           rotmat1(1,2)=-cp*ct*sk-sp*ck
           rotmat1(1,3)=cp*st
           rotmat1(2,1)=sp*ct*ck+cp*sk
           rotmat1(2,2)=-sp*ct*sk+cp*ck
           rotmat1(2,3)=sp*st
           rotmat1(3,1)=-st*ck
           rotmat1(3,2)=st*sk
           rotmat1(3,3)=ct

           cp=cos(ph2)
           sp=sin(ph2)
           ct=cos(th2)
           st=sin(th2)
           ck=cos(ch2)
           sk=sin(ch2)
               
           dn2(1)=cp*st
           dn2(2)=sp*st
           dn2(3)=ct

c           dn1norm=0.d0
c           dn2norm=0.d0
c           do i=1,3
c               dn1norm=dn1norm+dn1(i)*dn1(i)
c               dn2norm=dn2norm+dn2(i)*dn2(i)
c           enddo
c           print *, dn1norm,dn2norm
           dn1n2=0.d0
           do i=1,3
              dn1n2=dn1n2+dn1(i)*dn2(i)
           enddo
c           print *, (2.d0*(1.d0-dn1n2))
c
c (n1+n2)^2=n1^2+n2^2+2 n1*n2=2*(1+n1*n2)
c
           rotmat2(1,1)=cp*ct*ck-sp*sk
           rotmat2(1,2)=-cp*ct*sk-sp*ck
           rotmat2(1,3)=cp*st
           rotmat2(2,1)=sp*ct*ck+cp*sk
           rotmat2(2,2)=-sp*ct*sk+cp*ck
           rotmat2(2,3)=sp*st
           rotmat2(3,1)=-st*ck
           rotmat2(3,2)=st*sk
           rotmat2(3,3)=ct
                      
           call caleng(com1,com2,E2h2o,rotmat1,rotmat2,
     +          ROwf,R1wf,R2wf,RMwf)
           call caleng(com1,com2plus,E2h2oplus,rotmat1,
     +          rotmat2,ROwf,R1wf,R2wf,RMwf)
           call caleng(com1,com2minus,E2h2ominus,rotmat1,
     +          rotmat2,ROwf,R1wf,R2wf,RMwf)
           if (beta*(E2h2o-E0) .LE. 100.) then
              sum(ir)=sum(ir)+exp(-beta*(E2h2o-E0))
              dmu(ir)=dmu(ir)+exp(-beta*(E2h2o-E0))
     +               *dsqrt(2.d0*(1.d0+dn1n2))
              sum2(ir)=sum2(ir)+exp(-2.d0*beta*(E2h2o-E0))
c              sumgrad=sumgrad+exp(-beta*(E2h2o-E0))*
c     +             (E2h2oplus-E2h2ominus)/dr2
           endif
           if(E2h2o.lt.Emin(ir)) then
              Emin(ir)=E2h2o
              gradmin=(E2h2oplus-E2h2ominus)/dr2
           endif
           
        enddo
        
c     ... scale the integral by the chebyshev's weight. since there are 3 azimuthal angles, phi, ph, and ch,
c     ... we need to raise the power of the weight to the third
c        force=sumgrad/sum(i)
        sum(ir)=sum(ir)/real(mc_max)
        dmu(ir)=dmu(ir)/real(mc_max)
        sum2(ir)=sum2(ir)/real(mc_max)
        sumgrad=32.d0*(pi**3)*sumgrad/real(mc_max)
      
c     write(6,*)'T=',tempt,'beta=',beta,'r=',r,'rho=',sumphi,'Emin=',
c     +           Emin
c  open(8,file="pmf_nojacobian_kcalpermol")
c      open(9,file="pmf_nojacobian_kJpermol")
c      open(10,file="pmf_withjacobian_kcalpermol")
c      open(11,file="pmf_withjacobian_kJpermol")
c      open(12,file="density_nojacobian")
c      open(13,file="density_withjacobian")

c         write(6,*) r,sum(ir),(-log(sum(ir))*4.184d0/beta)+E0*4.184d0,
c     +              Emin(i)*4.184d0,force*4.184d0,
c     +              1./beta*sum(ir)*r*r/dadr/rstep*4.184d0,dadr/rstep
c         flush(6)


 10   continue
      z=0.d0 
      dadr=0.d0
      dmuavg=0.d0
      do 20 ir=1,nstep
         r=rmin+(ir-1.d0)*rstep
c error on rho below (absolute)
         rhoerror=dsqrt((sum2(ir)-sum(ir)*sum(ir))/real(mc_max))
c relative error below
         rhoerrorrel=rhoerror/sum(ir)
c absolute error of pmf
         write(8,*) r,(-log(sum(ir)*r*r)/beta),rhoerrorrel/beta
         write(9,*) r,(-4.184d0*log(sum(ir)*r*r)/beta),
     *              rhoerrorrel/beta*4.184d0
         write(10,*) r,(-log(sum(ir))/beta)+E0,rhoerrorrel/beta
         write(11,*) r,(-4.184d0*log(sum(ir))/beta)+E0*4.184d0,
     *              rhoerrorrel/beta*4.184d0
         write(12,*) r,(sum(ir)*r*r)
         write(13,*) r,(sum(ir))
         write(14,*) r,(Emin(ir)),Emin(ir)*4.184

         dadr=dadr+sum(ir)*r*r

         write(15,*) r,sum(ir)*r*r/dadr/rstep

         dmuavg=dmuavg+dmu(ir)*r*r
        z=z+sum(ir)*r*r
         write(16,*) r,dmuavg/z,dmu(ir)/sum(ir),dmu(ir)

c     write(6,*)(commin(i),i=1,3),(Eulmin(i),i=1,3)
         
 20   continue

      end
      subroutine caleng(com1,com2,E2H2O,rotmat1,rotmat2,
     +                 ROwf,R1wf,R2wf,RMwf)
c ... this subroutine calculates SPC/Fw potential betwee two rigid waters given the coordinates
c ... of their centres of mass and their respective Euler angles
      implicit double precision(a-h,o-z)

      parameter(zero=0.0)
      dimension ROwf(3),R1wf(3),R2wf(3),RMwf(3),
     +          com1(3),com2(3),Eulan1(3),
     +          Eulan2(3),RO1sf(3),RO2sf(3),R11sf(3),R12sf(3),
     +          R21sf(3),R22sf(3),RM1sf(3),RM2sf(3),vec(3),
     +          rotmat2(3,3),rotmat1(3,3),crossA(3),crossB(3)

c ... SPC/WF parameters
c      parameter(epsoo=0.1554253,sigoo=3.165492,qo=-0.82,qh=0.41,
c. ..... SPC
c      parameter(epsoo=.15535372848948374760,sigoo=3.15365,
c     +                qo=-0.82,qh=0.41,
c .. TIP4P below
      parameter(epsoo=0.154875717017208413d0,sigoo=3.15365d0,
c TIP4P dlpoly
c      parameter(epsoo=0.154875717017208413d0,sigoo=3.154d0,
     +          qo=-1.040d0,qh=0.520d0,
     +          br2ang=0.52917721092d0,hr2kcl=627.509469d0)
 
c ... WFF coordinates in Angs
c      data ROwf/zero,zero,0.06562d0/,R1wf/0.7557d0,zero,-0.5223d0/,
c     +     R2wf/-0.7557d0,zero,-0.5223d0/

c PN mod below
c ... WFF coordinates in Angs
c      data ROwf/zero,zero,0.0650555d0/,
c     +     R1wf/0.7569503d0,zero,-0.5210104d0/,
c     +     R2wf/-0.7569503d0,zero,-0.5210104d0/,
c     +     RMwf/zero,zero,-0.0849445d0/


c ... prepare rotational matrix for water 1
c ... obtain the SFF coordinates for H, H, and O of water 1
c .. why?
c      call rottrn(rotmat1,ROwf,RO1sf,com1)
c      call rottrn(rotmat1,R1wf,R11sf,com1)
c      call rottrn(rotmat1,R2wf,R21sf,com1)
c      call rottrn(rotmat1,RMwf,RM1sf,com1)
c      do i=1,3
c         RO1sf(i)=ROwf(i)
c         R11sf(i)=R1wf(i)
c         R21sf(i)=R2wf(i)
c         RM1sf(i)=RMwf(i)
c      enddo

      do i=1,3
         RO1sf(i)=0.d0
      enddo
      call DGEMV ('N', 3, 3, 1.d0, rotmat1, 3, ROwf, 1, 1.d0, RO1sf, 1 )

c      call rottrn(rotmat2,R1wf,R12sf,com2)
      do i=1,3
         R11sf(i)=0.d0
      enddo
      call DGEMV ('N', 3, 3, 1.d0, rotmat1, 3, R1wf, 1, 1.d0, R11sf, 1 )

c      call rottrn(rotmat2,R2wf,R22sf,com2)
      do i=1,3
         R21sf(i)=0.d0
      enddo
      call DGEMV ('N', 3, 3, 1.d0, rotmat1, 3, R2wf, 1, 1.d0, R21sf, 1 )

c      call rottrn(rotmat2,RMwf,RM2sf,com2)
      do i=1,3
         RM1sf(i)=0.d0
      enddo
      call DGEMV ('N', 3, 3, 1.d0, rotmat1, 3, RMwf, 1, 1.d0, RM1sf, 1 )


c ... prepare rotational matrix for water 2
c ... obtain the SFF coordinates for H, H, and O of water 2
c      call rottrn(rotmat2,ROwf,RO2sf,com2)
      do i=1,3
         RO2sf(i)=com2(i)
      enddo
      call DGEMV ('N', 3, 3, 1.d0, rotmat2, 3, ROwf, 1, 1.d0, RO2sf, 1 )

c      call rottrn(rotmat2,R1wf,R12sf,com2)
      do i=1,3
         R12sf(i)=com2(i)
      enddo
      call DGEMV ('N', 3, 3, 1.d0, rotmat2, 3, R1wf, 1, 1.d0, R12sf, 1 )

c      call rottrn(rotmat2,R2wf,R22sf,com2)
      do i=1,3
         R22sf(i)=com2(i)
      enddo
      call DGEMV ('N', 3, 3, 1.d0, rotmat2, 3, R2wf, 1, 1.d0, R22sf, 1 )


c      call rottrn(rotmat2,RMwf,RM2sf,com2)
      do i=1,3
         RM2sf(i)=com2(i)
      enddo
      call DGEMV ('N', 3, 3, 1.d0, rotmat2, 3, RMwf, 1, 1.d0, RM2sf, 1 )

c$$$      dotC2=0.d0
c$$$      crossA(1)=(RO1sf(2)-R11sf(2))*(RO1sf(3)-R21sf(3))
c$$$     +         -(RO1sf(2)-R21sf(2))*(RO1sf(3)-R11sf(3))
c$$$      crossA(2)=-(RO1sf(1)-R11sf(1))*(RO1sf(3)-R21sf(3))
c$$$     +         +(RO1sf(1)-R21sf(1))*(RO1sf(3)-R11sf(3))
c$$$      crossA(3)=(RO1sf(1)-R11sf(1))*(RO1sf(2)-R21sf(2))
c$$$     +         -(RO1sf(1)-R21sf(1))*(RO1sf(2)-R11sf(2))
c$$$
c$$$      crossB(1)=(RO2sf(2)-R12sf(2))*(RO2sf(3)-R22sf(3))
c$$$     +         -(RO2sf(2)-R22sf(2))*(RO2sf(3)-R12sf(3))
c$$$      crossB(2)=-(RO2sf(1)-R12sf(1))*(RO2sf(3)-R22sf(3))
c$$$     +         +(RO2sf(1)-R22sf(1))*(RO2sf(3)-R12sf(3))
c$$$      crossB(3)=(RO2sf(1)-R12sf(1))*(RO2sf(2)-R22sf(2))
c$$$     +         -(RO2sf(1)-R22sf(1))*(RO2sf(2)-R12sf(2))

      do i=1,3
c        dotC2=dotC2+crossA(i)*crossB(i)
c        dotC2=dotC2+crossA(i)*(RO2sf(i)-R12sf(i))

      enddo
c      print *, dotC2

c ... calculate water dimer energies through SPC/WF formula
      E2H2O=0.d0
c ... O-O interaction
      roo=0.0d0
      rMM=0.0d0
      do i=1,3
        roo=roo+(RO1sf(i)-RO2sf(i))*(RO1sf(i)-RO2sf(i))
        rMM=rMM+(RM1sf(i)-RM2sf(i))*(RM1sf(i)-RM2sf(i))
      enddo
c      roo=sqrt(roo)
      rMM=sqrt(rMM)
c     write(6,*)'roo=',roo
      roo4=roo*roo
      roo6=roo4*roo
      roo12=roo6*roo6
      roo=sqrt(roo)
c      sigoo12=967747.17204467749476871048d0
c      sigoo6=983.74141523302631297392d0
c      epsoo4=.619504d0
c      A=sigoo12*epsoo4
c      C=sigoo6*epsoo4

c      A=599523.24407036588671919521d0
c AMBER value below
      A=5.99896595E+05
c SPC A
c      A=630249.62860670026275013051d0
c      A=601154.85315725817743160592d0
c      B=609.4317417025207329925d0
c amber value below
      B=6.09865468E+02
c SPC below
c      B=625.81668141130045058174d0
c      B=610.68355286022045797912d0
c ... o2lj is the L-J term between two oxygens in the unit of Kcal/Mol
c      o2lj=(sigoo12/roo12-sig2oo6/roo6)*epsoo4

      o2lj=A/roo12-B/roo6
c      print *, o2lj
c ... H-O, H-H and O-O Columbic interaction
      rho1=0.0
      rho2=0.0
      rho3=0.0
      rho4=0.0
      rhh1=0.0
      rhh2=0.0
      rhh3=0.0
      rhh4=0.0
      do i=1,3
c        rho1=rho1+(RO1sf(i)-R12sf(i))*(RO1sf(i)-R12sf(i))
c        rho2=rho2+(RO1sf(i)-R22sf(i))*(RO1sf(i)-R22sf(i))
c        rho3=rho3+(RO2sf(i)-R11sf(i))*(RO2sf(i)-R11sf(i))
c        rho4=rho4+(RO2sf(i)-R21sf(i))*(RO2sf(i)-R21sf(i))
        rho1=rho1+(RM1sf(i)-R12sf(i))*(RM1sf(i)-R12sf(i))
        rho2=rho2+(RM1sf(i)-R22sf(i))*(RM1sf(i)-R22sf(i))
        rho3=rho3+(RM2sf(i)-R11sf(i))*(RM2sf(i)-R11sf(i))
        rho4=rho4+(RM2sf(i)-R21sf(i))*(RM2sf(i)-R21sf(i))
        rhh1=rhh1+(R11sf(i)-R12sf(i))*(R11sf(i)-R12sf(i))
        rhh2=rhh2+(R11sf(i)-R22sf(i))*(R11sf(i)-R22sf(i))
        rhh3=rhh3+(R21sf(i)-R12sf(i))*(R21sf(i)-R12sf(i))
        rhh4=rhh4+(R21sf(i)-R22sf(i))*(R21sf(i)-R22sf(i))
      enddo
      rho1=sqrt(rho1)
      rho2=sqrt(rho2)
      rho3=sqrt(rho3)
      rho4=sqrt(rho4)
      rhh1=sqrt(rhh1)
      rhh2=sqrt(rhh2)
      rhh3=sqrt(rhh3)
      rhh4=sqrt(rhh4)
c      write(6,*)rho1,rho2,rho3,rho4,rhh1,rhh2,rhh3,rhh4
 
c ... ohcolm is the coulumbic term between O and H from different H2O in the unit of Hartree
      ohcolm=qo*qh*(1./rho1+1./rho2+1./rho3+1./rho4)
c ... hhcolm is ... between H and H ...
      hhcolm=qh*qh*(1./rhh1+1./rhh2+1./rhh3+1./rhh4)
c ... oocolm is ... between O and O ...
      oocolm=qo*qo*(1./rMM)
c      oocolm=qo*qo*(1./roo)


c ... scale the columbic interactions to the unit of Kcal/mol
c      ohcolm=ohcolm*hr2kcl
c      hhcolm=hhcolm*hr2kcl
c      oocolm=oocolm*hr2kcl

c      E2H2O=o2lj
      E2H2O=o2lj+(ohcolm+oocolm+hhcolm)*hr2kcl*br2ang
c      print *, E2H2O
c      E2H2O=0.d0
c     write(6,'(13(1x,f10.5))')
c    +          (com2(i),i=1,3),
c    +          (Eulan2(i),i=1,3),
c    +          ohcolm,hhcolm,oocolm,E2H2O
      
      

      end
c-----------------------------------------------------------------------
      subroutine matpre(rotmatglobal,thetagrid,phigrid,nlgrid,
     +                  ncgrid)
      implicit double precision(a-h,o-z)
      integer index
      dimension rotmatglobal(3,3,nlgrid*ncgrid*ncgrid)
      dimension phigrid(ncgrid),thetagrid(nlgrid)
      
      do iph=1,ncgrid
         phi=phigrid(iph)
         do ith=1,nlgrid
            theta=thetagrid(ith)
            do ich=1,ncgrid
               chi=phigrid(ich)
               
               cp=cos(phi)
               sp=sin(phi)
               ct=cos(theta)
               st=sin(theta)
               ck=cos(chi)
               sk=sin(chi)
               
               index=(iph*(nlgrid-1)+ith)*(ncgrid-1)+ich

               rotmatglobal(1,1,index)=cp*ct*ck-sp*sk
               rotmatglobal(1,2,index)=-cp*ct*sk-sp*ck
               rotmatglobal(1,3,index)=cp*st
               rotmatglobal(2,1,index)=sp*ct*ck+cp*sk
               rotmatglobal(2,2,index)=-sp*ct*sk+cp*ck
               rotmatglobal(2,3,index)=sp*st
               rotmatglobal(3,1,index)=-st*ck
               rotmatglobal(3,2,index)=st*sk
               rotmatglobal(3,3,index)=ct

c               rotmatglobal(1,1,index)=ck*cp-ct*sp*sk
c               rotmatglobal(1,2,index)=ck*sp+ct*cp*sk 
c               rotmatglobal(1,3,index)=sk*st 
c               rotmatglobal(2,1,index)=-sk*cp-ct*sp*ck 
c               rotmatglobal(2,2,index)=-sk*sp+ct*cp*ck 
c               rotmatglobal(2,3,index)=ck*st
c               rotmatglobal(3,1,index)=st*sp
c               rotmatglobal(3,2,index)=-st*cp
c               rotmatglobal(3,3,index)=ct

            enddo
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rottrn(rotmat,rwf,rsf,rcom)
      implicit double precision(a-h,o-z)
      dimension rotmat(3,3),rwf(3),rsf(3),rcom(3)

      do i=1,3
        rsf(i)=rcom(i)
        do j=1,3
          rsf(i)=rsf(i)+rotmat(i,j)*rwf(j)
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------
      SUBROUTINE GAULEGF(X1,X2,X,W,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X1,X2,X(N),W(N)
      PARAMETER (EPS=3.D-14)
      parameter(Pi=3.14159265358979323846d+00)
      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO 12 I=1,M
        Z=COS(pi*(I-.25D0)/(N+.5D0))
1       CONTINUE
          P1=1.D0
          P2=0.D0
          DO 11 J=1,N
            P3=P2
            P2=P1
            P1=((2.D0*J-1.D0)*Z*P2-(J-1.D0)*P3)/J
11        CONTINUE
          PP=N*(Z*P1-P2)/(Z*Z-1.D0)
          Z1=Z
          Z=Z1-P1/PP
        IF(ABS(Z-Z1).GT.EPS)GO TO 1
        X(I)=XM-XL*Z
        X(N+1-I)=XM+XL*Z
        W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
        W(N+1-I)=W(I)
12    CONTINUE
      RETURN
      END

c-----------------------------------------------------------------
      SUBROUTINE GAUCHEB(X1,X2,X,N)
      implicit double precision(a-h,o-z)
      dimension x(N)
      parameter(Pi=3.14159265358979323846d+00)

      xl=(x2-x1)/2.d0
      do i=1,n
        phi=pi*dfloat(2*i-1)/dfloat(2*n)
        x(i)=x1+xl*(cos(phi)+1.d0)
      enddo

      return
      end
