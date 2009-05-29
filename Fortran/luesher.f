C Hi Al,

C If the variable XX is supposed
C to get a random value, just say 

C XX=rndnum() .

C For the initialisation, in the main program just after reading in the
C old configs ( the REWIND and READ for units 40,50) a statement
C  CALL rini(iseed)
C is needed.

C I believe that is all. In the routine random.f, just before the end
C statement, I have added a small positive number to the random number.
C This is because the radnoms appear as arguments of logs and
C my machine does not accept that. If yours takes proper action, you can
C remove that. Or else, you may have to change its size according to machine
C specification. As it is, it is something like the smalles 64 bit number
C on some Sun.

C Good luck with the implementation. If it does not work let me know and I
C can try.

C Cheers, Owe.


*  file random.f, Ulli Wolff, 15.9.95
*
*  random generator modified RCARRY, 
*  see M. Luescher, Comp. Phys. Comm. 79(1994)100
*  and literature quoted there
*
*  This generator is relatively slow, but of KNOWN good quality
*  The speed strongly depends on the value ithrow
*  One way to speed up would be vectoriztion over parallel copies
*
*  ALWAYS CHANGE PARAMETERS IN BOTH ROUTINES OF THIS FILE
*
*  values for ithrow: 0 for simple uses, 24 for normal, 199 for
*  practically no detectable correlations, 365 for tests in delicated cases;
*  the bigger ithrow, the slower the generator!
*
*  0 <= rndnum < 1  (rndnum=0 DOES occur occasionally)
*  
**************************************************************************

      real*8 function rndnum()
*
      integer ir,is,irs,ir1,ikeep,ithrow,ivl1,ivl,ibase,iv,icall
      real basein
      parameter (ir=24,is=10,irs=ir-is,ir1=ir-1)
      parameter (ikeep=ir,ithrow=24,ivl1=ikeep+ithrow,ivl=ivl1+ir1)
      parameter (ibase=2**24,basein=1.0/ibase)
*
*  state of the generator (to be saved and reloaded for continuation):
      common/rrand/ iv(0:ivl),icall
      save /rrand/
*
      integer i,j,ivn
*
      if (icall.eq.ikeep) then
*  disregard ithrow numbers:
        do j=1,ithrow
          ivn=iv(irs+icall)-iv(icall)
          icall=icall+1
*  carry bit:
          if (ivn.lt.0) then
            iv(icall)=iv(icall)+1
            ivn=ivn+ibase
          endif
*
          iv(icall+ir1)=ivn
        enddo
*  copy last ir numbers to the beginning:
        do i=0,ir1
          iv(i)=iv(ivl1+i)
        enddo
      icall=0
      endif
*
*  basic cycle:
*
      ivn=iv(irs+icall)-iv(icall)
      icall=icall+1
*  carry bit:
      if (ivn.lt.0) then
        iv(icall)=iv(icall)+1
        ivn=ivn+ibase
      endif
*
      iv(icall+ir1)=ivn
*  convert to floating point:
      rndnum=float(ivn)*basein
      if (rndnum.eq.0.0) rndnum = 4.5d-9
      end


C**************************************************************************

c     This gets num random numbers at one time, to save on overheads.
c     There is the pesky IF (icall .EQ. ikeep) line here.

c     SInce ikeep is only 24, this statement is unavoidable.

C**************************************************************************

      SUBROUTINE MULTI_RNDNUM(num,rndnum)
*
      INTEGER num
      REAL*8 rndnum(num)

      integer ir,is,irs,ir1,ikeep,ithrow,ivl1,ivl,ibase,iv,icall
      real basein
      parameter (ir=24,is=10,irs=ir-is,ir1=ir-1)
      parameter (ikeep=ir,ithrow=24,ivl1=ikeep+ithrow,ivl=ivl1+ir1)
      parameter (ibase=2**24,basein=1.0/ibase)
*
*  state of the generator (to be saved and reloaded for continuation):
      common/rrand/ iv(0:ivl),icall
      save /rrand/
*
      integer i,j,ivn

      INTEGER nn
*
      DO nn = 1,num

      if (icall.eq.ikeep) then
*  disregard ithrow numbers:
        do j=1,ithrow
          ivn=iv(irs+icall)-iv(icall)
          icall=icall+1
*  carry bit:
          if (ivn.lt.0) then
            iv(icall)=iv(icall)+1
            ivn=ivn+ibase
          endif
*
          iv(icall+ir1)=ivn
        enddo
*  copy last ir numbers to the beginning:
        do i=0,ir1
          iv(i)=iv(ivl1+i)
        enddo
      icall=0
      endif
*
*  basic cycle:
*
      ivn=iv(irs+icall)-iv(icall)
      icall=icall+1
*  carry bit:
      if (ivn.lt.0) then
        iv(icall)=iv(icall)+1
        ivn=ivn+ibase
      endif
*
      iv(icall+ir1)=ivn
*  convert to floating point:
      rndnum(nn)=float(ivn)*basein
      if (rndnum(nn).eq.0.0d0) rndnum(nn)=4.5d-9

       ENDDO

      end


**************************************************************************
*==========================================================
**************************************************************************

*  initialize from a random seed iran, 0 <= iran < 259200, 
*  with the help of a "quick and dirty generator
*
      subroutine rini(iran)
*
      integer ir,is,irs,ir1,ikeep,ithrow,ivl1,ivl,ibase,iv,icall
      real basein
      parameter (ir=24,is=10,irs=ir-is,ir1=ir-1)
      parameter (ikeep=ir,ithrow=24,ivl1=ikeep+ithrow,ivl=ivl1+ir1)
      parameter (ibase=2**24,basein=1.0/ibase)
*
*  state of the generator (to be saved and reloaded for continuation):
      common/rrand/ iv(0:ivl),icall
      save /rrand/
*
      integer im,ia,ic,iran,jran,i,ifac
*
*  parameters for auxiliary generator used for initialization:
      parameter(im=259200,ia=7141,ic=54773)
*
      if(iran.lt.0.or.iran.ge.im)stop
     X 'rini: iran out of range'
      jran=iran
*  warm up the auxiliary generator a little
      do i=1,10
*  cycle of auxiliary generator:
        jran=mod(jran*ia+ic,im)
      enddo
*
      ifac=(ibase-1)/im
*  initialize 0 <= iv < ibase:
      do i=0,ir1
        jran=mod(jran*ia+ic,im)
        iv(i)=ifac*jran
      enddo
      icall=0
      write(*,'(/)')
      write(*,*)' rini: random #s newly initialized, iran = ',iran
      end
*==========================================================
