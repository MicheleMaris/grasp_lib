__DESCRIPTION__="""
 grasp_lib.py V 0.6 - 3 Feb 2012 - 23 Mar 2012 (0.4) - 2012 Nov 27 (0.5) - 2013 Dec 12 -
   M.Maris, M.Sandri, F.Villa
   From a set  of routines created by M.Sandri e F.Villa

 This library allows to import a grasp file in GRASP format and to convert it in an healpix map, it also
 performs interpolation of the grasp map

 In GRASP convention theta and phi are colatitude and longitude, with 
    phi    meridian circle longitude, [0,180]deg
    theta  on meridian circle polar distance, [-180,180] deg
 a point with theta < 0 denotes a point located at the same polar distance of abs(theta) but with longitude pi+180deg

 The usual polar convention is 
    phi    meridian halfcircle longitude, [-180,180]deg
    theta  on meridian circle polar distance, [0,180] deg
"""

def thetaUVphiUV2UV(thetaUV,phiUV,deg=True) :
   """converts thetaUV and phiUV into U=x0=sin(thetaUV)*cos(phiUV), V=y0=sin(thetaUV)*sin(phiUV)"""
   from numpy import pi, sin, cos
   if deg : return sin(180./pi*thetaUV)*cos(180./pi*phiUV),sin(180./pi*thetaUV)*sin(180./pi*phiUV)
   return sin(thetaUV)*cos(phiUV),sin(thetaUV)*sin(phiUV)

def UV2thetaUVphiUV(U,V,deg=True) :
   from numpy import pi, arctan2, arccos,arcsin,sin,cos,array,mod
   f=(180./pi) if deg else 1.
   phiUV=arctan2(V,U)
   A=cos(phiUV)
   B=sin(phiUV)
   thetaUV=arcsin((A*U+B*V)/(A*A+B*B))
   if deg : return f*thetaUV,mod(f*phiUV,360.)
   return f*thetaUV,mod(phiUV,2.*pi)

def phitheta2longcolat(phi_grasp,theta_grasp) :
   """ converts GRASP (phi,theta) coordinates into standard (long,colat) 
   upon request returns a structure """
   import numpy as np
   _long = phi_grasp*1.
   colat=theta_grasp*1.
   idx = np.where(colat < 0)[0]
   if len(idx) > 0 :
      _long[idx]=_long[idx]+180.
      colat[idx]=np.abs(colat[idx])
   return _long,colat

def longcolat2phitheta(_long,colat) :
   """converts ususal polar (long,colat) coordinates into GRASP (phi,theta) returns a structure"""
   import numpy as np
   phi=_long*1.
   theta=colat*1.
   idx = np.where(phi >= 180.)[0]
   if len(idx) > 0 :
      phi[idx]=phi[idx]-180.
      theta[idx]=-theta[idx]
   return phi,theta

def longcolat2rowcol(_long,colat,phi0,dphi,theta0,dtheta) :
   """ converts (long,colat) into index of phi and of theta in the matrix """
   import numpy as np
   phi,theta=longcolat2phitheta(_long,colat)
   return (pt.phi-phi0)/dphi,(pt.theta-theta0)/dtheta

def ipix2rowcol(nside,ipix,phi0,dphi,theta0,dtheta,nest=False) :
   """ converts an healpix ipix (ring) into index of phi and of theta in the matrix"""
   from healpy import pix2ang_ring
   colat,_long=pix2ang(nside,ipix,nest=nest)
   return longcolat2rowcol(_long/np.pi*180.,colat/np.pi*180.,phi0,dphi,theta0,dtheta)

def nside2ipix(nside,Reversed=False) :
   """ converts nside into a list of pixels (ring) 
       reversed = True means the orderring is reversed
   """
   import numpy as np
   if not Reversed : return np.arange(12*int(nside)*int(nside))
   return np.arange(12*int(nside)*int(nside)-1,-1,-1)

from grid2d import *

#class GraspMapsCube :
   #"a list of grasp maps is used to integrate over a grasp map in band"
   #def _class_line_cube :
      #def __init__(self,listNames) :
         #self.Name=[]
         #for k in range(len(listNames)) :
            #try :
               #self.plane(open(listNames[k],'r').readlines())
               #self.Name.append(listNames)
            #except :
               #print k," impossible to read, skipped"
         #self.N=len(self.Name)
      #def __len__(self) :
         #return self.N
      #def __getitem__(self,i) :
         #import numpy as np
         #l=[]
         #for k in range(len(self)) :
            #l.append(self.plane[k][i]))
         #return l
   #def __init__(self,listNames) :
      #import numpy as np
      #self._line=-1
      #self._cube=_class_line_cube(listNames)
   #def __len__(self) :
      #return len(self._cube)
   #def _fetch_line(self) :
      #self._line=+1
      #return self._cube[self._line]

def components2cocross(r1,i1,r2,i2) :
   "given r1,i1,r2,i2 return Eco,Ecross"
   re=[r1,r2]
   im=[i1,i2]
   p=[]
   p.append(r1**2+i1**2)
   p.append(r2**2+i2**2)
   if p[0].max() > p[1].max() :
      ico=0
      icross=1
   else :
      ico=1
      icross=0
   Eco=re[ico]*complex(0,1.)*re[ico]
   Ecross=re[icross]*complex(0,1.)*re[icross]
   return Eco,Ecross
   
def cocross2rhclhc(Eco,Ecross) :
   "given Eco,Ecross return Erhc,Elhc"
   isqrt2=2.**(-0.5)
   Erhc=(Eco-complex(0,1.)*Ecross)*isqrt2
   Elhc=(Eco+complex(0,1.)*Ecross)*isqrt2
   return Erhc,Elhc

def components2rhclhc(r1,i1,r2,i2) :
   "given r1,i1,r2,i2 return Erhc,Elhc"
   Eco,Ecross=components2cocross(r1,i1,r2,i2)
   return cocross2rhclhc(Eco,Ecross)
   
def polarization_ellipse_from_fields(Erhc,Elhc) :
   "given Erhc,Elhc returns rmajor, rminor, directivity, psi_pol_rad"
   import numpy as np
   isqrt2=2.**(-0.5)
   rmajor=abs(abs(Erhc)+abs(Elhc))*isqrt2
   rminor=abs(abs(Erhc)-abs(Elhc))*isqrt2
   directivity=abs(Erhc)**2+abs(Elhc)**2
   aa=(Erhc/Elhc)**0.5
   psi_pol=np.arctan2(aa.imag,aa.real)
   return rmajor,rminor,directivity,psi_pol
   
def components2polarization_ellipse(r1,i1,r2,i2) :
   "given r1,i1,r2,i2 return rmajor, rminor, directivity, psi_pol_rad"
   Erhc,Elhc=components2rhclhc(r1,i1,r2,i2)
   return polarization_ellipse_from_fields(Erhc,Elhc)

class GraspMap(MapGrid) :
   def __init__(self,inputfile,skiplines,CounterPhi=[1e6,1.],silent=False,useCounterPhi=False,closeColumn=False,Pickle=False,periodicColumn=False,badConversionValue=0.) :
      """badConversionValue = Value to replace samples with problems in converting strings to numbers"""
      MapGrid.__init__(self)
      self._init_failed = True
      if Pickle :
         self.load(inputfile)
         self._init_failed = False
         return
      self.info['graps_file']=inputfile.strip()
      self.info['projection']='GRASP-CUT'
      self.info['ReferenceDocument']="LFI BEAMS DELIVERY: FORMAT SPECIFICATIONS\nM. Sandri\nPL-LFI-PST-TN-044, 1.0,July 2003"
      self.info['numbad']=-1
      if inputfile.strip() == '' :
         return
      self.get_cuts(inputfile,skiplines,badConversionValue=badConversionValue)#,CounterPhi=CounterPhi,silent=silent,useCounterPhi=useCounterPhi) 
      if periodicColumn : 
         self.set_col_periodic()
      if closeColumn : 
         print 'Closing rows at right'
         for k in self.M.keys() :
            if k!='_row_values' and k!='_col_values' and k!='_row_index' and k!='_col_index' :
               for r in range(self.R['n']) :
                  self.M[k][r][-1]=self.M[k][self.R['n']-1-r][0]
      #if closeColumn : self.right_close_col()
   def get_cuts(self,inputfile,skiplines,CounterPhi=[1e6,1.],silent=False,useCounterPhi=False,badConversionValue=0.) :
      """ get_cuts
      
 this program get the cuts file
 inputfile = name of the input file

 the output is a structure with an entry for each cut in phi


 COUNTER_PHI v.z. HEADER_PHI
   By default the field PHI in each cut is the one declared in the
   header of the block of the input grasp file.
   The value is also returned as HEADER_PHI in the structure
   There are some cases in which the header is bad formatted and 
   the PHI is not reliable.
   To solve this problem GET_CUTS provides an internal PHI calculator
   assuming PHI increments on constant steps, the value from
   the calculator is in COUNTER_PHI
   The counter is tuned by using the keyword
      CounterPhi = [phi0,step] (default [1d6,1d0])
   so the default is simply a counter of the number of cuts. 
   The value 1e6. as first value is to assure the COUNTER_PHI is not 
   confused with an angle.
   To make COUNTER_PHI as PHI it is sufficient to add the keyword
      /useCounterPhi
   so that the PHI will be forced to be COUNTER_PHI instead of
   HEADER_PHI.
   At last the tag PHI_TYPE in the structure specifies wether PHI 
   comes from the HEADER or from the COUNTER
   
   In self.info['numbad'] is the number of lines which can not be properly decoded into float,
   they are marked bad (1) in the flag_bad map.
   badConversionValue = Value to replace samples with problems in converting strings to numbers
   
      """
      import sys
      import numpy as np
      import copy
      CounterPhi0=CounterPhi[0]*1.
      CounterDeltaPhi=CounterPhi[1]*1.
      deltaphi=1.
      phi0=0.
      header='header'
      fileinput='fileinput'
      thetai=0.
      dtheta=0.
      ntheta=1
      phi=0.
      k1=1
      k2=1
      k3=1
      comp1r=0.
      comp1i=0.
      comp2r=0.
      comp2i=0.
   #
   #********************************
   #
      self.clean()
      self.info['hdr']={'file':inputfile}
      self.info['wrong_lines']=[]
      print "Reading ",inputfile
      try :
         h=open(inputfile,'r').readlines()
         self.mapname=inputfile
         self.info['inputfile']=inputfile
      except :
         print "File %s not found"%inputfile
         return
      # removes the new line
      for i in range(len(h)) :
         h[i] = h[i].split('\n')[0].split('\r')[0]
      # skips a given number of lines 
      self.info['hdr']['skipped']=[]
      self.info['hdr']['skiplines']=skiplines*1
      if skiplines > 0 :
         for line in h[0:(skiplines-1)] :
            self.info['hdr']['skipped'].append(line.split('\n')[0])
         h=h[skiplines:]
      # skips all the lines until it reaches an header
      notHeader = True
      icount=-1
      while notHeader :
         icount+=1
         ll=h[icount].split('\n')[0].strip().split()
         try : 
            lla = np.array(ll)
            notHeader = False
         except :
            notHeader = True
         if not notHeader :
            if len(ll) != 7 :
               notHeader = True
      icount-=1
      if icount > 0 :
         for k in h[0:icount] :
            self.info['hdr']['skipped'].append(k)
            self.info['hdr']['skiplines']+=1
         h=h[icount:]
      # the second line of the first block gives the number of lines per block
      currentline=1
      ll=h[currentline].split('\n')[0].strip().split()
      try :
         thetai = float(ll[0])
         dtheta = float(ll[1])
         ntheta = int(ll[2])
         header_phi = float(ll[3])
         k1 = int(ll[4])
         k2 = int(ll[5])
         k3 = int(ll[6])
      except :
         return h[currentline-1],'level 1',currentline,h
      self.info['nlines']=len(h)
      self.info['blocksize']=ntheta+2
      self.info['nblocks'] = len(h)/(ntheta+2)
      nblocks = self.info['nblocks']
      self.info['thetai']=np.zeros(nblocks)
      self.info['dtheta']=np.zeros(nblocks)
      self.info['ntheta']=np.zeros(nblocks,dtype='int')
      self.info['phi']=np.zeros(nblocks)
      self.info['k1']=np.zeros(nblocks,dtype='int')
      self.info['k2']=np.zeros(nblocks,dtype='int')
      self.info['k3']=np.zeros(nblocks,dtype='int')
      self.info['line']=np.zeros(nblocks,dtype='string')
      self.info['fail']=np.zeros(nblocks,dtype='int')+1
      self.info['iline']=np.zeros(nblocks,dtype='int')
      if (ntheta+2)*self.info['nblocks']-len(h) != 0 :
         print "Error: too much or too few lines to form the required number of blocks"
         print "Nblocks : ",self.info['nblocks']
         print "lines   : ",len(h)
         print "lines per block : ",ntheta+2
         print "lines in blocks : ",(ntheta+2)*nblocks
         print "residual lines  : ",(ntheta+2)*nblocks-len(h)
         return None
      print self.info['nblocks']," x ",ntheta," elements"
      # decomposes all the headers
      for i_block in range(nblocks) :
         ii = i_block*self.info['blocksize']
         self.info['iline'][i_block]=ii*1
         self.info['line'][i_block]=h[ii]+''
         ll=h[ii+1].split('\n')[0].strip().split()
         try :
            self.info['thetai'][i_block] = float(ll[0])
            self.info['dtheta'][i_block] = float(ll[1])
            self.info['ntheta'][i_block] = int(ll[2])
            self.info['phi'][i_block] = float(ll[3])
            self.info['k1'][i_block] = int(ll[4])
            self.info['k2'][i_block] = int(ll[5])
            self.info['k3'][i_block] = int(ll[6])
            self.info['fail'][i_block] = 0
         except :
            print "Fail to decode block %d, line %d\n'%s'\n"%(i_block,ii,h[ii+1])
      if self.info['fail'].sum() > 0 :
         print "fail to decode blocks"
         return
      # sets the phi along the x axis of the grid i.e. the columns
      self.set_col_scale('phi','deg',self.info['phi'])
      self.dphi=self.C['delta']
      self.phi0=self.C['min']
      # sets the theta along the y axis of the grid  i.e. the rows
      self.set_row_scale('theta','deg',np.arange(self.info['ntheta'][0])*self.info['dtheta'][0]+self.info['thetai'][0])
      self.dtheta=self.R['delta']
      self.theta0=self.R['min']
      #initialize private compoenents used for debug
      self.newmap('_line_index',dtype='int')
      # initializes the five component matrices
      self.newmap('r1')
      self.newmap('i1')
      self.newmap('r2')
      self.newmap('i2')
      self.newmap('power')
      self.newmap('flag_bad',dtype='int')
      #
      # fill the component matrices
      # each block is for a given phi, i.e. a given Column 
      # each value is for a given theta, i,e, a given Row
      self.newmap('phi')
      self.newmap('theta')
      self.info['numbad']=0
      for i_block in range(nblocks) :
         ii = i_block*self.info['blocksize']+2
         for i_raw in range(self.R['n']) :
            iline = ii+i_raw
            self.M['phi'][i_raw,i_block]=self.C['v'][i_block]*1
            self.M['theta'][i_raw,i_block]=self.R['v'][i_raw]*1
            self.M['_row_values'][i_raw,i_block]=self.R['v'][i_raw]*1
            self.M['_col_values'][i_raw,i_block]=self.C['v'][i_block]*1
            self.M['_row_index'][i_raw,i_block]=i_raw*1
            self.M['_col_index'][i_raw,i_block]=i_block*1
            self.M['_line_index'][i_raw,i_block]=iline*1
            ll=h[iline].split('\n')[0].strip().split()
            lla=np.zeros(4)
            nbad=0
            for icol in range(4) :
               try :
                  lla[icol]=np.array(ll[icol],dtype=float)
               except :
                  lla[icol]=badConversionValue
                  nbad+=1
            if nbad == 0 : 
               self.M['flag_bad'][i_raw,i_block]=0
            else :
               self.info['wrong_lines'].append("skiplines %d, line %d,  row %d, block %d\n%s"%(self.info['hdr']['skiplines'],iline,i_raw,i_block,h[iline]))
               print "Invalid litteral in file, skiplines %d, line %d,  row %d, block %d\n"%(self.info['hdr']['skiplines'],iline,i_raw,i_block)
               print ">"+" ".join(ll)+"< out >",lla,'<'
               self.M['flag_bad'][i_raw,i_block]=1
               self.info['numbad']+=1
            self.M['r1'][i_raw,i_block]=lla[0]*1
            self.M['i1'][i_raw,i_block]=lla[1]*1
            self.M['r2'][i_raw,i_block]=lla[2]*1
            self.M['i2'][i_raw,i_block]=lla[3]*1
      self.M['power']=self.M['r1']**2+self.M['i1']**2+self.M['r2']**2+self.M['i2']**2
      # stores longitudes and latitudes
      self.newmap('long')
      self.newmap('colat')
      self.M['long']=self.M['_row_values']*0
      self.M['colat']=self.M['_row_values']*0
      self.M['long'],self.M['colat']=phitheta2longcolat(self.M['phi'],self.M['theta'])
      self._init_failed = False
   def initFailed(self) :
      return self._init_failed
   def haveBadSamples(self) :
      if self.initFailed(): return False
      return self.info['numbad']>0
   def formatGrasp(self) : return {'float':' %17.10e','int':' %4d'}
   def recompose_header(self) :
      import copy
      hdr=[]
      if self.info['hdr']['skiplines']>0 :
         hdr=copy.deepcopy(self.info['hdr']['skipped'])
      return hdr
   def recompose_block_header(self,i_block) :
      import numpy as np
      fmtF=self.formatGrasp()['float']
      fmtI=self.formatGrasp()['int']
      if i_block < 0 : return
      if i_block >= self.info['nblocks'] : return
      ll=''
      ll+=fmtF%self.info['thetai'][i_block]
      ll+=fmtF%self.info['dtheta'][i_block]
      ll+=fmtI%self.info['ntheta'][i_block]
      ll+=fmtF%self.info['phi'][i_block]
      ll+=fmtI%self.info['k1'][i_block]
      ll+=fmtI%self.info['k2'][i_block]
      ll+=fmtI%self.info['k3'][i_block]
      return ['Planck,',ll.upper()]
   def recompose_block_data(self,i_block,tab0=None) :
      import numpy as np
      import copy
      fmtF=self.formatGrasp()['float']
      fmtI=self.formatGrasp()['int']
      if i_block < 0 : return
      if i_block >= self.info['nblocks'] : return
      tab=[]
      if type(tab0)==type([]) : tab=copy.deepcopy(tab0) 
      if type(tab0)==type('') : tab=[tab0]
      for i_raw in range(self.R['n']) :
         ll=''
         ll+=fmtF%self.M['r1'][i_raw,i_block]
         ll+=fmtF%self.M['i1'][i_raw,i_block]
         ll+=fmtF%self.M['r2'][i_raw,i_block]
         ll+=fmtF%self.M['i2'][i_raw,i_block]
         tab.append(ll.upper())
      return tab
   def recompose_block(self,i_block,tab0=None,fmt='%18.10e') :
      if i_block < 0 : return
      if i_block >= self.info['nblocks'] : return
      tab=self.recompose_block_header(i_block) 
      return self.recompose_block_data(i_block,tab0=tab)
   def recompose_map(self,tab0=None) :
      tab=[]
      for i_block in range(self.info['nblocks']) :
         for l in self.recompose_block(i_block) :
            tab.append(l)
      return tab
   def FourColumnsPower(self,power1Name='p1',power2Name='p2',powerName='power') :
      "a FourColumns map has r1=sqrt(p1), i1=0, r2=sqrt(p2), i2=0"
      new=self.copy()
      new.info['ktype']=1
      if self.M.has_key(power1Name) and self.M.has_key(power1Name) :
         new.info['ncomp']=2
         new.M['r1']=self[power1Name]**0.5
         new.M['r2']=self[power2Name]**0.5
         new.M['i1']=self[power1Name]*0
         new.M['i2']=self[power1Name]*0
      elif self.M.has_key(power) :
         new.info['ncomp']=1
         new.M['r1']=self[power]**0.5
         new.M['r2']=self[power]*0
         new.M['i1']=self[power]*0
         new.M['i2']=self[power]*0
      else :
         print "the map shall contain ",power1Name,power2Name," or ",powerName
         return
      return new
   def grasp2longcolat(self,phi_grasp,theta_grasp) :
      """ converts GRASP (phi,theta) coordinates into standard (long,colat) 
      upon request returns a structure """
      import numpy as np
      _long = phi_grasp*1.
      colat=theta_grasp*1.
      idx = np.where(colat < 0)[0]
      if len(idx) > 0 :
         _long[idx]=_long[idx]+180.
         colat[idx]=np.abs(colat[idx])
      return _long,colat
   def longcolat2grasp(self,_long,colat) :
      """converts ususal polar (long,colat) coordinates into GRASP (phi,theta) returns a structure"""
      import numpy as np
      phi=_long*1.
      theta=colat*1.
      idx = np.where(phi >= 180.)[0]
      if len(idx) > 0 :
         phi[idx]=phi[idx]-180.
         theta[idx]=-theta[idx]
      return phi,theta
   def longcolat2rowcol(self,_long,colat) :
      """ converts (long,colat) into index of phi and of theta in the matrix """
      import numpy as np
      phi,theta=self.longcolat2grasp(_long,colat)
      return (theta-self.theta0)/self.dtheta,(phi-self.phi0)/self.dphi
   def ipix2longcolat(self,nside,ipix,nest=False,deg=True) :
      """ converts an healpix ipix (ring) into index of phi and of theta in the matrix"""
      from healpy import pix2ang
      import numpy as np
      colat,_long=pix2ang(nside,ipix,nest=nest)
      if deg : return _long*180./np.pi,colat*180./np.pi
      return _long,colat
   def ipix2rowcol(self,nside,ipix,nest=False,deg=False) :
      """ converts an healpix ipix (ring) into index of phi and of theta in the matrix"""
      from healpy import pix2ang
      import numpy as np
      colat,_long=pix2ang(nside,ipix,nest=nest)
      if deg : return self.longcolat2rowcol(_long,colat)
      return self.longcolat2rowcol(_long/np.pi*180.,colat/np.pi*180.)
   def nside2ipix(self,nside,Reversed=False) :
      """ converts nside into a list of pixels (ring) 
         reversed = True means the orderring is reversed
      """
      return nside2ipix(nside,Reversed=Reversed)
   def parseColatRange(self,colatrange) :
      prs=(colatrange.strip()).split(',')
      left = [prs[0][0],float(prs[0][1:])]
      right = [prs[1][-1],float(prs[1][0:-1])]
      return left,right
   def healpix(self,nside,mapname='power',nest=False,Reversed=False,rot=[0.,0.],colatrange=None) :
      """converts to healpix or a stack of healpix maps of given nside
         colatrange=None , takes all the map
         colatrange=']a,b['
         colatrange='[a,b['
         colatrange=']a,b]'
         for gridmap
      """
      import numpy as np
      import healpy as H
      if colatrange==None :
         ipix=self.nside2ipix(nside,Reversed=Reversed)
         if rot[0]==0. and rot[1]==0 :
            _long,colat = self.ipix2longcolat(nside,ipix)
         
         phi,theta=self.longcolat2grasp(_long,colat)
         r1=self.bilinearXY(mapname,phi,theta)
         return r1 
      else :
         fact=180./np.pi
         prs=(colatrange.strip()).split(',')
         left = [prs[0][0],float(prs[0][1:])]
         right = [prs[1][-1],float(prs[1][0:-1])]
         NPEQ=12*nside/2
         print left,right
         ipixmin=H.ang2pix(nside,left[1]/fact,0)-NPEQ
         ipixmax=H.ang2pix(nside,right[1]/fact,0)+NPEQ
         if ipixmin < 0 : ipixmin=0
         if ipixmax > 12*nside**2-1 : ipixmax=12*nside**2
         ipix = np.arange(ipixmin,ipixmax)
         colat,Long = H.pix2ang(nside,ipix)
         fl=np.ones(len(colat))
         if left[1] == ']' :
            fl*=(left[1]/fact)<colat
         else :
            fl*=(left[1]/fact)<=colat
         if right[1] == '[' :
            fl*=colat<(right[1]/fact)
         else :
            fl*=colat<=(right[1]/fact)
         idx=np.where(fl)[0]
         ipix=ipix[idx]
         colat=colat[idx]*fact
         Long=Long[idx]*fact
         fl=None
         idx=None
         phi,theta=self.longcolat2grasp(Long,colat)
         r1=self.bilinearXY(mapname,phi,theta)
         return r1,ipix
   def polarplot(self,mapname,long_step=2,colat_step=10,log=None,cm=None,grayBad="#707070",adAxes=True,area=[12.,1.],cmap='hsv',vmin=-30,vmax=-0.1) :
      import numpy as np
      from matplotlib import colors as Colors
      import pylab
      from matplotlib import pyplot as plt
      from matplotlib import cm
      try :
         _cm=cm.__dict__[cmap]
      except :
         print "required cmap ",cmap," no found, replaced with 'hsv'"
         print "allowed values "
         print cm.__dict__.keys()
      if adAxes : ax = plt.subplot(111, polar=True)
      y = self.M[mapname]*1
      if log == 'e' or log == 'ln' :
         y=np.log(y)
      elif log == '10' or log == 'log10' :
         y=np.log(y)/np.log(10)
      elif log == '2' :
         y=np.log(y)/np.log(2)
      else :
         try :
            b=np.float(log) 
         except :
            b=None
         if b!= None : 
            y=np.log(y)/np.log(b)
      shape2=y.shape
      shape1=shape2[0]*shape2[1]
      idxLong=np.arange(0,shape2[1],colat_step)
      idxColat=np.arange(0,shape2[0],long_step)
      for i in idxColat : 
         tt=np.pi/180*self['long'][i][idxLong]
         cc=self['colat'][i][idxLong]
         aa=(area[0]-area[1])*(1-np.cos(np.pi/180.*cc))/2.+area[1]
         print i,cc.min(),cc.max(),aa.min(),aa.max()
         try :
            c=plt.scatter(tt,cc, c=y[i][idxLong], s=aa, cmap=_cm,edgecolor='none',vmin=vmin,vmax=vmax) 
         except :
            pass
      plt.axis([0,360,0,180])
      plt.title(self.mapname)
   def mercatore(self,reversed=False) :
      """converts a CUT map to a MERCATOR MAP"""
      import numpy as np
      import copy
      M=MercatoreMap()
      shape=self.shape
      halfrow=(self.R['n']-1)/2
      newrow=(self.R['n']-1)/2+1
      newcol=self.C['n']*2
      for k in self.M.keys() :
         M.M[k]=[]
         if reversed :
            for drow in np.arange(newrow-1,-1,-1) :
               M.M[k].append(np.concatenate((self.M[k][halfrow+drow],self.M[k][halfrow-drow])))
         else :
            for drow in np.arange(0,newrow) :
               M.M[k].append(np.concatenate((self.M[k][halfrow+drow],self.M[k][halfrow-drow])))
         M.M[k]=np.array(M.M[k])
      M.shape=[M.M['phi'].shape[0],M.M['phi'].shape[1]]
      M.M['long']=copy.deepcopy(M.M['phi'])
      M.M['colat']=copy.deepcopy(M.M['theta'])
      # long is phi where theta > 0 and phi+180 where theta < 0
      M.M['long']+=180*(M.M['colat']<0)
      # but for theta=0 the algorithm fails, so this is a patch
      idx2=np.where(np.sign(M.M['theta']).ptp(axis=1)==2)[0][0]
      for k in np.where(np.sign(M.M['theta']).ptp(axis=1)==0)[0] : 
         M.M['long'][k]=M.M['long'][idx2]
      # colat is abs(theta)
      M.M['colat']=abs(M.M['colat'])
      M.C=copy.deepcopy(self.C)
      M.R['name']='long'
      M.C['v']=M.M['long'][0]*1
      M.C['n']=len(M.C['v'])
      M.C['min']=M.C['v'].min()
      M.C['max']=M.C['v'].max()
      M.R=copy.deepcopy(self.R)
      M.R['name']='colat'
      M.R['v']=M.M['colat'][:,0]*1
      M.R['n']=len(M.R['v'])
      M.R['min']=M.R['v'].min()
      M.R['max']=M.R['v'].max()
      M.info['projection']='GRASP-CUT,Mercatore'
      return M

class MercatoreMap(MapGrid) :
   """A Mercator map: a map with x=column=longitude, y=row=colatitude"""
   def __init__(self,*args) :
      MapGrid.__init__(self)
      self.info['projection']='GRASP-Mercatore'
   def right_close_col(self,period=False,right_col_value=None) :
      MapGrid.right_close_col(self,period=period,right_col_value=None)
      if right_col_value != None :
         try :
            self.C['v'][-1]=float(right_col_value)
         except : 
            pass
         self.C['max'] =  self.C['v'].max()
         try :
            self.M['long'][:,-1]=float(right_col_value)
         except : 
            pass
   #def resample(self,Long,Colat) :
      #new=self.copy()
   def unitPixelArea(self) : 
      import numpy as np
      return np.deg2rad(self.R['delta'])*np.deg2rad(self.C['delta'])
   def radialIntegral(self,arg,returnJustIntegral=False,thetaRange=None,asStruct=False) :
      import numpy as np
      sinColat=np.sin(np.deg2rad(self.M['colat']))
      if type(arg) == type('') :
         try :
            field=self[arg]*sinColat*self.unitPixelArea()
         except :
            return None
      else :
         try :
            field=arg*sinColat*self.unitPixelArea()
         except :
            return None
      dIdtheta=np.zeros(self.R['n'])
      for itheta in range(len(dIdtheta)) : dIdtheta[itheta]=0.5*(field[itheta][1:]+field[itheta][0:-1]).sum()
      midH=0.5*(dIdtheta[1:]+dIdtheta[0:-1])
      if returnJustIntegral and thetaRange==None: return np.sort(midH).sum()
      Itheta=np.zeros(self.R['n'])
      Itheta[0]=dIdtheta[0]*1
      for itheta in range(1,len(Itheta)) : Itheta[itheta]=np.sort(midH[0:itheta+1]).sum()
      if asStruct : return {'colat':self.M['colat'].mean(axis=1),'dIdcolat':dIdtheta,'Icolat':Itheta,'method':'spherical,trapezoidal'}
      return self.M['colat'].mean(axis=1),dIdtheta,Itheta,'spherical,trapezoidal'

      
class DirectionalMapMoments_Base_Old:
   def __init__(self,method,maxorder) :
      """computs the directional moments on a map using healpix integration"""
      import numpy as np
      self.method=method
      self.TreatNaNasZero=True
      self.TreatInfAsZero=True
      self.TreatOutBoundsAsZero=True
      if Bounds == None :
         self.Bounds =np.array([-np.inf,np.inf])
      else :
         self.Bounds =Bounds
      if type(Map) == type('') :
         self.load(Map)
         return
      self.exclusionRadius = 180. 
      self.mu=None
      self.nside=None
      self.npix=None
      self.pixelArea = None
      self.n=None
      self.maxorder=maxorder
      self.Sum=np.zeros([self.maxorder,self.maxorder,self.maxorder])
      for xp in range(maxorder) :
         for yp in range(maxorder) :
            for zp in range(maxorder) :
               self.Sum[xp,yp,zp]=(m*x**xp*y**yp*z**zp).sum()
   def __getitem__(self,xp,yp,zp) :
      return self.Sum[xp,yp,zp]
   def calc_integral(self) :
      Pi4=4*np.pi
      Integral = self.Sum*self.pixelArea
      norm=Integral[0,0,0]
      Sx=Integral[1,0,0]/norm
      Sy=Integral[0,1,0]/norm
      Sz=Integral[0,0,1]/norm
      line=',%e,%e,%e,%e'%(Sx,Sy,Sz,norm)
      return Sx,Sy,Sz,norm,line
   def save(self,pickle_file) :
      import pickle
      if type(pickle_file)==type('') :
         self.filename=pickle_file
         try :
            pickle.dump(self.__dict__,open(pickle_file,'w'))
         except :
            return False
      else :
         try :
            pickle.dump(self.__dict__,pickle_file)
         except :
            return False
      return True
   def load(self,pickle_file) :
      import pickle
      if type(pickle_file)==type('') :
         self.filename=pickle_file
         try :
            self.__dict__=pickle.load(open(pickle_file,'r'))
         except :
            return False
      else :
         try :
            self.__dict__=pickle.load(pickle_file)
         except :
            return False
      return True

class DirectionalMapMoments_GRD(DirectionalMapMoments_Base_Old) :
   def __init__(self,GrdMap,exclusionRadius=None,reverse=False,excludeOut=False,NormalizedByBeam=False,asDict =True,Nside=None,ipixVec=None,maxorder=3,TreatNaNasZero=True,TreatInfAsZero=True,TreatOutBoundsAsZero=True,Bounds=None) :
      """computes the directional moments on a GRD map """
      import numpy as np
      DirectionalMapMoments_Base.__init__(self,'grd',maxorder)
   
class DirectionalMapMoments_CUT_Mercatore(DirectionalMapMoments_Base_Old) :
   def __init__(self,GrdMap,exclusionRadius=None,reverse=False,excludeOut=False,NormalizedByBeam=False,asDict =True,Nside=None,ipixVec=None,maxorder=3,TreatNaNasZero=True,TreatInfAsZero=True,TreatOutBoundsAsZero=True,Bounds=None) :
      """computes the directional moments on a CUT map managed as a Mercatore map"""
      import numpy as np
      DirectionalMapMoments_Base.__init__(self,'cut-mercatore',maxorder)

class DirectionalMapMoments_Healpix(DirectionalMapMoments_Base_Old) :
   def __init__(self,Map,exclusionRadius=None,reverse=False,excludeOut=False,NormalizedByBeam=False,asDict =True,Nside=None,ipixVec=None,maxorder=3,TreatNaNasZero=True,TreatInfAsZero=True,TreatOutBoundsAsZero=True,Bounds=None) :
      """computes the directional moments on a map using healpix integration"""
      import numpy as np
      import healpy as H
      DirectionalMapMoments_Base.__init__(self,'healpix',maxorder)
      self.TreatNaNasZero=TreatNaNasZero
      self.TreatInfAsZero=TreatNaNasZero
      self.TreatOutBoundsAsZero=TreatOutBoundsAsZero
      if Bounds == None :
         self.Bounds =np.array([-np.inf,np.inf])
      else :
         self.Bounds =Bounds
      if type(Map) == type('') :
         self.load(Map)
         return
      if exclusionRadius==None :
         self.exclusionRadius = 180. if excludeOut else 0.
      else :
         self.exclusionRadius = exclusionRadius
      self.mu=np.cos(self.exclusionRadius/180.*np.pi)
      self.nside=int(np.sqrt(len(Map)/12.)) if Nside == None else int(Nside) 
      self.npix=int(12.*self.nside**2)
      self.pixelArea = 4.*np.pi/float(12.*self.nside**2)
      if ipixVec == None :
         v=np.array(H.pix2vec(self.nside,nside2ipix(self.nside)))
         idxAllowed=np.where(v[2]<=self.mu)
      else :
         v=np.array(H.pix2vec(self.nside,ipixVec))
         idxAllowed=ipixVec
      m=Map[idxAllowed]
      x=v[0][idxAllowed]
      y=v[1][idxAllowed]
      z=v[2][idxAllowed]
      if self.TreatOutBoundsAsZero : 
         idx = np.where((m<self.Bounds[0])*(self.Bounds[1]<m))[0]
         if len(idx) > 0 : m[idx]=0.
      if self.TreatNaNasZero : 
         idx = np.where(np.isnan(m))[0]
         if len(idx) > 0 : m[idx]=0.
      if self.TreatInfAsZero : 
         idx = np.where(1-np.isfinite(m))[0]
         if len(idx) > 0 : m[idx]=0.
      self.n=len(m)
      self.maxorder=maxorder
      self.Sum=np.zeros([self.maxorder,self.maxorder,self.maxorder])
      for xp in range(maxorder) :
         for yp in range(maxorder) :
            for zp in range(maxorder) :
               self.Sum[xp,yp,zp]=(m*x**xp*y**yp*z**zp).sum()
   def calc_integral(self) :
      Pi4=4*np.pi
      Integral = iBSM.Sum*iBSM.pixelArea
      norm=Integral[0,0,0]
      Sx=Integral[1,0,0]/norm
      Sy=Integral[0,1,0]/norm
      Sz=Integral[0,0,1]/norm
      line=',%e,%e,%e,%e'%(Sx,Sy,Sz,norm)
      return Sx,Sy,Sz,norm,line
      
      

#from numpy import nan

#def cuts2matrix(self,No180Pad=No180Pad)
#;
#; converts a structure from get_cuts in a set of matrices, the output is a structure
#; theta = theta values converted from GRASP convention to usual polar convention
#; phi   = phi values from GRASP convention to usual polar convention
#; c1r, c1i, c2r, c2i = components 1 and 2 real (r) and imaginary (i) parts
#; GRASP = a structure containing the original grasp theta and phi
#;
#; NOTE = usually GRASP does not generates cuts for PHI=180deg since it is
#; nothing else than PHI=0deg, by default CUTS2MATRIX add a PHI=180deg
#; cut. To exclude this set the /No180Pad keyword
#;
#;
#if not keyword_set(No180Pad) then Pad180 = 1 else Pad180 = 0
#nphi=n_tags(cuts)
#m1=dblarr(nphi+Pad180,cuts.(0).ntheta)
#theta=m1
#phi=m1
#long=m1
#colat=m1
#c1r=m1
#c1i=m1
#c2r=m1
#c2i=m1
#vtheta=cuts.(0).theta
#vphi=dblarr(nphi+Pad180)
#vlong=dblarr(nphi+Pad180)
#for iphi=0,nphi-1 do begin
   #theta[iphi,*]=cuts.(iphi).theta
   #phi[iphi,*]=cuts.(iphi).phi
   #c1r[iphi,*]=cuts.(iphi).c1r
   #c1i[iphi,*]=cuts.(iphi).c1i
   #c2r[iphi,*]=cuts.(iphi).c2r
   #c2i[iphi,*]=cuts.(iphi).c2i
   #vphi[iphi]=cuts.(iphi).phi
#endfor
#if Pad180 ne 0 then begin
   #; performs padding of 180 deg
   #theta[nphi,*]=-reverse(cuts.(0).theta)
   #phi[nphi,*]=cuts.(0).phi+180.
   #c1r[nphi,*]=reverse(cuts.(0).c1r)
   #c1i[nphi,*]=reverse(cuts.(0).c1i)
   #c2r[nphi,*]=reverse(cuts.(0).c2r)
   #c2i[nphi,*]=reverse(cuts.(0).c2i)
   #vphi[nphi]=180.
#endif
#return,create_struct( $
   #'theta',theta $
  #,'phi',phi $
  #,'c1r',c1r $
  #,'c1i',c1i $
  #,'c2r',c2r $
  #,'c2i',c2i $
  #,'GRASP',create_struct('theta',vtheta,'phi',vphi) $
  #,'nphi',nphi $
  #,'Pad180',Pad180 $
  #,'phi0',vphi[0] $
  #,'Delta_phi',vphi[1]-vphi[0] $
  #,'theta0',vtheta[0] $
  #,'Delta_theta',vtheta[1]-vtheta[0] $
  #)
#end

#function cuts_grid,cuts
#; derives gridding parameters for a structure produced by get_cuts
   #phi0=cuts.(0).phi
   #theta0=cuts.(0).theta[0]
   #dphi = cuts.(1).phi-cuts.(0).phi
   #dtheta = cuts.(0).theta[1]-cuts.(0).theta[0]
   #xmax = -1e6
   #for kk = 0,n_tags(cuts)-1 do begin
      #_xmax=max(abs(cuts.(0).theta))
      #if _xmax gt xmax then xmax = _xmax
   #endfor
   #return,create_struct('phi0',phi0,'theta0',theta0,'dphi',dphi,'dtheta',dtheta,'thetamax',xmax)
#end

#function componentsTOstoke,component1_real,component1_imaginary,component2_real,component2_imaginary,direct=direct
#;
#; converts component maps to stokes
#;
  #E1 = component1_real^2+component1_imaginary^2
  #E2 = component2_real^2+component2_imaginary^2
  #SI = E1+E2
  #SQ = E1-E2
  #E1 = sqrt(E1)
  #E2 = sqrt(E2)
  #F1 = ATAN(component1_imaginary,component1_real)
  #F2 = ATAN(component2_imaginary,component2_real)
  #SU = 2*E1*E2*COS(F2 - F1)
  #SV = 2*E1*E2*SIN(F2 - F1)
  #return,create_struct('I',SI,'Q',SQ,'U',SU,'V',SV,'F1',F1,'F2',F2)
#end

#function cuts2cartesian,cuts,side=side,stokes=stokes
#; converts a cut in a matrix using cartesian polar coordinates
  #if not keyword_set(side) then side = 600
  #phi0=cuts.(0).phi
  #theta0=cuts.(0).theta[0]
  #dphi = cuts.(1).phi-cuts.(0).phi
  #dtheta = cuts.(0).theta[1]-cuts.(0).theta[0]
  #npix=side+1
  #xmin=-max(abs(cuts.(0).theta))
  #xmax = -1e6
  #for kk = 0,n_tags(cuts)-1 do begin
     #_xmax=max(abs(cuts.(0).theta))
     #if _xmax gt xmax then xmax = _xmax
  #endfor
  #ix0=(npix-1)/2
  #iy0=(npix-1)/2
  #xmap=dblarr(npix,npix)
  #ymap=dblarr(npix,npix)
  #for r = 0,npix-1 do xmap[r,*]=(double(indgen(npix))/double(npix-1)-0.5)*xmax*2
  #for c = 0,npix-1 do ymap[*,c]=(double(indgen(npix))/double(npix-1)-0.5)*xmax*2
  #colatmap=sqrt(xmap^2+ymap^2)
  #longmap=atan(ymap,xmap)/!dpi*180.
  #idx = where(longmap lt 0,count)
  #if count gt 0 then longmap[idx]=360+longmap[idx]
  #pt=longcolat2phitheta(longmap,colatmap)
  #rc=longcolat2rowcol(longmap,colatmap,dphi=dphi,dtheta=dtheta,phi0=phi0,theta0=theta0)
  #slm=cuts2matrix(cuts)
  #c1r=map_interpolate(rc.iphi,rc.itheta,slm.c1r)
  #c1i=map_interpolate(rc.iphi,rc.itheta,slm.c1i)
  #c2r=map_interpolate(rc.iphi,rc.itheta,slm.c2r)
  #c2i=map_interpolate(rc.iphi,rc.itheta,slm.c2i)
  #out=create_struct('x',xmap,'y',ymap,'ix0',ix0,'iy0',iy0,'colat',colatmap,'long',longmap,'phi',pt.phi,'theta',pt.theta,'iphi',rc.iphi,'itheta',rc.itheta)
  #if keyword_set(stokes) then begin
     #return,create_struct(out,'stokes',componentsTOstoke(c1r,c1i,c2r,c2i))
  #endif
  #return,create_struct(out,'power',c1r^2+c1i^2+c2r^2+c2i^2)
#end

#function cuts2healpix,nside,cuts,reversed=reversed,ipix=ipix,onlyPower=onlyPower,dbi=dbi,stokes=stokes
#; convert a cuts into an healpix max by bilinear interpolation
#; if /onlyPower returns just the power otherwise returns
#;    c1r = component 1 real part
#;    c1i = component 1 imaginary part
#;    c2r = component 2 real part
#;    c2i = component 2 imaginary part
#;    power
#; if /dbi power is converted in 10*alog10(power) (if /stokes this is not done)
#; if /stokes return stokes parameters and F1 and F2 instead of components c1, c2
#;
  #if not keyword_set(reversed) then ipix=nside2ipix(nside) else ipix=nside2ipix(nside,/reversed)
  #slm = cuts2matrix(cuts)
  #rc=ipix2rowcol(nside,ipix,slm.phi0,slm.delta_phi,slm.theta0,slm.delta_theta)
  #r1=map_interpolate(rc.iphi,rc.itheta,slm.c1r)
  #i1=map_interpolate(rc.iphi,rc.itheta,slm.c1i)
  #r2=map_interpolate(rc.iphi,rc.itheta,slm.c2r)
  #I2=map_interpolate(rc.iphi,rc.itheta,slm.c2i)
  #if keyword_set(stokes) then return,componentsTOstoke(r1,i1,r2,i2)
  #power = r1^2+i1^2+r2^2+i2^2
  #if keyword_set(dbi) then power=10d0*alog10(power)
  #if keyword_set(onlyPower) then return,power
  #return,create_struct('c1r',r1,'c2r',r2,'c1i',i1,'c2i',i2,'power',power)
#end

#function beamSumMap,map,exclusionRadius=exclusionRadius,reverse=reverse,v=v,excludeOut=excludeOut,notNormalizedByBeam=notNormalizedByBeam,asStruct=asStruct
#; computs the directional moments on a map
  #if not keyword_set(excludeOut) then excludeOut = 0
  #if not keyword_set(exclusionRadius) then $
     #if excludeOut then exclusionRadius=[180d0] else exclusionRadius=[0d0]
  #npix=n_elements(map)
  #nside=long(sqrt(npix/12.))
  #pix2vec_ring,nside,indgen(npix,/long),v
  #z = v[*,2]
  #for i = 0,2 do v[*,i]=v[*,i]*map
  #sss = dblarr(10,n_elements(exclusionRadius))
  #sss[9,*]=4d0*!dpi/double(npix)
  #sss[5,*]=npix
  #sss[6,*]=nside
  #for ir=0,n_elements(exclusionRadius)-1 do begin
     #xr=exclusionRadius[ir]
     #sss[7,ir]=xr
     #mu=cos(xr/180d0*!dpi)
     #imin = 0
     #imax = npix-1l
     #count=-1
     #if excludeOut then $
        #if xr eq 180. then imax=npix-1 else imax = min(where(mu ge z,count)) $
     #else $
        #if xr eq 0. then imin=0 else imin = min(where(mu ge z,count))
     #print,ir,xr,excludeOut,ir,count,imin,imax
     #sss[8,ir]=count
     #sss[4,ir]=total(map[imin:imax])
     #for ic = 0,2 do sss[ic,ir] = total(v[imin:imax,ic])
     #if not keyword_set(notNormalizedByBeam) then for ic = 0,2 do sss[ic,ir] = sss[ic,ir]/sss[4,ir]
     #sss[3,ir]=sqrt(sss[0,ir]^2+sss[1,ir]^2+sss[2,ir]^2)  
     #for ic = 0,2 do sss[ic,ir]=sss[ic,ir]/sss[3,ir]
  #endfor
  #if not keyword_set(asStruct) then return,sss
  #; computes the polar deflection
  #polar_deflection=dblarr(n_elements(exclusionRadius))
  #longitude_deflection=dblarr(n_elements(exclusionRadius))
  #for ir=0,n_elements(exclusionRadius)-1 do begin
     #polar_deflection=acos(sss[2,ir])*180d0/!dpi
     #normxy = sqrt(total(sss[0:1,ir]^2))
     #longitude_deflection=atan(sss[1,ir]/normxy,sss[0,ir]/normxy)*180d0/!dpi
  #endfor
  #return,create_struct($
     #'vSx',sss[0,*] $
     #,'vSy',sss[1,*] $
     #,'vSz',sss[2,*] $
     #,'S',sss[3,*]   $
     #,'beam_sum',sss[4,*] $
     #,'npix',long(sss[5,*]) $
     #,'nside',long(sss[6,*]) $
     #,'exclusionRadius',sss[7,*] $
     #,'excludedPixels',long(sss[8,*]) $
     #,'pixelArea',sss[9,*] $
     #,'deflection_polar_deg',polar_deflection $
     #,'deflection_longitude_deg',longitude_deflection $
     #) 
#end

#function beamSums,nside,cuts,exclude_angle=exclude_angle,map=map,asStruct=asStruct
#; computes Sx, Sy. Sz (directional integrals) for a beam map
  #map=cuts2healpix(nside,cuts,ipix=ipix,/onlyPower)
  #pix2vec_ring,nside,ipix,v
  #if keyword_set(exclude_angle) then begin
     #print,exclude_angle
     #ang = 180./dpi*acos(v[*,2])
     #idx = where(ang > exclude_angle,count) 
     #if count gt 0 then begin
        #v1=dblarr(n_elements(idx),3)
        #for i=0,2 do v1[*,i]=v[idx,i]
        #v=v1
     #endif else begin
       #print,"Error all pixels removed"
       #return,0
     #endelse
   #endif
  #sss = dblarr(7)
  #sss[6]=nside
  #sss[5]=12d0*double(nside)^2
  #sss[4]=total(map) 

  #; returns the versors
  #for i=0,2 do sss[i] = (total(v[*,i]*map))/sss[4]

  #; normalize
  #sss[3]=sqrt(total(sss[0:2]^2))
  #for i=0,2 do sss[i] = sss[i]/sss[3]

  #if not keyword_set(asStruct) then return,sss

  #; computes the polar deflection
  #polar_deflection_deg=acos(sss[2]/sss[3])*180d0/!dpi
  #normxy = sqrt(total(sss[0:1]^2))
  #longitude_deflection_deg=atan(sss[1]/normxy,sss[0]/normxy)*180d0/!dpi
  #return,create_struct( $
       #'vSx',sss[0] $
      #,'vSy',sss[1] $
      #,'vSz',sss[2] $
      #,'S',sss[3]   $
      #,'beam_sum',sss[4] $
      #,'npix',long(sss[5]) $
      #,'nside',long(sss[6]) $
      #,'deflection_polar_deg',polar_deflection_deg $
      #,'deflection_longitude_deg',longitude_deflection_deg $
     #)
#end

#function beamSumS2,lcuts,hcuts,nside=nside,map=map,returnFirst=returnFirst,returnSecond=returnSecond
#; computes Sx, Sy. Sz for a beam using two maps, 
#; lcuts = a lowress map of cuts
#; hcuts = an highres map of cuts
#; /returnFirst returns just the high resolution map (no summation)
#; /returnSecond returns just the second map (no summation)
#;
  #; high resolution integral
  #hpa=cuts_grid(hcuts)
  #if not keyword_set(nside) then nside = 1024l
  #map=dblarr(12l*nside*nside)
  #radius = -1e6
  #for kk = 0,n_tags(hcuts)-1 do begin
     #_xmax=max(abs(hcuts.(kk).theta))
     #if _xmax gt radius then radius = _xmax
  #endfor
  #query_disc,nside,[0.,0.,1.],radius,ipix,/deg,/inclusive
  #slm = cuts2matrix(hcuts)
  #rc=ipix2rowcol(nside,ipix,hpa.phi0,hpa.dphi,hpa.theta0,hpa.dtheta)
#;dphi=hpa.dphi,dtheta=hpa.dtheta,phi0=hpa.phi0,theta0=hpa.theta0)
  #r1=map_interpolate(rc.iphi,rc.itheta,slm.c1r)
  #i1=map_interpolate(rc.iphi,rc.itheta,slm.c1i)
  #r2=map_interpolate(rc.iphi,rc.itheta,slm.c2r)
  #I2=map_interpolate(rc.iphi,rc.itheta,slm.c2i)
  #map[ipix] = r1^2+r2^2+i1^2+i2^2
  #if keyword_set(returnFirst) then return,map
  #query_disc,nside,[0.,0.,-1.],180.-radius,ipix,/deg
#;ipix=nside2ipix(nside)
  #slm = cuts2matrix(lcuts)
  #lpa=cuts_grid(lcuts)
  #rc=ipix2rowcol(nside,ipix,lpa.phi0,lpa.dphi,lpa.theta0,lpa.dtheta)
  #r1=map_interpolate(rc.iphi,rc.itheta,slm.c1r)
  #i1=map_interpolate(rc.iphi,rc.itheta,slm.c1i)
  #r2=map_interpolate(rc.iphi,rc.itheta,slm.c2r)
  #I2=map_interpolate(rc.iphi,rc.itheta,slm.c2i)
  #map[ipix] = r1^2+r2^2+i1^2+i2^2
  #if keyword_set(returnSecond) then return,map
  #return,beamSumMap(map,/reverse,/asStruct)
#end

#function radialDependence,mapname,listRadiiDeg=listRadiiDeg
  #if not keyword_set(listRadiiDeg) then listRadiiDeg=[0.1,0.5,1.,1.5,2.,2.5,5.,7.5,10.,20.,30.,40.,50,60,70,80.,85.,90,100,110,120,130,140,150,160,170,180]
  #read_fits_map,mapname,mapX
  #xxx = beamSumMap(mapX,exclusionRadius=listRadiiDeg,/excludeOut)
  #sss=xxx
  #for ic=0,2 do sss[ic,*]=sss[ic,*]/sss[4,*]
  #sss[3,*]=sqrt(sss[0,*]^2+sss[1,*]^2+sss[2,*]^2)
  #for ic=0,2 do sss[ic,*]=sss[ic,*]/sss[3,*]
  #radius=sss[7,*]
  #dump=sss[3,*]
  #polar_deflection=acos(sss[2,*])/!dpi*180.*60.
  #longitudinal_deflection=atan(sss[1,*],sss[0,*])/!dpi*180.
  #return,create_struct('name',mapname,'long_def',longitudinal_deflection,'pol_def',polar_deflection,'dump',dump,'radius',radius)
#end


#function readgrd, fileinput
#; reads a grd file
#; (deprecated)
#xs = 0.d
#xs = 0.d
#ye = 0.d
#ye = 0.d
#str='!5 '
#ktype = 1l      ; --> data type format
#nset = 1l       ; --> number of beams in the file
#icomp = 1l      ; --> field component
#ncomp = 1l      ; --> number of components
#igrid = 1l      ; --> type of field grid
#ix = 1l         ; --> center of the beam
#iy = 1l         ;     (ix,iy)
#c1 = 0.d
#c2 = 0.d
#c3 = 0.d
#c4 = 0.d
#nx = 1l
#ny = 1l
#klimit = 1l
#openr,1,fileinput
#for line=0,100 do begin
  #if (strtrim(str,2) ne '++++') then begin
    #readf,1,str
    #print,str
  #endif else begin
    #goto, jump1
  #endelse
#endfor
#jump1: readf,1,ktype
#readf,1,nset,icomp,ncomp,igrid
#readf,1,ix,iy
#readf,1,xs,ys,xe,ye
#readf,1,nx,ny,klimit
#dx = (xe - xs)/(nx-1)
#x = findgen(nx)*dx + xs
#dy = (ye - ys)/(ny-1)
#y = findgen(ny)*dy + ys
#print,'Reading ', fileinput
#print,'grid of ', nx,' x ', ny,' points'
#c1r = dblarr(nx,ny)
#c1i = dblarr(nx,ny)
#c2r = dblarr(nx,ny)
#c2i = dblarr(nx,ny)
#for i=0,nx-1 do begin
  #for j=0,ny-1 do begin
    #readf,1,c1,c2,c3,c4
    #c1r(j,i) = c1
    #c1i(j,i) = c2
    #c2r(j,i) = c3
    #c2i(j,i) = c4
  #endfor
#endfor
#close,1
#power = c1r^2 + c1i^2 + c2r^2 + c2i^2
#res = { x : x , $
        #y : y , $
        #power : power $
      #}
#return, res
#end

def readgrd(inputfile) :
   """
   ; reads a grd file
   ; (deprecated)
   Reference document: 
   LFI BEAMS DELIVERY: FORMAT SPECIFICATIONS
   M. Sandri
   PL-LFI-PST-TN-044, 1.0,July 2003
   """
   #xs = 0.d
   #xs = 0.d
   #ye = 0.d
   #ye = 0.d
   #str='!5 '
   #ktype = 1l      ; --> data type format
   #nset = 1l       ; --> number of beams in the file
   #icomp = 1l      ; --> field component
   #ncomp = 1l      ; --> number of components
   #igrid = 1l      ; --> type of field grid
   #ix = 1l         ; --> center of the beam
   #iy = 1l         ;     (ix,iy)
   #c1 = 0.d
   #c2 = 0.d
   #c3 = 0.d
   #c4 = 0.d
   #nx = 1l
   #ny = 1l
   #klimit = 1l
   try :
      h=open(inputfile,'r').readlines()
   except :
      print "File %s not found"%inputfile
      return
   # removes the new line
   for i in range(len(h)) :
      h[i] = h[i].split('\n')[0].split('\r')[0]
   currentline=0
   while(h[currentline] != '++++') : 
      currentline+=1
      if currentline == len(h):
         print "Error marker ++++ not found"
         return h
   infos=h[0:currentline]
   currentline +=1
   print h[currentline]
   ktype = int(h[currentline])
   currentline +=1
   print h[currentline]
   ll = h[currentline].split()
   nset  = int(ll[0])
   icomp = int(ll[1]) 
   ncomp = int(ll[2])
   igrid = int(ll[3])
   currentline +=1
   print h[currentline]
   ll = h[currentline].split()
   ix = int(ll[0])
   iy = int(ll[1]) 
   currentline +=1
   print h[currentline]
   ll = h[currentline].split()
   xs    = float(ll[0])
   ys    = float(ll[1]) 
   xe    = float(ll[2])
   ye    = float(ll[3])
   currentline +=1
   print h[currentline]
   ll = h[currentline].split()
   nx = int(ll[0])
   ny = int(ll[1]) 
   klimit = int(ll[2]) 
   dx = (xe - xs)/float(nx-1)
   xcen=ix*dx
   x = np.arange(nx)*dx + xs+xcen
   dy = (ye - ys)/float(ny-1)
   ycen=iy*dy
   y = np.arange(ny)*dy + ys+ycen
   print 'Reading ', inputfile
   print 'grid of ', nx,' x ', ny,' points'
   print 'ix ', ix,' iy ', iy
   c1r = np.zeros([ny,nx])
   c1i = np.zeros([ny,nx])
   c2r = np.zeros([ny,nx])
   c2i = np.zeros([ny,nx])
   for j in range(ny) :
      for i in range(nx) :
         currentline +=1
         ll = h[currentline].split()
         c1r[j,i] = float(ll[0])
         c1i[j,i] = float(ll[1])
         c2r[j,i] = float(ll[2])
         c2i[j,i] = float(ll[3])
   return {'x':x,'y':y,'ri':c1r,'r2':c2r,'i1':c1i,'i2':c2i,'power':c1r**2 + c1i**2 + c2r**2 + c2i**2,'infos':infos}

class GridMap(MapGrid) :
   def __init__(self,inputfile,skiplines=0,silent=False,closeColumn=False,Pickle=False,nodata=False,addPolar=True,addUV=True,addP1P2=True) :
      MapGrid.__init__(self)
      if Pickle :
         self.load(inputfile)
         return
      self.info['grd_file']=inputfile.strip()
      self.info['projection']='GRASP-GRD'
      self.info['ReferenceDocument']="LFI BEAMS DELIVERY: FORMAT SPECIFICATIONS\nM. Sandri\nPL-LFI-PST-TN-044, 1.0,July 2003"
      if inputfile.strip() == '' : return
      self.get_grd(inputfile,nodata=nodata)#,CounterPhi=CounterPhi,silent=silent,useCounterPhi=useCounterPhi) 
      if closeColumn :
         self.right_close_col()
         for k in self.M.keys() :
            if k!='_row_values' and k!='_col_values' and k!='_row_index' and k!='_col_index' :
               for r in range(self.R['n']) :
                  self.M[k][r][-1]=self.M[k][self.R['n']-1-r][0]
      if addPolar :
         self.M['colat']=np.rad2deg(np.arcsin((self.M['_col_values']**2+self.M['_row_values']**2)**0.5))
         self.M['long']=np.mod(np.rad2deg(np.arctan2(self.M['_row_values'],self.M['_col_values'])),360.)
      if addUV :
         self.newmap(self.R['name'],unit='',value=self['_row_values'])
         self.newmap(self.C['name'],unit='',value=self['_col_values'])
      if addP1P2 :
         self.newmap('p1',unit='',value=self['r1']**2+self['i1']**2)
         self.newmap('p2',unit='',value=self['r2']**2+self['i2']**2)
   def end_of_header_marker(self) :
      """returns the string used as a marker of end of header"""
      return '++++'
   def formatGrasp(self) : return {'float':' %17.10e','int':' %11d'}
   def get_grd(self,inputfile,nodata=False) :
      """
      reads a grd file
      Reference document: 
      LFI BEAMS DELIVERY: FORMAT SPECIFICATIONS
      M. Sandri
      PL-LFI-PST-TN-044, 1.0,July 2003
      
      Beware: the faster moving index is the column,
      in IDL readgrd inverts rows with columns to take in
      columns reading columns as they would be the slowest index.
      """
      import sys
      import numpy as np
      import copy
      try :
         h=open(inputfile,'r').readlines()
         self.mapname=inputfile
         self.info['inputfile']=inputfile
      except :
         print "File %s not found"%inputfile
         return
      # removes the new line and other special characters
      for i in range(len(h)) :
         h[i] = h[i].split('\n')[0].split('\r')[0].strip()
      currentline=0
      while(h[currentline].strip() != self.end_of_header_marker()) : 
         currentline+=1
         if currentline == len(h):
            print "Error marker %s not found" % self.end_of_header_marker()
            return h
      self.info['header']=copy.deepcopy(h[0:currentline])
      currentline +=1
      self.info['ktype']=int(h[currentline])
      currentline +=1
      print h[currentline].strip()
      ll = h[currentline].split()
      self.info['nset']=int(ll[0])
      self.info['icomp']= int(ll[1]) 
      self.info['ncomp']= int(ll[2])
      self.info['igrid'] = int(ll[3])
      currentline +=1
      print h[currentline].strip()
      ll = h[currentline].split()
      self.info['ix'] = int(ll[0])
      self.info['iy'] = int(ll[1]) 
      currentline +=1
      print h[currentline].strip()
      ll = h[currentline].split()
      self.info['xs']    = float(ll[0])
      self.info['ys']    = float(ll[1]) 
      self.info['xe']    = float(ll[2])
      self.info['ye']    = float(ll[3])
      currentline +=1
      print h[currentline].strip()
      ll = h[currentline].split()
      self.info['nx'] = int(ll[0])
      self.info['ny'] = int(ll[1]) 
      self.info['klimit'] = int(ll[2]) 
      #computed parameters
      # X are columns
      self.info['dx'] = (self.info['xe']-self.info['xs'])/float(self.info['nx']-1)
      self.info['xcen'] = self.info['dx']*self.info['ix']
      self.set_col_scale('U','uv',np.arange(self.info['nx'])*self.info['dx']+ self.info['xs']+self.info['xcen'])
      # Y are rows
      self.info['dy'] = (self.info['ye']-self.info['ys'])/float(self.info['ny']-1)
      self.info['ycen'] = self.info['dy']*self.info['iy']
      self.info['grd_file']=inputfile.strip()
      self.set_row_scale('V','uv',np.arange(self.info['ny'])*self.info['dy']+ self.info['ys']+self.info['ycen'])
      print 'Reading ', inputfile
      print 'grid of ', self.info['nx'],' x ', self.info['ny'],' points'
      if nodata :
         return
      # maps used for debug
      self.newmap('_line_index','')
      # compoenent maps
      self.newmap('r1','')
      self.newmap('i1','')
      self.newmap('r2','')
      self.newmap('i2','')
      self.newmap('power','')
      #self.newmap('rho1','')
      #self.newmap('rho2','')
      #self.newmap('phi1','')
      #self.newmap('phi2','')
      for r in range(self.R['n']) :
         for c in range(self.C['n']) :
            currentline +=1
            self.M['_row_values'][r,c]=self.R['v'][r]*1.
            self.M['_row_index'][r,c]=r*1.
            self.M['_col_values'][r,c]=self.C['v'][c]*1.
            self.M['_col_index'][r,c]=c*1.
            self.M['_line_index'][r,c]=currentline*1.
            ll = h[currentline].split()
            self.M['r1'][r,c]=float(ll[0])
            self.M['i1'][r,c]=float(ll[1])
            self.M['r2'][r,c]=float(ll[2])
            self.M['i2'][r,c]=float(ll[3])
            self.M['power'][r,c]=float(ll[0])**2+float(ll[1])**2+float(ll[2])**2+float(ll[3])**2
            #self.M['rho1'][r,c]=(float(ll[0])**2+float(ll[1])**2)**0.5
            #self.M['rho2'][r,c]=(float(ll[2])**2+float(ll[3])**2)**0.5
            #self.M['phi1'][r,c]=np.arctan2(float(ll[1]),float(ll[0]))
            #self.M['phi2'][r,c]=np.arctan2(float(ll[3]),float(ll[2]))
   def UV(self,Vectors=False) :
      """returns the U, V matrix
         if Vectors=True the values of R and C are returned
      """
      if Vectors :
         return self.C['v'],self.R['v']
      return self.M['_col_values'],self.M['_row_values']
   def thetaUVphiUV(self,deg=True) :
      """returns the thetaUV, phiUV matrix
      """
      return UV2thetaUVphiUV(self.M['_col_values'],self.M['_row_values'],deg=deg)
   def cart3d(self) :
      """returns the x,y,z matrices
      """
      import numpy as np
      theta,phi=UV2thetaUVphiUV(self.M['_col_values'],self.M['_row_values'],deg=False)
      return np.cos(phi)*np.sin(theta),np.sin(phi)*np.sin(theta),np.cos(theta)
   def recompose_header(self,*arg,**karg) :
      "keywords: inhdr header in input"
      import copy
      fmtF=self.formatGrasp()['float']
      fmtI=self.formatGrasp()['int']
      hdr=None
      if karg.has_key('inhdr') :
         hdr=copy.deepcopy(karg['inhdr'])
      if hdr == None :
         hdr=copy.deepcopy(self.info['header'])
      if len(arg)>0 :
         if type(arg[0]) == type('') :
            hdr.append(arg[0])
         else :
            for k in arg[0] :
               hdr.append(k)
      hdr.append('In the lines after the header marker defined by 4 "+" characters')
      hdr.append('line 1 : ktype')
      hdr.append('line 2 : nset icomp ncomp igrid')
      hdr.append('line 3 : ix iy')
      hdr.append('line 4 : xs ys xe ye')
      hdr.append('line 5 : nx ny klimit')
      hdr.append(self.end_of_header_marker())
      ll=''
      ll+=fmtI%int(fmtI%int(self.info['ktype']))
      hdr.append(ll.upper())
      ll=''
      ll+=fmtI%int(self.info['nset'])
      ll+=fmtI%int(self.info['icomp'])
      ll+=fmtI%int(self.info['ncomp'])
      ll+=fmtI%int(self.info['igrid'])
      hdr.append(ll.upper())
      ll=''
      ll+=fmtI%int(self.info['ix'])
      ll+=fmtI%int(self.info['iy'])
      hdr.append(ll.upper())
      ll=''
      ll+=fmtF%float(self.info['xs'])
      ll+=fmtF%float(self.info['ys'])
      ll+=fmtF%float(self.info['xe'])
      ll+=fmtF%float(self.info['ye'])
      hdr.append(ll.upper())
      ll=''
      ll+=fmtI%int(self.info['nx'])
      ll+=fmtI%int(self.info['ny'])
      ll+=fmtI%int(self.info['klimit'])
      hdr.append(ll.upper())
      return hdr
   def recompose_map(self,*arg,**karg) :
      import copy
      fmtF=self.formatGrasp()['float']
      fmtI=self.formatGrasp()['int']
      if len(arg) > 0 :
         lst=copy.deepcopy(arg[0])
      else :
         lst=[]
      for r in range(self.R['n']) :
         for c in range(self.C['n']) :
            ll = fmtF%(self.M['r1'][r,c])
            ll += fmtF%(self.M['i1'][r,c])
            ll += fmtF%(self.M['r2'][r,c])
            ll += fmtF%(self.M['i2'][r,c])
            lst.append(ll.upper())
      return lst
   def FourColumnsPower(self,power1Name='p1',power2Name='p2',powerName='power') :
      "a FourColumns map has r1=sqrt(p1), i1=0, r2=sqrt(p2), i2=0"
      new=self.copy()
      new.info['ktype']=1
      if self.M.has_key(power1Name) and self.M.has_key(power1Name) :
         new.info['ncomp']=2
         new.M['r1']=self[power1Name]**0.5
         new.M['r2']=self[power2Name]**0.5
         new.M['i1']=self[power1Name]*0
         new.M['i2']=self[power1Name]*0
      elif self.M.has_key(power) :
         new.info['ncomp']=1
         new.M['r1']=self[power]**0.5
         new.M['r2']=self[power]*0
         new.M['i1']=self[power]*0
         new.M['i2']=self[power]*0
      else :
         print "the map shall contain ",power1Name,power2Name," or ",powerName
         return
      return new
   def ipix2longcolat(self,nside,ipix,nest=False,deg=True) :
      """ converts an healpix ipix (ring) into index of phi and of theta in the matrix"""
      from healpy import pix2ang
      import numpy as np
      colat,_long=pix2ang(nside,ipix,nest=nest)
      if deg : return _long*180./np.pi,colat*180./np.pi
      return _long,colat
   def nside2ipix(self,nside,Reversed=False) :
      """ converts nside into a list of pixels (ring) 
         reversed = True means the orderring is reversed
      """
      return nside2ipix(nside,Reversed=Reversed)
   def healpix(self,nside,mapname='power',nest=False,Reversed=False,colatrange=None,returnAll=False,usePeriodicalInterpolator=True) :
      """converts to healpix or a stack of healpix maps of given nside
         colatrange=None , takes all the map
         colatrange=']a,b['
         colatrange='[a,b['
         colatrange=']a,b]'
      """
      import numpy as np
      import healpy as H
      if colatrange==None :
         ipix=self.nside2ipix(nside,Reversed=Reversed)
         phiUV,thetaUV = self.ipix2longcolat(nside,ipix,deg=False)
      else :
         fact=180./np.pi
         prs=(colatrange.strip()).split(',')
         left = [prs[0][0],float(prs[0][1:])]
         right = [prs[1][-1],float(prs[1][0:-1])]
         NPEQ=12*nside/2
         print left,right
         ipixmin=H.ang2pix(nside,left[1]/fact,0)-NPEQ
         ipixmax=H.ang2pix(nside,right[1]/fact,0)+NPEQ
         if ipixmin < 0 : ipixmin=0
         if ipixmax > 12*nside**2-1 : ipixmax=12*nside**2
         ipix = np.arange(ipixmin,ipixmax)
         colat,Long = H.pix2ang(nside,ipix)
         fl=np.ones(len(colat))
         if left[1] == ']' :
            fl*=(left[1]/fact)<colat
         else :
            fl*=(left[1]/fact)<=colat
         if right[1] == '[' :
            fl*=colat<(right[1]/fact)
         else :
            fl*=colat<=(right[1]/fact)
         idx=np.where(fl)[0]
         ipix=ipix[idx]
         thetaUV=colat[idx]
         phiUV=Long[idx]
         fl=None
         idx=None
         colat=None
         Long=None
      U,V = thetaUVphiUV2UV(thetaUV,phiUV,deg=False)
      r1=self.bilinearXY(mapname,U,V)
      if returnAll : return r1,ipix,U,V,thetaUV,phiUV
      return r1,ipix
   def healpix_pixelArea(self,nside) :
      import numpy as np
      return 4*pi/(12.*65536.**2)
   def healpix_integral(self,nside,mapname,Reversed=False,colatrange=None) :
      import numpy as np
      h=self.healpix(nside,mapname=mapname,nest=False,Reversed=Reversed,colatrange=colatrange,returnAll=False)
      pixelaArea=4*np.pi/(12.*65536.**2)
      return sort(h[0]).sum()*pixelaArea
   def maximumRadius(self) :
      "returns the largest possible radius for an inscribed circle"
      a=self.C['v'].ptp()/2.
      b=self.R['v'].ptp()/2.
      return a if a < b else b
   def circularMask(self,*arg) :
      "returns a mask for the largest possible circle inscribed in the map"
      mask = np.array(self.radius()<=self.maximumRadius(),dtype='int')
      if len(arg) == 0 : return mask
      try :
         self.M[arg]*=mask
      except :
         print arg," not a valid name"
   def unitPixelArea(self) : 
      import numpy as np
      return self.R['delta']*self.C['delta']
   def radialIntegral(self,arg,method='planar,simpson',returnJustIntegral=False,asStruct=False,nRadii=51) :
      """
         returns a radial integral from 0 up to the maximum possible radius divided in nRadii steps
         integration method given by "method" keyword, default 'raster,planar,trapezoidal'
            raster,planar,trapezoidal : 
               raster over the rows of the grd map forcing to zero any sample outside the wanted ring
            planar,trapezoidal : 
               uses direct trapezoidal integration
            planar,simpson : 
               uses direct simpson integration
         some test shows the difference between planar,simpson and planar,trapezoidal is order of magnitudes larger
         than the difference between planar,trapezoidal and raster,planar,trapezoidal
         so default is planar,simpson
      """
      import numpy as np
      if returnJustIntegral :
         mm=self.circularMask()*self[arg] if type(arg) == type('') else self.circularMask()*arg
         return self.simpson2d(mm) if 'planar,simpson' else self.trapz2d(mm)
      oradius=np.arange(nRadii)/float(nRadii-1)*self.maximumRadius()
      oradius[nRadii-1]=self.maximumRadius()
      _r=self.radius()
      mm=self[arg] if type(arg) == type('') else arg
      Itheta=np.zeros(len(oradius))
      #
      if method=='raster,planar,trapezoidal' :
         for j in range(len(oradius)) :
            u=mm*1
            u[np.where(_r>oradius[j])]=0.
            acc=np.zeros(u.shape[0])
            for r in range(u.shape[0]) :
               x=u[r,1:]
               acc[r]=((x[1:]+x[0:-1])/2.).sum()
            Itheta[j]=(acc[1:]+acc[0:-1]).sum()/2.
         Itheta*=self.unitPixelArea()
      #
      elif method=='planar,trapezoidal' :
         for j in range(len(oradius)) :
            Itheta[j]=self.trapz2d(mm,outerCut=oradius[j])
      #
      elif method=='planar,simpson' :
         for j in range(len(oradius)) :
            Itheta[j]=self.simpson2d(mm,outerCut=oradius[j])
      #
      else :
         print "Unknown integration method %s"%method
         return None
      #
      if asStruct : return {'colat':oradius.mean(axis=1),'dIdcolat':Itheta[1:]-Itheta[0:-1],'Icolat':Itheta,'method':method}
      return oradius,Itheta[1:]-Itheta[0:-1],Itheta,method

if __name__=='__main__' :
   import numpy as np

   print "load"
   h={}
   h['30']=GraspMap('mb/FB_LFI27_SWE_X_FM_1-0.cut',12)
   h['27']=GraspMap('outband/FB_LFI27_SWE_X_FM_1-0_27.cut',0)
   h['33']=GraspMap('outband/FB_LFI27_SWE_X_FM_1-0_33.cut',0)

   lambda0=1/30.
   lambda1=1/27.


   WL=np.array([1/27.,1/30.])

   deltaWL=WL[1]-WL[0]

   print "parameters"
   PARS={}
   for k in ['r1','i1','r2','i2'] :
      A=h['27'].M[k]*1
      B=h['30'].M[k]*1
      A.shape=A.shape[0]*A.shape[1]
      B.shape=B.shape[0]*B.shape[1]
   #   PARS[k]=np.polyfit(WL,np.array([A,B]),1)
      PARS[k]=np.array([(B-A)/deltaWL,A])

   C=h['33'].M['power']*1
   C.shape=C.shape[0]*C.shape[1]

   print "lambda interpolate"
   ipt = {}
   for nu in [27,28,29,30,31,32,33] :
      ipt[str(nu)] ={}
      for k in ['r1','i1','r2','i2'] :
         #ipt[str(nu)][k]=np.polyval(PARS[k],1/float(nu))
         ipt[str(nu)][k]=PARS[k][0]*(1/float(nu)-WL[0])+PARS[k][1]
      ipt[str(nu)]['power']=ipt[str(nu)]['r1']**2+ipt[str(nu)]['i1']**2+ipt[str(nu)]['r2']**2+ipt[str(nu)]['i2']**2
      print nu,1/float(nu),10*np.log10(ipt[str(nu)]['power'].max())
   print 10*np.log10(C.max())


   R=h['27'].M['power']*1 ; R.shape=R.shape[0]*R.shape[1]
   G=h['30'].M['power']*1 ; G.shape=G.shape[0]*G.shape[1]
   B=h['33'].M['power']*1 ; B.shape=B.shape[0]*B.shape[1]
   I=np.arange(len(R))*3

   sys.exit()
   from matplotlib import pyplot as plt
   plt.close('all')
   plt.figure()
   plt.plot(I-3./(33.-27.),10*np.log10(R),'r.')
   plt.plot(I+0./(33.-27.),10*np.log10(G),'g.')
   plt.plot(I+3./(33.-27.),10*np.log10(B),'b.')
   plt.show()
