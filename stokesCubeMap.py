VERSION='V 1.4 - 2014 Jun 4 - '

from grid2d import MapGrid

class stokesCubeMap :
   def __init__(self,*Arg) :
      arm_alias={'x':'S','y':'M'}
      self._nameRIMO=None
      self._angularCut=None
      self._Nsamples=-1
      self.File = []
      self.Component = []
      self.Instrument = []
      self.Channel = []
      self.Horn = []
      self.Arm = []
      self.FreqCode = []
      self.Version = []
      self.FreqMHz = []
      self.Polarization = []
      self.Beamdata=[]
      self.BeamdataQ=[]
      self.BeamdataU=[]
      self.BeamdataV=[]
      self.header={}
      if len(Arg) < 2 : 
         self._arm=None
         self._LstFileName=None
         return
      self._arm=Arg[0]
      self._LstFileName=Arg[1]
      for k in open(self._LstFileName,'r') :
         kk=k.split('/')[-1].split('.')[0].split('_')
         if kk[4]==self._arm :
            self.File.append(k.split('.')[0]+'.stokes')
            self.Component.append(kk[0])
            self.Instrument.append(kk[1])
            self.Channel.append(kk[2])
            self.Horn.append(kk[3])
            self.Arm.append(kk[4])
            self.FreqCode.append(kk[5])
            self.Version.append(kk[6])
            self.FreqMHz.append(float(kk[5][1:]))
            self.Polarization.append(arm_alias[kk[4]])
      for k in self.keys() : self.__dict__[k]=np.array(self.__dict__[k])
   def fill_from_fits(self) :
      import numpy as np
      self.Beamdata=[]
      self.BeamdataQ=[]
      self.BeamdataU=[]
      self.BeamdataV=[]
      self.header={}
      isFirst=True
      for k in self.File :
         print k
         p,x,b,q,u,v = self.get_fits(k)
         self.Beamdata.append(b)
         self.BeamdataQ.append(q)
         self.BeamdataU.append(u)
         self.BeamdataV.append(v)
         if isFirst :
            self.header={'p':p,'x':x}
            isFirst=False
      for k in ['','Q','U','V'] :
         nn='Beamdata'+k
         self.__dict__[nn]=np.array(self.__dict__[nn]).transpose()
      self._Nsamples=self.Beamdata.shape[0]
   def fitsType(self) :
      if self.header['x'].has_key('HIERARCH Nx') :
         return 'grd'
      else :
         return 'cut'
   def getGeometry(self) :
      import numpy as np
      if self.fitsType()=='grd' :
         return
      else :
         geom={}
         for k in ['objType','Ntheta','Nphi','Mintheta','Maxtheta'] :
            geom[k]=self.header['x']['HIERARCH '+k]
         geom['colat']=np.rad2deg(np.arange(geom['Ntheta'])*(geom['Maxtheta']-geom['Mintheta'])/float(geom['Ntheta']-1)+geom['Mintheta'])
         geom['long']=np.arange(geom['Nphi'])/float(geom['Nphi']-1)*360.
         return geom
   def apply_angularCut(self,angularCut) :
      import numpy as np
      if self.fitsType()=='grd' : return
      if angularCut == None : return
      self._angularCut=angularCut
      gg=self.getGeometry()
      idx=np.where(gg['colat']<self._angularCut)[0].max()
      imax=idx*gg['Nphi']+gg['Nphi']
      self.Beamdata[:imax,:]=0
      self.BeamdataQ[:imax,:]=0
      self.BeamdataU[:imax,:]=0
      self.BeamdataV[:imax,:]=0
   def __len__(self) : 
      return len(self.FreqMHz)
   def __getitem__(self,this) :
      return self.__dict__[this]
   def keys(self) :
      l=[]
      for k in self.__dict__.keys() :
         if k[0]!='_' : l.append(k)
      return l
   def copy(self) :
      import copy
      return copy.deepcopy(self)
   def argslice(self,idx) :
      out=self.copy()
      for k in self.keys() :
         out.__dict__[k]=self.__dict__[k][idx]
   def get_fits(self,File,ihdu=1,fitsfile=None) :
      import pyfits
      if fitsfile == None :
         t=pyfits.open(File)
      else :
         t=fitsfile
      p=t[0].header
      x=t[hdu].header
      b=t[hdu].data.field('Beamdata')
      q=t[hdu].data.field('BeamdataQ')
      u=t[hdu].data.field('BeamdataU')
      v=t[hdu].data.field('BeamdataV')
      if fitsfile == None :
         t.close()
      return p,x,b,q,u,v
   def _template(self) :
      import numpy as np
      return self.Beamdata.shape
   def isGRD(self) :
      "Returns True if the file is a GRD"
      return self.header['x'].has_key('HIERARCH Nx')
   def Nelements(self) :
      return self.header['x']['NAXIS2']
   def NrowsNcols(self) :
      "returns the number of rows and cols"
      if self.isGRD() :
         return (self.header['x']['HIERARCH Ny'],self.header['x']['HIERARCH Nx'])
      return (self.header['x']['HIERARCH Ntheta'],self.header['x']['HIERARCH Nphi'])
   def reshape2d(self) :
      "reshape to 2d the matrices"
      for k in self.keys() :
         self.__dict__[k].shape=self.NrowsNcols()
   def reshape1d(self) :
      "reshape to 1d the matrices"
      for k in self.keys() :
         self.__dict__[k].shape=self.__dict__[k].size
   def rows_cols_idx(self) :
      "returns the maps of rows and cols index"
      import numpy as np
      nr,nc=self.NrowsNcols()
      row=np.zeros([nr,nc])
      for k in range(nr) : row[k,:]=k
      row.shape=self.Nelements()
      col=np.zeros([nr,nc])
      for k in range(nc) : row[:,k]=k
      col.shape=self.Nelements()
      return row,col
   def rows_cols_values(self) :
      "returns the maps of rows and cols values"
      import numpy as np
      nr,nc=self.NrowsNcols()
      row=np.zeros([nr,nc])
      for k in range(nr) : row[k,:]=k
      row.shape=self.Nelements()
      col=np.zeros([nr,nc])
      for k in range(nc) : row[:,k]=k
      col.shape=self.Nelements()
      return row,col
   def interp(self,idx,FreqMHz) :
      import numpy as np
      b=np.interp(FreqMHz,self.FreqMHz,self.Beamdata[idx])
      q=np.interp(FreqMHz,self.FreqMHz,self.BeamdataQ[idx])
      u=np.interp(FreqMHz,self.FreqMHz,self.BeamdataU[idx])
      v=np.interp(FreqMHz,self.FreqMHz,self.BeamdataV[idx])
      return b,q,u,v
   def resample(self,FreqMHz) :
      import copy
      import numpy as np
      out=stokesCubeMap()
      out.header=self.header
      out.FreqMHz=FreqMHz
      out._Nsamples=self._Nsamples
      out.Beamdata = np.zeros([self._Nsamples,len(out.FreqMHz)])
      out.BeamdataQ = np.zeros([self._Nsamples,len(out.FreqMHz)])
      out.BeamdataU = np.zeros([self._Nsamples,len(out.FreqMHz)])
      out.BeamdataV = np.zeros([self._Nsamples,len(out.FreqMHz)])
      for ii in range(self._Nsamples) : 
         b,q,u,v=self.interp(ii,out.FreqMHz)
         out.Beamdata[ii]=b
         out.BeamdataQ[ii]=q
         out.BeamdataU[ii]=u
         out.BeamdataV[ii]=v
      return out
   def average(self,FreqMHz,Weight,Method=None,nameRIMO=None) :
      import numpy as np
      import time
      import copy
      out=stokesCubeMap()
      out.File=copy.deepcopy(self.File)
      out._angularCut=self._angularCut
      out.header=self.header
      out.Beamdata = np.zeros([self._Nsamples])
      out.BeamdataQ = np.zeros([self._Nsamples])
      out.BeamdataU = np.zeros([self._Nsamples])
      out.BeamdataV = np.zeros([self._Nsamples])
      out._Nsamples=0
      dw=(FreqMHz[1:]-FreqMHz[:-1])*0.5
      out._Norm=((Weight[1:]+Weight[:-1])*dw).sum()
      out._Method=Method
      tic=time.time()
      for ii in range(self._Nsamples) : 
         b,q,u,v=self.interp(ii,FreqMHz)
         xx=b*Weight ; out.Beamdata[ii] = ((xx[1:]+xx[:-1])*dw).sum()/out._Norm
         xx=q*Weight ; out.BeamdataQ[ii] = ((xx[1:]+xx[:-1])*dw).sum()/out._Norm
         xx=u*Weight ; out.BeamdataU[ii] = ((xx[1:]+xx[:-1])*dw).sum()/out._Norm
         xx=v*Weight ; out.BeamdataV[ii] = ((xx[1:]+xx[:-1])*dw).sum()/out._Norm
      out._elapsed_time=time.time()-tic
      out._Method=Method
      out._nameRIMO=nameRIMO
      return out
   def tofits(self,fitsname,Author='M.Maris',creator='',version='',doNotWrite=False,clobber=True) :
      "saves the file in fits"
      import numpy as np
      from collections import OrderedDict
      import pyfits
      import dict2fits
      import time
      import copy
      #o=pyfits.open(fitsname)
      out = OrderedDict()
      out['Beamdata']=np.array(self['Beamdata'],dtype='float32')
      out['BeamdataQ']=np.array(self['BeamdataQ'],dtype='float32')
      out['BeamdataU']=np.array(self['BeamdataU'],dtype='float32')
      out['BeamdataV']=np.array(self['BeamdataV'],dtype='float32')
      T=dict2fits.Table(out)
      #T.header=copy.deepcopy(self.header['x'])
      T.header.update('TUNIT1',self.header['x']['TUNIT1'])
      T.header.update('TUNIT2',self.header['x']['TUNIT2'])
      T.header.update('TUNIT3',self.header['x']['TUNIT3'])
      T.header.update('TUNIT4',self.header['x']['TUNIT4'])
      #
      #test wether the original file was a GRD or a CUT file
      print ' copying Header'
      if self.header['x'].has_key('HIERARCH Nx') :
         # a grd file 
         for k in ['objType','Nx','Ny','Xcentre','Ycentre','Xdelta','Ydelta'] :
            T.header.update('HIERARCH '+k,self.header['x']['HIERARCH '+k])
            print "   %s : '%s' in '%s'"%(k,self.header['x']['HIERARCH '+k],T.header['HIERARCH '+k])
      else :
         # a CUT file
         for k in ['objType','Ntheta','Nphi','Mintheta','Maxtheta'] :
            T.header.update('HIERARCH '+k,self.header['x']['HIERARCH '+k])
            print "   %s : '%s' in '%s'"%(k,self.header['x']['HIERARCH '+k],T.header['HIERARCH '+k])
      print
      T.header.update('HIERARCH SUM BEAMDATA',self['Beamdata'].sum(),'sum of Beamdata')
      T.header.update('HIERARCH SUM BEAMDATAQ',self['BeamdataQ'].sum(),'sum of BeamdataQ')
      T.header.update('HIERARCH SUM BEAMDATAU',self['BeamdataU'].sum(),'sum of BeamdataU')
      T.header.update('HIERARCH SUM BEAMDATAV',self['BeamdataV'].sum(),'sum of BeamdataV')
      T.header.update('HIERARCH RIMO',self._nameRIMO if self._nameRIMO != None else '','')
      T.header.update('HIERARCH ANGULARCUT',self._angularCut if self._angularCut != None else 'None','angular cut [deg]')
      T.header.update('HIERARCH DATE',time.asctime(),'')
      T.header.update('HIERARCH CREATOR',creator,'')
      T.header.update('HIERARCH CONTACT',Author,'')
      if self._Method != None :
         for k in self._Method :
            T.header.update('HIERARCH '+k[0],k[1],k[2])
      T.header.add_comment('')
      T.header.add_comment('Beam band averaged')
      T.header.add_comment('Follows the list of input files used')
      for i in range(len(self.File)) :
         l='%s'%(self.File[i])
         T.header.add_comment(l)
      T.header.add_comment('')
      if not doNotWrite :
         print " Writing to:",fitsname
         T.writeto(fitsname,clobber=clobber)
      return T

class stokesMap(MapGrid) :
   def __init__(self,StokesFileName,mode='readonly') :
      from grid2d import GridAxis
      self.__new_info__()
      self.__info__['fits_readout_time_sec']=-1.
      if StokesFileName==None : return
      # connects to the fits file
      self.fits_connect(StokesFileName,mode=mode)
      # gets the header of the first hdu
      self.fits_load_hdu(1,justHeader=True)
      # formats the mapgrid according to the header content
      self.__info__['geometry']=self.getGeometry()
      if self.isGRD() :
         MapGrid.__init__(self,GridAxis('y','',self.__info__['geometry']['y']),GridAxis('x','',self.__info__['geometry']['x']))
      else :
         MapGrid.__init__(self,GridAxis('long','deg',self.__info__['geometry']['long']),GridAxis('colat','deg',self.__info__['geometry']['colat']))
      # gets the first hdu
      self.fits_load_hdu(1)
   def __new_info__(self) :
      from collections import OrderedDict
      self.__info__=OrderedDict()
   def fits_info(self) :
      "return infos on the fits file"
      if self.__info__['fitsfile']==None :
         return "File not connected"
      return self.__info__['fitsfile'].info()
   def fits_primary_header(self) : 
      "returns the current fits primary header content"
      return self.__info__['header']['p']
   def fits_extended_header(self) : 
      "returns the current fits extended header content"
      return self.__info__['header']['x']
   def fits_connect(self,FileName,mode='readonly') : 
      "connect (open and keep control) to a fits file"
      import pyfits
      self.__info__['StokesFileName']=FileName
      self.__info__['fitsfile']=pyfits.open(FileName,mode)
      self.__info__['ihdu']=-1
      self.__info__['header']={'p':None,'x':None}
      if self.__info__['fitsfile']==None :
         print "Error: file %s not found"%StokesFileName
         MapGrid.__init__(self)
         return 
      self.__info__['header']['p']=self.__info__['fitsfile'][0].header
   def fits_unconnect(self) : 
      "unconnect (close and left control) the fits file"
      if self.__info__['fitsfile']==None :
         print "File not connected"
      self.__info__['fitsfile'].close()
      self.__info__['fitsfile']=None
   def fits_load_hdu(self,ihdu,justHeader=False) :
      """fits_load_hdu(ihdu) 
      Load a fits hdu ihdu.
      If succesfull returns (ihdu, fits_readout_time_sec)"""
      from grid2d import GridAxis
      import time
      self.__info__['fits_readout_time_sec']=-1.
      if self.__info__['fitsfile']==None :
         print "File not connected"
      try :
         x=self.__info__['fitsfile'][ihdu].header
      except :
         print "hdu : ",ihdu," does not exists"
         return
      self.__info__['header']['x']=x
      self.__info__['ihdu']=ihdu
      if justHeader : return
      tic=time.time()
      b=self.__info__['fitsfile'][ihdu].data.field('Beamdata') ; b.shape=(self.R['n'],self.C['n']) ; self.newmap('Beamdata',value=b)
      q=self.__info__['fitsfile'][ihdu].data.field('BeamdataQ') ; q.shape=(self.R['n'],self.C['n']) ; self.newmap('BeamdataQ',value=q)
      u=self.__info__['fitsfile'][ihdu].data.field('BeamdataU') ; u.shape=(self.R['n'],self.C['n']) ; self.newmap('BeamdataU',value=u)
      v=self.__info__['fitsfile'][ihdu].data.field('BeamdataV') ; v.shape=(self.R['n'],self.C['n']) ; self.newmap('BeamdataV',value=v)
      self.__info__['fits_readout_time_sec']=time.time()-tic
      return self.__info__['ihdu'],self.__info__['fits_readout_time_sec']
   def copy(self,skipFields=None,skipFits=True) :
      "makes a copy, without fits informations"
      import copy
      out=stokesMap(None)
      for k in self.__dict__.keys() :
         if k !='__info__' :
            out.__dict__[k]=copy.deepcopy(self.__dict__[k])
         else :
            out.__new_info__()
            for k1 in self.__dict__[k].keys() :
               if (k1+'      ')[0:4].lower()!='fits' :
                  out.__dict__[k][k1]=copy.deepcopy(self.__dict__[k][k1])
      return out
   def fitsType(self) :
      if self.__info__['header']['x'].has_key('HIERARCH Nx') :
         return 'grd'
      else :
         return 'cut'
   def isGRD(self) :
      "Returns True if the file is a GRD"
      return self.__info__['header']['x'].has_key('HIERARCH Nx')
   def getGeometry(self) :
      import numpy as np
      geom={}
      if self.fitsType()=='grd' :
         for k in ['objType','Nx','Ny','Xcentre','Ycentre','Xdelta','Ydelta'] :
            geom[k]=self.__info__['header']['x']['HIERARCH '+k]
         geom['x']=(np.arange(geom['Nx'])-geom['Xcentre'])*geom['Xdelta']
         geom['y']=(np.arange(geom['Ny'])-geom['Ycentre'])*geom['Ydelta']
      else :
         for k in ['objType','Ntheta','Nphi','Mintheta','Maxtheta'] :
            geom[k]=self.__info__['header']['x']['HIERARCH '+k]
         geom['colat']=np.rad2deg(np.arange(geom['Ntheta'])*(geom['Maxtheta']-geom['Mintheta'])/float(geom['Ntheta']-1)+geom['Mintheta'])
         geom['long']=np.arange(geom['Nphi'])/float(geom['Nphi']-1)*360.
      return geom
   def coadd(self,that) :
      for k in ['Beamdata','BeamdataQ','BeamdataU','BeamdataV'] :
         self[k]+=that[k]
   def scale(self,that) :
      for k in ['Beamdata','BeamdataQ','BeamdataU','BeamdataV'] :
         self[k]*=that
   def __sub__(self,that) :
      new=self.copy()
      for k in ['Beamdata','BeamdataQ','BeamdataU','BeamdataV'] :
         try :
            new[k]-=that[k]
         except :
            new[k]-=that
      return new
   def __add__(self,that) :
      new=self.copy()
      for k in ['Beamdata','BeamdataQ','BeamdataU','BeamdataV'] :
         try :
            new[k]+=that[k]
         except :
            new[k]+=that
      return new
   def __pow__(self,that) :
      new=self.copy()
      for k in ['Beamdata','BeamdataQ','BeamdataU','BeamdataV'] :
         new[k]=new[k]**that
      return new
            
   
      
      