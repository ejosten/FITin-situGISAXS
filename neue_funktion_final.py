#-*- coding: utf-8 -*-

  # FITin-situGISAXS
 
#     In-situ GISAXS analysis and simulation/fit software. 
#     
#     The software enables rapid analysis of a large number of GISAXS data sets. In our case, this yields quantitative information 
#     on the time-resolved evolution of lattice constants, lattice constant distributions and tilt angles of the nanoparticle
#     superlattices with respect to the substrate. Moreover, by quantification of the coherent and incoherent intensities, it is
#     possible to quantify the relative ratio of nanoparticles in ordered assemblies/superlattices. 
#     
#     
# The software run under the plot.py software https://sourceforge.net/projects/plotpy/.
#     
#     
# Please read the accompanying LICENSE file before downloading the software. By downloading the software, you are agreeing 
#     to cite the corresponding paper: E. Josten et al.  Superlattice growth and rearrangement during evaporation-induced 
#     nanoparticle self-assembly.
#    
#     
#     
#     Copyright (C) {2017}  {A.Glavic, E.Josten}
#     
#     Contact: e.josten@hzdr.de
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or any 
#     later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.



'''
  Definition of mesocrystal GISAXS peak shape including instrumental
  resolution function and structural tilt and size distribution.
  To work the data has to be defined on a regular grid, otherwise the 
  FFT convolution algorithm won't work.
  
  The calculations can be performed with multiprocessing, splitting each
  peak calculation to the different processes.
'''

from numpy import sqrt, where, arctan2, abs, exp,pi, zeros_like, unique, minimum, sin, cos
from scipy.signal import fftconvolve
from plot_script.fit_data import FitFunction3D, register_class
from multiprocessing import Pool
import atexit

LAMBDA=4.51 # A

class PeakShapeResTilt(FitFunction3D):
  '''
    Fit a Gaussian or Lorentzian peak function with radial tilting distribution
    and a beam shape function convoluted with fft.
  '''

  # define class variables.
  name="Mesopeak + Resolution"
  # default parameters
  parameters=[0.058, 0.021, 0.0001, 0.0001,
              0., 0.001,
              0.,0.,0.,
              0.00531166, 0.000941805, 9233.57,
              0.3, 0.1,
              0.,
              ]
  # names of parameters
  parameter_names=['a*', 'c*', 'γ_y', 'γ_z', # peak parameters
                   'σ_tilt', 'σ_r', # tilt parameter in °
                   'Qy_off','Qz_off','Qphi_off', # instrumental offsets
                   'width' , 'height', 'decay', # instrumental parameters
                   'α_i', 'α_c',
                   'BG', # background
                   ]
  fit_function_text='Meso: a^*=[a*] c^*=[c*] γ_y=[γ_y|2] γ_z=[γ_z|2] σ_{tilt}=[σ_tilt|2]'
  parameter_description={'a*': 'Horizontal reciprocal lattice parameter',
                         'c*': 'Vertical reciprocal lattice parameter',
                         'Qy_off': 'Measurement Qy offset',
                         'Qz_off': 'Measurement Qz offset',
                         'Qphi_off': 'Measurement tilt offset',
                         'γ_y': 'Lorenz width in Qy-direction',
                         'γ_z': 'Lorenz width in Qz-direction',
                         'σ_tilt': 'Width of tilting angle',
                         'σ_r': 'Width of radial broughdening (distribution of lattice size)',
                         'σ_y': 'Width of Qy broughdening (distribution of a parameter)',
                         'σ_z': 'Width of Qz broughdening (distribution of c parameter)',
                         'width': 'Beam width horizontal (window function)',
                         'height': 'Beam width vertical (widnow function)',
                         'decay': 'Speed of exponential intensity decay outside of Beam window',
                         'BG': 'Background',
                         }
  
  use_rotation=True # default is an rotational average
  regions=None
  max_iter=50 # restrain the maximum number of iterations
  pool=None

  def __init__(self, initial_parameters=[], use_mp=False):
    FitFunction3D.__init__(self, initial_parameters)
    self.refine_parameters=[0,1,2,3,4,5]
    self.parameter_names=list(self.parameter_names)
    self.constrains={}
    if use_mp:
      # create a pool of worker processes
      self.activate_multiprocessing()
  
  def activate_multiprocessing(self):
    '''
      Create a worker pool.
    '''
    if self.pool is not None:
      self.deactivate_multiprocessing()
    self.pool=Pool()

  def deactivate_multiprocessing(self):
    '''
      Terminate the worker pool. Must be called before exiting
      the program, otherwise the processes stay alive and CTRL+c
      is needed to exit the program.
    '''
    if self.pool is not None:
      self.pool.close()
      self.pool.join()
      self.pool=None

  def add_region(self, xf, xt, yf, yt, Ii, Qyi, Qzi):
    '''
      Add an other peak with a defined fit region, intensity and position.
    '''
    if self.regions is None:
      self.regions=[]
    idx=len(self.regions)
    self.regions.append([xf, xt, yf, yt])
    self.parameters.append(Ii)
    self.parameters.append(Qyi)
    self.parameters.append(Qzi)
    self.parameter_names.append("I_%i"%idx)
    self.parameter_names.append("Qy_%i"%idx)
    self.parameter_names.append("Qz_%i"%idx)
    
  def fit_function(self, p, x, y):
    '''
      Combined intensity function for several peak regions including
      resolution and tilting.
      Each peak is a convolution of the resolution function a
      Gauusian/Lorentzian peak shape and a gaussian tilt and radial distribution.
    '''
    # get parameters by name
    astar=p[0]
    cstar=p[1]
    gamma_y=p[2]
    gamma_z=p[3]
    sigma_tilt=p[4]/180.*pi
    sigma_r=p[5]
    Qy_off=p[6]
    Qz_off=p[7]
    Qphi_off=p[8]/180.*pi
    width=p[9]
    height=p[10]
    decay=p[11]
    alpha_i=p[12]*pi/180.
    alpha_c=p[13]*pi/180.
    BG=p[14]
    if self.use_rotation:
      self.parameter_names[4:6]=['σ_tilt', 'σ_r']
    else:
      self.parameter_names[4:6]=['σ_y', 'σ_z']
    # if no region is defined the selected fit area is taken as region
    if self.regions is None:
      raise IndexError, "no peak defined"
    else:
      # each region has corresponding intensity, Qy- and Qz-position parameters
      regions=[[reg[0], reg[1], reg[2], reg[3],p[15+3*i], p[16+3*i], p[17+3*i]] \
                          for i,reg in enumerate(self.regions)]
    I=zeros_like(x)
    background=zeros_like(x)
    I_list=[]
    param_list=[]
    res_list=[]
    idx_list=[]
    pix_size=None
    # peaks for all regions
    for xf, xt, yf, yt, Ii, HK, L in regions:
      # get indices of points in the region
      idx=where((x>=xf)&(x<=xt)&(y>=yf)&(y<=yt))
      # calc position with phi offset
      Qyi=HK*astar
      Qzi=L*cstar
      
      # reshape x and y to two dimensional array
      ux=unique(x[idx])
      uy=unique(y[idx])
      xitems=len(ux)
      yitems=len(uy)
      if pix_size is None:
        pix_size=(ux[1]-ux[0])*(uy[1]-uy[0])
      # calculate pixel size for normaliztaion
      Qy=x[idx].reshape(yitems, xitems)
      Qz=y[idx].reshape(yitems, xitems)
      # correct instrumental missalignment
      # 0-offset
      Qy=Qy-Qy_off
      Qz=Qz-Qz_off
      # sample tilting
      Qytmp=Qy
      Qy= cos(Qphi_off)*Qytmp+sin(Qphi_off)*Qz
      Qz=-sin(Qphi_off)*Qytmp+cos(Qphi_off)*Qz
      
      in_list=[Qy, Qz, Qyi, Qzi, gamma_y, gamma_z,
                sigma_tilt, sigma_r, 
                width, height, decay, 
                alpha_i, alpha_c,
                Ii, self.use_rotation]
      idx_list.append(idx)
      if self.pool is None:
        # calculate one peak
        I_list.append(calc_I(in_list))
      else:
        # send calculation of one peak to the multiprocessing pool
        res_list.append(self.pool.apply_async(calc_I, args=(in_list,)))
    if self.pool is not None:
      # fetch the results of the multiprocessing calculations
      I_list=map(lambda item: item.get(), res_list)
    for idx, Ip in zip(idx_list, I_list):
      # add peak intensity to output
      I[idx]+=Ip.flatten()
      background[idx]=BG
    return I+background # add background
  
  def to_list(self):
    '''
      Get all important parameters as a list.
    '''
    return [self.parameters, self.parameter_names, 
            self.regions, [self.x_from, self.x_to, self.y_from, self.y_to],
            self.refine_parameters]
  
  def from_list(self, param_list):
    '''
      Set all important parameters from a list.
    '''
    self.parameters=list(param_list[0])
    self.parameter_names=list(param_list[1])
    self.regions=list(param_list[2])
    self.x_from,self.x_to,self.y_from,self.y_to=param_list[3]
    self.refine_parameters=list(param_list[4])
  
  def add_mirrors(self):
    '''
      For each region add the mirror peak on the other side of the specular line.
    '''
    regions=self.regions
    regparams=[[self.parameters[15+3*i],self.parameters[16+3*i],self.parameters[17+3*i]]
                 for i in range(len(regions))]
    pstart=15+3*len(regions)
    for i, regi, parami in zip(range(len(regions)), regions, regparams):
      self.add_region(-regi[1], -regi[0], regi[2], regi[3], parami[0], -parami[1],parami[2])
      # constrain all peak parameters to the original peak
      self.constrains[pstart+3*i]={'bounds': [None, None], 'tied': '[I_%i]'%i}
      self.constrains[pstart+1+3*i]={'bounds': [None, None], 'tied': '-[Qy_%i]'%i}
      self.constrains[pstart+2+3*i]={'bounds': [None, None], 'tied': '[Qz_%i]'%i}
      for j in range(3):
        # fit parameter if original was selected
        if (15+j+3*i) in self.refine_parameters:
          self.refine_parameters.append(pstart+j+3*i)
    self.x_from=-self.x_to

def Qz_obs(Qz, alpha_i, alpha_c):
  k0=2.*pi/LAMBDA
  return k0*(sin(alpha_i)+
          sqrt(sin(alpha_c)**2+
              (Qz*LAMBDA/(2.*pi)-sqrt(sin(alpha_i)**2-sin(alpha_c)**2))**2
              ))


def gaussian_tilt(Qy, Qz, Qy0, Qz0, sigma_tilt, sigma_r):
  '''
    Return a peak with gaussian tilt and radius components.
  '''
  Qr=sqrt(Qy**2+Qz**2)
  Qphi=arctan2(Qz,Qy)
  Qr0=sqrt(Qy0**2+Qz0**2)
  Qphi0=arctan2(Qz0,Qy0)
  #tilt_norm=1./(2.*pi*sigma_r*sigma_tilt)
  tilt_shape=exp(-0.5*(((Qr-Qr0)/(sigma_r*Qr0))**2+((Qphi-Qphi0)/sigma_tilt)**2))
  return tilt_shape#*tilt_norm

def beam(x, y, width, height, fall=8000., pix_size=1e-10):
  '''
    Calculate a beam shape function as a window horizontal and
    vertical widnow with exponentially decaying edges.
  '''
  px=exp(-fall*(abs(x)-width/2.))
  py=exp(-fall*(abs(y)-height/2.))
  if height<=pix_size:
    pos=(y<=pix_size)&(y>=-pix_size)
    py[pos]=minimum(1.-abs(y[pos])/pix_size,0.)
  if width<=pix_size:
    pos=(x<=pix_size)&(x>=-pix_size)
    px[pos]=minimum(1.-abs(y[pos])/pix_size, 0.)
  Peak=minimum(px, 1.)*minimum(py,1.)
  
  return Peak#/Peak.sum()

def gaussian(x,y,sigma_x, sigma_y):
  '''
    Not normalized Gaussian peak shape function.
  '''
  #norm=1./(2.*pi*sigma_x*gamma_y)
  G=exp(-0.5*((x/sigma_x)**2+(y/sigma_y)**2))
  return G#*norm

def lorentzian(x,y,gamma_x, gamma_y, pix_size):
  '''
    Not normalized lorentzian peak shape function.
  '''
  #norm=1./((pi*gamma_x*gamma_y))
  L=1./(1.+((x/gamma_x)**2+(y/gamma_y)**2))
  # make sure the maximum of the function is not reduced by too small gamma values
  if gamma_x<=pix_size and gamma_y<=pix_size:
    L[(x<=pix_size/2.)&(x>=-pix_size/2.)&(y<=pix_size/2.)&(y>=-pix_size/2.)]=1.
  elif gamma_x<=pix_size:
    pos=(x<=pix_size/2.)&(x>=-pix_size/2.)
    L[pos]=1./(1.+(y[pos]/gamma_y)**2)
  elif gamma_y<=pix_size:
    pos=(y<=pix_size/2.)&(y>=-pix_size/2.)
    L[pos]=1./(1.+(x[pos]/gamma_x)**2)
  return L#*norm
    
def calc_I(in_list):
  '''
    Calculate intensity on one region.
  '''
  Qy, Qz, Qyi, Qzi=in_list[:4]
  gamma_y, gamma_z=in_list[4:6]
  sigma_tilt, sigma_r=in_list[6:8]
  width, height, decay=in_list[8:11]
  alpha_i, alpha_c=in_list[11:13]
  Ii=in_list[13]
  use_rotation=in_list[14]
  # calculate approx. pixe size
  pix_size=(Qz.max()-Qz.min())/len(Qz)
  if sigma_tilt==0:
    # no tilting of mesos
    # calculate peak shape function
    Ip=lorentzian(Qy-Qyi,Qz-Qzi,gamma_y, gamma_z, pix_size)
  else:
    # tilting of mesos
    # calculate peak shape function
    peak_shape=lorentzian(Qy-Qy.mean(), Qz-Qz.mean(), gamma_y, gamma_z, pix_size)
    if use_rotation:
      # calculate tilt and radial distribution function
      gauss_shape=gaussian_tilt(Qy,Qz,Qyi,Qzi, sigma_tilt, sigma_r)
    else:
      # calculate gaussian distribution in Qy and Qz direction
      gauss_shape=gaussian(Qy-Qyi,Qz-Qzi, sigma_tilt*180./pi*Qyi, sigma_r*Qzi)
    # convolute tilt and radial distribution with peak shape
    Ip=fftconvolve(peak_shape, gauss_shape, 'same')
  if height!=0 and width!=0:
    # calculate beam shape function
    # the peak position is moved due to refraction effects
    Qzi_obs=Qz_obs(Qzi, alpha_i, alpha_c)
    beam_shape=beam(Qy-Qy.mean(), Qz-Qz.mean()-Qzi_obs+Qzi, width, height, decay, pix_size)
    # convolute peak shape with beam shape
    Ip=fftconvolve(beam_shape, Ip, 'same')
  return Ii*Ip/Ip.sum() # normalize to integrated intensity
  
class BeamResolution(FitFunction3D):
  '''
    A pure resolution function to refine the peak shape with
    e.g. silicon total reflection measurement.
  '''
  name="Resolution"
  parameters=[10000., 0., 0.177, 0.00531166, 0.000941805, 9233.57]
  parameter_names=['I0', 'x0','y0','width','height','decay']
  fit_function_text='Resolution'
  
  def fit_function(self, p, x, y):
    I0,x0,y0,w,h,d=p
    return I0*beam(x-x0,y-y0,w,h,d)
    
# Add class to the available fit functions in Plot.py
register_class(PeakShapeResTilt)
register_class(BeamResolution)

#[0.031860800000000002, 0.055555899999999998, 0.10974299999999999, 0.0018228000000000001, 0.0019904900000000001, 5.2483300000000002, 0.00013660899999999999, 0.0028664099999999998, 0.0, 0.020954199999999999, 0.058215099999999999, 0.089450399999999999]
#[[0.035, 0.08, 0.085, 0.13], [0.035, 0.08, 0.07, 0.115]]