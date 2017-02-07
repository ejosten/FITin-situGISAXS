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
  Load a set of measurements and fit a peak function for each dataset.
'''

import os,sys
import shutil
import time
from glob import glob
from copy import deepcopy
from plot_script.read_data.kws2 import read_edf_file
from form_factors import FitSphereQres
from neue_funktion_final import PeakShapeResTilt


stdout=None

print """
    {FITin-situGISAXS}  Copyright (C) {2017}  {A.Glavic, E.Josten}
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details. 
    
    By downloading the software, you are agreeing to cite the corresponding paper: 
    E. Josten et al.  Superlattice growth and rearrangement during evaporation-induced 
    nanoparticle self-assembly.
  
    Usage:
    -Open one dataset and make your correction:
        d=filter_dataset(dataset(), [0.,10.,0.,10.])
        session.active_file_data.append(d)
        remove_bg(d)
        remove_sphere(d)
        # <ctrl>+N
        
    -Fit the peak function with the desired xy-region and define the free parameters
    
    -Run the fits for all files in the active directory:
        z.B. plots_for_all(preset='s113')
          plots_for_all(fit_function=None,
                        name='Fits',
                        store_data=False,
                        correction=None, 
                        remove_background=False,
                        preset=None)"""
  
SPHERE_PARAMS=[
                0.000577067, # I0
                0.,  # BG
                50.0518, # R
                0.0638503, # δR
                0.0035, # ΔΘ vorher
                1e-10, # Δλ
                4.51, # λ
                ]
                
               # 0.00135482, # ΔΘ vorher
CUBES_PARAMS=[
                0.000577067, # I0
                0.,  # BG
                60.3707, # R
                0.0556317, # δR
                0.0035, # ΔΘ
                1e-10, # Δλ
                4.51, # λ
                ]
        #  x-from x-to  y-from y-to
#BG_REGION=[0.15,  0.16, 0.21,  0.229]   #vorher
BG_REGION=[0.151,  0.161, 0.22,  0.229] #P321, P227



def iteration_update(step_add=None,info=''):
  '''
    Define a function which updates the user on the current fit iteration.
  '''
  stdout.write(info+'\n')
  stdout.flush()
  return 0

def fit_data_list(datasets, # list of file names
                  initial_fit=None, # Fit function to be used, if None uses fit-0 from active dataset 
                  correction=None, # form factor correction function to be used
                  remove_background=False, # remove the background
                  plot_prefix='Result_', # name of session.file_data list to append datasets to
                  incregion=None, # region to extract incoherent scattering
                  ii_region=None, # region to calculate integrated intensity
                  param_file=None, # file to store the text data to
                  ): 
  '''
    Starting from a given fit (if None uses the first 
    fit of the active dataset) all datasets in a list of
    names are read and fit. The result parameters for each
    file are returned as a list.
  '''
  global stdout, do_plot
  if stdout is None:
    try:
      # open dialog for output
      stdout=plot_gui.status_dialog
      stdout.show_all()
      do_plot=True
      session.file_data['temp']=[None]
      session.active_file_data=session.file_data['temp']
      session.active_file_name='temp'
      plot_gui.measurement=session.file_data['temp']
      plot_gui.index_mess=0
    except:
      stdout=sys.stdout
      do_plot=False
  from config import gnuplot_preferences
  gnuplot_preferences.settings_3dmap+='set format cb "10^{%L}"\nset cblabel offset 1.5\n'
  gnuplot_preferences.plotting_parameters_3d='w points palette ps .1 pt 5'
  if do_plot:
    import gtk
    keep_running=gtk.CheckButton(label='Keep Running')
    keep_running.set_active(True)
    stdout.vbox.pack_end(keep_running, False)
    keep_running.show()
    stdout.present()
  if initial_fit is None:
    # get first fit from active dataset
    fit=dataset().fit_object.functions[0][0]
  else:
    fit=initial_fit
  # create multiprocessing framework
  fit.activate_multiprocessing()
  
  if param_file is not None:
    param_file=open(param_file+'.txt', 'w')
    param_file.write( "\t".join(fit.parameter_names)+"\t"
                      "\t".join(["δ"+name for name in fit.parameter_names])+
                      "\t#\tTime\n"
                     )
    param_file.flush()
  
  parameters=[]
  ff=None # formfactor intensity array
  incidxs=None
  ii_idxs=None
  first_time=None
  for i, ds in enumerate(datasets):
    # iterate through all datasets
    stdout.write("  Reading %s...\n" % ds)
    stdout.flush()
    # read dataset and extract number and time
    data=read_edf_file(ds)[0]
    info=map(lambda line: map(str.strip, line.split(':',1)), data.info.splitlines())
    info=dict(filter(lambda item: len(item)==2, info))
    file_number=int(ds.rsplit('_',1)[1].split("ccd")[0])
    gmtime=time.strptime(info['HMStartTime'].split('.')[0], '%Y-%m-%dT%H:%M:%S')
    float_time=time.mktime(gmtime)
    if remove_background:
      stdout.write("    ...background...\n")
      stdout.flush()
      remove_bg(data)      
    if correction is not None:
      stdout.write("    ...correcting...\n")
      stdout.flush()
      # filter the dataset for the fitted region
      dataf=filter_dataset_fitregion(data, fit)
      # clear all references to the data and delete the temporary file
      data.store_data()
      os.remove(data.tmp_export_file)
      del(data)
      data=dataf
      if ff is None:
        # calculate the form factor correction intensities on first call
        ff=correction(data)
      data.z/=ff
    if incregion is not None:
      stdout.write("    ...incoherent background...\n")
      stdout.flush()
      if incidxs is None:
        incidxs=np.where((data.x>=incregion[0])&(data.x<=incregion[1])&
                         (data.y>=incregion[2])&(data.y<=incregion[3]))
      incdata=float(data.z[incidxs].sum())/len(data.z[incidxs])
      bg_param_index=fit.parameter_names.index('BG')
      fit.parameters[bg_param_index]=incdata
      if bg_param_index in fit.refine_parameters:
        fit.refine_parameters.remove(bg_param_index)
        
      stdout.write("    ...=%g...\n"%incdata)
      stdout.flush()  
    
    
    if ii_region is not None:
      stdout.write("    ...integrated intensity")
      stdout.flush()
      if ii_idxs is None:
        ii_idxs=np.where((data.x>=ii_region[0])&(data.x<=ii_region[1])&
                         (data.y>=ii_region[2])&(data.y<=ii_region[3]))
      ii_data=data.z[ii_idxs].sum()/len(data.z[incidxs])
      if incregion is not None:
        ii_data-=incdata
      stdout.write(" %f\n"%(float(ii_data)))
      stdout.flush()
    stdout.write("    ...refineing...\n")
    stdout.flush()
    fitmsg, covar=fit.refine(data.x, data.y, data.z, None,#data.z.error, ##None=ignore errors, data.z.error=berücksichtigung
                              progress_bar_update=iteration_update) #Fit Routine
    errors=[]
    if fitmsg!='':
      stdout.write("Error in Fit: "+fitmsg)
      stdout.flush()
    errors=[np.sqrt(covar[i][i]) for i in range(len(fit.parameters))]
    if first_time is None:
      first_time=float_time
    float_time-=first_time
    if ii_region is None:
      fit_parameters=list(fit.parameters)+errors+[file_number, float_time]
    else:
      fit_parameters=list(fit.parameters)+[float(ii_data)]+\
                          errors+[float(ii_data.error)]+[file_number, float_time]
    parameters.append(fit_parameters)
    if param_file is not None:
      param_file.write( "\t".join(["%g"%param for param in  fit_parameters])+"\n")
      param_file.flush()

    stdout.write("    ...result:\n  \t%s\n  \t%s\n+/-\t%s\n\n" % (
                                            repr(fit.parameter_names),
                                            repr(fit.parameters),
                                            repr(errors)))
    stdout.flush()
    # create fit data for plotting
    data.is_matrix_data=False
    data.plot_options.xrange=[fit.x_from,fit.x_to]
    data.plot_options.yrange=[fit.y_from,fit.y_to]
    data.plot_options.zrange=[1.,1000.]
    data.plot_together_zindex=-1
    fit_data=mds.MeasurementData([],x=0,y=1,zdata=2)
    fit_data.data.append(data.x)
    fit_data.data.append(data.y)
    fit_data.data.append(mds.PhysicalProperty(data.z.dimension, data.z.unit,
                          fit(data.x,data.y)))
    data.plot_together.append(fit_data)
    data.sample_name=''
    data.short_info='#%i'%file_number
    session.picture_width='2400'
    session.pictire_height='1200'
    if do_plot:
      # plot an image and copy the result
      plot_gui.measurement[0]=data
      replot()
      shutil.copy(session.TEMP_DIR+'/plot_temp.png', plot_prefix+'%04i.png'%file_number)
    else:
      # directly plot as file
      mdp.gnuplot_plot_script(session,
                           [data],
                            'temp_plot',
                           '.png',
                           data.short_info,
                           [data.short_info],
                           False,
                           plot_prefix+'%04i.png'%file_number)
      
    # clear all references to the data and delete the temporary file
    data.store_data()
    os.remove(data.tmp_export_file)
    del(data)
    if do_plot and not keep_running.get_active():
      # stop iteration due to user input
      break
  # terminate other processes
  fit.deactivate_multiprocessing()
  if do_plot:
    stdout.vbox.remove(keep_running)
  return np.array(parameters).transpose()

def plots_from_paramlist(parameters, parameter_names):
  # create plots of parameters 
  idx=mds.PhysicalProperty('File No.', '', parameters[-2])
  tm=mds.PhysicalProperty('Time', 's', parameters[-1])
  for i,name in enumerate(parameter_names):
    ds=mds.MeasurementData(x=0,y=2)
    ds.data.append(tm)
    ds.data.append(idx)
    if np.all(parameters[len(parameter_names)+i]==0):
      ycol=mds.PhysicalProperty(name, '', parameters[i])
    else:
      ycol=mds.PhysicalProperty(name, '', parameters[i],
                                        parameters[len(parameter_names)+i])
    ds.data.append(ycol)
    session.active_file_data.append(ds)

def plots_for_all(fit_function=None,
                  name='Fits',
                  store_data=True,
                  correction=None, 
                  remove_background=False,
                  param_file=None,
                  preset=None):
  '''
    Perform analysis for all files in the active directory
    starting from the last file. Options can be given
    when called or using the preset name.
  '''
  global params, initial_fit
  files=glob('*.edf')
  files.sort()
  files.reverse()
  if preset is None:
    initial_fit=None
    incregion=None
    ii_region=None
  else:
    preset=fit_presets[preset]
    initial_fit=PeakShapeResTilt(list(preset['initial_fit']))
    initial_fit.use_rotation=preset['use_rotation']
    # add all peak regions and paramters
    min_x=1000.
    max_x=-1000.
    min_y=1000.
    max_y=-1000.
    initial_fit.refine_parameters=range(8)+[13,14] #all free peak parameters + background and α_c
    initial_fit.regions=[]
    for xf,xt,yf,yt,Ii,Qyi,Qzi in preset['peaks']:
      # set the fit region as maximum from all peak ranges
      min_x=min(min_x, xf)
      max_x=max(max_x, xt)
      min_y=min(min_y, yf)
      max_y=max(max_y, yt)
      # add the peak to the fit function
      i=len(initial_fit.regions)
      initial_fit.add_region(xf, xt, yf, yt, Ii, Qyi, Qzi)
      #if i in [0]:
      initial_fit.refine_parameters+=[15+3*i] # free intensity
      #else:
      #  initial_fit.refine_parameters+=[15+3*i, 17+3*i] # free intensity and Qz
    # constrain intensities to be >= 0
    initial_fit.constrains[14]={'bounds': [0., None], 'tied': ''} # α_c>=0
    for i in range(len(initial_fit.regions)):
      if (15+3*i) in initial_fit.constrains:
        initial_fit.constrains[15+3*i]['bounds'][0]=0.
      else:
        initial_fit.constrains[15+3*i]={'bounds': [0., None], 'tied': ''}
    initial_fit.x_from=min_x
    initial_fit.x_to=max_x
    initial_fit.y_from=min_y
    initial_fit.y_to=max_y
    if preset['mirror_peaks']:
      initial_fit.add_mirrors()
    correction=preset['correction']
    remove_background=preset['remove_background']
    incregion=preset['incregion']
    ii_region=preset['ii_region']
    if param_file is None:
      param_file=preset['param_file']
    if name == 'Fits':
      name=preset['name']
  # fit all datasets
  params=fit_data_list(files, 
                       correction=correction, 
                       remove_background=remove_background,
                       initial_fit=initial_fit,
                       plot_prefix=param_file+'_',
                       ii_region=ii_region,
                       incregion=incregion, param_file=param_file)
  #if param_file is not None:
  #  params[-1]-=params[-1].min()
  #  np.savetxt(param_file+'.txt', params.transpose(), fmt='% 12g', delimiter='\t')
  # create plots for parameters
  session.file_data[name]=[]
  session.active_file_data=session.file_data[name]
  session.active_file_name=name
  plot_gui.measurement=session.file_data[name]
  plots_from_paramlist(params, initial_fit.parameter_names)
  try:
    plot_gui.rebuild_menus()
    replot()
  except:
    pass

def remove_sphere(ds):
  '''
    Devide the measured intensity by the particle formfactor.
  '''
  ff=ff_sphere(ds)
  ds.z/=ff

def ff_sphere(ds):
  # sphere form factor
  FF=FitSphereQres(SPHERE_PARAMS)
  return FF(np.sqrt(ds.x**2+ds.y**2))
  
def remove_cube(ds):
  '''
    Devide the measured intensity by the particle formfactor.
  '''
  ff=ff_cube(ds)
  ds.z/=ff
  
  #Erhöhung der Minimas zur Vermedung der Ringe 
  #Q=np.sqrt(ds.x**2+ds.y**2)
  #ds.z/=ff+5*Q**(-3)+30000
  
def ff_cube(ds):
  # cube form factor
  FF=FitSphereQres(CUBES_PARAMS)
  return FF(np.sqrt(ds.x**2+ds.y**2))
  
def remove_bg(ds):
  '''
    Subtract background.
  '''
  rg=BG_REGION
  bg=ds.z[(ds.x>=rg[0])&(ds.x<=rg[1])&\
               (ds.y>=rg[2])&(ds.y<=rg[3])]
  bg=bg.sum()/len(bg)
  ds.z-=bg   #P231,P227
  #ds.z-=8  #P225

def filter_dataset(ds, region):
  '''
    Get a new dataset for a region defined by [x0, x1, y0, y1]
  '''
  ids=np.where((ds.x>=region[0])&(ds.x<=region[1])&\
               (ds.y>=region[2])&(ds.y<=region[3]))
  return ds[ids]

def filter_dataset_fitregion(ds, fit):
  '''
    Get a new dataset for a region defined by PeakShapeResTilt peak regions.
  '''
  ids=ds.x!=ds.x # all False
  for region in fit.regions:
    ids=ids|((ds.x>=region[0])&(ds.x<=region[1])&(ds.y>=region[2])&(ds.y<=region[3]))
  return ds[ids]


# Paramters for different Peaks
fit_presets={
  
  
   'P231ai03': {                 # General peak parameters
          'initial_fit': [  #'a*', 'c*', 'γ_y', 'γ_z'
                           0.058643, 0.0208204, 0.000276603, 0.000432451, 
                           #'σ_tilt', 'σ_r',
                           4.92965, 0.0179373, 
                           #'Qy_off','Qz_off','Qphi_off' is fixed
                           -0.00219964, 0.00659945, -0.577801, 
                           #'width' , 'height', 'decay' is fixed
                           0.00497968, 0.000754677, 3195.79, 
                           #'α_i', 'α_c'
                           0.3, 0.0000244498,
                           # BG
                           27.,
                         ],
          'peaks': [  # parameters for each single peak
                      # xfrom, xto, yfrom, yto,  I,   H/K     L
                      [0.026, 0.09, 0.154, 0.185, 5e5, 1, 8], # 1 0 8
                      [0.026, 0.087, 0.073, 0.11, 5e5, 1, 4],# 1 0 4
                      [0.026, 0.087, 0.089, 0.127, 5e5, 1, 5],# 1 0 5 
                      [0.07, 0.116, 0.049, 0.0888199, 5e5, np.sqrt(3.), 3], # 1 1 3
                      [0.09, 0.14, 0.067, 0.107, 5e5, 2, 4], # 2 0 4
                      [0.07, 0.116, 0.104658385093, 0.146, 5e5, np.sqrt(3.), 6] # 1 1 6
              
                    ],
                    #  x-from, x-to, y-from, y-to
           'mirror_peaks': True, # add peaks from the left side
           'use_rotation': True, 
           'correction': ff_sphere,
           'remove_background': True,
           'param_file': 'P231_Spheres_short_ai03_timedependence_final',
           'name': 'P231_Spheres_short_ai03',
           'ii_region': [0.035, 0.08, 0.096, 0.13], # Integrated inteisity around 0 1 5
           'incregion': [0.02 , 0.03, 0.095,   0.11],
          }, 
          
    'P226ai03': {                 # General peak parameters
          'initial_fit': [ #'a*', 'c*', 'γ_y', 'γ_z'
                           0.056734, 0.0200253, 0.0000322504, 0.000158157, 
                           #'σ_tilt', 'σ_r',
                           2.73872, 0.0162109, 
                           #'Qy_off','Qz_off','Qphi_off' is fixed
                           -0.00189297, 0.00698159, -0.564942, 
                           #'width' , 'height', 'decay' is fixed
                           0.00497968, 0.000754677, 3195.79, 
                           #'α_i', 'α_c'
                           0.3, 0.0000244498,
                           # BG
                           27.,
                         ],
          'peaks': [  # parameters for each single peak
                      # xfrom, xto, yfrom, yto,  I,   H/K     L
                      [0.026, 0.09, 0.154, 0.185, 5e5, 1, 8], # 1 0 8
                      [0.026, 0.087, 0.073, 0.11, 5e5, 1, 4],# 1 0 4
                      [0.026, 0.087, 0.089, 0.127, 5e5, 1, 5],# 1 0 5 
                      [0.07, 0.116, 0.049, 0.0888199, 5e5, np.sqrt(3.), 3], # 1 1 3
                      [0.09, 0.14, 0.067, 0.107, 5e5, 2, 4], # 2 0 4
                      [0.07, 0.116, 0.104658385093, 0.146, 5e5, np.sqrt(3.), 6] # 1 1 6
              
                    ],
                    #  x-from, x-to, y-from, y-to
           'mirror_peaks': True, # add peaks from the left side
           'use_rotation': True, 
           'correction': ff_sphere,
           'remove_background': True,
           'param_file': 'P226_Spheres_middle_ai03_timedependence_final',
           'name': 'P226_Spheres_midlle_ai03',
           'ii_region': [0.035, 0.08, 0.095, 0.13], # Integrated inteisity around 0 1 5
           'incregion': [0.02 , 0.03, 0.095,   0.11],
          },  
     
  
  
  
  #'P231ai03': {                 # General peak parameters
          #'initial_fit': [ #'a*', 'c*', 'γ_y', 'γ_z'
                           #0.0586, 0.02084, 0.0017, 0.0019, 
                           ##'σ_tilt', 'σ_r',
                           #5.2, 0.0085, 
                           ##'Qy_off','Qz_off','Qphi_off' is fixed
                           #-0.0021, 0.0054, -0.614, 
                           ##'width' , 'height', 'decay' is fixed
                           #0.00497968, 0.000754677, 3195.79, 
                           ##'α_i', 'α_c'
                           #0.3, 0.142,
                           ## BG
                           #3,
                         #],
          #'peaks': [  # parameters for each single peak
                      ## xfrom, xto, yfrom, yto,  I,   H/K     L
                      #[0.026, 0.09, 0.154, 0.19, 5e5, 1, 8], # 1 0 8
                      #[0.026, 0.087, 0.073, 0.11, 5e5, 1, 4],# 1 0 4
                      #[0.026, 0.087, 0.089, 0.127, 5e5, 1, 5],# 1 0 5 
                      #[0.07, 0.116, 0.049, 0.0888199, 5e5, np.sqrt(3.), 3], # 1 1 3
                      #[0.09, 0.14, 0.067, 0.107, 5e5, 2, 4], # 2 0 4
                      #[0.07, 0.116, 0.104658385093, 0.15, 5e5, np.sqrt(3.), 6] # 1 1 6
              
                    #],
                    ##  x-from, x-to, y-from, y-to
           #'mirror_peaks': True, # add peaks from the left side
           #'use_rotation': True, 
           #'correction': ff_sphere,
           #'remove_background': True,
           #'param_file': 'P231_Spheres_short_ai03_timedependence',
           #'name': 'P231_Spheres_short_ai03',
           #'ii_region': [0.035, 0.08, 0.096, 0.13], # Integrated inteisity around 0 1 5
           #'incregion': [0.02 , 0.03, 0.095,   0.11],
          #},  
  #'P231ai03test1': {                 # General peak parameters
          #'initial_fit': [ #'a*', 'c*', 'γ_y', 'γ_z'
                           #0.0583594, 0.0206211, 0.0012697, 0.0016383, 
                           ##'σ_tilt', 'σ_r',
                           #4.43711, 0.007965, 
                           ##'Qy_off','Qz_off','Qphi_off' is fixed
                           #-0.00219, 0.0047469, -0.622864, 
                           ##'width' , 'height', 'decay' is fixed
                           #0.00497968, 0.000754677, 3195.79, 
                           ##'α_i', 'α_c'
                           #0.3, 0.205425,
                           ## BG
                           #10,
                         #],
          #'peaks': [  # parameters for each single peak
                      ## xfrom, xto, yfrom, yto,  I,   H/K     L
                      #[0.026, 0.09, 0.154, 0.19, 5e5, 1, 8], # 1 0 8
                      #[0.026, 0.087, 0.073, 0.11, 5e5, 1, 4],# 1 0 4
                      #[0.026, 0.087, 0.089, 0.127, 5e5, 1, 5],# 1 0 5 
                      #[0.07, 0.116, 0.049, 0.0888199, 5e5, np.sqrt(3.), 3], # 1 1 3
                      #[0.09, 0.14, 0.067, 0.107, 5e5, 2, 4], # 2 0 4
                      #[0.07, 0.116, 0.104658385093, 0.15, 5e5, np.sqrt(3.), 6] # 1 1 6
              
                    #],
                    ##  x-from, x-to, y-from, y-to
           #'mirror_peaks': True, # add peaks from the left side
           #'use_rotation': True, 
           #'correction': ff_sphere,
           #'remove_background': True,
           #'param_file': 'P231_Spheres_short_ai03_timedependence_test1',
           #'name': 'P231_Spheres_short_ai03_test1',
           #'ii_region': [0.035, 0.08, 0.096, 0.13], # Integrated inteisity around 0 1 5
           #'incregion': [0.02 , 0.03, 0.095,   0.11],
          #},  
          
          
     #'P231ai03test1eng': {                 # General peak parameters
          #'initial_fit': [ #'a*', 'c*', 'γ_y', 'γ_z'
                           #0.0583594, 0.0206211, 0.0012697, 0.0016383, 
                           ##'σ_tilt', 'σ_r',
                           #4.43711, 0.007965, 
                           ##'Qy_off','Qz_off','Qphi_off' is fixed
                           #-0.00219, 0.0047469, -0.622864, 
                           ##'width' , 'height', 'decay' is fixed
                           #0.00497968, 0.000754677, 3195.79, 
                           ##'α_i', 'α_c'
                           #0.3, 0.205425,
                           ## BG
                           #10,
                         #],
          #'peaks': [  # parameters for each single peak
                      ## xfrom, xto, yfrom, yto,  I,   H/K     L
                      #[0.026, 0.09, 0.154, 0.185, 5e5, 1, 8], # 1 0 8
                      #[0.026, 0.087, 0.073, 0.11, 5e5, 1, 4],# 1 0 4
                      #[0.026, 0.087, 0.089, 0.127, 5e5, 1, 5],# 1 0 5 
                      #[0.07, 0.116, 0.049, 0.0888199, 5e5, np.sqrt(3.), 3], # 1 1 3
                      #[0.09, 0.14, 0.067, 0.107, 5e5, 2, 4], # 2 0 4
                      #[0.07, 0.116, 0.104658385093, 0.146, 5e5, np.sqrt(3.), 6] # 1 1 6
              
                    #],
                    ##  x-from, x-to, y-from, y-to
           #'mirror_peaks': True, # add peaks from the left side
           #'use_rotation': True, 
           #'correction': ff_sphere,
           #'remove_background': True,
           #'param_file': 'P231_Spheres_short_ai03_timedependence_test1_eng',
           #'name': 'P231_Spheres_short_ai03_test1_eng',
           #'ii_region': [0.035, 0.08, 0.096, 0.13], # Integrated inteisity around 0 1 5
           #'incregion': [0.02 , 0.03, 0.095,   0.11],
          #},  
    
    #'P226ai03': {                 # General peak parameters
          #'initial_fit': [ #'a*', 'c*', 'γ_y', 'γ_z'
                           #0.0586, 0.02084, 0.0017, 0.0019, 
                           ##'σ_tilt', 'σ_r',
                           #5.2, 0.0085, 
                           ##'Qy_off','Qz_off','Qphi_off' is fixed
                           #-0.0021, 0.0054, -0.614, 
                           ##'width' , 'height', 'decay' is fixed
                           #0.00497968, 0.000754677, 3195.79, 
                           ##'α_i', 'α_c'
                           #0.3, 0.142,
                           ## BG
                           #3,
                         #],
          #'peaks': [  # parameters for each single peak
                      ## xfrom, xto, yfrom, yto,  I,   H/K     L
                      #[0.026, 0.09, 0.154, 0.19, 5e5, 1, 8], # 1 0 8
                      #[0.026, 0.087, 0.073, 0.11, 5e5, 1, 4],# 1 0 4
                      #[0.026, 0.087, 0.089, 0.127, 5e5, 1, 5],# 1 0 5 
                      #[0.07, 0.116, 0.049, 0.0888199, 5e5, np.sqrt(3.), 3], # 1 1 3
                      #[0.09, 0.14, 0.067, 0.107, 5e5, 2, 4], # 2 0 4
                      #[0.07, 0.116, 0.104658385093, 0.15, 5e5, np.sqrt(3.), 6] # 1 1 6
              
                    #],
                    ##  x-from, x-to, y-from, y-to
           #'mirror_peaks': True, # add peaks from the left side
           #'use_rotation': True, 
           #'correction': ff_sphere,
           #'remove_background': True,
           #'param_file': 'P226_Spheres_middle_ai03_timedependence',
           #'name': 'P226_Spheres_midlle_ai03',
           #'ii_region': [0.035, 0.08, 0.096, 0.13], # Integrated inteisity around 0 1 5
           #'incregion': [0.02 , 0.03, 0.095,   0.11],
          #},  
     
             
    #'P226ai03test1': {         # General peak parameters
          #'initial_fit': [ #'a*', 'c*', 'γ_y', 'γ_z'
                              #0.0562635, 0.019694, 0.000299194, 0.000897285, 
                           ##'σ_tilt', 'σ_r',
                           #1.32479, 0.00616268, 
                           ##'Qy_off','Qz_off','Qphi_off' is fixed
                           #-0.00187526, 0.00686513, -0.750099, 
                           ##'width' , 'height', 'decay' is fixed
                           #0.00497968, 0.000754677, 3195.79, 
                           ##'α_i', 'α_c'
                           #0.3, 0.0439488,
                           ## BG
                           #6.7, 
                         #],
          #'peaks': [  # parameters for each single peak
                      ## xfrom, xto, yfrom, yto,  I,   H/K     L
                      #[0.026, 0.09, 0.154, 0.19, 5e5, 1, 8], # 1 0 8
                      #[0.026, 0.087, 0.073, 0.11, 5e5, 1, 4],# 1 0 4
                      #[0.026, 0.087, 0.089, 0.127, 5e5, 1, 5],# 1 0 5 
                      #[0.07, 0.116, 0.049, 0.0888199, 5e5, np.sqrt(3.), 3], # 1 1 3
                      #[0.09, 0.14, 0.067, 0.107, 5e5, 2, 4], # 2 0 4
                      #[0.07, 0.116, 0.104658385093, 0.15, 5e5, np.sqrt(3.), 6] # 1 1 6
              
                    #],
                    ##  x-from, x-to, y-from, y-to
           #'mirror_peaks': True, # add peaks from the left side
           #'use_rotation': True, 
           #'correction': ff_sphere,
           #'remove_background': True,
           #'param_file': 'P226_Spheres_middle_ai03_timedependence_test1',
           #'name': 'P226_Spheres_midlle_ai03_test1',
           #'ii_region': [0.035, 0.08, 0.096, 0.13], # Integrated inteisity around 0 1 5
           #'incregion': [0.02 , 0.03, 0.095,   0.11],
          #},  
  
          
    #'P227ai03': {                 # General peak parameters
          #'initial_fit': [# a*         c*         γ_y        γ_z       
                          #4.790e-02,  2.790e-02,  1.563e-05,  1.099e-02,
                          ##σ_y        σ_z        Qy_off     Qz_off    
                          #1.118e-01,  0.05, -2.740e-03,  5.526e-03,
                          ##Qphi_off   width      height     decay     
                          #-4.281e-01,  0.00497968, 0.000754677, 3195.79,
                          ##α_i        α_c        BG         I_0       
                          #3.000e-01,  1.488e-03,  5.000e+00,
                         #],
          #'peaks': [  # parameters for each single peak
                      ## xfrom, xto, yfrom, yto,  I,   H/K     L
                      #[0.02, 0.075, 0.063, 0.11, 5e5, 1, 3],#103
                      #[0.025, 0.08,  0.11, 0.18, 5e5, 1, 5],#105
                      #[0.07, 0.11, 0.044, 0.08, 5e5, 2, 2],#202
                      #[0.065, 0.12, 0.09, 0.15, 5e5, 2, 4],#204              
                    #],
                    ##  x-from, x-to, y-from, y-to
           #'mirror_peaks': True, # add peaks from the left side
           #'use_rotation': False, 
           #'correction': ff_cube,
           #'remove_background': True,
           #'param_file': 'P227_Cubes_short_ai03_timedependence',
           #'name': 'P227_Cubes_short_ai03',
           #'ii_region': [0.025, 0.07, 0.12, 0.17], # Integrated inteisity around 105
           #'incregion': [0.02 , 0.03, 0.11, 0.13],
          #},  
          
   #'P227ai03test1': {                 # General peak parameters
          #'initial_fit': [# a*         c*         γ_y        γ_z       
                          #4.790e-02,  2.790e-02,  1.563e-05,  1.099e-02,
                          ##σ_y        σ_z        Qy_off     Qz_off    
                          #1.118e-01,  0.05, -2.740e-03,  5.526e-03,
                          ##Qphi_off   width      height     decay     
                          #-4.281e-01,  0.00497968, 0.000754677, 3195.79,
                          ##α_i        α_c        BG         I_0       
                          #3.000e-01,  1.488e-03,  5.000e+00,
                         #],
          #'peaks': [  # parameters for each single peak
                      ## xfrom, xto, yfrom, yto,  I,   H/K     L
                      #[0.02, 0.075, 0.063, 0.11, 5e5, 1, 3],#103
                      #[0.025, 0.08,  0.11, 0.18, 5e5, 1, 5],#105
                      #[0.07, 0.11, 0.044, 0.08, 5e5, 2, 2],#202
                      #[0.065, 0.12, 0.09, 0.15, 5e5, 2, 4],#204              
                    #],
                    ##  x-from, x-to, y-from, y-to
           #'mirror_peaks': True, # add peaks from the left side
           #'use_rotation': False, 
           #'correction': ff_cube,
           #'remove_background': True,
           #'param_file': 'P227_Cubes_short_ai03_timedependence_test1',
           #'name': 'P227_Cubes_short_ai03_test1',
           #'ii_region': [0.025, 0.07, 0.12, 0.17], # Integrated inteisity around 105
           #'incregion': [0.02 , 0.03, 0.11, 0.13],
          #},   
          
   #'P232ai03': {                 # General peak parameters
          #'initial_fit': [ # a*         c*         γ_y        γ_z 
                          #4.790e-02,  2.790e-02,  1.563e-05,  1.099e-02,
                          ##σ_y        σ_z        Qy_off     Qz_off    
                          #1.118e-01,  0.05, -2.740e-03,  5.526e-03,
                          ##Qphi_off   width      height     decay     
                          #-4.281e-01,  0.00497968, 0.000754677, 3195.79,
                          ##α_i        α_c        BG         I_0       
                          #3.000e-01,  1.488e-03,  5.000e+00,
                         #],
          #'peaks': [  # parameters for each single peak
                      ## xfrom, xto, yfrom, yto,  I,   H/K     L
                      #[0.02, 0.08, 0.063, 0.11, 5e5, 1, 3],#103
                      #[0.02, 0.08,  0.11, 0.18, 5e5, 1, 5],#105
                      #[0.07, 0.11, 0.044, 0.08, 5e5, 2, 2],#202
                      #[0.065, 0.12, 0.09, 0.15, 5e5, 2, 4],#204    
                      
                    #],
                    ##  x-from, x-to, y-from, y-to
           #'mirror_peaks': True, # add peaks from the left side
           #'use_rotation': False, 
           #'correction': ff_cube,
           #'remove_background': True,
           #'param_file': 'P232_Cubes_short_wout_ai03_timedependence',
           #'name': 'P232_Cubes_short_wout_ai03',
           #'ii_region': [0.025, 0.07, 0.12, 0.17], # Integrated inteisity around 105
           #'incregion': [0.02 , 0.03, 0.11, 0.13],
          #},  
  #'P225ai03': {                 # General peak parameters
          #'initial_fit': [ # a*         c*         γ_y        γ_z 
                          #4.475e-02, 2.644e-02,  4.114e-05,  3.340e-03,
                          ##σ_y        σ_z        Qy_off     Qz_off    
                          #7.325e-03, -5.161e-04, -1.230e-03,  3.028e-03,
                          ##Qphi_off   width      height     decay     
                          #-6.366e-02,  4.980e-03,  7.547e-04,  3.196e+03,
                          ##α_i        α_c        BG      
                          #3.000e-01, 8.997e-02,  1.000e+01,
                   #],
          #'peaks': [  # parameters for each single peak
                      ## xfrom, xto, yfrom, yto,  I,   H/K     L
                      #[0.024, 0.062, 0.065, 0.1, 5e5, 1, 3],#103
                      #[0.025, 0.06,  0.12, 0.15, 5e5, 1, 5],#105
                      #[0.07, 0.095, 0.044, 0.08, 5e5, 2, 2],#202
                      #[0.069, 0.1, 0.09, 0.125, 5e5, 2, 4],#204              
                      #[0.05, 0.075, 0.09, 0.13, 5e5,  np.sqrt(2), 4],#114
                      #[0.04, 0.09,  0.145, 0.18, 5e5, np.sqrt(2), 6],#116
                    #],
                    ##  x-from, x-to, y-from, y-to
           #'mirror_peaks': True, # add peaks from the left side
           #'use_rotation': False, 
           #'correction': ff_cube,
           #'remove_background': True,
           #'param_file': 'P225_Cubes_long_wout_ai03_timedependence',
           #'name': 'P225_Cubes_long_ai03',
           #'ii_region': [0.025, 0.057, 0.12, 0.15], # Integrated inteisity around 105
           #'incregion': [0.02 , 0.03, 0.11, 0.13],
          #}, 
          
   #'P225ai03test1': {                 # General peak parameters
          #'initial_fit': [ # a*         c*         γ_y        γ_z 
                          #4.475e-02, 2.644e-02,  4.114e-05,  3.340e-03,
                          ##σ_y        σ_z        Qy_off     Qz_off    
                          #7.325e-03, -5.161e-04, -1.230e-03,  3.028e-03,
                          ##Qphi_off   width      height     decay     
                          #-6.366e-02,  4.980e-03,  7.547e-04,  3.196e+03,
                          ##α_i        α_c        BG      
                          #3.000e-01, 8.997e-02,  1.000e+01,
                   #],
          #'peaks': [  # parameters for each single peak
                      ## xfrom, xto, yfrom, yto,  I,   H/K     L
                      #[0.024, 0.062, 0.065, 0.1, 5e5, 1, 3],#103
                      #[0.025, 0.06,  0.12, 0.15, 5e5, 1, 5],#105
                      #[0.07, 0.095, 0.044, 0.08, 5e5, 2, 2],#202
                      #[0.069, 0.1, 0.09, 0.125, 5e5, 2, 4],#204              
                      #[0.05, 0.075, 0.09, 0.13, 5e5,  np.sqrt(2), 4],#114
                      #[0.04, 0.09,  0.145, 0.18, 5e5, np.sqrt(2), 6],#116
                    #],
                    ##  x-from, x-to, y-from, y-to
           #'mirror_peaks': True, # add peaks from the left side
           #'use_rotation': False, 
           #'correction': ff_cube,
           #'remove_background': True,
           #'param_file': 'P225_Cubes_long_wout_ai03_timedependence_test1',
           #'name': 'P225_Cubes_long_ai03_test1',
           #'ii_region': [0.025, 0.057, 0.12, 0.15], # Integrated inteisity around 105
           #'incregion': [0.02 , 0.03, 0.11, 0.13],
          #},         
          
  }
