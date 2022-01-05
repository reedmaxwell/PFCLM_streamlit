# Import the ParFlow package
#
from parflow import Run
import os
import shutil 
import streamlit as st
from parflow.tools.fs import get_absolute_path
from parflowio.pyParflowio import PFData
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio

# set up directories
base_dir = get_absolute_path(".")
print(base_dir)

st.set_page_config(page_title="Interactive PF-SC", page_icon=None, layout='centered', initial_sidebar_state='auto')
### run ParFlow
st.write(" Set up and Run ParFlow Single Column")


# Set our Run Name 
PFCLM_SC = Run("PFCLM_SC")


stopt = 8760

#-----------------------------------------------------------------------------
# File input version number
#-----------------------------------------------------------------------------
PFCLM_SC.FileVersion = 4

#-----------------------------------------------------------------------------
# Process Topology
#-----------------------------------------------------------------------------

PFCLM_SC.Process.Topology.P = 1
PFCLM_SC.Process.Topology.Q = 1
PFCLM_SC.Process.Topology.R = 1

#-----------------------------------------------------------------------------
# Computational Grid
#-----------------------------------------------------------------------------
PFCLM_SC.ComputationalGrid.Lower.X = 0.0
PFCLM_SC.ComputationalGrid.Lower.Y = 0.0
PFCLM_SC.ComputationalGrid.Lower.Z = 0.0

PFCLM_SC.ComputationalGrid.DX      = 2.0
PFCLM_SC.ComputationalGrid.DY      = 2.0
PFCLM_SC.ComputationalGrid.DZ      = 0.1

PFCLM_SC.ComputationalGrid.NX      = 1
PFCLM_SC.ComputationalGrid.NY      = 1
PFCLM_SC.ComputationalGrid.NZ      = 20

#-----------------------------------------------------------------------------
# The Names of the GeomInputs
#-----------------------------------------------------------------------------
PFCLM_SC.GeomInput.Names = 'domain_input'

#-----------------------------------------------------------------------------
# Domain Geometry Input
#-----------------------------------------------------------------------------
PFCLM_SC.GeomInput.domain_input.InputType = 'Box'
PFCLM_SC.GeomInput.domain_input.GeomName  = 'domain'

#-----------------------------------------------------------------------------
# Domain Geometry
#-----------------------------------------------------------------------------
PFCLM_SC.Geom.domain.Lower.X = 0.0
PFCLM_SC.Geom.domain.Lower.Y = 0.0
PFCLM_SC.Geom.domain.Lower.Z = 0.0

PFCLM_SC.Geom.domain.Upper.X = 2.0
PFCLM_SC.Geom.domain.Upper.Y = 2.0
PFCLM_SC.Geom.domain.Upper.Z = 2.0

PFCLM_SC.Geom.domain.Patches = 'x_lower x_upper y_lower y_upper z_lower z_upper'


#--------------------------------------------
# variable dz assignments
#------------------------------------------

PFCLM_SC.Solver.Nonlinear.VariableDz = True
PFCLM_SC.dzScale.GeomNames           = 'domain'
PFCLM_SC.dzScale.Type                = 'nzList'
PFCLM_SC.dzScale.nzListNumber        = 20

# cells start at the bottom (0) and moves up to the top
# domain is 3.21 m thick, root zone is down to 19 cells 
# so the root zone is 2.21 m thick
PFCLM_SC.Cell._0.dzScale.Value  = 10.0   # first cell is 10*0.1 1m thick
PFCLM_SC.Cell._1.dzScale.Value  = 5.0    # next cell is 5*0.1 50 cm thick
PFCLM_SC.Cell._2.dzScale.Value  = 1.0   
PFCLM_SC.Cell._3.dzScale.Value  = 1.0
PFCLM_SC.Cell._4.dzScale.Value  = 1.0
PFCLM_SC.Cell._5.dzScale.Value  = 1.0
PFCLM_SC.Cell._6.dzScale.Value  = 1.0
PFCLM_SC.Cell._7.dzScale.Value  = 1.0
PFCLM_SC.Cell._8.dzScale.Value  = 1.0
PFCLM_SC.Cell._9.dzScale.Value  = 1.0
PFCLM_SC.Cell._10.dzScale.Value = 1.0
PFCLM_SC.Cell._11.dzScale.Value = 1.0
PFCLM_SC.Cell._12.dzScale.Value = 1.0
PFCLM_SC.Cell._13.dzScale.Value = 1.0
PFCLM_SC.Cell._14.dzScale.Value = 1.0
PFCLM_SC.Cell._15.dzScale.Value = 1.0
PFCLM_SC.Cell._16.dzScale.Value = 1.0
PFCLM_SC.Cell._17.dzScale.Value = 1.0
PFCLM_SC.Cell._18.dzScale.Value = 1.0
PFCLM_SC.Cell._19.dzScale.Value = 0.1   #0.1* 0.1 = 0.01  1 cm top layer

#-----------------------------------------------------------------------------
# Perm
#-----------------------------------------------------------------------------
PFCLM_SC.Geom.Perm.Names              = 'domain'
PFCLM_SC.Geom.domain.Perm.Type        = 'Constant'
PFCLM_SC.Geom.domain.Perm.Value       = 0.001465
PFCLM_SC.Perm.TensorType              = 'TensorByGeom'
PFCLM_SC.Geom.Perm.TensorByGeom.Names = 'domain'
PFCLM_SC.Geom.domain.Perm.TensorValX  = 1.0
PFCLM_SC.Geom.domain.Perm.TensorValY  = 1.0
PFCLM_SC.Geom.domain.Perm.TensorValZ  = 1.0

#-----------------------------------------------------------------------------
# Specific Storage
#-----------------------------------------------------------------------------

PFCLM_SC.SpecificStorage.Type              = 'Constant'
PFCLM_SC.SpecificStorage.GeomNames         = 'domain'
PFCLM_SC.Geom.domain.SpecificStorage.Value = 1.0e-4

#-----------------------------------------------------------------------------
# Phases
#-----------------------------------------------------------------------------

PFCLM_SC.Phase.Names = 'water'

PFCLM_SC.Phase.water.Density.Type     = 'Constant'
PFCLM_SC.Phase.water.Density.Value    = 1.0

PFCLM_SC.Phase.water.Viscosity.Type   = 'Constant'
PFCLM_SC.Phase.water.Viscosity.Value  = 1.0

#-----------------------------------------------------------------------------
# Contaminants
#-----------------------------------------------------------------------------
PFCLM_SC.Contaminants.Names = ''


#-----------------------------------------------------------------------------
# Gravity
#-----------------------------------------------------------------------------

PFCLM_SC.Gravity = 1.0

#-----------------------------------------------------------------------------
# Setup timing info
#-----------------------------------------------------------------------------

PFCLM_SC.TimingInfo.BaseUnit     = 1.0
PFCLM_SC.TimingInfo.StartCount   = 0
PFCLM_SC.TimingInfo.StartTime    = 0.0
PFCLM_SC.TimingInfo.StopTime     = stopt
PFCLM_SC.TimingInfo.DumpInterval = 1.0
PFCLM_SC.TimeStep.Type           = 'Constant'
PFCLM_SC.TimeStep.Value          = 1.0


#-----------------------------------------------------------------------------
# Porosity
#-----------------------------------------------------------------------------

PFCLM_SC.Geom.Porosity.GeomNames    = 'domain'

PFCLM_SC.Geom.domain.Porosity.Type  = 'Constant'
PFCLM_SC.Geom.domain.Porosity.Value = 0.3

#-----------------------------------------------------------------------------
# Domain
#-----------------------------------------------------------------------------
PFCLM_SC.Domain.GeomName = 'domain'

#-----------------------------------------------------------------------------
# Mobility
#-----------------------------------------------------------------------------
PFCLM_SC.Phase.water.Mobility.Type  = 'Constant'
PFCLM_SC.Phase.water.Mobility.Value = 1.0

#-----------------------------------------------------------------------------
# Relative Permeability
#-----------------------------------------------------------------------------

PFCLM_SC.Phase.RelPerm.Type        = 'VanGenuchten'
PFCLM_SC.Phase.RelPerm.GeomNames   = 'domain'

PFCLM_SC.Geom.domain.RelPerm.Alpha = 2.0
PFCLM_SC.Geom.domain.RelPerm.N     = 2.0

#---------------------------------------------------------
# Saturation
#---------------------------------------------------------

PFCLM_SC.Phase.Saturation.Type        = 'VanGenuchten'
PFCLM_SC.Phase.Saturation.GeomNames   = 'domain'

PFCLM_SC.Geom.domain.Saturation.Alpha = 2.0
PFCLM_SC.Geom.domain.Saturation.N     = 3.0
PFCLM_SC.Geom.domain.Saturation.SRes  = 0.2
PFCLM_SC.Geom.domain.Saturation.SSat  = 1.0

#-----------------------------------------------------------------------------
# Wells
#-----------------------------------------------------------------------------
PFCLM_SC.Wells.Names = ''


#-----------------------------------------------------------------------------
# Time Cycles
#-----------------------------------------------------------------------------
PFCLM_SC.Cycle.Names = 'constant'
PFCLM_SC.Cycle.constant.Names = 'alltime'
PFCLM_SC.Cycle.constant.alltime.Length = 1
PFCLM_SC.Cycle.constant.Repeat = -1

#-----------------------------------------------------------------------------
# Boundary Conditions: Pressure
#-----------------------------------------------------------------------------
PFCLM_SC.BCPressure.PatchNames = 'x_lower x_upper y_lower y_upper z_lower z_upper'

PFCLM_SC.Patch.x_lower.BCPressure.Type          = 'FluxConst'
PFCLM_SC.Patch.x_lower.BCPressure.Cycle         = 'constant'
PFCLM_SC.Patch.x_lower.BCPressure.alltime.Value = 0.0

PFCLM_SC.Patch.y_lower.BCPressure.Type          = 'FluxConst'
PFCLM_SC.Patch.y_lower.BCPressure.Cycle         = 'constant'
PFCLM_SC.Patch.y_lower.BCPressure.alltime.Value = 0.0

#PFCLM_SC.Patch.z_lower.BCPressure.Type = 'FluxConst'
PFCLM_SC.Patch.z_lower.BCPressure.Type          = 'DirEquilRefPatch'
PFCLM_SC.Patch.z_lower.BCPressure.RefGeom       = 'domain'
PFCLM_SC.Patch.z_lower.BCPressure.RefPatch      = 'z_lower'
PFCLM_SC.Patch.z_lower.BCPressure.Cycle         = 'constant'
PFCLM_SC.Patch.z_lower.BCPressure.alltime.Value = 0.0

PFCLM_SC.Patch.x_upper.BCPressure.Type          = 'FluxConst'
PFCLM_SC.Patch.x_upper.BCPressure.Cycle         = 'constant'
PFCLM_SC.Patch.x_upper.BCPressure.alltime.Value = 0.0

PFCLM_SC.Patch.y_upper.BCPressure.Type          = 'FluxConst'
PFCLM_SC.Patch.y_upper.BCPressure.Cycle         = 'constant'
PFCLM_SC.Patch.y_upper.BCPressure.alltime.Value = 0.0

PFCLM_SC.Patch.z_upper.BCPressure.Type          = 'OverlandFlow'
PFCLM_SC.Patch.z_upper.BCPressure.Cycle         = 'constant'
PFCLM_SC.Patch.z_upper.BCPressure.alltime.Value = 0.0

#---------------------------------------------------------
# Topo slopes in x-direction
#---------------------------------------------------------

PFCLM_SC.TopoSlopesX.Type              = 'Constant'
PFCLM_SC.TopoSlopesX.GeomNames         = 'domain'
PFCLM_SC.TopoSlopesX.Geom.domain.Value = 0.05

#---------------------------------------------------------
# Topo slopes in y-direction
#---------------------------------------------------------

PFCLM_SC.TopoSlopesY.Type              = 'Constant'
PFCLM_SC.TopoSlopesY.GeomNames         = 'domain'
PFCLM_SC.TopoSlopesY.Geom.domain.Value = 0.00

#---------------------------------------------------------
# Mannings coefficient
#---------------------------------------------------------

PFCLM_SC.Mannings.Type               = 'Constant'
PFCLM_SC.Mannings.GeomNames          = 'domain'
PFCLM_SC.Mannings.Geom.domain.Value  = 2.e-6

#-----------------------------------------------------------------------------
# Phase sources:
#-----------------------------------------------------------------------------

PFCLM_SC.PhaseSources.water.Type              = 'Constant'
PFCLM_SC.PhaseSources.water.GeomNames         = 'domain'
PFCLM_SC.PhaseSources.water.Geom.domain.Value = 0.0

#-----------------------------------------------------------------------------
# Exact solution specification for error calculations
#-----------------------------------------------------------------------------

PFCLM_SC.KnownSolution = 'NoKnownSolution'

#-----------------------------------------------------------------------------
# Set solver parameters
#-----------------------------------------------------------------------------

PFCLM_SC.Solver         = 'Richards'
PFCLM_SC.Solver.MaxIter = 17600

PFCLM_SC.Solver.Nonlinear.MaxIter           = 100
PFCLM_SC.Solver.Nonlinear.ResidualTol       = 1e-5
PFCLM_SC.Solver.Nonlinear.EtaChoice         = 'Walker1'
PFCLM_SC.Solver.Nonlinear.EtaValue          = 0.01
PFCLM_SC.Solver.Nonlinear.UseJacobian       = False
PFCLM_SC.Solver.Nonlinear.DerivativeEpsilon = 1e-12
PFCLM_SC.Solver.Nonlinear.StepTol           = 1e-30
PFCLM_SC.Solver.Nonlinear.Globalization     = 'LineSearch'
PFCLM_SC.Solver.Linear.KrylovDimension      = 100
PFCLM_SC.Solver.Linear.MaxRestarts          = 5
PFCLM_SC.Solver.Linear.Preconditioner       = 'PFMG'
PFCLM_SC.Solver.PrintSubsurf                = False
PFCLM_SC.Solver.Drop                        = 1E-20
PFCLM_SC.Solver.AbsTol                      = 1E-9

#Writing output options for ParFlow
#  PFB only no SILO
PFCLM_SC.Solver.PrintSubsurfData         = True
PFCLM_SC.Solver.PrintPressure            = True
PFCLM_SC.Solver.PrintSaturation          = True
PFCLM_SC.Solver.PrintCLM                 = True
PFCLM_SC.Solver.PrintMask                = True
PFCLM_SC.Solver.PrintSpecificStorage     = True
PFCLM_SC.Solver.WriteSiloMannings        = False
PFCLM_SC.Solver.WriteSiloMask            = False
PFCLM_SC.Solver.WriteSiloSlopes          = False
PFCLM_SC.Solver.WriteSiloSaturation      = False

#---------------------------------------------------
# LSM / CLM options
#---------------------------------------------------

# set LSM options to CLM
PFCLM_SC.Solver.LSM              = 'CLM'
# specify type of forcing, file name and location
PFCLM_SC.Solver.CLM.MetForcing   = '1D'
#PFCLM_SC.Solver.CLM.MetFileName = 'forcing_1.txt'
#PFCLM_SC.Solver.CLM.MetFileName = 'WY2017_ph_forcing.txt'
PFCLM_SC.Solver.CLM.MetFileName = 'LocalForcing_LowSWE.txt'
#PFCLM_SC.Solver.CLM.MetFileName  = 'narr_1hr.txt'
PFCLM_SC.Solver.CLM.MetFilePath  = '/Users/reed/Projects/HydroFrame-ML/PFCLM_SC/forcing'

# Set CLM Plant Water Use Parameters
PFCLM_SC.Solver.CLM.EvapBeta       = 'Linear'
PFCLM_SC.Solver.CLM.VegWaterStress = 'Saturation'
PFCLM_SC.Solver.CLM.ResSat         = 0.2
PFCLM_SC.Solver.CLM.WiltingPoint   = 0.2
PFCLM_SC.Solver.CLM.FieldCapacity  = 1.00
PFCLM_SC.Solver.CLM.IrrigationType = 'none'
PFCLM_SC.Solver.CLM.RootZoneNZ     =  19
PFCLM_SC.Solver.CLM.SoiLayer       =  15

#Writing output options for CLM
#  PFB only, no SILO, no native CLM logs
PFCLM_SC.Solver.PrintLSMSink        = False
PFCLM_SC.Solver.CLM.CLMDumpInterval = 1
PFCLM_SC.Solver.CLM.CLMFileDir      = 'output/'
PFCLM_SC.Solver.CLM.BinaryOutDir    = False
PFCLM_SC.Solver.CLM.IstepStart      = 1
PFCLM_SC.Solver.WriteCLMBinary      = False
PFCLM_SC.Solver.WriteSiloCLM        = False
PFCLM_SC.Solver.CLM.WriteLogs       = False
PFCLM_SC.Solver.CLM.WriteLastRST    = True
PFCLM_SC.Solver.CLM.DailyRST        = False
PFCLM_SC.Solver.CLM.SingleFile      = True

#---------------------------------------------------
# Initial conditions: water pressure
#---------------------------------------------------

PFCLM_SC.ICPressure.Type                 = 'HydroStaticPatch'
PFCLM_SC.ICPressure.GeomNames            = 'domain'
PFCLM_SC.Geom.domain.ICPressure.Value    = -1.0
PFCLM_SC.Geom.domain.ICPressure.RefGeom  = 'domain'
PFCLM_SC.Geom.domain.ICPressure.RefPatch = 'z_upper'

os.chdir(base_dir)

# User can adjust parameters
PFCLM_SC.Geom.domain.Perm.Value = st.number_input("Hydraulic Conductivity Value [m/h]", min_value=1e-4, max_value=1.0, value=PFCLM_SC.Geom.domain.Perm.Value, step=None, format=None, key=None, help=None)
PFCLM_SC.Geom.domain.Porosity.Value = st.number_input("Porosity Value [-]", min_value=.1, max_value=.7, value=PFCLM_SC.Geom.domain.Porosity.Value, step=None, format=None, key=None, help=None)

#-----------------------------------------------------------------------------
# Run ParFlow 
#-----------------------------------------------------------------------------
control1 = st.button(" (re)Run ParFlow ")
if control1 == True:
    with st.spinner('Running ParFlow'):
        # copy CLM files
        os.chdir(base_dir+"/output")
        shutil.copyfile(base_dir+'/inputs/drv_clmin.dat', 'drv_clmin.dat')
        shutil.copyfile(base_dir+'/inputs/drv_vegm.dat', 'drv_vegm.dat')
        shutil.copyfile(base_dir+'/inputs/drv_vegp.dat', 'drv_vegp.dat')
        os.chdir(base_dir)
        # run PF
        PFCLM_SC.run(base_dir+"/output")

        st.success('Done!')




####
## Plot PF CLM results
####
control2 = st.button(" Plot Results ")
if control2 == True:
    # intialize data and time arrays
 # intialize data and time arrays
    data    = np.zeros([8,8760])  # an array where we store the PF output as columns
    time    = np.zeros([8760])    # time array, we will probably want to swap with a date
    forcing = np.zeros([8,8760])  # array to load in the meterological forcing

    # load forcing file as numpy array, note should fix so count is 1-8760 instead of 0-8759
    # variables are:
    # 0 DSWR: Downward Visible or Short-Wave radiation [W/m2].
    # 1 DLWR: Downward Infa-Red or Long-Wave radiation [W/m2]
    # 2 APCP: Precipitation rate [mm/s]
    # 3 Temp: Air temperature [K]
    # 4 UGRD: West-to-East or U-component of wind [m/s] 
    # 5 VGRD: South-to-North or V-component of wind [m/s]
    # 6 Press: Atmospheric Pressure [pa]
    # 7 SPFH: Water-vapor specific humidity [kg/kg]
    #
    #ffname = 'forcing/forcing_1.txt'
    ffname = (base_dir+'/forcing/LocalForcing_LowSWE.txt')
    forcing_ls = np.loadtxt(ffname,max_rows=8760)
    print(np.shape(forcing_ls))


    # loop over a year of files (8760 hours) and load in the CLM output
    # then map specific variables to the data array which holds things for analysis
    # and plotting
    # for icount in range(0, 8759):
    #     base = (base_dir+"/output/PFCLM_SC.out.clm_output.{:05d}.C.pfb")
    #     filename = base.format(icount+1)
    #     data_obj = PFData(filename)
    #     data_obj.loadHeader()
    #     data_obj.loadData()
    #     data_arr = data_obj.getDataAsArray()
    #     data_obj.close()
    #     data[0,1,icount] = data_arr[0,0,0]  #total (really, it is net) latent heat flux (Wm-2)
    #     data[0,2,icount] = data_arr[4,0,0]  #net veg. evaporation and transpiration and soil evaporation (mms-1)
    #     data[0,3,icount] = data_arr[10,0,0] #SWE (mm)
    #     # base = "output/PFCLM_SC.LS.out.press.{:05d}.pfb"
    #     # filename = base.format(icount)
    #     # data_obj = PFData(filename)
    #     # data_obj.loadHeader()
    #     # data_obj.loadData()
    #     # data_arr = data_obj.getDataAsArray()
    #     # data_obj.close()
    #     # data[4,icount] = (np.sqrt(slope)/mannings) * np.maximum(data_arr[19,0,0],0.0)**(5.0/3.0)
    #     time[icount] = icount
    slope    = 0.05
    mannings = 2.e-6

    for icount in range(1, 8760):
        base = (base_dir+"/output/PFCLM_SC.out.clm_output.{:05d}.C.pfb")
        filename = base.format(icount)
        data_obj = PFData(filename)
        data_obj.loadHeader()
        data_obj.loadData()
        data_arr = data_obj.getDataAsArray()
        data_obj.close()
        data[1,icount] = data_arr[0,0,0]  #net latent heat flux (Wm-2)
        data[2,icount] = data_arr[4,0,0]  #net veg. evaporation and transpiration and soil evaporation (mms-1)
        data[3,icount] = data_arr[10,0,0] #SWE (mm)
        data[5,icount] = data_arr[2,0,0]  #net sensible heat flux (Wm-2)
        data[6,icount] = data_arr[3,0,0]  #net sensible heat flux (Wm-2)
        base = (base_dir+"/output/PFCLM_SC.out.press.{:05d}.pfb")
        filename = base.format(icount)
        data_obj = PFData(filename)
        data_obj.loadHeader()
        data_obj.loadData()
        data_arr = data_obj.getDataAsArray()
        data_obj.close()
        data[4,icount] = (np.sqrt(slope)/mannings) * np.maximum(data_arr[19,0,0],0.0)**(5.0/3.0)
        time[icount] = icount

    # # Plot LH Flux, SWE and Runoff
    # fig, ax = plt.subplots()
    # ax2 = ax.twinx()
    # #ax.plot(time[1:17520],forcing[0:17519,2], color='g')
    # #ax2.plot(time[1:17520],data[3,1:17520], color='b')
    # #ax2.plot(time[1:17520], forcing[0:17519,3], color='r')
    # ax.plot(time,forcing_ls[:,2], color='g')
    # ax2.plot(time,data[0,3,:], color='b')
    # ax2.plot(time, forcing_ls[:,3], color='r')
    # ax.set_xlabel('Time, WY [hr]')
    # ax.set_ylabel('Precip [mm/s]')
    # ax2.set_ylabel('Temp [K], SWE[mm')
    # st.pyplot(fig)
    fig = make_subplots(rows=2, cols=2,column_widths=[0.5, 0.5],row_heights=[0.3, 0.7],horizontal_spacing=0.1, vertical_spacing=0.1,specs=[[{"type": "scatter","rowspan": 2}, {"type": "scattergeo"}],[None, {"type": "scatter","l":0,"r":0, "t":0, "b":0}]],subplot_titles=("P, ET, R","Site Location","LH, SH, G for Single Column Simulation"))
    #fig = go.Figure()
    fig.add_trace(go.Scatter(x=time[1:8760], y=(data[2,1:8760]*3.6), name='ET [m/h]', line=dict(color='red', width=1)), 1, 1)
    fig.add_trace(go.Scatter(x=time[1:8760], y=(forcing[1:8760,2]*3.6), name='precip [m/h]',line=dict(color='green', width=1)),1,1)
    fig.append_trace(go.Scatter(x=time[1:8760], y=data[4,1:8760], name = 'Runoff [m/h]',line=dict(color='blue', width=1)),1,1)
    fig.update_xaxes(title_text='Time [h] WY', row=1,col=1)
    fig.update_yaxes(title_text='Water Flux [m/h]', row=1,col=1)

    # domain location, should read in from metadata
    lat = (34.9604,) 
    lon = (-97.9789,)

    fig.add_trace(go.Scattergeo(lon=lon, lat=lat,name='Site Location',
        mode = 'markers', marker=dict(color='black')),1,2)
    fig.update_geos(
        resolution=50,
        showcoastlines=True, coastlinecolor="RebeccaPurple",
        showland=True, landcolor="LightGray",
        showocean=True, oceancolor="LightBlue",
        showlakes=True, lakecolor="Blue",
        showrivers=True, rivercolor="Blue")
    #fitbounds="locations")
    fig.update_layout(geo_scope='usa')  #,margin={"r":0,"t":0,"l":0,"b":0})

    #LH, SH, G plot
    #
    fig.add_trace(go.Scatter(x=time[1:8760], y=data[1,1:8760], name='LH [W/m2]',
                            line=dict(color='red', width=1)),2,2)
    fig.append_trace(go.Scatter(x=time[1:8760], y=data[5,1:8760], name = 'SH [W/m2]',
                            line=dict(color='blue', width=1)),2,2)
    fig.append_trace(go.Scatter(x=time[1:8760], y=data[6,1:8760], name='G [W/m2]',
                            line=dict(color='green', width=1) 
    ),2,2)

    fig.update_xaxes(title_text='Time [h] WY',row=2, col=2)
    fig.update_yaxes(title_text='Energy Flux [W/m2]', row=2, col=2)

    # Set theme, margin, and annotation in layout
    fig.update_layout(
        width=1000,
        height=500)
    #    margin=dict(
    #        l=10,
    #        r=10,
    #        b=10,
    #        t=30,
    #        pad=0
    #    ))
    st.plotly_chart(fig)
