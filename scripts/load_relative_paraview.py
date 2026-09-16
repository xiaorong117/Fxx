"""pvpython load_paraview.py OUTPUT [--render]; or execute in ParaView Python shell."""
from pathlib import Path
import sys,json
from paraview.simple import *
out=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
data=PVDReader(FileName=str(out/'states.pvd'));data.UpdatePipeline()
scene=GetAnimationScene();scene.UpdateAnimationUsingDataTimeSteps()
view=GetActiveViewOrCreate('RenderView');view.ViewSize=[1400,1000]
network=Show(data,view);network.DiffuseColor=[0.5,0.5,0.5];network.Opacity=0.15;ColorBy(network,None)
for name,color in [('inlet_boundary.vtp',[1.,0.,0.]),('outlet_boundary.vtp',[0.,0.,1.])]:
    reader=XMLPolyDataReader(FileName=[str(out/name)])
    g=Glyph(Input=reader,GlyphType='Sphere');g.GlyphType.Radius=1.;g.ScaleArray=['POINTS','radius'];g.ScaleFactor=1.;g.GlyphMode='All Points'
    show=Show(g,view);ColorBy(show,None);show.DiffuseColor=color
connections=XMLPolyDataReader(FileName=[str(out/'boundary_connections.vtp')]);Show(connections,view).LineWidth=2
active=Threshold(Input=data);active.Scalars=['POINTS','active_gas'];active.LowerThreshold=0.5;active.UpperThreshold=1
g=Glyph(Input=active,GlyphType='Sphere');g.GlyphType.Radius=1.;g.ScaleArray=['POINTS','bubble_equivalent_radius'];g.ScaleFactor=1.;g.GlyphMode='All Points'
display=Show(g,view);ColorBy(display,('POINTS','Sg'));GetColorTransferFunction('Sg').RescaleTransferFunction(0,1)
display.SetScalarBarVisibility(view,True)
annotation=AnnotateTimeFilter(Input=data);Show(annotation,view)
ResetCamera(view);camera=(list(view.CameraPosition),list(view.CameraFocalPoint),list(view.CameraViewUp))
if '--render' in sys.argv:
    for index,t in enumerate([data.TimestepValues[0],data.TimestepValues[-1]]):
        scene.AnimationTime=t;view.ViewTime=t
        view.CameraPosition,view.CameraFocalPoint,view.CameraViewUp=camera
        SaveScreenshot(str(out/f'bubbles_{index}.png'),view)
    Hide(g,view);network.Opacity=0.6;ColorBy(network,('POINTS','C'));GetColorTransferFunction('C').RescaleTransferFunction(0,0.8)
    SaveScreenshot(str(out/'concentration.png'),view)
    ColorBy(network,('POINTS','Pl_relative_Pa'));GetColorTransferFunction('Pl_relative_Pa').RescaleTransferFunction(-10,110)
    SaveScreenshot(str(out/'pressure_boundaries.png'),view)
    Hide(data,view);local=PVDReader(FileName=str(out/'local_states.pvd'));local.UpdatePipeline();ld=Show(local,view)
    ColorBy(ld,('CELLS','advective_component_flux_mol_s'));GetColorTransferFunction('advective_component_flux_mol_s').RescaleTransferFunction(-1e-9,1e-9)
    ResetCamera(view);SaveScreenshot(str(out/'local_advection.png'),view)
    ColorBy(ld,('CELLS','diffusive_component_flux_mol_s'));GetColorTransferFunction('diffusive_component_flux_mol_s').RescaleTransferFunction(-1e-11,1e-11)
    SaveScreenshot(str(out/'local_diffusion.png'),view)
SaveState(str(out/'relative_pressure.pvsm'))
