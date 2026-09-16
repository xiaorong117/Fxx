#!/usr/bin/env python3
"""Convert legacy PolyData to verified XML/PVD and export fixed-scale views."""
import csv
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'relative-vis-deps'))
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk


def read(path):
    if path.suffix=='.vtp': r=vtk.vtkXMLPolyDataReader()
    else: r=vtk.vtkPolyDataReader();r.ReadAllScalarsOn();r.ReadAllVectorsOn()
    r.SetFileName(str(path));r.Update()
    if r.GetErrorCode() or r.GetOutput().GetNumberOfPoints()==0: raise RuntimeError(f'VTK read failed {path}')
    result=vtk.vtkPolyData();result.DeepCopy(r.GetOutput());return result


def save(path,data):
    w=vtk.vtkXMLPolyDataWriter();w.SetFileName(str(path));w.SetInputData(data)
    w.SetDataModeToAppended();w.SetCompressorTypeToZLib()
    if not w.Write(): raise RuntimeError(f'VTK write failed {path}')
    verified=read(path)
    if (verified.GetNumberOfPoints()!=data.GetNumberOfPoints() or
            verified.GetNumberOfCells()!=data.GetNumberOfCells()):
        raise RuntimeError(f'XML topology readback mismatch: {path}')
    for label,getter in [('point',vtk.vtkPolyData.GetPointData),
                         ('cell',vtk.vtkPolyData.GetCellData)]:
        expected=getter(data);actual=getter(verified)
        expected_names={expected.GetArrayName(i) for i in range(expected.GetNumberOfArrays())}
        actual_names={actual.GetArrayName(i) for i in range(actual.GetNumberOfArrays())}
        if expected_names!=actual_names:
            raise RuntimeError(f'XML {label}-field readback mismatch: {path}')
        for name in expected_names:
            a,b=expected.GetArray(name),actual.GetArray(name)
            if (a.GetNumberOfTuples()!=b.GetNumberOfTuples() or
                    a.GetNumberOfComponents()!=b.GetNumberOfComponents()):
                raise RuntimeError(f'XML {label}-array shape mismatch: {path}: {name}')


def point_subset(data,indices):
    out=vtk.vtkPolyData();pts=vtk.vtkPoints();cells=vtk.vtkCellArray()
    for k,i in enumerate(indices): pts.InsertNextPoint(data.GetPoint(int(i)));cells.InsertNextCell(1);cells.InsertCellPoint(k)
    out.SetPoints(pts);out.SetVerts(cells)
    for j in range(data.GetPointData().GetNumberOfArrays()):
        a=data.GetPointData().GetArray(j)
        b=numpy_to_vtk(vtk_to_numpy(a)[indices].copy(),deep=True);b.SetName(a.GetName());out.GetPointData().AddArray(b)
    return out


def line_subset(data,indices):
    if len(indices)==0:
        raise RuntimeError('local region contains no internal edge')
    used=sorted({data.GetCell(int(i)).GetPointId(j) for i in indices
                 for j in range(data.GetCell(int(i)).GetNumberOfPoints())})
    remap={old:new for new,old in enumerate(used)}
    out=vtk.vtkPolyData();points=vtk.vtkPoints();lines=vtk.vtkCellArray()
    for old in used: points.InsertNextPoint(data.GetPoint(old))
    out.SetPoints(points)
    for i in indices:
        cell=data.GetCell(int(i));lines.InsertNextCell(cell.GetNumberOfPoints())
        for j in range(cell.GetNumberOfPoints()): lines.InsertCellPoint(remap[cell.GetPointId(j)])
    out.SetLines(lines)
    used_array=np.asarray(used,dtype=int)
    for j in range(data.GetPointData().GetNumberOfArrays()):
        a=data.GetPointData().GetArray(j)
        b=numpy_to_vtk(vtk_to_numpy(a)[used_array].copy(),deep=True)
        b.SetName(a.GetName());out.GetPointData().AddArray(b)
    for j in range(data.GetCellData().GetNumberOfArrays()):
        a=data.GetCellData().GetArray(j);b=numpy_to_vtk(vtk_to_numpy(a)[indices].copy(),deep=True)
        b.SetName(a.GetName());out.GetCellData().AddArray(b)
    return out


def choose_local_region(data,virtual,component,main_component):
    pd=data.GetPointData();coords=vtk_to_numpy(data.GetPoints().GetData())
    main_real=np.flatnonzero((~virtual)&(component==main_component))
    if not len(main_real): raise RuntimeError('main-flow component contains no real node')
    bubbles=np.flatnonzero((vtk_to_numpy(pd.GetArray('active_gas'))>0)&
                           (component==main_component))
    selected=int(bubbles[len(bubbles)//2]) if len(bubbles) else int(main_real[len(main_real)//2])
    reason='active_main_flow_bubble' if len(bubbles) else 'fallback_main_flow_real_node_no_active_bubble'
    centre=coords[selected];radius=max(np.ptp(coords[:,0])*0.08,1e-5)
    local_ids=np.asarray([],dtype=int);local_edges=[]
    for _ in range(12):
        local_ids=np.flatnonzero((np.linalg.norm(coords-centre,axis=1)<=radius)&
                                 (~virtual)&(component==main_component))
        local_set=set(int(x) for x in local_ids)
        local_edges=[i for i in range(data.GetNumberOfCells())
                     if data.GetCell(i).GetNumberOfPoints()==2 and
                     all(data.GetCell(i).GetPointId(j) in local_set for j in range(2))]
        if local_edges: return selected,reason,centre,local_ids,local_edges
        radius*=2
    raise RuntimeError('unable to construct a nonempty main-flow local edge region')


def main(output):
    output=Path(output);series=ET.parse(output/'states.pvd').getroot().find('Collection')
    entries=[(float(x.attrib['timestep']),x.attrib['file']) for x in series]
    # Solver may append legacy files after previous XML conversion.
    entries=sorted(set(entries));xml=ET.Element('VTKFile',type='Collection',version='0.1',byte_order='LittleEndian')
    collection=ET.SubElement(xml,'Collection');local_xml=ET.Element('VTKFile',type='Collection',version='0.1',byte_order='LittleEndian');lc=ET.SubElement(local_xml,'Collection')
    first=read(output/entries[0][1]);pd=first.GetPointData()
    ids=vtk_to_numpy(pd.GetArray('node_id')).astype(int);flags=vtk_to_numpy(pd.GetArray('boundary_flag'))
    virtual=vtk_to_numpy(pd.GetArray('is_virtual_boundary')).astype(bool)
    component=vtk_to_numpy(pd.GetArray('component_id')).astype(int)
    comp=json.loads((output/'component_boundary_and_makeup_audit.json').read_text())['main_flow_component_id']
    selected,selection_reason,centre,local_ids,local_edges=choose_local_region(
        first,virtual,component,comp)
    boundary_edges=np.flatnonzero(vtk_to_numpy(first.GetCellData().GetArray('is_boundary_edge'))>0)
    manifest=json.loads((output/'run_manifest.json').read_text())
    pore_path=Path(manifest['configuration']['io']['pore_file'])
    if not pore_path.is_absolute(): pore_path=PROJECT.parent/pore_path
    pores=np.loadtxt(pore_path)
    scale=manifest['configuration']['geometry']['input_length_to_m']
    for name,flag in [('inlet_boundary.vtp',1),('outlet_boundary.vtp',2)]:
        ix=np.flatnonzero(virtual&(flags==flag));data=point_subset(first,ix)
        for field,values in [('radius',pores[ix,6]*scale),('connected_real_node_id',pores[ix,2])]:
            arr=numpy_to_vtk(values.copy(),deep=True);arr.SetName(field);data.GetPointData().AddArray(arr)
        save(output/name,data)
    save(output/'boundary_connections.vtp',line_subset(first,boundary_edges))
    for t,filename in entries:
        src=output/filename;dst=src.with_suffix('.vtp');data=read(src)
        if not dst.exists() or dst.stat().st_mtime<src.stat().st_mtime: save(dst,data)
        ET.SubElement(collection,'DataSet',timestep=str(t),group='',part='0',file=dst.name)
        local_path=output/('local_'+dst.name);save(local_path,line_subset(data,local_edges))
        ET.SubElement(lc,'DataSet',timestep=str(t),group='',part='0',file=local_path.name)
        # Every frame: verify derived concentration and each boundary flux, not just their sum.
        a=data.GetPointData();c=vtk_to_numpy(a.GetArray('C'));q=vtk_to_numpy(data.GetCellData().GetArray('liquid_flux_m3_s'))
        fa=vtk_to_numpy(data.GetCellData().GetArray('advective_component_flux_mol_s'));fd=vtk_to_numpy(data.GetCellData().GetArray('diffusive_component_flux_mol_s'))
        for i in boundary_edges:
            cell=data.GetCell(int(i));x,y=cell.GetPointId(0),cell.GetPointId(1);r=y if virtual[x] else x
            if c[x]!=c[y]:
                raise RuntimeError('per-edge copied concentration audit failed')
            if t>0:
                if fd[i]!=0 or not np.isclose(fa[i],q[i]*c[r],rtol=1e-12,atol=1e-30):
                    raise RuntimeError('per-edge copy concentration/flux audit failed')
    ET.ElementTree(xml).write(output/'states.pvd',encoding='utf-8',xml_declaration=True)
    ET.ElementTree(local_xml).write(output/'local_states.pvd',encoding='utf-8',xml_declaration=True)
    written=[(float(x.attrib['timestep']),x.attrib['file']) for x in ET.parse(output/'states.pvd').getroot().find('Collection')]
    if written!=[(t,Path(f).with_suffix('.vtp').name) for t,f in entries]:
        raise RuntimeError('PVD physical-time readback mismatch')
    local_verified=read(output/('local_'+Path(entries[0][1]).with_suffix('.vtp').name))
    local_node_ids=set(vtk_to_numpy(local_verified.GetPointData().GetArray('node_id')).astype(int).tolist())
    expected_local_node_ids={int(ids[p]) for edge in local_edges for p in
                             [first.GetCell(edge).GetPointId(0),first.GetCell(edge).GetPointId(1)]}
    if local_node_ids!=expected_local_node_ids:
        raise RuntimeError('local VTP node-ID/connectivity readback mismatch')
    (output/'visualization_audit.json').write_text(json.dumps({'xml_readback_passed':True,'vtk_version':vtk.vtkVersion.GetVTKVersion(),
        'local_node_ids':sorted(local_node_ids),'representative_node_id':int(ids[selected]),
        'representative_selection':selection_reason,'local_edge_count':len(local_edges),
        'main_flow_component':int(comp),'local_bounds_centre_m':centre.tolist(),
        'bubble_glyph':'Sphere radius=1; scale factor=1; scale array=bubble_equivalent_radius; physical radius, no magnification',
        'paraview_render_validated':False,'frame_count':len(entries)},indent=2))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=list(csv.DictReader((output/'diagnostics.tsv').open(),delimiter='\t'));ts=[float(r['time_s']) for r in rows]
    fig,axes=plt.subplots(3,2,figsize=(12,10))
    groups=[['free_mol_main_flow_component','dissolved_mol_total'],['gas_volume_main_flow_component'],
            ['active_main_flow'],['cumulative_component_inflow_mol','cumulative_component_outflow_mol'],
            ['cumulative_liquid_volume_change','cumulative_net_liquid_inflow'],['boundary_advective_component_outflow','boundary_diffusive_component_outflow']]
    for ax,keys in zip(axes.flat,groups):
        for k in keys: ax.plot(ts,[float(r[k]) for r in rows],label=k)
        ax.set_xlabel('physical time [s]');ax.legend(fontsize=7);ax.grid(True)
    fig.tight_layout();fig.savefig(output/'history.png',dpi=160);fig.savefig(output/'history.pdf');plt.close(fig)
    import shutil,subprocess
    if shutil.which('pvpython') and shutil.which('xvfb-run'):
        with (output/'paraview_read_render.log').open('w') as log:
            rc=subprocess.run(['xvfb-run','-a','pvpython',str(PROJECT/'scripts/load_relative_paraview.py'),str(output),'--render'],stdout=log,stderr=subprocess.STDOUT,timeout=180).returncode
        expected_images=[output/name for name in ('bubbles_0.png','bubbles_1.png','concentration.png',
            'pressure_boundaries.png','local_advection.png','local_diffusion.png')]
        images_valid=rc==0
        if images_valid:
            from PIL import Image
            for path in expected_images:
                try:
                    with Image.open(path) as picture: picture.verify()
                    images_valid=images_valid and path.stat().st_size>0
                except Exception: images_valid=False
        audit=json.loads((output/'visualization_audit.json').read_text())
        audit['paraview_render_validated']=images_valid
        audit['paraview_render_exit_code']=rc
        audit['paraview_expected_images']=[str(p) for p in expected_images]
        (output/'visualization_audit.json').write_text(json.dumps(audit,indent=2))
    shutil.copyfile(PROJECT/'scripts/load_relative_paraview.py',output/'load_paraview.py')
    import hashlib
    for p in output.iterdir():
        if p.is_file() and p.name!='run_manifest.json' and (p.suffix in ('.vtp','.pvd','.png','.pdf','.pvsm') or p.name in ('load_paraview.py','visualization_audit.json','paraview_read_render.log')):
            manifest['outputs'][str(p.resolve())]=hashlib.sha256(p.read_bytes()).hexdigest()
    (output/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


PROJECT=Path(__file__).resolve().parents[1]
if __name__=='__main__': main(sys.argv[1])
