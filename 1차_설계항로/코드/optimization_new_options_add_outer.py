"""Outer IAF -> IF ring -> FAF, and opposite-FAF Departure, with new-options outputs.

Run: python -B 수정본_20261103/코드/optimization_new_options_add_outer.py

Reuses optimization_new_options.py for independent crossing-height analysis and
vertical max-min selection; optimization_path.py supplies geometry primitives,
candidate selection, MILP, and Landing_Area evaluation. It does NOT reconstruct
Outer routes using the Inner RSR helper. Outer paths are new CCW LSL curves.

An IF is the unique intersection of the complete LSL and the supplied inner
ring polyline, not a separate clockwise-ring join. Thus no extra turn or
heading discontinuity is introduced at IF. Source geometry heights are used
for Outer/IF, FAF `alt` for endpoints. User authorized climbing after IF when
FAF is above IF; phase reversals are reported, not silently treated as descent.
Absolute segment gradients, bank <=30 deg, the 180 deg shape screen, one IF
crossing, and confinement to the fitted Outer ring are retained.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, asdict, replace
import hashlib
import json
import math
from pathlib import Path
import time

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import LineString, Point, Polygon

import optimization_path as core
import optimization_new_options as new

BASE = Path(__file__).resolve().parents[1]


@dataclass
class Context:
    config: object
    traffic: object
    nodes: np.ndarray
    headings: np.ndarray
    center: np.ndarray
    outer_radius: float
    outer_z: float
    if_line: object
    if_polygon: object
    if_radius: float
    if_z: float
    fafs: list
    areas: list
    inputs: list
    audit: dict
    if_speed_kt: float = 70.


def fit_circle(xy):
    xy = np.unique(np.asarray(xy,dtype=float),axis=0)
    origin = xy.mean(axis=0)
    local = xy-origin
    fit,*_ = np.linalg.lstsq(np.column_stack((2*local,np.ones(len(local)))),
                            np.sum(local*local,axis=1),rcond=None)
    center = origin+fit[:2]
    radius = float(math.sqrt(fit[2]+np.dot(fit[:2],fit[:2])))
    residual = float(np.abs(np.linalg.norm(xy-center,axis=1)-radius).max())
    if residual>.05:
        raise ValueError(f'Ring points are not circular: max residual {residual:g}m.')
    return center,radius,residual


def prepare(args):
    qgis = args.qgis.resolve()
    outer_path,if_path = args.outer_iaf.resolve(),args.if_ring.resolve()
    outer = core.load(outer_path)
    points = [core.single_point(g) for g in outer.geometry]
    if not all(p.has_z for p in points):
        raise ValueError('Outer IAF points require geometry Z.')
    xyz = np.array([p.coords[0][:3] for p in points],dtype=float)
    if not np.isfinite(xyz).all() or np.ptp(xyz[:,2])>1e-7:
        raise ValueError('Outer ring must have finite, uniform geometry Z.')
    center,radius,residual = fit_circle(xyz[:,:2])
    outer_z = float(xyz[0,2])
    inner = core.load(if_path)
    if len(inner)!=1 or inner.geometry.iloc[0].geom_type!='LineString':
        raise ValueError('IF ring must be one closed LineString Z.')
    geometry = inner.geometry.iloc[0]
    if not geometry.has_z:
        raise ValueError('IF ring geometry must have Z.')
    inner_xyz = np.asarray(geometry.coords,dtype=float)
    if not np.isfinite(inner_xyz).all() or np.ptp(inner_xyz[:,2])>1e-7:
        raise ValueError('IF ring must have finite, uniform geometry Z.')
    if_line = LineString(inner_xyz[:,:2])
    if not if_line.is_ring:
        raise ValueError('IF ring must be closed and simple.')
    if_center,if_radius,if_residual = fit_circle(inner_xyz[:,:2])
    if np.linalg.norm(if_center-center)>.05 or if_radius>=radius:
        raise ValueError('IF and Outer rings must be concentric, with IF inside Outer.')
    if_z = float(inner_xyz[0,2])
    config = core.Config(iaf_count=args.iaf_count,ring_z=outer_z,start_kt=90.,end_kt=60.)
    traffic = core.TrafficConfig(multistarts=args.multistarts,candidate_limit=args.candidate_limit)
    _,_,fafs,areas,base_audit,inputs = core.prepare(qgis,config)
    polygon = Polygon(inner_xyz[:,:2])
    if any(not polygon.contains(Point(*f['xyz'][:2])) for f in fafs):
        raise ValueError('All FAF points must be strictly inside the IF ring for a one-crossing model.')
    theta0 = math.atan2(xyz[0,1]-center[1],xyz[0,0]-center[0])
    theta = theta0+np.arange(args.iaf_count)*core.TAU/args.iaf_count
    nodes = center+radius*np.column_stack((np.cos(theta),np.sin(theta)))
    headings = theta+math.pi/2  # Outer CCW tangent
    audit = dict(base_audit,iaf_source_count=len(outer),iaf_unique_count=len(np.unique(xyz[:,:2],axis=0)),
                 ring_center=center.tolist(),ring_radius_m=radius,ring_fit_residual_m=residual,
                 generated_iaf_count=len(nodes),outer_source_points=len(outer),outer_center=center.tolist(),
                 outer_radius_m=radius,outer_geometry_z=outer_z,outer_fit_residual_m=residual,
                 if_radius_m=if_radius,if_geometry_z=if_z,if_fit_residual_m=if_residual,
                 if_altitude_attribute=inner['altitude'].tolist() if 'altitude' in inner else None,
                 actual_if_selection='unique LSL intersection with source IF ring; continuous XY',
                 outer_direction='CCW / LSL',if_direction='continuous route tangent; not Inner CW merge',
                 faf_above_if=[int(f['faf_id']) for f in fafs if f['xyz'][2]>if_z],
                 ascent_after_if_authorized=True,if_reference_points=args.iaf_count)
    return Context(config,traffic,nodes,headings,center,radius,outer_z,if_line,polygon,if_radius,if_z,
                   fafs,areas,sorted(set(inputs+[outer_path,if_path])),audit)


def mirror(xy):
    return np.asarray(xy)*np.array([1.,-1.])


def lsl_lengths(start,h0,end,h1,radius):
    return core.rsr(mirror(start),-h0,mirror(end),-h1,radius)


def sample_lsl(start,h0,end,start_z,radius,lengths,config):
    target = np.array([end[0],-end[1],end[2]])
    xyz,s = core.sample_rsr(mirror(start),-h0,target,radius,lengths,replace(config,ring_z=start_z))
    xyz[:,1] *= -1
    return xyz,s


def endpoints(operation,ni,fi,ctx):
    faf = ctx.fafs[fi]
    if operation=='STAR':
        return ctx.nodes[ni],ctx.headings[ni],faf['xyz'],faf['heading'],ctx.outer_z
    return faf['xyz'][:2],faf['heading']+math.pi,np.r_[ctx.nodes[ni],ctx.outer_z],ctx.headings[ni],float(faf['xyz'][2])


def insert_if(xyz,s,point,s_if):
    nearest = int(np.argmin(np.abs(s-s_if)))
    if abs(s[nearest]-s_if)<1e-7:
        xyz = xyz.copy()
        xyz[nearest,:2] = point
        return xyz,s.copy(),nearest
    index = int(np.searchsorted(s,s_if))
    return np.insert(xyz,index,[point[0],point[1],0.],axis=0),np.insert(s,index,s_if),index


def phase_z(xyz,s,index,z0,z_if,z1):
    if not 0<index<len(xyz)-1:
        raise ValueError('IF cannot coincide with a route endpoint.')
    xyz = xyz.copy()
    xyz[:index+1,2] = z0+(z_if-z0)*s[:index+1]/s[index]
    xyz[index:,2] = z_if+(z1-z_if)*(s[index:]-s[index])/(s[-1]-s[index])
    xyz[0,2],xyz[index,2],xyz[-1,2] = z0,z_if,z1
    return xyz


def find_if(xyz,s,ctx):
    line = LineString(xyz[:,:2])
    hit = line.intersection(ctx.if_line)
    if hit.geom_type=='Point':
        points = [np.array(hit.coords[0])]
    elif hit.geom_type=='MultiPoint':
        points = [np.array(p.coords[0]) for p in hit.geoms]
    else:
        return None
    unique = []
    for p in points:
        if not any(np.linalg.norm(p-q)<1e-6 for q in unique):
            unique.append(p)
    if len(unique)!=1:
        return None
    point = unique[0]
    along = float(line.project(Point(*point)))
    chord = np.r_[0.,np.cumsum(np.linalg.norm(np.diff(xyz[:,:2],axis=0),axis=1))]
    i = min(max(int(np.searchsorted(chord,along,side='right'))-1,0),len(xyz)-2)
    ratio = np.clip((along-chord[i])/(chord[i+1]-chord[i]),0,1)
    s_if = float(s[i]+ratio*(s[i+1]-s[i]))
    return point,s_if


def outer_xyz(row,ctx):
    ni,fi = int(row.iaf_id)-1,int(row.faf_index)
    start,h0,end,_,z0 = endpoints(row.operation,ni,fi,ctx)
    xyz,s = sample_lsl(start,h0,end,z0,float(row.radius_m),
        np.array([row.arc1_m,row.straight_m,row.arc2_m]),ctx.config)
    xyz,s,index = insert_if(xyz,s,np.array([row.if_x,row.if_y]),float(row.if_s_m))
    return phase_z(xyz,s,index,z0,ctx.if_z,float(end[2])),s,index


def z_mode(z0,zm,z1):
    def word(a,b):
        return 'CLIMB' if b>a+1e-8 else ('DESCEND' if b<a-1e-8 else 'LEVEL')
    return word(z0,zm)+'_'+word(zm,z1)


def generate(ctx,out,hashes):
    config,traffic = ctx.config,ctx.traffic
    speed = config.start_kt*core.KNOT
    minimum_radius = speed**2/(core.G*math.tan(math.radians(config.max_bank)))
    radii = [minimum_radius*1.000001,500.,750.,1000.]
    signature = json.loads(json.dumps(dict(version=1,inputs=hashes,config=asdict(config),
        max_turn_deg=traffic.max_total_turn_deg,max_arc_deg=traffic.max_arc_deg,radii=radii,
        outer_z=ctx.outer_z,if_z=ctx.if_z,if_speed_kt=ctx.if_speed_kt)))
    meta,csv = out/'candidate_run.json',out/'outer_candidate_metrics.csv'
    if meta.exists() and csv.exists():
        old = json.loads(meta.read_text(encoding='utf-8'))
        if old['signature']==signature:
            print('Reusing validated Outer candidate inventory.',flush=True)
            return pd.read_csv(csv),old['generation']
    limit = min(math.tan(math.radians(config.max_slope_deg)),config.max_ft_nm*.3048/1852)
    rows,rejected = [],dict(degenerate=0,turn=0,outer_boundary=0,if_crossings=0,gradient=0,self_crossing=0)
    started = time.monotonic()
    for operation in ('STAR','SID'):
        for fi,faf in enumerate(ctx.fafs):
            for ni in range(len(ctx.nodes)):
                start,h0,end,h1,z0 = endpoints(operation,ni,fi,ctx)
                for radius in radii:
                    lengths = lsl_lengths(start,h0,end[:2],h1,radius)
                    if lengths is None:
                        rejected['degenerate']+=1
                        continue
                    turn = math.degrees((lengths[0]+lengths[2])/radius)
                    if turn>traffic.max_total_turn_deg+1e-8 or math.degrees(max(lengths[0],lengths[2])/radius)>traffic.max_arc_deg+1e-8:
                        rejected['turn']+=1
                        continue
                    max_rad = core.maximum_ring_radius(mirror(start),-h0,radius,lengths,mirror(ctx.center))
                    if max_rad>ctx.outer_radius+1e-7:
                        rejected['outer_boundary']+=1
                        continue
                    xyz,s = sample_lsl(start,h0,end,z0,radius,lengths,config)
                    if not LineString(xyz[:,:2]).is_simple:
                        rejected['self_crossing']+=1
                        continue
                    found = find_if(xyz,s,ctx)
                    if found is None:
                        rejected['if_crossings']+=1
                        continue
                    point,s_if = found
                    xyz,s,index = insert_if(xyz,s,point,s_if)
                    xyz = phase_z(xyz,s,index,z0,ctx.if_z,float(end[2]))
                    gradients = np.abs(np.diff(xyz[:,2]))/np.linalg.norm(np.diff(xyz[:,:2],axis=0),axis=1)
                    if gradients.max()>limit+1e-10:
                        rejected['gradient']+=1
                        continue
                    count,overlap,gaps,_ = core.landing_metrics(xyz,ctx.areas)
                    bank = math.degrees(math.atan(speed**2/(core.G*radius)))
                    length = float(s[-1])
                    base = config.length_weight*length/1000+config.bank_weight*bank/config.max_bank
                    score = (base+config.crossing_weight*count+config.overlap_weight*overlap
                             +config.clearance_weight*sum(math.exp(-g/config.clearance_scale_m) for g in gaps))
                    rows.append(dict(cand_id=len(rows)+1,operation=operation,pair=faf['pair'],faf_index=fi,
                        faf_id=faf['faf_id'],iaf_id=ni+1,radius_m=radius,length_m=length,
                        arc1_m=float(lengths[0]),straight_m=float(lengths[1]),arc2_m=float(lengths[2]),
                        turn_deg=turn,bank_deg=bank,slope_deg=math.degrees(math.atan(float(gradients.max()))),
                        ft_nm=float(gradients.max()*1852/.3048),if_x=float(point[0]),if_y=float(point[1]),
                        if_z=ctx.if_z,if_s_m=s_if,if_vert=index,leg1_m=s_if,leg2_m=length-s_if,
                        leg1_deg=math.degrees(math.atan((ctx.if_z-z0)/s_if)),
                        leg2_deg=math.degrees(math.atan((float(end[2])-ctx.if_z)/(length-s_if))),
                        zmode=z_mode(z0,ctx.if_z,float(end[2])),cross_n=count,overlap_m=overlap,
                        min_gap_m=min(gaps) if gaps else np.nan,base_cost=base,score=score,
                        max_rad_m=max_rad,inside=True,shape_ok=True))
            print(f'Outer {operation} FAF {faf["faf_id"]}: accepted={len(rows)}, {time.monotonic()-started:.1f}s',flush=True)
    columns = ['cand_id','operation','pair','faf_index','faf_id','iaf_id','radius_m','length_m',
        'arc1_m','straight_m','arc2_m','turn_deg','bank_deg','slope_deg','ft_nm','if_x','if_y','if_z','if_s_m',
        'if_vert','leg1_m','leg2_m','leg1_deg','leg2_deg','zmode','cross_n','overlap_m','min_gap_m',
        'base_cost','score','max_rad_m','inside','shape_ok']
    frame = pd.DataFrame(rows,columns=columns)
    frame.to_csv(csv,index=False,encoding='utf-8-sig')
    generation = dict(total_combinations=2*len(ctx.nodes)*len(ctx.fafs)*len(radii),feasible=len(frame),
        rejected=rejected,radii_m=radii,minimum_radius_m=minimum_radius,seconds=time.monotonic()-started,
        if_selection='one actual crossing of complete LSL with source IF ring')
    meta.write_text(json.dumps(dict(signature=signature,generation=generation),ensure_ascii=False,indent=2),encoding='utf-8')
    return frame,generation


def candidate_pool(frame,ctx):
    """Retain all FAF identities even when an operation has no feasible path."""
    if frame.empty:
        return frame.copy(),{p:[] for p in range(1,ctx.audit['pair_count']+1)},dict(
            strict_pairs=[],penalized_pairs=[],missing_pairs=list(range(1,ctx.audit['pair_count']+1)),
            eligible_candidates=0,shortlisted_candidates=0)
    # Core's helper expects every FAF identity to appear in the inventory.
    # Missing families get disabled placeholders (never eligible for selection).
    extras = []
    for fi,faf in enumerate(ctx.fafs):
        if not frame.faf_id.eq(faf['faf_id']).any():
            extra = {c:np.nan for c in frame.columns}
            extra.update(cand_id=-faf['faf_id'],operation='STAR',pair=faf['pair'],faf_id=faf['faf_id'],
                         faf_index=fi,iaf_id=1,shape_ok=False,inside=False,cross_n=0,score=math.inf)
            extras.append(extra)
    if extras:
        data = pd.concat([frame,pd.DataFrame(extras)],ignore_index=True)
    else:
        data = frame
    pool,options,policy = core.traffic_pool(data,'inside_only',ctx.config,ctx.traffic,ctx.audit)
    # Some whole pairs can have zero candidates and be absent from data.
    for p in range(1,ctx.audit['pair_count']+1):
        options.setdefault(p,[])
    policy['missing_pairs'] = [int(p) for p,v in options.items() if not v]
    return pool,options,policy


def write_coverage(ctx,frame,out):
    limit=min(math.tan(math.radians(ctx.config.max_slope_deg)),ctx.config.max_ft_nm*.3048/1852)
    rows=[]
    for faf in ctx.fafs:
        data=frame[frame.faf_id==faf['faf_id']]
        rows.append(dict(pair=faf['pair'],faf_id=faf['faf_id'],faf_alt=float(faf['xyz'][2]),
            if_alt=ctx.if_z,outer_alt=ctx.outer_z,approach_candidates=int((data.operation=='STAR').sum()),
            departure_candidates=int((data.operation=='SID').sum()),
            min_outer_if_length_m=abs(ctx.outer_z-ctx.if_z)/limit,
            min_if_faf_length_m=abs(ctx.if_z-faf['xyz'][2])/limit,
            approach_if_faf_climbs=bool(faf['xyz'][2]>ctx.if_z)))
    pd.DataFrame(rows).to_csv(out/'candidate_coverage.csv',index=False,encoding='utf-8-sig')


def xy_matrices(pool,models,out,key,traffic):
    cache = out/'outer_xy_matrices.npz'
    if cache.exists():
        with np.load(cache,allow_pickle=False) as old:
            if old['key'].item()==key:
                d,c = old['distance'].copy(),old['crossing'].copy()
                g = np.maximum(0.,1-d/traffic.spacing_target_m)**2
                np.fill_diagonal(g,0.)
                return d,c,g
    n = len(pool)
    lines = np.empty(n,dtype=object)
    lines[:] = [m['line'].simplify(1e-8,preserve_topology=False) for m in models]
    d,c = np.zeros((n,n)),np.zeros((n,n),dtype=bool)
    def block(first):
        result=[]
        for i in range(first,min(first+32,n)):
            hits = shapely.intersects(lines[i],lines[i+1:])
            values = np.zeros(len(hits))
            if (~hits).any():
                values[~hits] = shapely.distance(lines[i],lines[i+1:][~hits])
            result.append((i,values,hits))
        return result
    completed=0
    with ThreadPoolExecutor(max_workers=4) as executor:
        tasks = [executor.submit(block,i) for i in range(0,n,32)]
        for task in as_completed(tasks):
            values=task.result()
            for i,dd,cc in values:
                d[i,i+1:],c[i,i+1:]=dd,cc
            completed+=len(values)
            if completed%256<32:
                print(f'Outer XY matrix {completed}/{n}',flush=True)
    d+=d.T
    c|=c.T
    g = np.maximum(0.,1-d/traffic.spacing_target_m)**2
    np.fill_diagonal(g,0.)
    np.savez_compressed(cache,key=key,distance=d,crossing=c)
    return d,c,g


def samples_output(ctx,out):
    gpd.GeoDataFrame(dict(iaf_id=np.arange(1,len(ctx.nodes)+1),alt=ctx.outer_z),
        geometry=[Point(x,y,ctx.outer_z) for x,y in ctx.nodes],crs=5186).to_file(
            out/'Outer_IAF_candidates_z.shp',engine='pyogrio',encoding='UTF-8')
    theta = np.arange(len(ctx.nodes))*core.TAU/len(ctx.nodes)
    xy = ctx.center+ctx.if_radius*np.column_stack((np.cos(theta),np.sin(theta)))
    gpd.GeoDataFrame(dict(ref_id=np.arange(1,len(xy)+1),alt=ctx.if_z),
        geometry=[Point(x,y,ctx.if_z) for x,y in xy],crs=5186).to_file(
            out/'IF_reference_points_z.shp',engine='pyogrio',encoding='UTF-8')


def export_case(pool,ids,out,ctx):
    out.mkdir(parents=True,exist_ok=True)
    records,geometries,points,point_attrs,vertices,legs,landings = [],[],[],[],[],[],[]
    limit = min(math.tan(math.radians(ctx.config.max_slope_deg)),ctx.config.max_ft_nm*.3048/1852)
    for index in ids:
        row=pool.loc[index]
        xyz,s,ii=outer_xyz(row,ctx)
        ni,fi=int(row.iaf_id)-1,int(row.faf_index)
        start,h0,end,h1,z0=endpoints(row.operation,ni,fi,ctx)
        assert np.allclose(xyz[0],np.r_[start,z0],rtol=0,atol=1e-7)
        assert np.allclose(xyz[-1],end,rtol=0,atol=1e-7)
        assert abs((h0+(row.arc1_m+row.arc2_m)/row.radius_m-h1+math.pi)%core.TAU-math.pi)<1e-8
        assert ctx.if_line.distance(Point(*xyz[ii,:2]))<1e-6 and xyz[ii,2]==ctx.if_z
        assert LineString(xyz[:,:2]).is_simple
        assert row.turn_deg<=ctx.traffic.max_total_turn_deg+1e-7
        assert row.max_rad_m<=ctx.outer_radius+1e-7 and row.bank_deg<=ctx.config.max_bank+1e-7
        assert (np.abs(np.diff(xyz[:,2]))/np.linalg.norm(np.diff(xyz[:,:2],axis=0),axis=1)).max()<=limit+1e-8
        assert find_if(xyz,s,ctx) is not None
        app=row.operation=='STAR'
        rec={k:row[k] for k in ['cand_id','operation','pair','faf_id','iaf_id','radius_m','length_m',
            'turn_deg','bank_deg','slope_deg','cross_n','overlap_m','leg1_m','leg2_m','leg1_deg','leg2_deg','zmode']}
        rec.update(path_type='LSL',if_x=xyz[ii,0],if_y=xyz[ii,1],if_z=ctx.if_z,
                   z_start=xyz[0,2],z_end=xyz[-1,2],max_rad_m=row.max_rad_m,
                   v_start=90. if app else 60.,v_if=ctx.if_speed_kt,v_end=60. if app else 90.)
        records.append(rec)
        geometries.append(LineString(xyz))
        positions=[('IAF',0 if app else -1),('IF',ii),('FAF',-1 if app else 0)]
        for role,vi in positions:
            p=xyz[vi]
            point_attrs.append(dict(route_id=int(row.cand_id),operation=row.operation,pair=int(row.pair),
                iaf_id=int(row.iaf_id),faf_id=int(row.faf_id),role=role,alt=float(p[2]),
                speed_kt=90. if role=='IAF' else (ctx.if_speed_kt if role=='IF' else 60.)))
            points.append(Point(*p))
        phase_names=['IAF_IF','IF_FAF'] if app else ['FAF_IF','IF_IAF']
        for phase,(lo,hi) in enumerate([(0,ii),(ii,len(xyz)-1)]):
            legs.append(dict(**rec,phase=phase_names[phase],geometry=LineString(xyz[lo:hi+1])))
        elapsed=0.
        speeds=(90.,ctx.if_speed_kt,60.) if app else (60.,ctx.if_speed_kt,90.)
        for phase,(lo,hi) in enumerate([(0,ii),(ii,len(xyz)-1)]):
            distance=s[lo:hi+1]-s[lo]
            d3=distance*math.hypot(1.,(xyz[hi,2]-xyz[lo,2])/distance[-1])
            v0,v1=speeds[phase]*core.KNOT,speeds[phase+1]*core.KNOT
            velocity=np.sqrt(v0*v0+(v1*v1-v0*v0)*distance/distance[-1])
            times=2*d3/(v0+velocity)
            for k in range(0 if phase==0 else 1,len(distance)):
                vi=lo+k
                vertices.append(dict(route_id=int(row.cand_id),operation=row.operation,pair=int(row.pair),
                    vertex=vi,x=xyz[vi,0],y=xyz[vi,1],z=xyz[vi,2],s_m=s[vi],
                    speed_kt=velocity[k]/core.KNOT,time_s=elapsed+times[k],is_if=vi==ii))
            elapsed+=times[-1]
        count,overlap,gaps,detail=core.landing_metrics(xyz,ctx.areas,details=True)
        assert count==row.cross_n and abs(overlap-row.overlap_m)<1e-5
        landings.extend(dict(operation=row.operation,pair=int(row.pair),**r) for r in detail)
    columns=['cand_id','operation','pair','faf_id','iaf_id','radius_m','length_m','turn_deg','bank_deg',
        'slope_deg','cross_n','overlap_m','leg1_m','leg2_m','leg1_deg','leg2_deg','zmode','path_type',
        'if_x','if_y','if_z','z_start','z_end','max_rad_m','v_start','v_if','v_end']
    result=gpd.GeoDataFrame(pd.DataFrame(records,columns=columns),geometry=geometries,crs=5186)
    for filename,data in [('Approach_z.shp',result[result.operation=='STAR']),
                          ('Departure_z.shp',result[result.operation=='SID']),('combined_routes_z.shp',result)]:
        path=out/filename
        data.to_file(path,engine='pyogrio',encoding='UTF-8',geometry_type='LineString Z')
        saved=gpd.read_file(path)
        assert len(saved)==len(data) and saved.has_z.all() and saved.crs.to_epsg()==5186
        for a,b in zip(saved.geometry,data.geometry):
            assert np.allclose(np.asarray(a.coords),np.asarray(b.coords),rtol=0,atol=1e-7)
    attrs=pd.DataFrame(point_attrs,columns=['route_id','operation','pair','iaf_id','faf_id','role','alt','speed_kt'])
    point_gdf=gpd.GeoDataFrame(attrs,geometry=points,crs=5186)
    for name,subset in [('connected_points_z',point_gdf)]+[(f'{role}_points_z',point_gdf[point_gdf.role==role]) for role in ['IAF','IF','FAF']]:
        subset.to_file(out/f'{name}.shp',engine='pyogrio',encoding='UTF-8',geometry_type='Point Z')
    if legs:
        leg_gdf=gpd.GeoDataFrame(legs,crs=5186)
        for op,phase in [('STAR','IAF_IF'),('STAR','IF_FAF'),('SID','FAF_IF'),('SID','IF_IAF')]:
            leg_gdf[(leg_gdf.operation==op)&(leg_gdf.phase==phase)].to_file(
                out/f'{"Approach" if op=="STAR" else "Departure"}_{phase}_z.shp',engine='pyogrio',encoding='UTF-8')
    result.drop(columns='geometry').to_csv(out/'selected_routes.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(vertices,columns=['route_id','operation','pair','vertex','x','y','z','s_m','speed_kt','time_s','is_if']).to_csv(
        out/'vertices.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(landings,columns=['operation','pair','area_pair','area_fid','overlap_m','min_abs_z','min_signed_z']).to_csv(
        out/'landing_crossings.csv',index=False,encoding='utf-8-sig')
    pair_rows=[]
    for i in range(len(result)):
        for j in range(i+1,len(result)):
            a,b=result.iloc[i],result.iloc[j]
            la,lb=LineString(np.asarray(a.geometry.coords)[:,:2]),LineString(np.asarray(b.geometry.coords)[:,:2])
            pair_rows.append(dict(route_a=int(a.cand_id),operation_a=a.operation,pair_a=int(a.pair),
                route_b=int(b.cand_id),operation_b=b.operation,pair_b=int(b.pair),min_xy_m=la.distance(lb),
                intersects=la.intersects(lb)))
    pd.DataFrame(pair_rows,columns=['route_a','operation_a','pair_a','route_b','operation_b','pair_b','min_xy_m','intersects']).to_csv(
        out/'route_pair_spacing.csv',index=False,encoding='utf-8-sig')
    crossing=new.crossing_audit(result,out)
    for pair,rows in result.groupby('pair'):
        assert len(rows)==2 and set(rows.operation)=={'STAR','SID'} and rows.faf_id.nunique()==2
    info=dict(selected_pairs=int(result.pair.nunique()),approach_count=int((result.operation=='STAR').sum()),
        departure_count=int((result.operation=='SID').sum()),missing_pairs=sorted(set(range(1,ctx.audit['pair_count']+1))-set(result.pair.astype(int))),
        length_m=float(result.length_m.sum()),landing_area_crossings=int(result.cross_n.sum()),
        minimum_xy_m=min((r['min_xy_m'] for r in pair_rows),default=None),**crossing)
    (out/'summary.json').write_text(json.dumps(info,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    plot(result,ctx,out)
    return result,info


def plot(result,ctx,out):
    colors=core.plt.get_cmap('tab10')
    theta=np.linspace(0,core.TAU,2001)
    for mode in ['combined','STAR','SID']:
        fig,ax=core.plt.subplots(figsize=(13,12))
        for area in ctx.areas:
            xy=np.asarray(area['geometry'].exterior.coords)
            ax.fill(*xy.T,color='tomato',alpha=.10)
        for radius,label in [(ctx.outer_radius,f'Outer CCW / {ctx.outer_z:g}m'),(ctx.if_radius,f'IF crossing ring / {ctx.if_z:g}m')]:
            xy=ctx.center+radius*np.column_stack((np.cos(theta),np.sin(theta)))
            ax.plot(*xy.T,'--',lw=1,color='gray',label=label)
        for _,row in result.iterrows():
            if mode!='combined' and row.operation!=mode:
                continue
            xyz=np.asarray(row.geometry.coords)
            color=colors((int(row.pair)-1)%10)
            ax.plot(*xyz[:,:2].T,ls='-' if row.operation=='STAR' else '--',color=color,
                label=f'{"AP" if row.operation=="STAR" else "DP"} P{int(row.pair)} F{int(row.faf_id)}')
            iaf=xyz[0] if row.operation=='STAR' else xyz[-1]
            faf=xyz[-1] if row.operation=='STAR' else xyz[0]
            ax.scatter(*iaf[:2],marker='o',s=25,color=color)
            ax.scatter(*faf[:2],marker='*' if row.operation=='STAR' else '^',s=65,color=color)
            ax.scatter(row.if_x,row.if_y,marker='D',s=25,color=color)
            k=min(len(xyz)-2,len(xyz)//2)
            ax.annotate('',xy=xyz[k+1,:2],xytext=xyz[k,:2],arrowprops=dict(arrowstyle='->',color=color))
        ax.set(aspect='equal',title=f'Outer IAF / IF / FAF: {mode}, CCW LSL',xlabel='EPSG:5186 X (m)',ylabel='EPSG:5186 Y (m)')
        ax.ticklabel_format(style='plain',useOffset=False)
        ax.legend(fontsize=7,ncol=2)
        ax.grid(alpha=.2)
        fig.tight_layout()
        fig.savefig(out/f'map_{"combined" if mode=="combined" else ("approach" if mode=="STAR" else "departure")}.png',dpi=160)
        core.plt.close(fig)
    fig,axes=core.plt.subplots(4,2,figsize=(14,13))
    for pair,rows in result.groupby('pair'):
        ax=axes.flat[int(pair)-1]
        for _,row in rows.iterrows():
            xyz=np.asarray(row.geometry.coords)
            ss=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(xyz[:,:2],axis=0),axis=1))]
            ax.plot(ss,xyz[:,2],ls='-' if row.operation=='STAR' else '--',label=row.operation)
            i=np.argmin(np.linalg.norm(xyz[:,:2]-np.array([row.if_x,row.if_y]),axis=1))
            ax.scatter(ss[i],ctx.if_z,marker='D',s=20)
        ax.set(title=f'Pair {int(pair)} / IF {ctx.if_z:g}m',xlabel='Horizontal distance (m)',ylabel='Z (m)')
        ax.legend(fontsize=8)
        ax.grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(out/'profiles.png',dpi=150)
    core.plt.close(fig)
    fig=core.plt.figure(figsize=(13,10))
    ax=fig.add_subplot(111,projection='3d')
    for _,row in result.iterrows():
        xyz=np.asarray(row.geometry.coords)
        ax.plot(*xyz.T,ls='-' if row.operation=='STAR' else '--',color=colors((int(row.pair)-1)%10))
        ax.scatter(row.if_x,row.if_y,ctx.if_z,marker='D',s=20)
    ax.set(title='Outer LSL paths / IF passage',xlabel='X (m)',ylabel='Y (m)',zlabel='Z (m)')
    fig.savefig(out/'view3d.png',dpi=160)
    core.plt.close(fig)


def report(out,summary):
    sections=[]
    labels=dict(no_xy_crossing='무교차 우선안',strict_no_xy_crossing='Landing_Area도 관통 금지',
                vertical_clearance_alternative='교차점 수직 이격 대안')
    for name,info in summary['results'].items():
        gap=info['minimum_crossing_z_m']
        text='교차점 없음' if gap is None else f'{gap:.3f}m'
        sections.append(f'''<h2>{labels[name]}</h2><p>FAF {info['selected_pairs']}쌍,
Approach {info['approach_count']}개 / Departure {info['departure_count']}개.
미연결 쌍 {info['missing_pairs']}. XY 교차 {info['planar_crossing_pairs']}쌍,
교점 최소 수직 간격 {text}. Landing_Area 관통 {info['landing_area_crossings']}건.</p>
<p><a href="{name}/Approach_z.shp">Approach SHP</a> · <a href="{name}/Departure_z.shp">Departure SHP</a> ·
<a href="{name}/IAF_points_z.shp">IAF</a> · <a href="{name}/IF_points_z.shp">IF</a> ·
<a href="{name}/FAF_points_z.shp">FAF</a> · <a href="{name}/selected_routes.csv">경로·구간별 구배</a></p>
<img src="{name}/map_combined.png"><details><summary>고도·속도 경로 / 3D</summary>
<img src="{name}/profiles.png"><img src="{name}/view3d.png"></details>''')
        if gap is not None:
            sections.append(f'<img src="{name}/crossing_z_spacing.png"><img src="{name}/crossing_vertical_links.png">')
    audit=summary['input_audit']
    html=f'''<!doctype html><html lang="ko"><meta charset="utf-8"><title>Outer IAF–IF–FAF</title>
<style>body{{font-family:Malgun Gothic,sans-serif;max-width:1250px;margin:32px auto;line-height:1.7}}img{{max-width:100%}}.note{{background:#fff1d7;padding:16px}}</style>
<h1>Outer IAF → IF → FAF / 반대편 Departure</h1>
<p>외곽 {audit['outer_geometry_z']:g}m, IF 기하 {audit['if_geometry_z']:g}m, FAF alt 사용.
Approach 속도 90→70→60kt, Departure 60→70→90kt. 두 방향 모두 새 CCW LSL.</p>
<p class="note">사용자 지시에 따라 IF보다 높은 FAF도 연결했습니다. 해당 Approach는 IF 이후 상승하고,
Departure는 FAF→IF에서 하강할 수 있습니다. 이를 계속 하강/계속 상승하는 절차라고 표현하지 않습니다.
실제 운항 승인·안전성 검증이 아닌 기하 후보 분석입니다.</p>
<p>전체 경로를 하나의 LSL로 생성한 뒤 원본 IF 링과의 유일한 교점을 IF로 잡습니다.
IF에서 내곽 시계방향 접선에 합류하지 않으며, 동일 경로의 연속 접선으로 통과합니다.
외곽 IAF 후보 {summary['config']['iaf_count']}개, 참고 IF 점 {audit['if_reference_points']}개를 저장하지만,
실제 선정 IF는 이 격자 중 하나가 아니라 원본 링 위의 연속 교점입니다.</p>
<p>Outer 경계 밖 우회·여러 번 IF 링 통과·총 선회 180° 초과·자기교차는 배제합니다.
각 구간의 절대 구배는 2.36°와 250ft/NM 중 엄격한 값, 뱅크 ≤30°를 유지합니다.
90kt 기반 최소반경 {summary['generation']['minimum_radius_m']:.2f}m 이상을 탐색합니다.</p>
<p>IF에서 Z를 250m에 맞춘 구간별 선형 보간을 사용하므로 종단 경사가 바뀔 수 있습니다.
수평 접선은 연속이지만 실제 롤·종단 전이·기체 성능·풍향·시간 분리는 포함하지 않습니다.
70m 고정 수직 기준이나 기존 Inner 선정 항로와의 공동 충돌 제약을 추가하지 않았습니다.</p>
<p><a href="candidate_coverage.csv">FAF별 후보 수와 필요한 최소 구간 길이</a>를 함께 저장했습니다.
부족한 쌍은 임의의 원형 우회·구배 완화로 채우지 않습니다. 미연결은 현재 격자·반경·형태 필터 안의 결과이며
모든 가능한 연속 항로가 불가능하다는 증명은 아닙니다.</p>
{''.join(sections)}
<p>Landing_Area는 가능한 FAF 쌍에서 관통 금지, 불가능한 쌍만 페널티로 처리합니다.
SHP는 저장 후 XYZ·IF 통과·양쪽 구배·LSL 방향을 검증합니다. 최적성은 선별 후보 집합의 범위입니다.</p>
<h2>재실행</h2><pre>python -B 수정본_20261103/코드/optimization_new_options_add_outer.py
python -B 수정본_20261103/코드/optimization_new_options_add_outer.py --self-check</pre></html>'''
    (out/'report.html').write_text(html,encoding='utf-8')
    (out/'README.md').write_text('''# Outer IAF–IF–FAF

report.html에 모든 결과 링크가 있습니다. Approach/Departure는 별도 Z SHP입니다.
IF_points_z.shp는 실제 선정된 연속 교점이며, IF_reference_points_z.shp는 참고 격자입니다.
단계별 SHP와 vertices.csv의 is_if·speed_kt·time_s로 IF 높이와 속도를 확인할 수 있습니다.
IF 250m보다 높은 FAF로의 상승 연결은 사용자 지시에 따른 예비 기하 모델입니다.
기존 Inner 결과와 원본 SHP는 수정하지 않습니다. 입력 DBF 수정 백업은 input_backups에 별도 보존합니다.
''',encoding='utf-8')


def self_check():
    config=core.Config()
    xyz,s=sample_lsl(np.array([0.,0.]),0.,np.array([1.,2.,350.]),350.,1.,np.array([math.pi/2,1.,0.]),config)
    assert np.allclose(xyz[-1],[1,2,350],rtol=0,atol=1e-8)
    a,b=xyz[1,:2]-xyz[0,:2],xyz[2,:2]-xyz[1,:2]
    assert a[0]*b[1]-a[1]*b[0]>0
    rng=np.random.default_rng(90)
    for _ in range(100):
        a,b=rng.normal(size=(2,2))*2000
        h0,h1=rng.uniform(-math.pi,math.pi,2)
        radius=float(rng.uniform(380,1000))
        lens=lsl_lengths(a,h0,b,h1,radius)
        xyz,s=sample_lsl(a,h0,np.r_[b,280],350.,radius,lens,config)
        assert np.allclose(xyz[-1],np.r_[b,280],rtol=0,atol=1e-7)
        assert abs((h0+(lens[0]+lens[2])/radius-h1+math.pi)%core.TAU-math.pi)<1e-9
    xyz=np.array([[0,0,350],[100,0,0],[200,0,300]],dtype=float)
    result=phase_z(xyz,np.array([0.,100.,200.]),1,350,250,300)
    assert np.array_equal(result[:,2],[350,250,300])
    from types import SimpleNamespace
    ring=Point(0,0).buffer(1.,quad_segs=128).exterior
    test_ctx=SimpleNamespace(if_line=LineString(ring.coords))
    xyz=np.array([[-2.,0.,350.],[0.,0.,214.]])
    s=np.array([0.,2.])
    point,s_if=find_if(xyz,s,test_ctx)
    assert np.allclose(point,[-1.,0.],atol=1e-8) and abs(s_if-1.)<1e-8
    split,ss,ii=insert_if(xyz,s,point,s_if)
    split=phase_z(split,ss,ii,350.,250.,214.)
    assert np.allclose(split[ii],[-1.,0.,250.],atol=1e-8)
    new.self_check()
    print('Outer self-check passed: CCW left arc, 100 LSL endpoints/headings, IF piecewise Z with permitted climb.')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--qgis',type=Path,default=BASE/'QGIS')
    parser.add_argument('--outer-iaf',type=Path,default=BASE/'QGIS/IAF_후보군/외곽_꼭짓점.shp')
    parser.add_argument('--if-ring',type=Path,default=BASE/'QGIS/내곽_2000_250.shp')
    parser.add_argument('--output',type=Path,default=BASE/'결과물/Outer')
    parser.add_argument('--iaf-count',type=int,default=1000)
    parser.add_argument('--candidate-limit',type=int,default=36)
    parser.add_argument('--multistarts',type=int,default=80)
    parser.add_argument('--search-seconds',type=float,default=45.)
    parser.add_argument('--z-tolerance',type=float,default=.5)
    parser.add_argument('--self-check',action='store_true')
    args=parser.parse_args()
    if args.self_check:
        self_check()
        return
    if args.iaf_count<2 or args.candidate_limit<1 or args.multistarts<1 or not all(
        math.isfinite(v) and v>0 for v in [args.search_seconds,args.z_tolerance]):
        raise ValueError('Invalid positive search/count settings.')
    out=args.output.resolve()
    out.mkdir(parents=True,exist_ok=True)
    ctx=prepare(args)
    hashes=core.fingerprint(ctx.inputs)
    frame,generation=generate(ctx,out,hashes)
    write_coverage(ctx,frame,out)
    samples_output(ctx,out)
    pool,options,policy=candidate_pool(frame,ctx)
    pool.to_csv(out/'candidate_shortlist.csv',index=False,encoding='utf-8-sig')
    summary=dict(config=asdict(ctx.config),traffic_config=asdict(ctx.traffic),input_audit=ctx.audit,
        input_sha256=hashes,generation=generation,policy=policy,results={},
        if_actual_geometry_z=ctx.if_z,if_climb_authorized=True,hard_vertical_threshold_m=None)
    models=[new.make_route(outer_xyz(row,ctx)[0]) for _,row in pool.iterrows()]
    key=hashlib.sha256(json.dumps(dict(version=1,inputs=hashes,config=asdict(ctx.config),
        candidate_csv_sha256=hashlib.sha256((out/'outer_candidate_metrics.csv').read_bytes()).hexdigest(),
        candidate_ids=pool.cand_id.astype(int).tolist()),sort_keys=True).encode()).hexdigest()
    d,c,g=xy_matrices(pool,models,out,key,ctx.traffic)
    ids,coverage=core.traffic_noncrossing_subset(pool,options,c,args.search_seconds)
    if ids:
        ids,spread=core.traffic_maxmin_spacing(pool,options,c,d,ctx.fafs,ids,seconds=args.search_seconds)
        active=set(pool.loc[ids,'pair'].astype(int))
        ids,_=core.traffic_select(pool,{p:v for p,v in options.items() if p in active},d,c,g,ctx.nodes,
            ctx.traffic,ctx.config,initial_indices=ids,forbid_crossings=True,
            minimum_spacing_m=max(0.,spread.get('achieved_minimum_m',0.)-1e-6))
    else:
        spread=dict(status='no_eligible_pairs')
    _,stats=export_case(pool,ids,out/'no_xy_crossing',ctx)
    assert stats['planar_crossing_pairs']==0
    summary['results']['no_xy_crossing']=dict(**stats,coverage_search=coverage,spacing_search=spread)
    strict_options={p:v for p,v in options.items() if p in policy['strict_pairs']}
    ids,coverage=core.traffic_noncrossing_subset(pool,strict_options,c,args.search_seconds)
    if ids:
        ids,spread=core.traffic_maxmin_spacing(pool,strict_options,c,d,ctx.fafs,ids,seconds=args.search_seconds)
        active=set(pool.loc[ids,'pair'].astype(int))
        ids,_=core.traffic_select(pool,{p:v for p,v in strict_options.items() if p in active},d,c,g,ctx.nodes,
            ctx.traffic,ctx.config,initial_indices=ids,forbid_crossings=True,
            minimum_spacing_m=max(0.,spread.get('achieved_minimum_m',0.)-1e-6))
    else:
        spread=dict(status='no_eligible_pairs')
    _,stats=export_case(pool,ids,out/'strict_no_xy_crossing',ctx)
    assert stats['landing_area_crossings']==0 and stats['planar_crossing_pairs']==0
    summary['results']['strict_no_xy_crossing']=dict(**stats,coverage_search=coverage,spacing_search=spread)
    active_count=sum(bool(v) for v in options.values())
    if active_count:
        z=new.build_vertical_matrix(pool,models,c,out,hashlib.sha256((key+':outer-piecewise-z-v1').encode()).hexdigest())
        seed,_=core.traffic_select(pool,options,d,c,g,ctx.nodes,ctx.traffic,ctx.config)
        ids,search=new.maximize_crossing_height(pool,options,c,z,seed,ctx.config,ctx.fafs,args.search_seconds,args.z_tolerance)
        floor=search.get('floor_m')
        forbidden=c if floor is None else c&(z<floor-1e-6)
        ids,_=core.traffic_select(pool,options,d,c,g,ctx.nodes,ctx.traffic,ctx.config,
            initial_indices=ids,forbid_crossings=True,forbidden_pairs=forbidden)
        _,stats=export_case(pool,ids,out/'vertical_clearance_alternative',ctx)
        assert floor is None or stats['minimum_crossing_z_m'] is None or stats['minimum_crossing_z_m']>=floor-1e-5
        summary['results']['vertical_clearance_alternative']=dict(**stats,vertical_search=search)
    else:
        _,stats=export_case(pool,[],out/'vertical_clearance_alternative',ctx)
        summary['results']['vertical_clearance_alternative']=dict(**stats,status='no_complete_opposite_faf_pair')
    assert core.fingerprint(ctx.inputs)==hashes
    summary['original_inputs_unchanged_during_generation']=True
    (out/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    report(out,summary)
    print(json.dumps({name:{k:v for k,v in info.items() if k not in ['coverage_search','spacing_search','vertical_search']}
        for name,info in summary['results'].items()},ensure_ascii=False,indent=2),flush=True)
    print(f'Complete: {out}',flush=True)


if __name__=='__main__':
    main()
