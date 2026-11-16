"""Editable assemblies, coordinate lifting and bounded map remapping.

Source intervals never change when fragments are reordered. Matrix remapping
allocates counts by fractional overlap at cuts inside a source bin; aligned
cuts and whole-bin permutations preserve the original counts exactly.
"""
from __future__ import annotations

import copy
import math
from pathlib import Path

import numpy as np


class Assembly:
    def __init__(self,fragments,groups):
        self.fragments=copy.deepcopy(fragments)
        self.groups=copy.deepcopy(groups)
        self.undo_stack=[]
        self.redo_stack=[]
        self.saved=None
        self.validate()
        self.mark_saved()

    def data(self):
        return dict(version=1,fragments=copy.deepcopy(self.fragments),groups=copy.deepcopy(self.groups))

    @classmethod
    def restore(cls,data):
        if data.get('version')!=1:
            raise ValueError('Unsupported assembly session')
        return cls(data['fragments'],data['groups'])

    @classmethod
    def read(cls,path,source_chrom='assembly',scale=1,by_name=False):
        if not isinstance(scale,int) or scale<1:
            raise ValueError('Coordinate scale must be a positive integer')
        fragments,groups,offset={},[],0
        for number,line in enumerate(Path(path).read_text(encoding='utf8').splitlines(),1):
            if not line.strip():
                continue
            try:
                if line.startswith('>'):
                    name,key_text,length_text=line[1:].split()
                    key,length=int(key_text),int(length_text)
                    if key in fragments:
                        raise ValueError('Duplicate fragment ID')
                    fragments[key]=dict(name=name,length=length,source_name=name,source_start=0,
                                        map_chrom=name if by_name else source_chrom,
                                        map_start=0 if by_name else offset,scale=scale)
                    offset+=length
                else:
                    groups.append([int(v) for v in line.split()])
            except (ValueError,TypeError) as exc:
                raise ValueError(f'Invalid .assembly line {number}: {exc}') from exc
        return cls(fragments,groups)

    @classmethod
    def from_chromsizes(cls,chromsizes):
        fragments={i:dict(name=name,length=length,source_name=name,source_start=0,
                         map_chrom=name,map_start=0,scale=1) for i,(name,length) in enumerate(chromsizes,1)}
        return cls(fragments,[[i] for i in fragments])

    def validate(self):
        # JSON object keys become strings.
        self.fragments={int(k):v for k,v in self.fragments.items()}
        if not self.fragments or not self.groups or any(not group for group in self.groups):
            raise ValueError('An assembly needs fragments and nonempty scaffold groups')
        seen=[]
        for key,f in self.fragments.items():
            if key<=0 or not isinstance(f['length'],int) or f['length']<=0:
                raise ValueError('Fragment IDs and lengths must be positive integers')
            if not f['name'] or any(c.isspace() for c in f['name']):
                raise ValueError('Fragment names cannot contain whitespace')
            if any(not isinstance(f[k],int) for k in ('source_start','map_start','scale')) or f['source_start']<0 or f['map_start']<0 or f['scale']<1:
                raise ValueError('Invalid source interval or coordinate scale')
            if any(not isinstance(f[k],str) or not f[k] or any(c.isspace() for c in f[k]) for k in ('source_name','map_chrom')):
                raise ValueError('Source names must be nonempty and cannot contain whitespace')
        for group in self.groups:
            for key in group:
                if not isinstance(key,int) or abs(key) not in self.fragments:
                    raise ValueError('Scaffold refers to an unknown fragment')
                seen.append(abs(key))
        if len(set(seen))!=len(seen) or set(seen)!=set(self.fragments):
            raise ValueError('Each fragment must occur exactly once')

    @property
    def dirty(self):
        return self.saved!=self.data()

    def mark_saved(self):
        self.saved=self.data()

    def edit(self,operation):
        before=self.data()
        try:
            operation()
            self.groups=[g for g in self.groups if g]
            self.validate()
        except Exception:
            self.fragments,self.groups=before['fragments'],before['groups']
            raise
        if self.data()!=before:
            self.undo_stack.append(before)
            self.undo_stack=self.undo_stack[-100:]
            self.redo_stack=[]

    def undo(self):
        if self.undo_stack:
            self.redo_stack.append(self.data())
            state=self.undo_stack.pop()
            self.fragments,self.groups=state['fragments'],state['groups']

    def redo(self):
        if self.redo_stack:
            self.undo_stack.append(self.data())
            state=self.redo_stack.pop()
            self.fragments,self.groups=state['fragments'],state['groups']

    def invert(self,group,first,last):
        if not 0<=first<=last<len(self.groups[group]):
            raise ValueError('Invalid fragment selection')
        self.edit(lambda:self.groups[group].__setitem__(slice(first,last+1),[-v for v in reversed(self.groups[group][first:last+1])]))

    def move(self,group,first,last,target,position):
        if not 0<=first<=last<len(self.groups[group]) or not 0<=position<=len(self.groups[target]):
            raise ValueError('Invalid move destination')
        if group==target and first<=position<=last+1:
            return
        def operation():
            selected=self.groups[group][first:last+1]
            del self.groups[group][first:last+1]
            at=position-(last-first+1) if group==target and position>last else position
            self.groups[target][at:at]=selected
        self.edit(operation)

    def split_group(self,group,position):
        if not 0<position<len(self.groups[group]):
            raise ValueError('Select a boundary inside the scaffold')
        def operation():
            tail=self.groups[group][position:]
            del self.groups[group][position:]
            self.groups.insert(group+1,tail)
        self.edit(operation)

    def join(self,group):
        if group+1>=len(self.groups):
            raise ValueError('There is no following scaffold to join')
        def operation():
            self.groups[group].extend(self.groups.pop(group+1))
        self.edit(operation)

    def split_fragment(self,group,position,offset):
        signed=self.groups[group][position]
        key=abs(signed)
        f=self.fragments[key]
        if not isinstance(offset,int) or not 0<offset<f['length']:
            raise ValueError('Cut must lie strictly inside the fragment')
        def operation():
            # offset is measured from the left in the current orientation.
            cut=offset if signed>0 else f['length']-offset
            a,b=max(self.fragments)+1,max(self.fragments)+2
            left,right=copy.deepcopy(f),copy.deepcopy(f)
            left.update(name=f['name']+':::fragment_1',length=cut)
            right.update(name=f['name']+':::fragment_2',length=f['length']-cut,
                         source_start=f['source_start']+cut,map_start=f['map_start']+cut)
            del self.fragments[key]
            self.fragments[a],self.fragments[b]=left,right
            self.groups[group][position:position+1]=[a,b] if signed>0 else [-b,-a]
        self.edit(operation)

    def layout(self):
        result=[]
        offset=0
        for group_id,group in enumerate(self.groups):
            for position,signed in enumerate(group):
                f=self.fragments[abs(signed)]
                result.append(dict(f,group=group_id,position=position,id=abs(signed),reverse=signed<0,start=offset,end=offset+f['length']))
                offset+=f['length']
        return result

    @property
    def length(self):
        return sum(f['length'] for f in self.fragments.values())

    def lift(self,position):
        if not 0<=position<self.length:
            raise ValueError('Coordinate outside assembly')
        for f in self.layout():
            if f['start']<=position<f['end']:
                offset=position-f['start']
                if f['reverse']:
                    offset=f['length']-1-offset
                return dict(chrom=f['map_chrom'],position=(f['map_start']+offset)//f['scale'],
                            source_name=f['source_name'],source_position=f['source_start']+offset,
                            reverse=f['reverse'],fragment=f['name'])

    def write(self,path):
        # Definitions stay in original source order, so consumers treat their
        # cumulative lengths as the original assembly coordinates.
        ordered=sorted(self.fragments,key=lambda k:(self.fragments[k]['map_chrom'],self.fragments[k]['map_start']))
        ids={key:i for i,key in enumerate(ordered,1)}
        lines=[f">{self.fragments[k]['name']} {ids[k]} {self.fragments[k]['length']}" for k in ordered]
        lines += [' '.join(str(ids[abs(v)]*(1 if v>0 else -1)) for v in group) for group in self.groups]
        Path(path).write_text('\n'.join(lines)+'\n',encoding='utf8')

    def write_agp(self,path,gap=0):
        if not isinstance(gap,int) or gap<0:
            raise ValueError('Gap length must be a nonnegative integer')
        lines=['##agp-version 2.1']
        for i,group in enumerate(self.groups,1):
            offset,part=0,1
            for j,signed in enumerate(group):
                f=self.fragments[abs(signed)]
                end=offset+f['length']
                lines.append(f"scaffold_{i}\t{offset+1}\t{end}\t{part}\tW\t{f['source_name']}\t{f['source_start']+1}\t{f['source_start']+f['length']}\t{'-' if signed<0 else '+'}")
                offset,part=end,part+1
                if gap and j+1<len(group):
                    lines.append(f'scaffold_{i}\t{offset+1}\t{offset+gap}\t{part}\tN\t{gap}\tscaffold\tyes\tproximity_ligation')
                    offset,part=offset+gap,part+1
        Path(path).write_text('\n'.join(lines)+'\n',encoding='utf8')

    def write_fasta(self,reference,path,gap=0):
        from .model import Region
        from .sequenceio import chromosome, fasta_index, reference_sequence, version
        # Validate all source coordinates before opening the destination.
        lengths={r[0]:r[1] for r in fasta_index(version(reference))}
        for f in self.fragments.values():
            name=chromosome(lengths,f['source_name'])
            if f['source_start']+f['length']>lengths[name]:
                raise ValueError('Assembly fragment exceeds FASTA source: '+f['source_name'])
        complement=str.maketrans('ACGTRYMKBDHVNacgtrymkbdhvn','TGCAYRKMVHDBNtgcayrkmvhdbn')
        with open(path,'w',encoding='ascii') as out:
            for i,group in enumerate(self.groups,1):
                out.write(f'>scaffold_{i}\n')
                buffer=''
                def emit(seq):
                    nonlocal buffer
                    buffer+=seq
                    n=(len(buffer)//80)*80
                    for at in range(0,n,80):
                        out.write(buffer[at:at+80]+'\n')
                    buffer=buffer[n:]
                for j,signed in enumerate(group):
                    f=self.fragments[abs(signed)]
                    lo,hi=f['source_start'],f['source_start']+f['length']
                    if signed>0:
                        for start in range(lo,hi,1_000_000):
                            emit(reference_sequence(reference,Region(f['source_name'],start,min(hi,start+1_000_000))))
                    else:
                        for end in range(hi,lo,-1_000_000):
                            seq=reference_sequence(reference,Region(f['source_name'],max(lo,end-1_000_000),end))
                            emit(seq.translate(complement)[::-1])
                    if gap and j+1<len(group):
                        for start in range(0,gap,1_000_000):
                            emit('N'*min(1_000_000,gap-start))
                if buffer:
                    out.write(buffer+'\n')


def assembly_matrix(assembly,path,start,end,resolution,state=None,max_bins=1500,edges=None):
    """Lift destination bins into bounded source tiles; never query a dense genome."""
    from .cache import RenderCache, source_chromosome
    from .model import Region
    if not 0<=start<end<=assembly.length or resolution<=0:
        raise ValueError('Invalid assembly interval or resolution')
    edges=np.asarray(edges if edges is not None else [*range(start,end,resolution),end],dtype=np.int64)
    if edges[0]!=start or edges[-1]!=end or np.any(np.diff(edges)<=0):
        raise ValueError("Invalid destination bin edges")
    n=len(edges)-1
    if n>max_bins:
        raise ValueError(f'Assembly preview exceeds {max_bins} bins; use a coarser resolution')
    owned=state is None
    state=state or RenderCache()
    try:
        reader=state.reader(path)
        lengths=dict(reader.chromsizes())
        levels=reader.resolutions()
        scales={f['scale'] for f in assembly.fragments.values()}
        if len(scales)!=1:
            raise ValueError('All assembly fragments must use the same coordinate scale')
        scale=next(iter(scales))
        source_resolution=resolution//scale
        cool=str(path).lower().endswith(('.cool','.mcool'))
        if resolution%scale or source_resolution<=0:
            raise ValueError('Preview resolution must be a positive multiple of the coordinate scale')
        if levels and source_resolution not in levels:
            choices=[v for v in levels if v<=source_resolution]
            if not cool or not choices:
                raise ValueError('Preview resolution must use a stored source resolution or a coarser Cooler resolution')
            return sparse_cool_assembly(assembly,path,edges,max(choices))
        groups: dict = {}
        for f in assembly.layout():
            lo,hi=max(start,f['start']),min(end,f['end'])
            if hi<=lo:
                continue
            name=source_chromosome(reader,f['map_chrom'])
            if name not in lengths or (f['map_start']+f['length'])/scale>lengths[name]+1:
                raise ValueError('Assembly source interval does not fit contact map: '+f['name'])
            first=max(0,int(np.searchsorted(edges,lo,side='right'))-1)
            last=min(n,int(np.searchsorted(edges,hi,side='left')))
            for dest in range(first,last):
                a,b=max(lo,int(edges[dest]))-f['start'],min(hi,int(edges[dest+1]))-f['start']
                if f['reverse']:
                    a,b=f['length']-b,f['length']-a
                a,b=(f['map_start']+a)/scale,(f['map_start']+b)/scale
                for source_bin in range(int(a//source_resolution),math.ceil(b/source_resolution)):
                    width=min(source_resolution,lengths[name]-source_bin*source_resolution)
                    overlap=min(b,(source_bin+1)*source_resolution)-max(a,source_bin*source_resolution)
                    if overlap<=0:
                        continue
                    key=(name,source_bin//128)
                    entries=groups.setdefault(key,{})
                    pair=(dest,source_bin%128)
                    entries[pair]=entries.get(pair,0)+overlap/width
        if len(groups)>128 and cool:
            return sparse_cool_assembly(assembly,path,edges,source_resolution)
        if len(groups)>128:
            raise ValueError('Preview spans more than 128 source tiles; zoom to a smaller assembly region')
        tiles=[]
        for (chrom,tile),entries in groups.items():
            rows=sorted({dest for dest,_ in entries})
            positions={dest:i for i,dest in enumerate(rows)}
            a=tile*128*source_resolution
            b=min(lengths[chrom],a+128*source_resolution)
            weights=np.zeros((len(rows),math.ceil((b-a)/source_resolution)),dtype=np.float32)
            for (dest,col),weight in entries.items():
                weights[positions[dest],col]=weight
            tiles.append((Region(chrom,a,b),rows,weights))
        result=np.zeros((n,n),dtype=np.float32)
        for i,(x,xrows,wx) in enumerate(tiles):
            for j in range(i,len(tiles)):
                y,yrows,wy=tiles[j]
                block=state.rectangle(path,x,y,source_resolution,None)
                if x==y:
                    block=np.array(block,copy=True)
                    diagonal=np.diag_indices_from(block)
                    block[diagonal]*=2
                mapped=wy @ block @ wx.T
                result[np.ix_(yrows,xrows)]+=mapped
                if i!=j:
                    result[np.ix_(xrows,yrows)]+=mapped.T
        result[np.diag_indices_from(result)]*=.5
        return result
    finally:
        if owned:
            state.close()


def render_assembly(job,state=None):
    import matplotlib

    from .contactview import ContactPanel
    from .model import Region
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from poreckit.layout.style import register_hic_colormaps
    register_hic_colormaps()
    from .mapdata import save_map_data, svg_plot_bounds
    assembly=Assembly.restore(job['assembly'])
    start,end,resolution=job['start'],job['end'],job['resolution']
    matrix=assembly_matrix(assembly,job['input'],start,end,resolution,state)
    region=Region('assembly',start,end)
    fig,ax=plt.subplots(figsize=(9,8))
    try:
        panel=ContactPanel(x=region,y=region,matrix=matrix,resolution=resolution,name='Edited assembly · raw contacts',log=job.get('log',True))
        panel.draw(ax)
        for f in assembly.layout():
            if start<f['start']<end:
                color='#177343' if f['position']==0 else '#444444'
                ax.axvline(f['start'],color=color,lw=.8 if f['position']==0 else .3)
                ax.axhline(f['start'],color=color,lw=.8 if f['position']==0 else .3)
            if f['end']>start and f['start']<end and (f['length']/(end-start)>.03):
                mid=(max(start,f['start'])+min(end,f['end']))/2
                ax.text(mid,start,f['name']+(' ←' if f['reverse'] else ' →'),rotation=90,fontsize=6,va='bottom')
        ax.patch.set_gid('contact-map-bounds')
        output=job['output']
        fig.savefig(output,bbox_inches='tight',dpi=150)
        spec=save_map_data(output,matrix,matrix,None,region,region,resolution)
        spec.update(svg_plot_bounds(output) or {})
        return dict(output=output,map_data=spec,resolution=resolution,start=start,end=end,assembly=assembly.data(),
                    fractional_cuts=any((f['map_start']%resolution or f['length']%resolution) for f in assembly.fragments.values()))
    finally:
        plt.close(fig)


def export_cool(assembly,source,path,resolution,state=None):
    """Write a bounded, corrected Cooler with real superscaffold boundaries."""
    import json

    import h5py
    lengths=[sum(assembly.fragments[abs(v)]['length'] for v in group) for group in assembly.groups]
    edges=[0]
    chrom,start,end=[],[],[]
    offset=0
    for i,length in enumerate(lengths):
        for pos in range(0,length,resolution):
            chrom.append(i)
            start.append(pos)
            end.append(min(pos+resolution,length))
            edges.append(offset+min(pos+resolution,length))
        offset+=length
    if len(chrom)>2500:
        raise ValueError('Corrected map export exceeds 2500 bins; choose a coarser stored resolution')
    matrix=assembly_matrix(assembly,source,0,assembly.length,resolution,state,max_bins=2500,edges=edges)
    a,b=np.triu_indices(len(chrom))
    values=matrix[a,b]
    keep=np.isfinite(values) & (values!=0)
    a,b,values=a[keep],b[keep],values[keep]
    with h5py.File(path,'w') as f:
        f.attrs.update({'format':'HDF5::Cooler','format-version':3,'bin-type':'fixed','bin-size':resolution,
                        'storage-mode':'symmetric-upper','nchroms':len(lengths),'nbins':len(chrom),
                        'nnz':len(values),'sum':float(values.sum()),'metadata':json.dumps(dict(
                            generator='poreckit-gui assembly editor',source=str(source),
                            remapping='fractional source-bin overlap; aligned cuts are exact',assembly=assembly.data()))})
        c=f.create_group('chroms')
        c['name']=np.array([f'scaffold_{i+1}'.encode() for i in range(len(lengths))])
        c['length']=np.asarray(lengths,dtype='int64')
        bins=f.create_group('bins')
        bins['chrom']=np.asarray(chrom,dtype='int32')
        bins['start']=np.asarray(start,dtype='int64')
        bins['end']=np.asarray(end,dtype='int64')
        pixels=f.create_group('pixels')
        pixels.create_dataset('bin1_id',data=a.astype('int64'),compression='gzip')
        pixels.create_dataset('bin2_id',data=b.astype('int64'),compression='gzip')
        pixels.create_dataset('count',data=values,compression='gzip')
        indexes=f.create_group('indexes')
        indexes['chrom_offset']=np.r_[0,np.cumsum([math.ceil(n/resolution) for n in lengths])].astype('int64')
        indexes['bin1_offset']=np.searchsorted(a,np.arange(len(chrom)+1)).astype('int64')


def export_assembly(job,state=None):
    import os
    import tempfile
    document=Assembly.restore(job['assembly'])
    destination=Path(job['output']).resolve()
    for source in [job['input'],job.get('reference','')]:
        if source and destination==Path(source).resolve():
            raise ValueError('The export destination must differ from source data')
    fd,temp=tempfile.mkstemp(prefix='.poreckit-export-',suffix=destination.suffix,dir=destination.parent)
    os.close(fd)
    try:
        if job['format']=='fasta':
            document.write_fasta(job['reference'],temp)
        elif job['format']=='cool':
            export_cool(document,job['input'],temp,int(job['resolution']),state)
        else:
            raise ValueError('Unsupported assembly export')
        os.replace(temp,destination)
    finally:
        Path(temp).unlink(missing_ok=True)
    return dict(output=str(destination))


def sparse_cool_assembly(assembly,path,edges,source_resolution):
    """Stream sparse Cooler pixels for overviews spanning many source contigs."""
    import h5py

    from .sequenceio import chromosome
    n=len(edges)-1
    src_ids,dest_ids,weights=[],[],[]
    with h5py.File(path,'r') as file:
        g=file[f'resolutions/{source_resolution}'] if 'resolutions' in file else file
        names=[v.decode() if isinstance(v,bytes) else str(v) for v in g['chroms/name'][:]]
        lengths=dict(zip(names,map(int,g['chroms/length'][:]),strict=False))
        offsets=g['indexes/chrom_offset'][:]
        for f in assembly.layout():
            lo,hi=max(int(edges[0]),f['start']),min(int(edges[-1]),f['end'])
            if hi<=lo:
                continue
            chrom=chromosome(lengths,f['map_chrom'])
            cid=names.index(chrom)
            if (f['map_start']+f['length'])/f['scale']>lengths[chrom]+1:
                raise ValueError('Assembly source interval exceeds contact map')
            first=max(0,int(np.searchsorted(edges,lo,side='right'))-1)
            last=min(n,int(np.searchsorted(edges,hi,side='left')))
            for dest in range(first,last):
                a,b=max(lo,int(edges[dest]))-f['start'],min(hi,int(edges[dest+1]))-f['start']
                if f['reverse']:
                    a,b=f['length']-b,f['length']-a
                a,b=(f['map_start']+a)/f['scale'],(f['map_start']+b)/f['scale']
                for source_bin in range(int(a//source_resolution),math.ceil(b/source_resolution)):
                    width=min(source_resolution,lengths[chrom]-source_bin*source_resolution)
                    overlap=min(b,(source_bin+1)*source_resolution)-max(a,source_bin*source_resolution)
                    if overlap>0:
                        src_ids.append(int(offsets[cid])+source_bin)
                        dest_ids.append(dest)
                        weights.append(overlap/width)
                    if len(src_ids)>2_000_000:
                        raise ValueError('Assembly overview exceeds two million source overlaps; use a coarser source map')
        order=np.argsort(src_ids,kind='stable')
        src_arr=np.asarray(src_ids,dtype=np.int64)[order]
        dest_arr=np.asarray(dest_ids,dtype=np.int64)[order]
        weights_arr=np.asarray(weights,dtype=np.float32)[order]
        result=np.zeros((n,n),dtype=np.float32)
        storage=g.attrs.get('storage-mode','symmetric-upper')
        if isinstance(storage,bytes):
            storage=storage.decode()
        if storage not in {'square','symmetric-upper'}:
            raise ValueError('Unsupported Cooler storage mode')
        for at in range(0,len(g['pixels/count']),262144):
            a=g['pixels/bin1_id'][at:at+262144]
            b=g['pixels/bin2_id'][at:at+262144]
            values=np.asarray(g['pixels/count'][at:at+262144],dtype=np.float32)
            values[a==b]*=2
            left_a,right_a=np.searchsorted(src_arr,a,'left'),np.searchsorted(src_arr,a,'right')
            left_b,right_b=np.searchsorted(src_arr,b,'left'),np.searchsorted(src_arr,b,'right')
            valid=(right_a>left_a)&(right_b>left_b)&np.isfinite(values)
            a,b,values=a[valid],b[valid],values[valid]
            left_a,right_a,left_b,right_b=left_a[valid],right_a[valid],left_b[valid],right_b[valid]
            for i in range(int(np.max(right_a-left_a,initial=0))):
                for j in range(int(np.max(right_b-left_b,initial=0))):
                    selected=(left_a+i<right_a)&(left_b+j<right_b)
                    ai,bi=left_a[selected]+i,left_b[selected]+j
                    x,y=dest_arr[ai],dest_arr[bi]
                    counts=values[selected]*weights_arr[ai]*weights_arr[bi]
                    np.add.at(result,(y,x),counts)
                    if storage=='symmetric-upper':
                        off=a[selected]!=b[selected]
                        np.add.at(result,(x[off],y[off]),counts[off])
        result[np.diag_indices_from(result)]*=.5
        return result
