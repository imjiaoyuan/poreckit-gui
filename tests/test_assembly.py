import json
from pathlib import Path

import h5py
import numpy as np
import pytest

from poreckit_gui.assembly import Assembly, assembly_matrix, export_assembly, export_cool, render_assembly
from poreckit_gui.demo import create_demo
from poreckit_gui.model import Region
from poreckit_gui.sequenceio import reference_sequence


def document():
    return Assembly.from_chromsizes([('a',4),('b',4),('c',4)])


def test_assembly_edits_lifting_and_atomic_undo():
    a=document()
    a.join(0)
    assert a.groups==[[1,2],[3]]
    a.invert(0,0,1)
    assert a.groups==[[-2,-1],[3]]
    assert a.lift(0)['chrom']=='b' and a.lift(0)['position']==3
    a.split_fragment(0,0,1)
    assert [a.fragments[abs(v)]['length'] for v in a.groups[0]]==[1,3,4]
    assert [a.lift(i)['position'] for i in range(4)]==[3,2,1,0]
    a.move(0,0,0,1,1)
    assert a.groups[1][-1]<0
    a.undo()
    assert len(a.groups[0])==3
    a.redo()
    assert len(a.groups[0])==2
    before=a.data()
    with pytest.raises(ValueError):
        a.edit(lambda:a.groups[0].append(999))
    assert a.data()==before
    assert a.dirty
    restored=Assembly.restore(json.loads(json.dumps(a.data())))
    assert restored.data()==a.data()


def test_assembly_agp_fasta_and_layout_roundtrip(tmp_path):
    source=tmp_path/'source.fa'
    source.write_text('>a\nAACC\n>b\nATGC\n>c\nTTAA\n')
    a=document()
    a.join(0)
    a.invert(0,1,1)
    a.split_fragment(0,1,2)
    a.write_fasta(str(source),tmp_path/'out.fa')
    assert reference_sequence(tmp_path/'out.fa',Region('scaffold_1',0,8))=='AACCGCAT'
    a.write_agp(tmp_path/'out.agp')
    agp=(tmp_path/'out.agp').read_text()
    assert '\tb\t3\t4\t-' in agp and '\tb\t1\t2\t-' in agp
    a.write(tmp_path/'out.assembly')
    loaded=Assembly.read(tmp_path/'out.assembly')
    assert loaded.length==a.length
    assert [[loaded.fragments[abs(v)]['length'] for v in g] for g in loaded.groups]==[[a.fragments[abs(v)]['length'] for v in g] for g in a.groups]
    assert [[v<0 for v in g] for g in loaded.groups]==[[v<0 for v in g] for g in a.groups]


def test_remapped_map_is_exact_permutation_and_preserves_counts(tmp_path):
    from poreckit import make_reader
    path=create_demo(tmp_path,[('chr1',40000),('chr2',40000)])
    a=Assembly.from_chromsizes([('chr1',40000),('chr2',40000)])
    original=assembly_matrix(a,path,0,80000,10000)
    a.join(0)
    a.invert(0,1,1)
    a.move(0,1,1,0,0)
    changed=assembly_matrix(a,path,0,80000,10000)
    order=[7,6,5,4,0,1,2,3]
    np.testing.assert_allclose(changed,original[np.ix_(order,order)])
    assert changed.sum()==pytest.approx(original.sum())
    # An interior cut allocates source bins by overlap, conserving total counts.
    a.split_fragment(0,0,15000)
    partial=assembly_matrix(a,path,0,80000,10000)
    np.testing.assert_allclose(partial,changed)
    result=render_assembly(dict(assembly=a.data(),input=path,start=0,end=80000,resolution=10000,output=str(tmp_path/'map.svg')))
    assert Path(result['output']).stat().st_size>1000
    assert Path(result['map_data']['path']).is_file()
    out=tmp_path/'corrected.cool'
    export_cool(a,path,str(out),10000)
    reader=make_reader(str(out))
    try:
        readback=reader.matrix('scaffold_1',0,80000,10000,None)
        np.testing.assert_allclose(readback,changed)
    finally:
        reader.close()


def test_export_failure_keeps_destination_and_source_intact(tmp_path):
    path=create_demo(tmp_path,[('chr1',40000)])
    a=Assembly.from_chromsizes([('chr1',40000)])
    dest=tmp_path/'out.cool'
    dest.write_bytes(b'previous')
    with pytest.raises(ValueError):
        export_assembly(dict(assembly=a.data(),input=path,output=str(dest),format='cool',resolution=123))
    assert dest.read_bytes()==b'previous'
    assert not list(tmp_path.glob('.poreckit-export-*'))
    with pytest.raises(ValueError,match='differ'):
        export_assembly(dict(assembly=a.data(),input=path,output=path,format='cool',resolution=10000))


def test_corrected_cooler_partial_scaffold_bins(tmp_path):
    path=create_demo(tmp_path,[('chr1',40000)])
    a=Assembly.from_chromsizes([('chr1',40000)])
    a.split_fragment(0,0,15000)
    a.split_group(0,1)
    export_cool(a,path,str(tmp_path/'split.cool'),10000)
    with h5py.File(tmp_path/'split.cool') as f:
        assert list(f['chroms/length'])==[15000,25000]
        assert list(f['bins/start'])==[0,10000,0,10000,20000]
        assert list(f['bins/end'])==[10000,15000,10000,20000,25000]
        assert list(f['indexes/chrom_offset'])==[0,2,5]


def test_assembly_editor_edits_live_map_and_restores_project(tmp_path,monkeypatch):
    import time

    from PySide6.QtWidgets import QApplication, QFileDialog, QMessageBox

    from poreckit_gui.app import Browser
    app=QApplication.instance() or QApplication([])
    b=Browser()
    editor=b.assembly_editor
    def wait(predicate):
        deadline=time.monotonic()+15
        while not predicate():
            app.processEvents()
            time.sleep(.01)
            assert time.monotonic()<deadline,editor.status.text()
    try:
        b.load_file(create_demo(tmp_path,[('chr1',40000),('chr2',40000)]))
        wait(lambda:not b.busy() and b.last_view is not None)
        b.tabs.setCurrentWidget(editor)
        editor.from_map()
        wait(lambda:not editor.busy() and editor.preview_result is not None)
        initial=editor.canvas.map_data.observed.copy()
        item=editor.tree.topLevelItem(1).child(0)
        item.setSelected(True)
        editor.invert()
        wait(lambda:not editor.busy())
        assert editor.document.groups==[[1],[-2]]
        assert not np.array_equal(editor.canvas.map_data.observed,initial)
        project=tmp_path/'project.json'
        monkeypatch.setattr(QFileDialog,'getSaveFileName',lambda *a,**kw:(str(project),'JSON (*.json)'))
        editor.save_session()
        assert not editor.document.dirty
        data=json.loads(project.read_text())
        assert data['input']==Path(b.input_path).name
        editor.document.invert(1,0,0)
        assert editor.document.dirty
        monkeypatch.setattr(QMessageBox,'question',lambda *a,**kw:QMessageBox.StandardButton.Discard)
        monkeypatch.setattr(QFileDialog,'getOpenFileName',lambda *a,**kw:(str(project),'JSON (*.json)'))
        editor.load_session()
        wait(lambda:not editor.busy())
        assert editor.document.groups==[[1],[-2]]
    finally:
        if editor.document:
            editor.document.mark_saved()
        b.close()


def test_coarse_sparse_assembly_preserves_contact_pairs(tmp_path):
    from poreckit_gui.assembly import sparse_cool_assembly
    path=create_demo(tmp_path,[('chr1',40000),('chr2',40000)])
    a=Assembly.from_chromsizes([('chr1',40000),('chr2',40000)])
    fine=assembly_matrix(a,path,0,80000,10000)
    coarse=assembly_matrix(a,path,0,80000,20000)
    # Contact counts on the diagonal occur once; off-diagonal source pixels
    # that merge into one coarse diagonal bin must not be counted twice.
    expected=fine.copy()
    expected[np.diag_indices_from(expected)]*=2
    expected=expected.reshape(4,2,4,2).sum(axis=(1,3))
    expected[np.diag_indices_from(expected)]*=.5
    np.testing.assert_allclose(coarse,expected)
    assert np.triu(coarse).sum()==pytest.approx(np.triu(fine).sum())
    edges=list(range(0,80001,10000))
    np.testing.assert_allclose(sparse_cool_assembly(a,path,edges,10000),fine)
    a.split_fragment(0,0,15000)
    a.move(0,0,0,1,1)
    edited=assembly_matrix(a,path,0,80000,10000)
    assert np.triu(edited).sum()==pytest.approx(np.triu(fine).sum())


def test_save_project_refuses_temporary_source_and_keeps_dirty_state(tmp_path,monkeypatch):
    import tempfile

    from PySide6.QtWidgets import QApplication, QFileDialog

    from poreckit_gui.assemblyview import AssemblyEditor
    _app=QApplication.instance() or QApplication([])
    class StubBrowser:
        def __init__(self):
            self.temp=tempfile.TemporaryDirectory(prefix='poreckit-stub-browser-')
    browser=StubBrowser()
    editor=AssemblyEditor(browser)
    try:
        # Source lives in the browser's temporary directory, like a synthetic demo.
        source=Path(browser.temp.name)/'synthetic-demo.cool'
        source.write_bytes(b'cool')
        a=document()
        a.join(0)
        editor.document=a
        editor.input_path=str(source)
        if not editor.resolution.currentText():
            editor.resolution.addItem('10000')
        project=tmp_path/'project.json'
        monkeypatch.setattr(QFileDialog,'getSaveFileName',lambda *a,**kw:(str(project),'JSON (*.json)'))
        editor.save_session()
        assert a.dirty
        assert 'temporary' in editor.status.text()
        assert not project.exists()
        # A persistent source still saves normally and clears the dirty flag.
        persistent=tmp_path/'my-map.cool'
        persistent.write_bytes(b'cool')
        editor.input_path=str(persistent)
        editor.save_session()
        assert not a.dirty
        assert json.loads(project.read_text())['input']==persistent.name
    finally:
        browser.temp.cleanup()
