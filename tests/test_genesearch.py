import gzip
import json
import time

from poreckit_gui.demo import create_demo
from poreckit_gui.genesearch import search_genes


def test_gtf_names_aliases_coordinate_conversion_and_index_invalidation(tmp_path):
    path = tmp_path / "genes.gtf"
    path.write_text(
        'chr1\tsrc\tgene\t101\t400\t.\t+\t.\tgene_id "G1"; gene_name "TP53";\n'
        'chr1\tsrc\texon\t151\t350\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "TP53";\n'
        'chr2\tsrc\texon\t201\t500\t.\t-\t.\tgene_id "G2"; gene_name "TP53BP1";\n'
        'chr3\tsrc\tgene\t1\t100\t.\t+\t.\tgene_id "G3"; gene_name "TP53";\n')
    chroms = [("1", 1000), ("2", 1000)]
    result = search_genes([str(path), str(path)], "tp53", chroms)
    assert result["total"] == 2
    assert result["results"][0]["name"] == "TP53"
    assert result["results"][0]["chrom"] == "1"
    assert result["results"][0]["start"] == 100
    assert result["results"][0]["end"] == 400
    assert search_genes([str(path)], "T1", chroms)["results"][0]["name"] == "TP53"
    path.write_text(path.read_text().replace("TP53", "NEWGENE"))
    assert search_genes([str(path)], "tp53", chroms)["total"] == 0
    assert search_genes([str(path)], "NEWGENE", chroms)["total"] == 2


def test_gff_and_bed_search_and_result_limits(tmp_path):
    path = tmp_path / "genes.gff3.gz"
    with gzip.open(path, "wt") as stream:
        stream.write('##gff-version 3\nchr1\tsrc\tgene\t11\t50\t.\t+\t.\tID=gene:G1;Name=Gene%20One;Alias=ABC,Alias2\n'
                     'chr1\tsrc\tmRNA\t11\t50\t.\t+\t.\tID=tx:T1;Parent=gene:G1\n'
                     '##FASTA\n>chr1\nACGT\n')
    gene = search_genes([str(path)], "G1", [("chr1", 1000)])["results"][0]
    assert gene["name"] == "Gene One" and gene["start"] == 10
    assert search_genes([str(path)], "Alias2", [("chr1", 1000)])["total"] == 1
    bed = tmp_path / "features.bed"
    bed.write_text("chr1\t10\t20\tABC1\nchr1\t30\t40\tABC2\n")
    result = search_genes([str(path), str(bed)], "ABC", [("chr1", 1000)], limit=2)
    assert result["total"] == 3 and len(result["results"]) == 2
    assert result["results"][0]["name"] == "Gene One"


def test_gene_search_gui_axis_selection_and_session_sources(tmp_path, monkeypatch):
    from PySide6.QtWidgets import QApplication, QFileDialog

    from poreckit_gui.app import Browser
    app = QApplication.instance() or QApplication([])
    b = Browser()
    b.show()
    b.auto_update.setChecked(False)
    errors = []
    b.error = lambda message: errors.append(str(message))
    def wait():
        deadline = time.monotonic()+30
        while b.busy() and time.monotonic() < deadline:
            app.processEvents()
            time.sleep(.01)
        app.processEvents()
        assert not b.busy(), b.log_view.toPlainText()
        assert not errors, errors
    try:
        b.load_file(create_demo(tmp_path / "data", [("chr1", 1200000), ("chr2", 800000)]))
        wait()
        genes = tmp_path / "genes.gtf"
        genes.write_text('2\tsrc\tgene\t300001\t350000\t.\t+\t.\tgene_id "G1"; gene_name "TARGET";\n')
        dialog = b.gene_search
        dialog.paths.append(str(genes))
        b.locus.setText("target")
        b.go()
        wait()
        assert dialog.isVisible()
        assert dialog.table.rowCount() == 1
        assert dialog.table.item(0, 1).text() == "chr2"
        assert b.locus.text().startswith("chr1:")
        y = b.current_viewport().y
        dialog.axis.setCurrentIndex(1)
        dialog.navigate()
        wait()
        assert b.current_viewport().x.chrom == "chr2"
        assert b.current_viewport().y == y
        assert b.current_viewport().x.start <= 300000 < b.current_viewport().x.end
        dialog.open_search()
        dialog.axis.setCurrentIndex(0)
        dialog.navigate()
        wait()
        assert b.current_viewport().diagonal
        session = tmp_path / "search-session.json"
        monkeypatch.setattr(QFileDialog, "getSaveFileName", lambda *a, **k: (str(session), ""))
        b.save_session()
        assert json.loads(session.read_text())["gene_annotations"] == ["genes.gtf"]
        dialog.paths.clear()
        monkeypatch.setattr(QFileDialog, "getOpenFileName", lambda *a, **k: (str(session), ""))
        b.load_session()
        wait()
        assert dialog.paths == [str(genes)]
    finally:
        b.close()
