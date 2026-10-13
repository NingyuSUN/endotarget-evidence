"""Boundary fixtures are synthetic and never scientific analysis inputs."""
import hashlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from scripts import m0_metadata as m

SCRATCH = Path(__file__).resolve().parents[1]/"work"

def setUpModule():
    SCRATCH.mkdir(parents=True, exist_ok=True)

D = "Characteristics[individual]"
A = "Characteristics[10x technology]"
C = "Characteristics[disease]"
P = "Characteristics[menstrual_cycle_stage]"

def table(rows, fields=None):
    fields = fields or ["Source Name", D, A, C, P]
    return ("\t".join(fields) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n").encode()

class DonorBoundaries(unittest.TestCase):
    def audit(self, rows):
        return m.audit_sdrf(table(rows))
    def test_technical_libraries_one_donor(self):
        x = self.audit([["lib1", "P1", "snRNAseq", "endometriosis", "Proliferative"],
                        ["lib2", "P1", "snRNAseq", "endometriosis", "Proliferative"]])
        self.assertEqual(x["library_rows"], 2)
        self.assertEqual(x["donor_units"], 1)
        self.assertEqual(x["candidate_proliferative_counts"][0]["donors"], 1)
    def test_same_id_across_assays_not_merged(self):
        x = self.audit([["l1", "P1", "snRNAseq", "normal", "Proliferative"],
                        ["l2", "P1", "scRNAseq", "normal", "Proliferative"]])
        self.assertEqual(x["donor_units"], 2)
        self.assertEqual(len(x["donors_by_assay"]), 2)
    def test_conflicting_clinical_labels_gate(self):
        x = self.audit([["l1", "P1", "snRNAseq", "normal", "Proliferative"],
                        ["l2", "P1", "snRNAseq", "endometriosis", "Proliferative"]])
        self.assertTrue(x["within_donor_conflicts"])
        self.assertEqual(x["candidate_proliferative_counts"], [])
        self.assertFalse(x["readiness"]["stage0_passed"])
    def test_missing_required_column_fails(self):
        with self.assertRaises(ValueError):
            m.audit_sdrf(table([["l1", "P1"]], ["Source Name", D]))
    def test_html_is_not_table(self):
        with self.assertRaises(ValueError):
            m.audit_sdrf(b"<html><body>404 not found</body></html>")
    def test_missing_donor_not_inferred(self):
        x = self.audit([["file_P1", "", "snRNAseq", "normal", "Proliferative"]])
        self.assertEqual(x["donor_units"], 0)
        self.assertEqual(len(x["rows_missing_identity"]), 1)
    def test_metadata_cannot_pass_stage0(self):
        x = self.audit([["l1", "P1", "snRNAseq", "normal", "Proliferative"]])
        self.assertEqual(x["readiness"]["final_eligibility"], "NOT_VERIFIED")
        self.assertFalse(x["readiness"]["stage0_passed"])
    def test_blank_label_is_conflict_and_missing(self):
        x = self.audit([["l1", "P1", "snRNAseq", "normal", "Proliferative"],
                        ["l2", "P1", "snRNAseq", "", "Proliferative"]])
        self.assertTrue(x["within_donor_conflicts"])
        self.assertTrue(x["missing_clinical_fields"])
        self.assertEqual(x["candidate_proliferative_counts"], [])
    def test_ragged_table_fails(self):
        with self.assertRaises(ValueError):
            m.audit_sdrf(table([["l1", "P1", "snRNAseq", "normal", "Proliferative", "extra"]]))

class FetchBoundaries(unittest.TestCase):
    def test_html_rejected_as_json(self):
        with self.assertRaises(ValueError):
            m.validate_payload(b"<html>Error</html>", "biostudies_json")
    def test_wrong_accession_rejected(self):
        with self.assertRaises(ValueError):
            m.validate_payload(b'{"accno":"E-MTAB-999","section":{}}', "biostudies_json")
    def test_over_limit_stream_fails(self):
        with tempfile.TemporaryDirectory(dir=SCRATCH) as td:
            out = Path(td)/"x.part"
            with self.assertRaises(ValueError):
                m.stream_to_part(io.BytesIO(b"x"*20), out, 10)
            self.assertFalse(out.exists())
    def test_stream_success_hash_and_bytes(self):
        with tempfile.TemporaryDirectory(dir=SCRATCH) as td:
            out = Path(td)/"x.part"
            size, sha = m.stream_to_part(io.BytesIO(b"hello"), out, 10)
            self.assertEqual((size,sha),(5,hashlib.sha256(b"hello").hexdigest()))
            self.assertEqual(out.read_bytes(), b"hello")
    def test_non_https_url_rejected(self):
        with tempfile.TemporaryDirectory(dir=SCRATCH) as td:
            with self.assertRaises(ValueError):
                m.fetch_one({"id":"x","url":"http://example.com","kind":"biostudies_json","filename":"x.json"}, Path(td))
    def test_cache_hash_mismatch_fails(self):
        with tempfile.TemporaryDirectory(dir=SCRATCH) as td:
            root = Path(td)
            (root/"data/metadata").mkdir(parents=True)
            (root/"receipts").mkdir()
            (root/"data/metadata/x.json").write_bytes(b'{"accno":"E-MTAB-14039","section":{}}')
            (root/"receipts/m0_metadata_x.json").write_text(json.dumps({"sha256":"wrong","bytes":42}))
            with self.assertRaises(ValueError):
                m.fetch_one({"id":"x","url":"https://example.com","kind":"biostudies_json","filename":"x.json"}, root)
    def test_valid_cache_reused_without_network(self):
        with tempfile.TemporaryDirectory(dir=SCRATCH) as td:
            root = Path(td)
            (root/"data/metadata").mkdir(parents=True)
            (root/"receipts").mkdir()
            data = b'{"accno":"E-MTAB-14039","section":{}}'
            (root/"data/metadata/x.json").write_bytes(data)
            receipt = {"sha256":hashlib.sha256(data).hexdigest(),"bytes":len(data),
                       "url":"https://example.com","kind":"biostudies_json","status":"downloaded"}
            (root/"receipts/m0_metadata_x.json").write_text(json.dumps(receipt))
            with patch("urllib.request.OpenerDirector.open", side_effect=AssertionError("network")):
                self.assertEqual(m.fetch_one({"id":"x","url":"https://example.com","kind":"biostudies_json","filename":"x.json"}, root)["status"], "reused_validated")
    def test_refresh_changed_source_preserves_original(self):
        with tempfile.TemporaryDirectory(dir=SCRATCH) as td:
            root = Path(td)
            (root/"data/metadata").mkdir(parents=True)
            (root/"receipts").mkdir()
            old = b'{"accno":"E-MTAB-14039","section":{}}'
            new = b'{"accno":"E-MTAB-14039","section":{"type":"Study"}}'
            (root/"data/metadata/x.json").write_bytes(old)
            receipt = {"sha256":hashlib.sha256(old).hexdigest(),"bytes":len(old),
                       "url":"https://example.com","kind":"biostudies_json","status":"downloaded"}
            (root/"receipts/m0_metadata_x.json").write_text(json.dumps(receipt))
            response = io.BytesIO(new)
            response.headers = {}
            response.geturl = lambda: "https://example.com"
            with patch("urllib.request.OpenerDirector.open", return_value=response):
                result = m.fetch_one({"id":"x","url":"https://example.com","kind":"biostudies_json","filename":"x.json"}, root, refresh=True)
            self.assertEqual(result["status"], "source_changed")
            self.assertEqual((root/"data/metadata/x.json").read_bytes(), old)
            self.assertEqual(Path(result["changed_path"]).read_bytes(), new)


class SourceRegression(unittest.TestCase):
    def test_true_source_phase_header(self):
        payload = table([["l1","ES1","snRNAseq","normal","Proliferative"]],
                        ["Source Name",D,A,C,"Characteristics[menstrual_cycle_stage]"])
        self.assertEqual(m.audit_sdrf(payload)["donor_units"],1)
    def test_library_id_difference_is_not_clinical_conflict(self):
        payload = table([["l1","ES1","snRNAseq","normal","Proliferative","L1"],
                         ["l2","ES1","snRNAseq","normal","Proliferative","L2"]],
                        ["Source Name",D,A,C,P,"Characteristics[library_id]"])
        audit = m.audit_sdrf(payload)
        self.assertEqual(audit["within_donor_conflicts"],[])
        self.assertEqual(audit["candidate_proliferative_counts"][0]["donors"],1)

class SoftBoundaries(unittest.TestCase):
    def test_valid_family_soft(self):
        import gzip
        payload = b"^SERIES = GSE266265\n!Series_title = Study\n^SAMPLE = GSM123\n!Sample_title = sample\n"
        self.assertEqual(m.validate_payload(gzip.compress(payload),"geo_soft_gz")["sample_records"],1)
    def test_compressed_bomb_rejected(self):
        import gzip
        with self.assertRaises(ValueError):
            m.validate_payload(gzip.compress(b"A"*(m.MAX_BYTES+1)),"geo_soft_gz")
    def test_wrong_soft_accession_rejected(self):
        import gzip
        with self.assertRaises(ValueError):
            m.validate_payload(gzip.compress(b"^SERIES = GSE1\n!Series_title = Wrong\n"),"geo_soft_gz")


class DeadlineRegression(unittest.TestCase):
    def test_trickling_http_read_is_interrupted(self):
        import http.client
        import socket
        import threading
        import time
        reader, writer = socket.socketpair()
        writer.sendall(b"HTTP/1.1 200 OK\r\nContent-Length: 30\r\n\r\n")
        stop = threading.Event()
        def trickle():
            try:
                for _ in range(30):
                    if stop.wait(0.01):
                        break
                    writer.sendall(b"x")
            except OSError:
                pass
        thread = threading.Thread(target=trickle)
        thread.start()
        response = http.client.HTTPResponse(reader)
        response.begin()
        try:
            with tempfile.TemporaryDirectory(dir=SCRATCH) as td:
                part = Path(td)/"slow.part"
                start = time.monotonic()
                with patch.object(m, "TOTAL_TIMEOUT", 0.05):
                    with self.assertRaises(TimeoutError):
                        m.stream_to_part(response, part, 100)
                elapsed = time.monotonic()-start
                self.assertLess(elapsed, 0.20, "A single trickling read exceeded the overall deadline")
                self.assertFalse(part.exists())
        finally:
            stop.set()
            response.close()
            reader.close()
            writer.close()
            thread.join(timeout=1)

    def test_overall_deadline_includes_open(self):
        import time
        def slow_open(*args, **kwargs):
            time.sleep(0.30)
            response = io.BytesIO(b'{"accno":"E-MTAB-14039","section":{}}')
            response.headers = {}
            response.geturl = lambda: "https://example.com"
            return response
        with tempfile.TemporaryDirectory(dir=SCRATCH) as td:
            start = time.monotonic()
            with patch.object(m, "TOTAL_TIMEOUT", 0.05):
                with patch("urllib.request.OpenerDirector.open", side_effect=slow_open):
                    with self.assertRaises(TimeoutError):
                        m.fetch_one({"id":"slow","url":"https://example.com","kind":"biostudies_json","filename":"x.json"},Path(td))
            self.assertLess(time.monotonic()-start,0.20)
            self.assertFalse((Path(td)/"data/metadata/x.json").exists())

if __name__ == "__main__":
    unittest.main()
