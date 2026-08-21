from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from solvelec.provenance import sha256_file, write_manifest

ROOT = Path(__file__).resolve().parents[1]


class ProvenanceTests(unittest.TestCase):
    def test_sha256_and_manifest(self) -> None:
        source = ROOT / "configs" / "campaign.yaml"
        self.assertEqual(len(sha256_file(source)), 64)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "manifest.json"
            write_manifest(output, ROOT, [source], campaign="smoke")
            data = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(data["campaign"], "smoke")
            self.assertEqual(data["inputs"][0]["sha256"], sha256_file(source))
            self.assertIn("capabilities", data)


if __name__ == "__main__":
    unittest.main()
