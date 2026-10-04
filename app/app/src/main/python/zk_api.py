"""Chaquopy entry point: exposes the ZK core to the WebView UI.

Called from Java (ZkBridge) via callAttr. All functions take/return
plain strings; results are JSON: {"ok": true, ...} or
{"ok": false, "error": "..."}.

The core/ package is copied next to this file by CI before the build
(see .github/workflows/build.yml); locally, symlink ../core here.
"""

import json
import traceback

from core import api as zk


def _err(msg):
    return json.dumps({"ok": False, "error": msg}, ensure_ascii=False)


def prove(equation, witness):
    """Prove knowledge of `witness` (int, as string) for `equation`."""
    try:
        w = int(witness.strip())
    except (ValueError, AttributeError):
        return _err("witness must be an integer")
    try:
        proof = zk.prove_equation(equation.strip(), w)
    except ValueError as e:
        return _err(str(e))
    except Exception:
        return _err("prove failed: " + traceback.format_exc(limit=1).strip())
    return json.dumps({"ok": True, "proof": proof}, ensure_ascii=False)


def verify(equation, proof_json):
    """Verify a proof (JSON string) against `equation`."""
    try:
        proof = json.loads(proof_json)
    except (ValueError, TypeError):
        return _err("proof is not valid JSON")
    try:
        valid = zk.verify_equation(equation.strip(), proof)
    except Exception:
        return _err("verify failed")
    return json.dumps({"ok": True, "valid": bool(valid)}, ensure_ascii=False)
