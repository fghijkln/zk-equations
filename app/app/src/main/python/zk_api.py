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


def prove(equation, witness, precision=""):
    """Prove knowledge of `witness` for `equation`.

    witness: integer string (exact mode) or decimal string (tolerance).
    precision: k as string; required in tolerance mode (final precision
    1/k), ignored in exact mode.
    """
    try:
        proof = zk.prove_equation(equation.strip(), witness.strip(),
                                  (precision or "").strip())
    except ValueError as e:
        return _err(str(e))
    except Exception:
        return _err("prove failed: " + traceback.format_exc(limit=1).strip())
    return json.dumps({"ok": True, "proof": proof}, ensure_ascii=False)


def verify(equation, proof_json, precision=""):
    """Verify a proof (JSON string) against `equation`.

    precision: k as string; if empty, taken from the proof itself.
    """
    try:
        proof = json.loads(proof_json)
    except (ValueError, TypeError):
        return _err("proof is not valid JSON")
    try:
        valid = zk.verify_equation(equation.strip(), proof,
                                   (precision or "").strip())
    except Exception:
        return _err("verify failed")
    return json.dumps({"ok": True, "valid": bool(valid)}, ensure_ascii=False)


def solve(equation, decimals="8"):
    """Numerically solve `equation`; return up to 5 real roots.

    decimals: digits after the point per root (default 8).
    Returns {"ok": true, "roots": ["0.52359878", ...]}.
    """
    try:
        roots = zk.solve_equation(equation.strip(), (decimals or "8").strip())
    except ValueError as e:
        return _err(str(e))
    except Exception:
        return _err("solve failed: " + traceback.format_exc(limit=1).strip())
    return json.dumps({"ok": True, "roots": roots}, ensure_ascii=False)


def certificate(proof_json):
    """Render a proof as a human-readable certificate (math-proof style).

    Returns {"ok": true, "cert": {...}} with bilingual prose.
    """
    try:
        proof = json.loads(proof_json)
    except (ValueError, TypeError):
        return _err("proof is not valid JSON")
    try:
        cert = zk.describe_proof(proof)
    except ValueError as e:
        return _err(str(e))
    except Exception:
        return _err("certificate failed: " + traceback.format_exc(limit=1).strip())
    return json.dumps({"ok": True, "cert": cert}, ensure_ascii=False)


# ---- committed inputs (debug/experimental) ----

def commit(equation, witness, precision=""):
    """Create a Pedersen commitment V to the equation's input wire.

    Returns {"ok": true, "V": hex, "gamma": hex}. gamma must be kept
    secret and passed to prove_committed.
    """
    try:
        r = zk.commit_input(equation.strip(), witness.strip(),
                            (precision or "").strip())
    except ValueError as e:
        return _err(str(e))
    except Exception:
        return _err("commit failed: " + traceback.format_exc(limit=1).strip())
    return json.dumps({"ok": True, "V": r["V"], "gamma": r["gamma"]},
                      ensure_ascii=False)


def prove_committed(equation, witness, gamma, precision=""):
    """Prove with the input wire publicly committed (see commit).

    Returns {"ok": true, "proof": {...}}; the proof carries "V".
    """
    try:
        proof = zk.prove_committed_equation(equation.strip(), witness.strip(),
                                            (gamma or "").strip(),
                                            (precision or "").strip())
    except ValueError as e:
        return _err(str(e))
    except Exception:
        return _err("prove failed: " + traceback.format_exc(limit=1).strip())
    return json.dumps({"ok": True, "proof": proof}, ensure_ascii=False)


def verify_committed(equation, V, proof_json, precision=""):
    """Verify a committed-input proof against equation and V.

    Returns {"ok": true, "valid": bool}.
    """
    try:
        proof = json.loads(proof_json)
    except (ValueError, TypeError):
        return _err("proof is not valid JSON")
    try:
        valid = zk.verify_committed_equation(equation.strip(), (V or "").strip(),
                                             proof, (precision or "").strip())
    except Exception:
        return _err("verify failed")
    return json.dumps({"ok": True, "valid": bool(valid)}, ensure_ascii=False)
