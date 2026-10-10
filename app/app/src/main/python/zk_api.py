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

    witness: decimal string (up to 12 places); an int is also accepted.
    precision: k as string, MANDATORY for every equation (final precision
    1/k: witness must be within 1/k of an exact root).
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

    precision: k as string, MANDATORY (must match the k used at prove
    time). Never taken from the proof itself.
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


def ode_coefficients(ode_body, deg):
    """Get Taylor coefficients for an ODE.

    ode_body: e.g. "y' = y, y(0) = 1" (without the ode() wrapper).
    deg: degree as string, e.g. "5".
    Returns {"ok": true, "coeffs": ["1", "1", "1/2", ...]}.
    """
    try:
        eq_str = "ode(%s, deg = %s)" % (ode_body.strip(), deg.strip())
        coeffs = zk.ode_coefficients(eq_str)
        strs = [str(c.numerator) if c.denominator == 1
                else "%d/%d" % (c.numerator, c.denominator)
                for c in coeffs]
    except ValueError as e:
        return _err(str(e))
    except Exception:
        return _err("ode_coefficients failed: " + traceback.format_exc(limit=1).strip())
    return json.dumps({"ok": True, "coeffs": strs}, ensure_ascii=False)


def prove_ode(ode_body, deg):
    """Prove an ODE: ode_body like "y' = y, y(0) = 1", deg like "5".
    Witness and k are auto-handled (not needed).
    """
    try:
        eq_str = "ode(%s, deg = %s)" % (ode_body.strip(), deg.strip())
        proof = zk.prove_equation(eq_str, "", "")
    except ValueError as e:
        return _err(str(e))
    except Exception:
        return _err("prove failed: " + traceback.format_exc(limit=1).strip())
    return json.dumps({"ok": True, "proof": proof}, ensure_ascii=False)


def verify_ode(ode_body, deg, proof_json):
    """Verify an ODE proof. deg may be coarser (smaller) than prove-time."""
    try:
        proof = json.loads(proof_json)
    except (ValueError, TypeError):
        return _err("proof is not valid JSON")
    try:
        eq_str = "ode(%s, deg = %s)" % (ode_body.strip(), deg.strip())
        valid = zk.verify_equation(eq_str, proof, "")
    except Exception:
        return _err("verify failed")
    return json.dumps({"ok": True, "valid": bool(valid)}, ensure_ascii=False)
