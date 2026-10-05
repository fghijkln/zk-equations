"""Fiat-Shamir transcript: turns the interactive Bulletproofs protocol
non-interactive by deriving every verifier challenge as SHA-256 of the
transcript so far (random-oracle model, cf. paper Section 4.4).

Domain separation prefix for every hash: b"zkeq-bp-v1".
"""

import hashlib

from . import curve_c as curve

DOMAIN = b"zkeq-bp-v1"


class Transcript:
    def __init__(self):
        self._buf = bytearray()

    def _absorb(self, label, data):
        self._buf += DOMAIN + b"|" + label + b"|" + data + b";"

    def append_point(self, label, Pt):
        """Append a curve point in SEC1 compressed form."""
        self._absorb(label.encode(), curve.compress(Pt))

    def append_scalar(self, label, s):
        """Append a scalar as 32-byte big-endian."""
        self._absorb(label.encode(), (s % curve.N).to_bytes(32, "big"))

    def append_str(self, label, text):
        self._absorb(label.encode(), text.encode())

    def append_int(self, label, v):
        self._absorb(label.encode(), str(v).encode())

    def challenge(self, label):
        """Squeeze a nonzero scalar challenge in Z_n^*.

        Rejection-resample with a counter if the digest is 0 mod n
        (probability ~2^-256, but let's be correct).
        """
        ctr = 0
        while True:
            d = hashlib.sha256(
                self._buf + DOMAIN + b"|chal|" + label.encode() + b"|" + bytes([ctr])
            ).digest()
            c = int.from_bytes(d, "big") % curve.N
            if c != 0:
                # Bind the challenge into the transcript so later challenges
                # depend on it.
                self._absorb(b"chal:" + label.encode(), d)
                return c
            ctr += 1
