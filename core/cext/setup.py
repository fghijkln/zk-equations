"""Build the constant-time C extension (debug/experimental).

    cd core/cext && python3 setup.py build_ext --inplace
Produces core/cext/curve_ext.<abi>.so, imported by core/curve_c.py.
Pure C99, no dependencies beyond Python headers.
"""
from setuptools import setup, Extension

setup(
    name="curve_ext",
    version="0.1.0",
    ext_modules=[
        Extension(
            "curve_ext",
            sources=["curve_ext.c"],
            extra_compile_args=["-O2", "-std=c99", "-Wall"],
        )
    ],
    # C extension only; don't package the test file or __pycache__
    py_modules=[],
)
