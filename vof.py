"""Convenience Python wrapper exposing the built extension and helpers.

This file re-exports the API from python/vof.py so tests and examples can
`import vof` from the repo root without installing a package.
"""
from python.vof import *  # noqa: F401,F403
