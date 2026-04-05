"""
Cache operations — re-exported from db.py for backwards compatibility.

All cache logic now lives in db.py alongside the shared connection management.
Existing imports like `from text2gene2.cache import cache_get` continue to work.
"""
from text2gene2.db import cache_get, cache_set, cache_delete  # noqa: F401
