"""
Shared database connections for text2gene2.

Provides a single lazy psycopg2 connection to the medgen Postgres DB on loki,
shared across all modules that need it (enrich, expand, and future local sources
like LOVD, local ClinVar, PubTator, etc.).

Also provides a shared async Redis client with connection reuse.

Usage:
    from text2gene2.db import get_medgen_conn, get_redis

    # Postgres (synchronous — wrap in asyncio.to_thread for async callers)
    conn = get_medgen_conn()
    if conn:
        with conn.cursor() as cur:
            cur.execute("SELECT ...")

    # Redis (async)
    redis = await get_redis()
    if redis:
        await redis.get("key")
"""
import json
import logging
from typing import Any

from text2gene2.config import settings

log = logging.getLogger(__name__)

# ── PostgreSQL ───────────────────────────────────────────────────────────────

_pg_conn = None


def get_medgen_conn():
    """
    Return a shared psycopg2 connection to the medgen DB, or None if unavailable.

    The connection is readonly with autocommit — safe for concurrent reads from
    multiple modules. Created on first call, reused thereafter. If a query fails
    and the caller sets the connection to a broken state, the next call will
    reconnect automatically.
    """
    global _pg_conn
    if not settings.medgen_db_url:
        return None
    try:
        import psycopg2
        if _pg_conn is None or _pg_conn.closed:
            _pg_conn = psycopg2.connect(settings.medgen_db_url)
            _pg_conn.set_session(readonly=True, autocommit=True)
            log.info("Connected to medgen DB at %s", settings.medgen_db_url.split("@")[-1])
        return _pg_conn
    except Exception as e:
        log.warning("medgen DB connect failed: %s", e)
        _pg_conn = None
        return None


def reset_medgen_conn():
    """Force reconnect on next get_medgen_conn() call. Use after a query error."""
    global _pg_conn
    if _pg_conn is not None:
        try:
            _pg_conn.close()
        except Exception:
            pass
    _pg_conn = None


# ── Redis ────────────────────────────────────────────────────────────────────

_redis_client = None
_redis_unavailable = False


async def get_redis():
    """
    Return a shared async Redis client, or None if Redis is unavailable.

    Connection is created and verified on first call, then reused. If Redis
    was previously unreachable, retries are suppressed for the lifetime of the
    process to avoid repeated connection timeouts on every cache operation.
    """
    global _redis_client, _redis_unavailable
    if _redis_unavailable:
        return None
    if _redis_client is not None:
        return _redis_client
    try:
        import redis.asyncio as aioredis
        client = aioredis.from_url(settings.redis_url, decode_responses=True)
        await client.ping()
        _redis_client = client
        return _redis_client
    except Exception:
        _redis_unavailable = True
        return None


# ── Cache helpers (convenience wrappers) ─────────────────────────────────────

# In-memory fallback when Redis is unavailable
_mem_cache: dict[str, Any] = {}


async def cache_get(key: str) -> Any | None:
    redis = await get_redis()
    if redis:
        try:
            val = await redis.get(key)
            return json.loads(val) if val else None
        except Exception as e:
            log.warning("Redis get failed: %s", e)
    return _mem_cache.get(key)


async def cache_set(key: str, value: Any, ttl: int = 3600) -> None:
    redis = await get_redis()
    if redis:
        try:
            await redis.set(key, json.dumps(value), ex=ttl)
            return
        except Exception as e:
            log.warning("Redis set failed: %s", e)
    _mem_cache[key] = value


async def cache_delete(key: str) -> None:
    redis = await get_redis()
    if redis:
        try:
            await redis.delete(key)
        except Exception as e:
            log.warning("Redis delete failed: %s", e)
    _mem_cache.pop(key, None)
