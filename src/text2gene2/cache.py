"""Redis-backed cache with in-memory fallback for dev (no Redis required)."""
import json
import logging
from typing import Any

log = logging.getLogger(__name__)

# Module-level in-memory fallback used when Redis is unavailable
_mem_cache: dict[str, Any] = {}


async def _get_redis():
    """Return a Redis client, or None if Redis is unavailable."""
    try:
        import redis.asyncio as aioredis
        from text2gene2.config import settings
        client = aioredis.from_url(settings.redis_url, decode_responses=True)
        await client.ping()
        return client
    except Exception:
        return None


async def cache_get(key: str) -> Any | None:
    redis = await _get_redis()
    if redis:
        try:
            val = await redis.get(key)
            return json.loads(val) if val else None
        except Exception as e:
            log.warning("Redis get failed: %s", e)
    return _mem_cache.get(key)


async def cache_set(key: str, value: Any, ttl: int = 3600) -> None:
    redis = await _get_redis()
    if redis:
        try:
            await redis.set(key, json.dumps(value), ex=ttl)
            return
        except Exception as e:
            log.warning("Redis set failed: %s", e)
    _mem_cache[key] = value


async def cache_delete(key: str) -> None:
    redis = await _get_redis()
    if redis:
        try:
            await redis.delete(key)
        except Exception as e:
            log.warning("Redis delete failed: %s", e)
    _mem_cache.pop(key, None)
