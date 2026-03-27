"""Simple async token-bucket rate limiters, one per external service."""
import asyncio
import time


class RateLimiter:
    """Token-bucket rate limiter for async code."""

    def __init__(self, rate: float):
        """
        Args:
            rate: maximum requests per second
        """
        self._rate = rate
        self._min_interval = 1.0 / rate
        self._last = 0.0
        self._lock = asyncio.Lock()

    async def acquire(self) -> None:
        async with self._lock:
            now = time.monotonic()
            wait = self._min_interval - (now - self._last)
            if wait > 0:
                await asyncio.sleep(wait)
            self._last = time.monotonic()


# One limiter per external service.  Instantiated at import time so they're
# shared across all requests within a process.
variant_validator = RateLimiter(rate=1.5)   # public API asks ≤ 2/sec
litvar2 = RateLimiter(rate=5.0)             # no documented limit; be polite
clinvar = RateLimiter(rate=3.0)             # 3/sec without NCBI key, 10 with
google_cse = RateLimiter(rate=1.0)          # 100/day free = ~0.001/sec; burst ok
europepmc = RateLimiter(rate=5.0)           # generous public API
