"""Abstract base class for all PMID sources."""
from abc import ABC, abstractmethod
from text2gene2.models import LVGResult, Source, SourceResult


class PMIDSource(ABC):
    source: Source

    @abstractmethod
    async def query(self, lvg: LVGResult) -> SourceResult:
        """Return PMIDs relevant to this variant from this source."""
        ...
