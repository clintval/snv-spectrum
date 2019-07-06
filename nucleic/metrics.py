import csv
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Dict, Iterable, List, TextIO, Tuple, Type

import attr

__all__ = ['Metric', 'MetricReader', 'MetricWriter', 'FgbioMetricReader', 'FgbioMetricWriter']


@attr.s(auto_attribs=True)
class Metric(ABC):
    """Abstract base class for any metric containing class.

    Examples:
        >>> import attr
        >>> @attr.s(auto_attribs=True)
        ... class CustomMetric(Metric):
        ...     count: int
        ...     rate: float
        >>> metric = CustomMetric(count = 2, rate = 0.2)

        >>> metric
        CustomMetric(count=2, rate=0.2)

        >>> metric.as_dict()
        {'count': 2, 'rate': 0.2}

        >>> metric.as_tuple()
        (2, 0.2)

        >>> CustomMetric.fields()
        ['count', 'rate']

    """

    def as_dict(self) -> Dict:
        """Transform this class to a dictionaary."""
        return attr.asdict(self)

    def as_tuple(self) -> Tuple:
        """Write the values of this metric to a line."""
        return attr.astuple(self)

    @classmethod
    def fields(cls) -> List[str]:
        """List the fields of this class."""
        fields: List[str] = [attribute.name for attribute in attr.fields(cls)]
        return fields


class MetricReader(csv.DictReader):
    """Read a MUT file."""

    def __init__(self, handle: TextIO, metric: Metric, delimiter = '\t') -> None:
        super().__init__(handle, delimiter = delimiter)
        self.metric: Metric = metric

    def __next__(self) -> Any:
        """Iterate through rows as dictionaries, then unpack into a :class:`Metric`."""
        item = super().__next__()
        record = self.metric(**item)
        return record

    @classmethod
    def read_path(cls: Type['MetricReader'], path: Path, metric: Metric) -> Dict[str, List[Metric]]:
        """Create a map of samples to :class:`Metric` from a file path."""
        with open(path, 'r') as handle:
            return list(cls(handle, metric))  # type: ignore  # error: "MetricReader" not callable


class MetricWriter(ABC):
    """Write metrics."""

    def __init__(self, handle: TextIO, metric: Metric) -> None:
        self.handle = handle
        self.metric = metric
        self.header = metric.fields()
        self._header_written: bool = False

    @property
    @abstractmethod
    def delimiter(self) -> str:
        """Text delimiter for output metric line."""
        return ','

    def write(self, record: Metric) -> None:
        """Write this record to the handle."""
        if self._header_written is False:
            self.handle.write(self.delimiter.join(self.header) + '\n')
            self._header_written = True
        self.handle.write(self.delimiter.join(map(str, record.as_tuple())) + '\n')
        return None

    def write_all(self, records: Iterable[Metric]) -> None:
        """Write all records to the handle."""
        for record in records:
            self.write(record)
        return None


class FgbioMetricWriter(MetricWriter):
    """Fulcrum Genomics `fgbio` metrics writer."""

    @property
    def delimiter(self) -> str:
        """Text delimiter for output metric line."""
        return '\t'
