import csv
from pathlib import Path
from typing import Any, Iterable, List, Optional, TextIO, Union

import attr

from pyfaidx import Fasta

from nucleic import Dna
from nucleic.io.fasta import subseq


@attr.s(auto_attribs=True, slots=True)
class Interval(object):
    """An interval which can be placed on a contig within a reference."""

    contig: Optional[str] = None
    start: Optional[int] = None
    end: Optional[int] = None
    name: Optional[str] = None

    @classmethod
    def fields(cls) -> List[str]:
        """List the fields of this class."""
        fields: List[str] = [attribute.name for attribute in attr.fields(cls)]
        return fields

    @classmethod
    def from_iterable(cls, iterable: Iterable) -> 'Interval':
        """Build this interval from a list."""
        vector = list(iterable)
        if len(vector) >= 4:
            contig, start, end, name, *_ = iterable
            return cls(str(contig), int(start), int(end), str(name))
        elif len(vector) == 3:
            contig, start, end = iterable
            return cls(str(contig), int(start), int(end))

    def is_closed(self) -> bool:
        """Test if this is a closed interval."""
        return self.start is not None and self.end is not None

    def is_placed(self) -> bool:
        """Test if this interval is located on a contig."""
        return self.contig is not None

    def reference_seq(self, reference: Union[Fasta, Path]) -> Dna:
        """Get the reference sequence of this interval from a FASTA."""
        if not self.is_closed():
            raise ValueError(f"Interval must have start and end coordinates: {self}")
        if not self.is_placed():
            raise ValueError(f"Interval have `contig` defined: {self}")
        self.contig
        return subseq(
            reference, self.contig, self.start, self.end
        )  # type: ignore  # MyPy has no clue.


class IntervalList(list):
    """A container of intervals."""

    def __init__(self, *args: Any) -> None:
        super().__init__(args)

    @classmethod
    def read_bed(cls, handle: TextIO) -> 'IntervalList':
        """Read a BED filehandle into an interval list."""
        reader = csv.reader(handle, delimiter='\t')
        return cls(*(Interval.from_iterable(line) for line in reader))

    def __repr__(self) -> str:
        items = ', '.join([str(item) for item in self])
        return f'{self.__class__.__qualname__}({items})'
