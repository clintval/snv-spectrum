import urllib.request as request
from tempfile import NamedTemporaryFile
from typing import Any

__all__ = ['fetch_gff_from_uniprot']


def fetch_gff_from_uniprot(uniprot_id: str) -> Any:
    """Fetch GFF records from the Uniprot server."""
    from BCBio import GFF

    url = f'http://www.uniprot.org/uniprot/{uniprot_id}.gff'

    with request.urlopen(url) as handle:
        content = handle.read().decode('utf8')

    def remove_tabs_from_line_ends(content: str) -> str:
        formatted = []
        for line in content.split('\n'):
            line = line.strip()
            formatted.append(line)
        return '\n'.join(formatted)

    handle = NamedTemporaryFile(mode='w+', delete=False)  # type: ignore
    handle.write(remove_tabs_from_line_ends(content))  # type: ignore
    handle.close()  # type: ignore
    records = list(GFF.parse(handle.name))
    return records
