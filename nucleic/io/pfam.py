import json
import urllib.request as request
from typing import Dict, List, Mapping, Optional, Tuple

__all__ = ['PfamGraphicFeature', 'PfamGraphicResponse']


class PfamGraphicFeature(object):
    """An object representing a Pfam graphic.
    Args:
        feature: A dictionary of graphic feature options.
    """

    def __init__(self, feature: Mapping) -> None:
        self.color: str = feature.get('colour', 'grey')
        self.display: Optional[str] = feature.get('display', None)
        self.end: Optional[int] = feature.get('end', None)
        self.endstyle: Optional[str] = feature.get('endStyle', None)
        self.link: Optional[str] = feature.get('href', None)
        self.start: Optional[int] = feature.get('start', None)
        self.startstyle: Optional[str] = feature.get('startStyle', None)
        self.text: Optional[str] = feature.get('text', '')
        self.type: Optional[str] = feature.get('type', None)
        self.metadata: Dict = feature.get('metadata', dict())

    def __repr__(self) -> str:
        return (
            f'PfamGraphicFeature('
            f'start={self.start} '
            f'end={self.end} '
            f'color="{self.color}" '
            f'link="{self.link}")'
        )


class PfamGraphicResponse(object):
    """An object representing a Pfam graphic server response."""

    def __init__(self, content: Mapping) -> None:
        self.regions: List[PfamGraphicFeature] = []

        self.length: Optional[int] = content.get('length', None)
        self.markups: Tuple[str] = content.get('markups', tuple())
        self.metadata: Dict = content.get('metadata', dict())
        self.motifs: Tuple[str] = content.get('motifs', tuple())

        for feature in content.get('regions', tuple()):
            self.regions.append(PfamGraphicFeature(feature))

    def __repr__(self) -> str:
        return (
            f'{self.__class__.__name__}(\n'
            f'  accession:    "{self.metadata.get("accession", "N/A")}"\n'
            f'  identifier:   "{self.metadata.get("identifier", "N/A")}"\n'
            f'  organism:     "{self.metadata.get("organism", "N/A")}")\n'
            f'  description:  "{self.metadata.get("description", "N/A")}"\n'
            f'  number of motifs:  {len(self.motifs)}\n'
            f'  number of regions: {len(self.regions)})'
        )


def fetch_pfam_graphics(pfam_id: str) -> List[PfamGraphicResponse]:
    url = f'http://pfam.xfam.org/protein/{pfam_id}/graphic'
    with request.urlopen(url) as handle:
        responses = list(map(PfamGraphicResponse, json.load(handle)))
    return responses
