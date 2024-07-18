#  MacSyFinder - Detection of macromolecular systems in protein dataset  #
#                using systems modelling and similarity search.          #
#  Authors: Sophie Abby, Bertrand Neron                                  #
#  Copyright (c) 2014-2024  Institut Pasteur (Paris) and CNRS.           #
#  See the COPYRIGHT file for details                                    #
#                                                                        #
#  This file is part of MacSyFinder package.                             #
#                                                                        #
#  MacSyFinder is free software: you can redistribute it and/or modify   #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  (at your option) any later version.                                   #
#                                                                        #
#  MacSyFinder is distributed in the hope that it will be useful,        #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
#  GNU General Public License for more details .                         #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with MacSyFinder (COPYING).                                     #
#  If not, see <https://www.gnu.org/licenses/>.                          #
##########################################################################
from __future__ import  annotations

import re

import yaml


class Maintainer:

    def __init__(self, name: str, email: str) -> None:
        self.name = name
        self.email = email


    def to_dict(self) -> dict[str, str]:
        return {'name': self.name , 'email': self.email}

    def __eq__(self, other):
        return self.name == other.name and self.email == other.email


class Metadata:
    """
    Handle package metadata
    """

    name = 'metadata.yml'

    def __init__(self, maintainer: Maintainer, short_desc: str) -> None:
        self.maintainer = maintainer
        self.short_desc = short_desc
        self.cite = []
        self.vers = None
        self.doc = ''
        self.license = ''
        self.copyright_date = ''
        self.copyright_holder = ''


    @staticmethod
    def load(path: str) -> Metadata:
        """
        Create a :class:`Metadata` object from a metadata file

        :param path: the path to the metadatafile in yaml format
        """
        with open(path) as raw_metadata:
            data = yaml.safe_load(raw_metadata)
            msgs = []
            try:
                maintainer = Maintainer(**data['maintainer'])
                new_meta = Metadata(maintainer,
                                    data['short_desc'])
            except (TypeError, KeyError):
                if 'short_desc' not in data:
                    msgs.append(f"The metadata file '{path}' is not valid: the element 'short_desc' "
                                "is required.")
                if 'maintainer' not in data:
                    msgs.append(f"The metadata file '{path}' is not valid: the element 'maintainer' "
                                "is required.")
                elif 'name' not in data['maintainer'] or 'email' not in data['maintainer']:
                    msgs.append(f"The metadata file '{path}' is not valid: "
                                f"the element 'maintainer' must have fields 'name' and 'email'.")
                raise ValueError('\n- ' + '\n- '.join(msgs)) from None

            for attr in 'cite', 'vers', 'doc', 'license', 'copyright':
                if attr in data:
                    if attr == 'copyright':
                        m = re.match(r'(\d+)([-, ]*)(\d*)(.*)', data['copyright'].strip())
                        if m:
                            date_from, date_sep, date_to, holder = m.groups()
                            if date_to:
                                new_meta.copyright_date = f"{date_from}{date_sep}{date_to}"
                            elif date_from:
                                new_meta.copyright_date = f"{date_from}"
                            new_meta.copyright_holder = holder.strip(' ,;.')
                    else:
                        setattr(new_meta, attr, data[attr])
        return new_meta


    def save(self, path: str) -> None:
        """
        Serialize this object in a metadata file

        :param path: The path to save the metadata
        """
        yaml_dict = dict()
        yaml_dict['maintainer'] = self.maintainer.to_dict()
        yaml_dict['short_desc'] = self.short_desc
        if self.cite:
            yaml_dict['cite'] = self.cite
        if self.doc:
            yaml_dict['doc'] = self.doc
        if self.license:
            yaml_dict['license'] = self.license
        if self._copyright_holder:
            yaml_dict['copyright'] = f"{self.copyright_date}, {self.copyright_holder}"
        if self.vers:
            yaml_dict['vers'] = self.vers

        with open(path, 'w') as metafile:
            yaml.dump(yaml_dict, metafile, allow_unicode=True, indent=2)

    @property
    def maintainer(self) -> Maintainer:
        """
        :return: dict with 2 keys
        :rtype:
        """
        return self._maintainer


    @maintainer.setter
    def maintainer(self, maintainer: Maintainer) -> None:
        self._maintainer = maintainer


    @property
    def short_desc(self) -> str:
        return self._short_desc


    @short_desc.setter
    def short_desc(self, desc: str) -> None:
        if not desc:
            raise ValueError("The field 'short_desc' is mandatory.")
        self._short_desc = desc.replace('\n', ' ')


    @property
    def cite(self) -> list[str]:
        return self._cite

    @cite.setter
    def cite(self, citations: list[str]) -> None:
        if citations:
            self._cite = citations
        else:
            self._cite = []

    @property
    def doc(self) -> str:
        return self._doc

    @doc.setter
    def doc(self, doc_link: str) -> None:
        self._doc = doc_link


    @property
    def vers(self) -> str | None:
        return self._vers


    @vers.setter
    def vers(self, vers: str | None) -> None:
        if vers:
            self._vers = str(vers)
        else:
            self._vers = None


    @property
    def license(self) -> str:
        return self._license


    @license.setter
    def license(self, license_val: str) -> None:
        self._license = license_val


    @property
    def copyright_date(self) -> str:
        return self._copyright_date


    @copyright_date.setter
    def copyright_date(self, value: str):
        self._copyright_date = str(value)


    @property
    def copyright_holder(self) -> str:
        return self._copyright_holder


    @copyright_holder.setter
    def copyright_holder(self, value: str):
        self._copyright_holder = value

    @property
    def copyright(self) -> str:
        if self.copyright_holder:
            return f"{self.copyright_date}, {self.copyright_holder}"
        else:
            return ''
