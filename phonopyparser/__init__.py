# Copyright 2015-2018 Lauri Himanen, Fawzi Mohamed, Ankit Kariryaa
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

from phonopyparser.phonopy_calculators import PhonopyCalculatorInterface

from nomad.parsing.parser import FairdiParser
from .metainfo import m_env


class PhonopyParser(FairdiParser):
    def __init__(self):
        super().__init__(
            name='parsers/phonopy', code_name='Phonopy', code_homepage='https://phonopy.github.io/phonopy/',
            mainfile_name_re=(r'(.*/phonopy-FHI-aims-displacement-0*1/control.in$)|(.*/phonon.yaml)')
        )

    def parse(self, filepath, archive, logger=None):
        self._metainfo_env = m_env

        interface = PhonopyCalculatorInterface(filepath, archive, logger)
        interface.parse()
