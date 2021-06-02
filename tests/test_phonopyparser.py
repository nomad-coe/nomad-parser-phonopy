#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import pytest
import numpy as np

from nomad.datamodel import EntryArchive
from phonopyparser import PhonopyParser


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


@pytest.fixture(scope='module')
def parser():
    return PhonopyParser()


def test_basic(parser):
    archive = EntryArchive()
    parser.parse('tests/data/Ge/phonopy-FHI-aims-displacement-01/control.in', archive, None)

    # need to assert values, no unbiased reference
    sec_thermo = archive.section_run[0].section_frame_sequence[0].section_thermodynamical_properties[0]
    assert len(sec_thermo.thermodynamical_property_temperature) == 11
    assert len(sec_thermo.thermodynamical_property_heat_capacity_C_v) == 11
    assert len(sec_thermo.vibrational_free_energy_at_constant_volume) == 11

    assert archive.section_run[0].section_method[0].x_phonopy_displacement.magnitude == 1e-12
    sec_scc = archive.section_run[0].section_single_configuration_calculation[0]
    assert np.shape(sec_scc.hessian_matrix) == (8, 8, 3, 3)
    assert np.shape(sec_scc.dos_phonon[0].total[0].value) == (201,)
    assert len(sec_scc.band_structure_phonon[0].band_structure_segment) == 10
