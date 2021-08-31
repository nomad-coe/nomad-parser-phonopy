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
    sec_thermo = archive.run[0].calculation[0].thermodynamics
    assert len(sec_thermo) == 11
    assert sec_thermo[0].temperature is not None
    assert sec_thermo[0].heat_capacity_c_v is not None
    assert sec_thermo[0].vibrational_free_energy_at_constant_volume is not None

    assert archive.run[0].method[0].x_phonopy_displacement.magnitude == 1e-12
    sec_scc = archive.run[0].calculation[0]
    assert np.shape(sec_scc.hessian_matrix) == (8, 8, 3, 3)
    assert np.shape(sec_scc.dos_phonon[0].total[0].value) == (201,)
    assert len(sec_scc.band_structure_phonon[0].segment) == 10

    assert archive.run[0].calculation[0].system_ref is not None
    assert archive.run[0].calculation[0].method_ref is not None


def test_vasp(parser):
    archive = EntryArchive()
    parser.parse('tests/data/vasp/phonopy.yaml', archive, None)