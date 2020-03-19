import numpy as np            # pylint: disable=unused-import
import typing                 # pylint: disable=unused-import
from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference
)
from nomad.metainfo.legacy import LegacyDefinition

from nomad.datamodel.metainfo import public

m_package = Package(
    name='phonopy_nomadmetainfo_json',
    description='None',
    a_legacy=LegacyDefinition(name='phonopy.nomadmetainfo.json'))


class x_phonopy_input(MCategory):
    '''
    Information about properties that concern phonopy calculations.
    '''

    m_def = Category(
        a_legacy=LegacyDefinition(name='x_phonopy_input'))


class section_method(public.section_method):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_method'))

    x_phonopy_displacement = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter',
        description='''
        Amplitude of the atom diplacement for the phonopy supercell
        ''',
        categories=[x_phonopy_input],
        a_legacy=LegacyDefinition(name='x_phonopy_displacement'))

    x_phonopy_symprec = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='meter',
        description='''
        Symmetry threshold for the space group identification of the crystal for which the
        vibrational properties are to be calculated
        ''',
        categories=[x_phonopy_input],
        a_legacy=LegacyDefinition(name='x_phonopy_symprec'))


class section_system(public.section_system):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_system'))

    x_phonopy_original_system_ref = Quantity(
        type=public.section_system,
        shape=[],
        description='''
        Original cell from which the supercell for the DFT calculations was constructed
        ''',
        a_legacy=LegacyDefinition(name='x_phonopy_original_system_ref'))


m_package.__init_metainfo__()
