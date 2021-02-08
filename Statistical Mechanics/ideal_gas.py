# Evaluation of Gas chemical potential by the assumption of ideal gas and rigid rotor

import numpy as np

h = 6.626 * 10 ** (-34)
haba = h / 2 / np.pi
kb = 1.38 * 10 ** (-23)


def eval_chem_pot_gas(mass, coordination, frequency, sym, temp):
    """
    Evaluate chemical potential of non-linear gas molecule.

    :param mass: '(the number of atom in molecule) size' matrix for atom mass (unit: kg)
    :param coordination: '(the number of atom in molecule, 3) size' matrix for atom coordination (unit: m)
    :param frequency: normal frequency (unit: Hz)
    :param sym: the symmetry number of molecule
    :param temp: temperature (unit: K) matrix where chemical potential is evaluated
    :return: contributions (unit: eV) from translation, vibration, zero-point-energy, and rotation
    """
    # evaluate total mass of ideal gas
    m_total = 0
    for mass_atom in mass:
        m_total = m_total + mass_atom

    # evaluate principal moments of inertia
    x_cm = 0
    y_cm = 0
    z_cm = 0
    for ind_atom, mass_atom in enumerate(mass):
        x_cm = x_cm + mass_atom * coordination(ind_atom, 0) / m_total
        y_cm = y_cm + mass_atom * coordination(ind_atom, 1) / m_total
        z_cm = z_cm + mass_atom * coordination(ind_atom, 2) / m_total

    inertia = np.zeros(3, 3)
    for ind_atom, mass_atom in enumerate(mass):
        inertia[0, 0] = inertia[0, 0] + mass_atom * ((coordination(ind_atom, 1) - y_cm) ** 2
                                                     + (coordination(ind_atom, 2) - z_cm) ** 2)
        inertia[1, 1] = inertia[1, 1] + mass_atom * ((coordination(ind_atom, 0) - x_cm) ** 2
                                                     + (coordination(ind_atom, 2) - z_cm) ** 2)
        inertia[2, 2] = inertia[2, 2] + mass_atom * ((coordination(ind_atom, 0) - x_cm) ** 2
                                                     + (coordination(ind_atom, 1) - y_cm) ** 2)
        inertia[0, 1] = inertia[0, 1] \
                        + mass_atom * (coordination(ind_atom, 0) - x_cm) * (coordination(ind_atom, 1) - y_cm)
        inertia[0, 2] = inertia[0, 2] \
                        + mass_atom * (coordination(ind_atom, 0) - x_cm) * (coordination(ind_atom, 2) - z_cm)
        inertia[1, 2] = inertia[1, 2] \
                        + mass_atom * (coordination(ind_atom, 1) - y_cm) * (coordination(ind_atom, 2) - z_cm)

    inertia[1, 0] = inertia[0, 1]
    inertia[2, 0] = inertia[0, 2]
    inertia[2, 1] = inertia[1, 2]
    principal_inertia, _ = np.linalg.eig(inertia)

    # evaluate chemical potential of translational motion of ideal gas
    trans = -kb / 1.6 * 10 ** 19 * temp * np.log((2 * np.pi * m_total / h ** 2) ** 1.5 * (kb * temp) ** 2.5 / 101325)

    # evaluate chemical potential of vibrational motion of ideal gas
    vib = 0
    zpe = 0
    for freq in frequency:
        vib = vib + kb / 1.6 * 10 ** 19 * temp * np.log((1 - np.exp(-h * freq / kb / temp)))
        zpe = zpe + 1 / 2 * h * freq / 1.6 * 10 ** 19

    # evaluate chemical potential of rotational motion of ideal gas
    rot = -kb / 1.6 * 10 ** 19 * temp * np.log(np.pi ** 0.5 / sym
                                               * np.sqrt(2 * principal_inertia[0] * kb * temp / haba ** 2)
                                               * np.sqrt(2 * principal_inertia[1] * kb * temp / haba ** 2)
                                               * np.sqrt(2 * principal_inertia[2] * kb * temp / haba ** 2))

    return trans, vib, zpe, rot


def eval_chem_pot_gas_linear(mass, coordination, frequency, sym, temp):
    """
    Evaluate chemical potential of linear gas molecule.

    :param mass: '(the number of atom in molecule) size' matrix for atom mass (unit: kg)
    :param coordination: '(the number of atom in molecule, 3) size' matrix for atom coordination (unit: m)
    :param frequency: normal frequency (unit: Hz)
    :param sym: the symmetry number of molecule
    :param temp: temperature (unit: K) matrix where chemical potential is evaluated
    :return: contributions (unit: eV) from translation, vibration, zero-point-energy, and rotation
    """
    # evaluate total mass of ideal gas
    m_total = 0
    for mass_atom in mass:
        m_total = m_total + mass_atom

    # evaluate principal moments of inertia
    x_cm = 0
    y_cm = 0
    z_cm = 0
    for ind_atom, mass_atom in enumerate(mass):
        x_cm = x_cm + mass_atom * coordination(ind_atom, 0) / m_total
        y_cm = y_cm + mass_atom * coordination(ind_atom, 1) / m_total
        z_cm = z_cm + mass_atom * coordination(ind_atom, 2) / m_total

    inertia = np.zeros(3, 3)
    for ind_atom, mass_atom in enumerate(mass):
        inertia[0, 0] = inertia[0, 0] + mass_atom * ((coordination(ind_atom, 1) - y_cm) ** 2
                                                     + (coordination(ind_atom, 2) - z_cm) ** 2)
        inertia[1, 1] = inertia[1, 1] + mass_atom * ((coordination(ind_atom, 0) - x_cm) ** 2
                                                     + (coordination(ind_atom, 2) - z_cm) ** 2)
        inertia[2, 2] = inertia[2, 2] + mass_atom * ((coordination(ind_atom, 0) - x_cm) ** 2
                                                     + (coordination(ind_atom, 1) - y_cm) ** 2)
        inertia[0, 1] = inertia[0, 1] \
                        + mass_atom * (coordination(ind_atom, 0) - x_cm) * (coordination(ind_atom, 1) - y_cm)
        inertia[0, 2] = inertia[0, 2] \
                        + mass_atom * (coordination(ind_atom, 0) - x_cm) * (coordination(ind_atom, 2) - z_cm)
        inertia[1, 2] = inertia[1, 2] \
                        + mass_atom * (coordination(ind_atom, 1) - y_cm) * (coordination(ind_atom, 2) - z_cm)

    inertia[1, 0] = inertia[0, 1]
    inertia[2, 0] = inertia[0, 2]
    inertia[2, 1] = inertia[1, 2]
    principal_inertia, _ = np.linalg.eig(inertia)

    # evaluate chemical potential of translational motion of ideal gas
    trans = -kb / 1.6 * 10 ** 19 * temp * np.log((2 * np.pi * m_total / h ** 2) ** 1.5 * (kb * temp) ** 2.5 / 101325)

    # evaluate chemical potential of vibrational motion of ideal gas
    vib = 0
    zpe = 0
    for freq in frequency:
        vib = vib + kb / 1.6 * 10 ** 19 * temp * np.log((1 - np.exp(-h * freq / kb / temp)))
        zpe = zpe + 1 / 2 * h * freq / 1.6 * 10 ** 19

    # evaluate chemical potential of rotational motion of ideal gas
    rot = -kb / 1.6 * 10 ** 19 * temp * np.log(1 / sym * 2 * principal_inertia[2] * kb * temp / haba ** 2)

    return trans, vib, zpe, rot
