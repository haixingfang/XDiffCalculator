# Author: Haixing Fang
# Email: haixingfang868@gmail.com
"""Generate bundled tube references using optional build-only SpekPy 2.5.4.
Install spekpy==2.5.4 separately; never import this module from the application.
"""
import csv
import os
import spekpy


def main():
    """Generate binned photon densities at 100 cm per mAs, with fixed filtration."""
    if spekpy.__version__ != '2.5.4':
        raise RuntimeError('Use SpekPy 2.5.4 to reproduce these reference spectra.')
    output = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'spectra')
    os.makedirs(output, exist_ok=True)
    for target, voltages in [('W', (40, 80, 120)), ('Cr', (40, 50)), ('Cu', (40, 50))]:
        for voltage in voltages:
            model = spekpy.Spek(kvp=voltage, th=12, targ=target, physics='casim', dk=0.05, z=100, mas=1)
            model.filter('Be', 1.0)
            model.filter('Al', 1.0)
            energies, density = model.get_spectrum()
            path = os.path.join(output, '{}_{}kV.csv'.format(target, voltage))
            with open(path, 'w', newline='', encoding='utf-8') as stream:
                stream.write('# SpekPy 2.5.4, physics=casim, reflection target {} at {} kV\n'.format(target, voltage))
                stream.write('# Anode angle 12 deg; filtration 1 mm Be + 1 mm Al; on-axis distance 100 cm; 1 mAs\n')
                stream.write('# Photon density: photons/cm^2/mAs/keV; energy-bin width 0.05 keV\n')
                stream.write('# Line strengths are bin-averaged; interpolation does not model intrinsic line width.\n')
                stream.write('# https://pypi.org/project/spekpy/2.5.4/\n')
                writer = csv.writer(stream)
                writer.writerow(['energy_keV', 'photon_density'])
                # Half-bin edge ramps prevent extrapolation beyond the tube endpoint.
                writer.writerow([max(0.001, float(energies[0]) - 0.025), 0])
                writer.writerows(('{:.8g}'.format(e), '{:.10g}'.format(v)) for e, v in zip(energies, density))
                writer.writerow([voltage, 0])
            print(path)


if __name__ == '__main__':
    main()
