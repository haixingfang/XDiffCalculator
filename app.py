# -*- coding: utf-8 -*-
# Author: Haixing Fang
# Email: haixingfang868@gmail.com
"""Lightweight offline X-ray powder and polycrystal diffraction calculator."""
from __future__ import division

import argparse
import csv
import json
import math
import os
import heapq
import queue
import threading
import tempfile
import traceback
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from snapshot import canvas_snapshot

from diffraction import cif_blocks, crystal_metadata, load_crystal, reflections, powder_lines, simulate_spots, polarization_factor
from geometry import HC_KEV_ANGSTROM, detector_geometry, maximum_two_theta, pixel_lab, ray_detector_hit, clip_detector_segment
from orientations import random_orientations, load_orientations, save_orientations
from orientations import euler_to_matrix, matrix_to_euler
from laue import laue_reflections, simulate_laue_spots, weight_laue_spots, group_laue_spots
from spectrum import source_spectrum, evaluate_spectrum, TUBE_VOLTAGES

ROOT = os.path.dirname(os.path.abspath(__file__))
DEFAULTS = dict(cif=os.path.join(ROOT, 'examples', 'Ni.cif'), block='Ni', mode='Powder rings',
                orientation_file='', number='200', seed='200', energy='40', distance='150', polarization='0.95',
                preset='EIGER2 X CdTe 4M', columns='2068', rows='2162', pixel_y='75', pixel_z='75',
                center_col='1033.5', center_row='1080.5', tilt_x='0', tilt_y='0', tilt_z='0', max_tth='40', start='0', end='2', step='0.25', minimum='0.1')
DEFAULTS.update(beam_mode='Monochromatic', energy_min='5', energy_max='80', omega='0', source='W', tube_voltage='120',
                spectrum_file='', laue_polarization='0', colour='Relative weight', laue_weight='Spectral',
                single_euler='269.178546842,27.3683882123,93.4990822005')


def example_cifs(directory=None):
    """List local example CIF filenames, case-insensitively, without parsing them."""
    directory = os.path.join(ROOT, 'examples') if directory is None else directory
    if not os.path.isdir(directory):
        return []
    return sorted([name for name in os.listdir(directory) if name.lower().endswith('.cif')
                   and os.path.isfile(os.path.join(directory, name))], key=str.lower)


def settings_for_disk(settings):
    """Preserve bundled example references across single-file extraction directories."""
    result = dict(settings)
    path = os.path.abspath(result['cif'])
    if os.path.dirname(path) == os.path.join(ROOT, 'examples'):
        result['cif'] = '@example/' + os.path.basename(path)
    return result


def settings_from_disk(settings):
    """Resolve explicitly marked bundled examples without changing external paths."""
    result = dict(settings)
    if result.get('mode') == 'Single orientation':
        result['mode'] = 'Multiple orientations'
    if 'single_u' in result:
        if 'single_euler' not in result:
            values = [float(value) for value in str(result['single_u']).replace(',', ' ').split()]
            if len(values) != 9:
                raise ValueError('Legacy single_u must contain nine matrix elements.')
            angles = matrix_to_euler([values[0:3], values[3:6], values[6:9]])
            result['single_euler'] = ','.join('{:.15g}'.format(value) for value in angles)
        del result['single_u']
    if 'single_euler' in result:
        multiple_orientations(result['single_euler'])
    value = result.get('cif', '')
    if isinstance(value, str) and value.startswith('@example/'):
        name = value[len('@example/'):]
        if name != os.path.basename(name) or name in ('', '.', '..'):
            raise ValueError('Invalid bundled-example reference.')
        result['cif'] = os.path.join(ROOT, 'examples', name)
    return result


def calculate(settings, cancel=None, progress=None):
    """Calculate one immutable input snapshot, without touching Tk widgets.
    Args:
        settings: String-valued input dictionary, following DEFAULTS.
        cancel: Optional threading.Event.
        progress: Optional completed/total grain callback.
    Returns:
        result: Detector, crystal summary, reflections, powder lines and spot events.
    """
    settings = dict(DEFAULTS, **settings_from_disk(settings))
    if settings['beam_mode'] == 'Polychromatic (Laue)':
        return calculate_laue(settings, cancel, progress)
    if settings['beam_mode'] != 'Monochromatic':
        raise ValueError('Unknown beam mode.')
    detector = detector_geometry(settings['distance'], settings['columns'], settings['rows'], settings['pixel_y'],
                                 settings['pixel_z'], settings['center_col'], settings['center_row'],
                                 settings.get('tilt_x', 0), settings.get('tilt_y', 0), settings.get('tilt_z', 0))
    energy = float(settings['energy'])
    maximum = float(settings['max_tth'])
    if not math.isfinite(maximum) or not 0 < maximum < 90:
        raise ValueError('Maximum 2theta must lie between 0 and 90 degrees.')
    minimum = float(settings['minimum'])
    if not math.isfinite(minimum) or not 0 <= minimum <= 100:
        raise ValueError('Minimum relative intensity must lie in [0,100].')
    structure = load_crystal(settings['cif'], settings['block'] or None)
    records = reflections(structure, energy, min(maximum, maximum_two_theta(detector)), cancel)
    value = str(settings.get('polarization', '0.95')).strip()
    polarization = None if not value else float(value)
    polarization_factor(0, polarization=polarization)
    lines = powder_lines(records, polarization)
    matrices, events = [], []
    if settings['mode'] == 'Random grains':
        matrices = random_orientations(int(settings['number']), int(settings['seed']))
    elif settings['mode'] == 'Imported orientations':
        matrices = load_orientations(settings['orientation_file'])
    elif settings['mode'] == 'Multiple orientations':
        matrices = multiple_orientations(settings['single_euler'])
    elif settings['mode'] != 'Powder rings':
        raise ValueError('Unknown sample mode.')
    if matrices:
        events = simulate_spots(records, matrices, energy, detector, settings['start'], settings['end'], settings['step'], cancel, progress, polarization)
    if cancel is not None and cancel.is_set():
        raise RuntimeError('Calculation cancelled.')
    if matrices:
        intensity_max = max([event['intensity'] for event in events] or [1])
    else:
        intensity_max = max([line['sum_f2']*line['lorentz']*max(polarization_factor(line['two_theta'], 0, polarization),
                            polarization_factor(line['two_theta'], math.pi/2, polarization)) for line in lines] or [1])
    return dict(settings=dict(settings), detector=detector, reflections=records, lines=lines, matrices=matrices, events=events,
                polarization=polarization, intensity_max=intensity_max,
                cell=list(structure.cell.parameters), spacegroup=structure.spacegroup.xhm(), wavelength=HC_KEV_ANGSTROM/energy)


def single_orientation(text):
    """Read Bunge phi1/Phi/phi2 in degrees and return crystal-to-sample U."""
    values = [float(value) for value in text.replace(',', ' ').split()]
    return euler_to_matrix(values)


def multiple_orientations(text):
    """Read one to ten semicolon-separated Bunge Euler triples in grain-ID order."""
    rows = text.split(';')
    if not 1 <= len(rows) <= 10:
        raise ValueError('Enter between 1 and 10 orientations.')
    matrices = []
    for index, row in enumerate(rows):
        try:
            matrices.append(single_orientation(row))
        except (ValueError, TypeError) as error:
            raise ValueError('Orientation {}: {}'.format(index, error))
    return matrices


def calculate_laue(settings, cancel=None, progress=None):
    """Calculate an origin-centred fixed-omega Laue snapshot with approximate weights."""
    minimum, low, high, omega, angle = [float(settings[key]) for key in ('minimum', 'energy_min', 'energy_max', 'omega', 'max_tth')]
    if not all(math.isfinite(v) for v in (minimum, low, high, omega, angle)) or not 0 <= minimum <= 100 or not 0 < low < high or not 0 < angle < 90:
        raise ValueError('Require 0 < Emin < Emax, finite omega, 0 < max 2theta < 90, and display threshold in [0,100].')
    if settings['mode'] == 'Powder rings':
        raise ValueError('Laue mode requires Multiple orientations, Random grains or Imported orientations; white-beam powder is not implemented.')
    if settings['mode'] in ('Multiple orientations', 'Single orientation'):
        matrices = multiple_orientations(settings['single_euler'])
    elif settings['mode'] == 'Random grains':
        matrices = random_orientations(int(settings['number']), int(settings['seed']))
    elif settings['mode'] == 'Imported orientations':
        matrices = load_orientations(settings['orientation_file'])
    else:
        raise ValueError('Unknown sample mode.')
    if not matrices:
        raise ValueError('At least one grain orientation is required.')
    detector = detector_geometry(settings['distance'], settings['columns'], settings['rows'], settings['pixel_y'], settings['pixel_z'],
                                 settings['center_col'], settings['center_row'], settings['tilt_x'], settings['tilt_y'], settings['tilt_z'])
    spectrum = source_spectrum(settings)
    polar = str(settings['laue_polarization']).strip()
    polar = None if not polar else float(polar)
    polarization_factor(0, polarization=polar)
    if settings['laue_weight'] not in ('Spectral', 'Uniform') or settings['colour'] not in ('Relative weight', 'Energy'):
        raise ValueError('Unknown Laue weighting or colour mode.')
    structure = load_crystal(settings['cif'], settings['block'] or None)
    records = laue_reflections(settings['cif'], settings['block'] or None, high, min(angle, maximum_two_theta(detector)))
    events = simulate_laue_spots(records, matrices, omega, detector, low, high, angle, cancel, progress)
    events = weight_laue_spots(events, spectrum, polar)
    groups = group_laue_spots(events)
    if settings['laue_weight'] == 'Uniform':
        for group in groups:
            group['intensity'] = 1.0
    maximum = max([g['intensity'] for g in groups] or [0]) or 1.0
    return dict(settings=dict(settings), detector=detector, reflections=records, lines=[], matrices=matrices, events=events, groups=groups,
                polarization=polar, intensity_max=maximum, cell=list(structure.cell.parameters), spacegroup=structure.spacegroup.xhm(),
                wavelength=None, spectrum=spectrum, laue=True)


def export_laue(result, base):
    """Write every Laue contribution and a self-contained spectrum/settings provenance file."""
    keys = ('grain', 'omega', 'energy_keV', 'wavelength_A', 'd', 'two_theta', 'azimuth', 'col', 'row', 'f2', 'polarization', 'spectral_density', 'spectral_weight')
    with open(base + '_laue.csv', 'w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['h', 'k', 'l'] + list(keys))
        for event in result['events']:
            writer.writerow(list(event['hkl']) + [event[key] for key in keys])
    with open(base + '_laue_settings.json', 'w', encoding='utf-8') as stream:
        json.dump(dict(settings=settings_for_disk(result['settings']), spectrum=result['spectrum'], matrices=result['matrices'],
                       euler_bunge_deg=[matrix_to_euler(matrix) for matrix in result['matrices']],
                       orientation_convention='Bunge ZXZ degrees; U = Rz(phi1) Rx(Phi) Rz(phi2), crystal to sample',
                       weight_model='S_E * F_squared * P; no Laue Jacobian, absorption, anomalous scattering or detector response'), stream, indent=2)


def intensity_color(fraction):
    """Map a linear intensity fraction to the same blue-cyan-yellow scale everywhere."""
    value = max(0, min(1, fraction))
    stops = [(30, 55, 145), (20, 185, 200), (255, 240, 70)]
    index = min(1, int(value*2))
    weight = value*2 - index
    rgb = [int(round(a + weight*(b - a))) for a, b in zip(stops[index], stops[index + 1])]
    return '#{:02x}{:02x}{:02x}'.format(*rgb)


def ring_segments(line, detector, polarization, maximum, guide=False):
    """Intersect a diffraction cone with the tilted detector and yield coloured segments."""
    angle = math.radians(line['two_theta'])
    cosine, sine = math.cos(angle), math.sin(angle)
    first = ray_detector_hit((cosine, sine, 0), detector, clip=False)
    for index in range(360):
        angle0, angle1 = math.radians(index), math.radians(index + 1)
        last = ray_detector_hit((cosine, sine*math.cos(angle1), sine*math.sin(angle1)), detector, clip=False)
        segment = clip_detector_segment(first, last, detector) if first is not None and last is not None else None
        first = last
        if segment is None:
            continue
        intensity = line['sum_f2']*line['lorentz']*polarization_factor(line['two_theta'], (angle0 + angle1)/2, polarization)
        color = '#46798e' if guide else intensity_color(intensity/maximum)
        yield segment + (color,)


class Application:
    """One-window desktop interface; expensive calculations run in a worker thread."""

    def __init__(self, root):
        """Create controls, detector canvas and reflection table."""
        self.root, self.result = root, None
        self.messages, self.cancel_event = queue.Queue(), threading.Event()
        self.worker, self.selected_line, self.drawn_spots = None, None, []
        self.hover_items, self.hover_xy = {}, None
        self.busy = False
        self.zoom = 1.0
        self.view_center = (0.5, 0.5)
        self.pan_anchor = None
        self.drag_start = None
        self.drag_moved = False
        self.rectangle_mode = tk.BooleanVar(value=False)
        self.laue_labels = tk.BooleanVar(value=False)
        self.metadata_job = None
        self.vars = dict((key, tk.StringVar(value=value)) for key, value in DEFAULTS.items())
        self.status = tk.StringVar(value='Choose a CIF and calculate. All grains are at the origin.')
        self.integrated = tk.BooleanVar(value=True)
        self.overlay = tk.BooleanVar(value=True)
        self.frame = tk.IntVar(value=0)
        root.title('XDiffCalculator | powder rings & polycrystal spots')
        self.logo = tk.PhotoImage(file=os.path.join(ROOT, 'assets', 'icon.png'))
        root.iconphoto(True, self.logo)
        root.geometry('1260x{}'.format(max(680, min(1000, root.winfo_screenheight() - 80))))
        root.minsize(1000, 680)
        root.protocol('WM_DELETE_WINDOW', self.close)
        with open(os.path.join(ROOT, 'detectors.json'), encoding='utf-8') as stream:
            self.presets = json.load(stream)

        background = ttk.Style(root).lookup('TFrame', 'background') or root.cget('background')
        self.workspace_panes = tk.PanedWindow(root, orient='horizontal', sashwidth=7, sashrelief='raised', borderwidth=0, background=background, opaqueresize=True)
        self.workspace_panes.pack(fill='both', expand=True)
        left = ttk.Frame(self.workspace_panes, padding=10)
        self.workspace_panes.add(left, minsize=360, width=380, stretch='never')
        ttk.Label(left, text=' XDiffCalculator', image=self.logo, compound='left', font=('Segoe UI', 17, 'bold')).pack(anchor='w')
        ttk.Label(left, text='Parallel beam / point grains / offline').pack(anchor='w', pady=(0, 10))
        panel = tk.Canvas(left, width=350, highlightthickness=0, background=background)
        scrollbar = ttk.Scrollbar(left, orient='vertical', command=panel.yview)
        scrollbar.pack(side='right', fill='y')
        panel.pack(fill='both', expand=True)
        panel.configure(yscrollcommand=scrollbar.set)
        form = ttk.Frame(panel)
        form_window = panel.create_window((0, 0), window=form, anchor='nw', width=342)
        panel.bind('<Configure>', lambda event: panel.itemconfigure(form_window, width=event.width))
        form.columnconfigure(1, weight=1)
        form.bind('<Configure>', lambda event: panel.configure(scrollregion=panel.bbox('all')))
        self.form, self.row = form, 0
        self.section('SAMPLE')
        self.example_name = tk.StringVar(value=os.path.basename(DEFAULTS['cif']))
        ttk.Label(form, text='Example crystal').grid(row=self.row, column=0, sticky='w', pady=3)
        self.examples = ttk.Combobox(form, textvariable=self.example_name, values=example_cifs(), state='readonly', width=20, height=12)
        self.examples.grid(row=self.row, column=1, columnspan=2, sticky='ew', padx=4)
        self.examples.configure(postcommand=lambda: self.examples.configure(values=example_cifs()))
        self.examples.bind('<<ComboboxSelected>>', self.select_example)
        self.row += 1
        self.field('CIF file', 'cif', button=('Browse', self.browse_cif))
        self.blocks = self.field('Structure block', 'block', choices=cif_blocks(DEFAULTS['cif']))
        self.crystal_info = tk.StringVar()
        self.crystal_grid = ttk.Frame(form)
        self.crystal_grid.grid(row=self.row, column=0, columnspan=3, sticky='ew', pady=(8, 4))
        self.crystal_values = [tk.StringVar(value='—') for index in range(7)]
        for index, heading in enumerate(('a (Å)', 'b (Å)', 'c (Å)', 'α (°)', 'β (°)', 'γ (°)', 'sgno')):
            self.crystal_grid.columnconfigure(index, weight=1, uniform='cell')
            ttk.Label(self.crystal_grid, text=heading, anchor='center', font=('Segoe UI', 9)).grid(row=0, column=index, sticky='ew', pady=(0, 3))
            ttk.Label(self.crystal_grid, textvariable=self.crystal_values[index], anchor='center', font=('Segoe UI', 9)).grid(row=1, column=index, sticky='ew')
        self.row += 1
        ttk.Label(form, textvariable=self.crystal_info, wraplength=330, justify='left').grid(row=self.row, column=0, columnspan=3, sticky='w')
        self.row += 1
        for key in ('cif', 'block'):
            self.vars[key].trace_add('write', self.schedule_crystal_info)
        self.update_crystal_info()
        self.sample_mode = self.field('Sample mode', 'mode', choices=['Powder rings', 'Multiple orientations', 'Random grains', 'Imported orientations'])
        euler_heading = ttk.Frame(form)
        euler_heading.grid(row=self.row, column=0, columnspan=3, sticky='ew', pady=(5, 0))
        ttk.Label(euler_heading, text='EulerZXZ (Bunge)').pack(side='left')
        self.euler_count = tk.StringVar()
        ttk.Label(euler_heading, textvariable=self.euler_count).pack(side='left', padx=8)
        self.add_orientation_button = ttk.Button(euler_heading, text='+', width=3, command=self.add_orientation)
        self.add_orientation_button.pack(side='right')
        self.row += 1
        self.euler_frame = ttk.Frame(form)
        self.euler_frame.grid(row=self.row, column=0, columnspan=3, sticky='ew')
        self.row += 1
        self.updating_euler = False
        self.update_euler_entries()
        self.vars['single_euler'].trace_add('write', self.update_euler_entries)
        ttk.Label(form, text='phi1, Phi, phi2 in degrees; crystal → sample.\nDefault: 27.37° from [001] in cubic IPF-Z.', wraplength=330).grid(row=self.row, column=0, columnspan=3, sticky='w')
        self.row += 1
        self.field('Number of grains', 'number')
        self.field('Random seed', 'seed')
        self.field('EulerZXZ file [deg]', 'orientation_file', button=('Browse', self.browse_orientations))
        self.section('BEAM & DETECTOR')
        self.field('Beam mode', 'beam_mode', choices=['Monochromatic', 'Polychromatic (Laue)'])
        self.energy_field = self.field('Mono energy [keV]', 'energy')
        laue_row_start = self.row
        self.laue_fields = []
        for title, key in [('Minimum energy [keV]', 'energy_min'), ('Maximum energy [keV]', 'energy_max'), ('Fixed omega [deg]', 'omega'), ('Laue polarization', 'laue_polarization')]:
            self.laue_fields.append(self.field(title, key))
        self.source_field = self.field('Source spectrum', 'source', choices=['Flat', 'W', 'Cr', 'Cu', 'CSV'])
        self.tube_field = self.field('Reference tube [kV]', 'tube_voltage', choices=['40', '80', '120'])
        self.spectrum_field = self.field('Spectrum CSV', 'spectrum_file', button=('Browse', self.browse_spectrum))
        self.weight_field = self.field('Laue markers', 'laue_weight', choices=['Spectral', 'Uniform'])
        self.colour_field = self.field('Laue colour', 'colour', choices=['Relative weight', 'Energy'])
        ttk.Button(form, text='Preview source spectrum', command=self.preview_spectrum).grid(row=self.row, column=0, columnspan=3, sticky='ew')
        self.row += 1
        ttk.Label(form, text='Tube references: 12° anode, 1 mm Be, 1 mm Al.\nLaue spectral weight = S(E) × F² × P; not calibrated intensity.\nUniform: geometry only. Imported CSV: photons/keV (relative allowed).', wraplength=330).grid(row=self.row, column=0, columnspan=3, sticky='w')
        self.row += 1
        self.laue_controls = [widget for widget in form.winfo_children() if widget.grid_info() and int(widget.grid_info()['row']) >= laue_row_start]
        self.mono_polarization = self.field('Mono polarization [-1,1]', 'polarization')
        ttk.Label(form, text='0.95: nearly horizontal; 0: unpolarized.\nBlank: disable polarization factor.', wraplength=330).grid(row=self.row, column=0, columnspan=3, sticky='w')
        self.row += 1
        self.field('Distance [mm]', 'distance')
        for title, key in [('Tilt X [deg]', 'tilt_x'), ('Tilt Y [deg]', 'tilt_y'), ('Tilt Z [deg]', 'tilt_z')]:
            self.field(title, key)
        ttk.Label(form, text='Tilts: right-handed about fixed lab X, Y, Z,\napplied X → Y → Z at the direct-beam point.\nDistance is along +X; beam centre stays fixed.', wraplength=330).grid(row=self.row, column=0, columnspan=3, sticky='w')
        self.row += 1
        preset = self.field('Detector preset', 'preset', choices=list(self.presets) + ['Custom'])
        preset.bind('<<ComboboxSelected>>', self.apply_preset)
        for title, key in [('Columns', 'columns'), ('Rows', 'rows'), ('Pixel Y [um]', 'pixel_y'), ('Pixel Z [um]', 'pixel_z'),
                           ('Beam centre column', 'center_col'), ('Beam centre row', 'center_row'), ('Max 2theta [deg]', 'max_tth')]:
            widget = self.field(title, key)
            if key in ('columns', 'rows', 'pixel_y', 'pixel_z'):
                widget.bind('<KeyRelease>', lambda event: self.vars['preset'].set('Custom'))
        self.section('ROTATION SCAN (GRAIN MODES)')
        self.scan_fields = []
        for title, key in [('Omega start [deg]', 'start'), ('Omega end (excluded) [deg]', 'end'), ('Frame width [deg]', 'step')]:
            self.scan_fields.append(self.field(title, key))
        self.scan_help = ttk.Label(form, wraplength=330, justify='left')
        self.scan_help.grid(row=self.row, column=0, columnspan=3, sticky='w', pady=5)
        self.row += 1
        self.vars['mode'].trace_add('write', self.update_scan_inputs)
        self.vars['beam_mode'].trace_add('write', self.update_beam_inputs)
        self.vars['source'].trace_add('write', self.update_beam_inputs)
        self.update_beam_inputs()
        self.update_scan_inputs()
        self.field('Display minimum [%]', 'minimum')
        buttons = ttk.Frame(form)
        buttons.grid(row=self.row, column=0, columnspan=3, pady=10, sticky='ew')
        self.run_button = ttk.Button(buttons, text='Calculate', command=self.run)
        self.run_button.pack(side='left')
        ttk.Button(buttons, text='Cancel', command=self.cancel_event.set).pack(side='left', padx=5)
        ttk.Button(buttons, text='Help', command=self.help).pack(side='left')
        self.row += 1
        for label, command in [('Save settings', self.save_settings), ('Load settings', self.load_settings),
                               ('Export CSV tables', self.export_csv), ('Export orientations (EulerZXZ)', self.export_orientations), ('Export snapshot as png', self.export_png)]:
            ttk.Button(form, text=label, command=command).grid(row=self.row, column=0, columnspan=3, sticky='ew', pady=2)
            self.row += 1

        right = ttk.Frame(self.workspace_panes, padding=10)
        self.workspace_panes.add(right, minsize=600, stretch='always')
        self.plot_panes = tk.PanedWindow(right, orient='vertical', sashwidth=7, sashrelief='raised', borderwidth=0, background=background, opaqueresize=True)
        self.plot_panes.pack(fill='both', expand=True)
        plot_area = ttk.Frame(self.plot_panes)
        table_area = ttk.Frame(self.plot_panes)
        self.plot_panes.add(plot_area, minsize=340, height=700, stretch='always')
        self.plot_panes.add(table_area, minsize=170, height=190, stretch='never')
        self.info = ttk.Label(plot_area, text='X beam | Z up | +omega CCW viewed from +Z', wraplength=820)
        self.info.pack(anchor='w', pady=(0, 6))
        self.canvas = tk.Canvas(plot_area, background='#101c2c', highlightthickness=0, height=120)
        self.canvas.pack(fill='both', expand=True)
        self.canvas.bind('<Configure>', lambda event: self.draw())
        self.canvas.bind('<Button-1>', self.start_drag)
        self.canvas.bind('<B1-Motion>', self.drag_view)
        self.canvas.bind('<ButtonRelease-1>', self.end_drag)
        self.canvas.bind('<MouseWheel>', self.zoom_wheel)
        self.canvas.bind('<Motion>', self.hover_info)
        self.canvas.bind('<Leave>', self.clear_hover)
        self.canvas.bind('<Double-Button-1>', self.spot_info)
        self.cursor_info = tk.StringVar(value='Move over the detector to inspect pixel coordinates and relative signal intensity.')
        cursor_label = tk.Label(plot_area, textvariable=self.cursor_info, anchor='w', justify='left', wraplength=740, height=3, width=1,
                                background=background, foreground=ttk.Style(root).lookup('TLabel', 'foreground') or '#000000', font='TkDefaultFont', borderwidth=0)
        cursor_label.pack(fill='x', pady=(3, 0))
        zoom_bar = ttk.Frame(plot_area)
        zoom_bar.pack(fill='x', pady=3)
        ttk.Checkbutton(zoom_bar, text='Rectangle zoom', variable=self.rectangle_mode, command=self.change_drag_mode).pack(side='left', padx=(0, 4))
        ttk.Button(zoom_bar, text='Restore view', command=self.restore_view).pack(side='left')
        self.zoom_info = tk.StringVar(value='1.00x | Wheel: zoom at cursor; left-drag: pan')
        ttk.Label(zoom_bar, textvariable=self.zoom_info).pack(side='left', padx=8)
        toolbar = ttk.Frame(plot_area)
        toolbar.pack(fill='x', pady=5)
        self.laue_labels_box = ttk.Checkbutton(plot_area, text='Label top 10 strongest Laue spots', variable=self.laue_labels, command=self.draw)
        self.laue_labels_box.pack(anchor='w')
        self.integrated_box = ttk.Checkbutton(toolbar, text='All frames (integrated scan)', variable=self.integrated, command=self.draw)
        self.integrated_box.pack(side='left')
        ttk.Checkbutton(toolbar, text='Ring overlay', variable=self.overlay, command=self.draw).pack(side='left', padx=8)
        ttk.Label(toolbar, text='Frame (zero-based):').pack(side='left')
        self.frame_box = ttk.Spinbox(toolbar, from_=0, to=0, width=7, textvariable=self.frame, command=self.draw)
        self.frame_box.pack(side='left')
        self.frame_box.bind('<Return>', lambda event: self.draw())
        self.view_help = ttk.Label(plot_area, wraplength=740, justify='left')
        self.view_help.pack(anchor='w', pady=(0, 5))
        self.update_scan_view()
        table_area.rowconfigure(0, weight=1)
        table_area.columnconfigure(0, weight=1)
        columns = ('index', 'hkl', 'd', 'inv_d', 'tth', 'mult', 'intensity', 'relative')
        self.table = ttk.Treeview(table_area, columns=columns, show='headings', height=8, selectmode='browse')
        for key, title, width in zip(columns, ['Index', 'HKL representative', 'd [Å]', '1/d [Å⁻¹]', '2theta [deg]', 'Multiplicities', 'Intensity [a.u.]', 'Relative [%]'], [45, 125, 80, 90, 90, 95, 110, 85]):
            self.table.heading(key, text=title)
            self.table.column(key, width=width, minwidth=width, anchor='center')
        table_scroll = ttk.Scrollbar(table_area, orient='vertical', command=self.table.yview)
        table_horizontal = ttk.Scrollbar(table_area, orient='horizontal', command=self.table.xview)
        self.table.configure(yscrollcommand=table_scroll.set, xscrollcommand=table_horizontal.set)
        table_scroll.grid(row=0, column=1, sticky='ns')
        table_horizontal.grid(row=1, column=0, sticky='ew')
        self.table.grid(row=0, column=0, sticky='nsew')
        self.table.bind('<<TreeviewSelect>>', self.highlight)
        table_note = ttk.Label(table_area, text='Table: powder-profile intensities (azimuth-averaged P). Plot: relative F² × P × L; grey rings in grain mode are guides. No absolute photon counts or absorption/detector response.', wraplength=740)
        self.table_note = table_note
        table_note.grid(row=2, column=0, columnspan=2, sticky='ew', pady=5)
        status_label = ttk.Label(table_area, textvariable=self.status, wraplength=820)
        status_label.grid(row=3, column=0, columnspan=2, sticky='ew')
        table_area.bind('<Configure>', lambda event: [label.configure(wraplength=max(200, event.width - 10)) for label in (table_note, status_label)])
        plot_area.bind('<Configure>', lambda event: [label.configure(wraplength=max(200, event.width - 10)) for label in (self.info, self.view_help, cursor_label)])
        root.after(100, self.poll)

    def update_euler_setting(self, *args):
        """Store all rows; retain the legacy settings key for single-row compatibility."""
        if not self.updating_euler:
            self.updating_euler = True
            self.vars['single_euler'].set(';'.join(','.join(variable.get() for variable in row) for row in self.euler_rows))
            self.updating_euler = False

    def update_euler_entries(self, *args):
        """Rebuild editable rows when loading settings or adding/removing a grain."""
        if not self.updating_euler:
            rows = self.vars['single_euler'].get().split(';')
            self.updating_euler = True
            for widget in self.euler_frame.winfo_children():
                widget.destroy()
            self.euler_rows, self.euler_widgets = [], []
            for column, label in enumerate(('ID', 'phi1 [deg]', 'Phi [deg]', 'phi2 [deg]')):
                ttk.Label(self.euler_frame, text=label).grid(row=0, column=column, sticky='w')
            for index, text in enumerate(rows):
                values = [value.strip() for value in text.split(',')] if ',' in text else text.split()
                values = (values + ['', '', ''])[:3]
                variables = [tk.StringVar(value=value) for value in values]
                self.euler_rows.append(variables)
                ttk.Label(self.euler_frame, text=str(index)).grid(row=index+1, column=0, padx=(0, 4))
                for column, variable in enumerate(variables, 1):
                    self.euler_frame.columnconfigure(column, weight=1)
                    entry = ttk.Entry(self.euler_frame, textvariable=variable, width=9)
                    entry.grid(row=index+1, column=column, sticky='ew', padx=(0, 4), pady=2)
                    self.euler_widgets.append(entry)
                    variable.trace_add('write', self.update_euler_setting)
                if len(rows) > 1:
                    button = ttk.Button(self.euler_frame, text='−', width=2, command=lambda index=index: self.remove_orientation(index))
                    button.grid(row=index+1, column=4)
                    self.euler_widgets.append(button)
            self.euler_values = self.euler_rows[0]
            self.euler_count.set('{} / 10'.format(len(rows)))
            self.updating_euler = False
            self.update_orientation_controls()

    def update_orientation_controls(self):
        """Enable manual rows only in the corresponding sample mode."""
        enabled = self.vars['mode'].get() == 'Multiple orientations'
        for widget in self.euler_widgets:
            widget.configure(state='normal' if enabled else 'disabled')
        self.add_orientation_button.configure(state='normal' if enabled and len(self.euler_rows) < 10 else 'disabled')

    def add_orientation(self):
        """Append an editable copy of the current default, up to ten grains."""
        if len(self.euler_rows) < 10:
            self.vars['single_euler'].set(self.vars['single_euler'].get() + ';' + DEFAULTS['single_euler'])

    def remove_orientation(self, index):
        """Remove a row, retaining at least one orientation and consecutive grain IDs."""
        rows = self.vars['single_euler'].get().split(';')
        if len(rows) > 1:
            del rows[index]
            self.vars['single_euler'].set(';'.join(rows))

    def schedule_crystal_info(self, *args):
        """Clear old metadata immediately and coalesce edits to the selected CIF/block."""
        self.crystal_info.set('Reading crystal metadata...')
        for variable in self.crystal_values:
            variable.set('—')
        if self.metadata_job is not None:
            self.root.after_cancel(self.metadata_job)
        self.metadata_job = self.root.after(150, self.update_crystal_info)

    def update_crystal_info(self):
        """Display lattice parameters, crystal system and short space-group name."""
        if self.metadata_job is not None:
            self.root.after_cancel(self.metadata_job)
            self.metadata_job = None
        try:
            metadata = crystal_metadata(self.vars['cif'].get(), self.vars['block'].get() or None)
            values = ['{:.5g}'.format(value) for value in metadata['cell']] + [metadata['sgno'] if metadata['sgno'] is not None else '—']
            for variable, value in zip(self.crystal_values, values):
                variable.set(value)
            text = 'Crystal system: {}\nShort name: {}'.format(metadata['crystal_system'], metadata['short_name']) if metadata['sgno'] is not None else 'Space group could not be resolved.'
        except Exception as error:
            for variable in self.crystal_values:
                variable.set('—')
            text = 'Crystal metadata unavailable: {}'.format(error)
        self.crystal_info.set(text)

    def view_transform(self):
        """Return detector-to-canvas scale/origin and a fixed plotting viewport."""
        d = self.result['detector']
        width, height = self.canvas.winfo_width(), self.canvas.winfo_height()
        left, top, right, bottom = 30, 35, max(31, width - 130), max(36, height - 35)
        scale = min((right - left)/d['columns'], (bottom - top)/d['rows'])*self.zoom
        x0 = (left + right)/2 - self.view_center[0]*d['columns']*scale
        y0 = (top + bottom)/2 - self.view_center[1]*d['rows']*scale
        return scale, x0, y0, (left, top, right, bottom)

    def zoom_view(self, factor, x=None, y=None):
        """Zoom about the viewport centre or a supplied cursor, preserving that point."""
        if self.result is None:
            return
        scale, x0, y0, viewport = self.view_transform()
        left, top, right, bottom = viewport
        x, y = ((left + right)/2, (top + bottom)/2) if x is None else (x, y)
        if not left <= x <= right or not top <= y <= bottom:
            return
        d = self.result['detector']
        old_zoom = self.zoom
        self.zoom = max(1, min(32, self.zoom*factor))
        new_scale = scale*self.zoom/old_zoom
        self.view_center = (((x - x0)/scale - (x - (left + right)/2)/new_scale)/d['columns'],
                            ((y - y0)/scale - (y - (top + bottom)/2)/new_scale)/d['rows'])
        self.view_center = tuple(max(0, min(1, value)) for value in self.view_center)
        if self.zoom == 1:
            self.view_center = (0.5, 0.5)
        self.draw()

    def zoom_wheel(self, event):
        """Zoom toward the mouse pointer on Windows."""
        if event.delta:
            self.cancel_drag()
            self.zoom_view(1.25 if event.delta > 0 else 0.8, event.x, event.y)
            self.hover_info(event)
        return 'break'

    def restore_view(self):
        """Fit the full detector without altering calculated data or intensity scale."""
        self.zoom, self.view_center = 1.0, (0.5, 0.5)
        self.rectangle_mode.set(False)
        self.change_drag_mode()
        self.draw()

    def cancel_drag(self):
        """Discard an unfinished pan or rectangle without changing the view."""
        self.pan_anchor, self.drag_start, self.drag_moved = None, None, False
        self.canvas.delete('zoom_rectangle')

    def change_drag_mode(self):
        """Switch between default left-drag pan and one-shot rectangle zoom."""
        self.cancel_drag()
        self.canvas.configure(cursor='crosshair' if self.rectangle_mode.get() else '')
        if hasattr(self, 'zoom_info'):
            self.zoom_info.set('{:.2f}x | Wheel: zoom at cursor; left-drag: {}'.format(self.zoom, 'select rectangle' if self.rectangle_mode.get() else 'pan'))

    def start_drag(self, event):
        """Begin left-button navigation; simple clicks still inspect spots."""
        self.cancel_drag()
        if self.result is None:
            return
        scale, x0, y0, (left, top, right, bottom) = self.view_transform()
        if left <= event.x <= right and top <= event.y <= bottom:
            self.drag_start = (event.x, event.y)
            if not self.rectangle_mode.get():
                self.start_pan(event)

    def drag_view(self, event):
        """Preview a selection rectangle or pan after a small drag threshold."""
        if self.drag_start is None:
            return
        x, y = self.drag_start
        if math.hypot(event.x - x, event.y - y) < 4 and not self.drag_moved:
            return
        self.drag_moved = True
        if self.rectangle_mode.get():
            scale, x0, y0, (left, top, right, bottom) = self.view_transform()
            self.canvas.delete('zoom_rectangle')
            self.canvas.create_rectangle(x, y, max(left, min(right, event.x)), max(top, min(bottom, event.y)), outline='#ffffff', dash=(4, 3), width=2, tags='zoom_rectangle')
        else:
            self.pan_view(event)
        self.hover_info(event)

    def end_drag(self, event):
        """Fit a selected rectangle, or inspect a spot only when no drag occurred."""
        if self.drag_start is None:
            return
        x, y = self.drag_start
        moved = self.drag_moved or math.hypot(event.x - x, event.y - y) >= 4
        if self.rectangle_mode.get() and moved:
            scale, x0, y0, (left, top, right, bottom) = self.view_transform()
            ex, ey = max(left, min(right, event.x)), max(top, min(bottom, event.y))
            if abs(ex - x) >= 4 and abs(ey - y) >= 4:
                d = self.result['detector']
                self.view_center = (max(0, min(1, ((x + ex)/2 - x0)/(scale*d['columns']))),
                                    max(0, min(1, ((y + ey)/2 - y0)/(scale*d['rows']))))
                self.zoom = min(32, self.zoom*min((right - left)/abs(ex - x), (bottom - top)/abs(ey - y)))
                self.rectangle_mode.set(False)
                self.change_drag_mode()
                self.draw()
        elif not moved and not self.rectangle_mode.get():
            self.spot_info(event)
        self.cancel_drag()

    def start_pan(self, event):
        """Start left-button panning inside the plotting viewport."""
        if self.result is None or self.zoom == 1:
            return
        scale, x0, y0, (left, top, right, bottom) = self.view_transform()
        if left <= event.x <= right and top <= event.y <= bottom:
            self.pan_anchor = (event.x, event.y, self.view_center, scale)

    def pan_view(self, event):
        """Translate the zoomed detector while keeping legend and labels fixed."""
        if self.pan_anchor is None:
            return
        x, y, center, scale = self.pan_anchor
        d = self.result['detector']
        self.view_center = (max(0, min(1, center[0] - (event.x - x)/(scale*d['columns']))),
                            max(0, min(1, center[1] - (event.y - y)/(scale*d['rows']))))
        self.draw()

    def update_scan_inputs(self, *args):
        """Explain scan inputs and disable them for orientation-averaged powder."""
        self.update_orientation_controls()
        powder = self.vars['mode'].get() == 'Powder rings' or self.vars['beam_mode'].get() == 'Polychromatic (Laue)'
        for widget in self.scan_fields:
            widget.configure(state='disabled' if powder else 'normal')
        text = 'Powder rings average over orientations; rotation scan settings do not apply.' if powder else (
            'Rotate about +Z; positive omega is CCW viewed from +Z toward the origin.\n'
            'Example: 0 to 2, width 0.25 gives 8 frames: [0,0.25), ... [1.75,2).\n'
            'Width is angular exposure, not calculation accuracy. Each ideal reflection belongs to one frame; no rocking width.\n'
            'Press Calculate after changing these inputs.')
        self.scan_help.configure(text=text)
        if self.vars['beam_mode'].get() == 'Polychromatic (Laue)':
            self.scan_help.configure(text='Laue uses the fixed omega above. No rotation integration or frame width is applied.')

    def update_beam_inputs(self, *args):
        """Enable beam-specific inputs without changing a completed calculation."""
        laue = self.vars['beam_mode'].get() == 'Polychromatic (Laue)'
        modes = ['Multiple orientations', 'Random grains', 'Imported orientations']
        self.sample_mode.configure(values=modes if laue else ['Powder rings'] + modes)
        if laue and self.vars['mode'].get() == 'Powder rings':
            self.vars['mode'].set('Multiple orientations')
        self.energy_field.configure(state='disabled' if laue else 'normal')
        self.mono_polarization.configure(state='disabled' if laue else 'normal')
        for widget in self.laue_controls:
            if laue:
                widget.grid()
            else:
                widget.grid_remove()
        for widget in self.laue_fields:
            widget.configure(state='normal' if laue else 'disabled')
        for widget in (self.source_field, self.weight_field, self.colour_field):
            widget.configure(state='readonly' if laue else 'disabled')
        self.tube_field.configure(state='readonly' if laue and self.vars['source'].get() in ('W', 'Cr', 'Cu') else 'disabled')
        voltages = TUBE_VOLTAGES.get(self.vars['source'].get(), ('40', '80', '120'))
        self.tube_field.configure(values=voltages)
        if self.vars['tube_voltage'].get() not in voltages:
            self.vars['tube_voltage'].set(voltages[0])
        self.spectrum_field.configure(state='normal' if laue and self.vars['source'].get() == 'CSV' else 'disabled')
        self.update_scan_inputs()

    def browse_spectrum(self):
        """Select a CSV photon spectral-density file."""
        path = filedialog.askopenfilename(filetypes=[('Spectrum CSV', '*.csv')])
        if path:
            self.vars['spectrum_file'].set(path)
            self.vars['source'].set('CSV')

    def preview_spectrum(self):
        """Plot current source settings on a small Tk canvas without extra dependencies."""
        try:
            settings = self.snapshot()
            spectrum = source_spectrum(settings)
            low, high = float(settings['energy_min']), float(settings['energy_max'])
            if not 0 < low < high or not math.isfinite(high):
                raise ValueError('Require 0 < Emin < Emax.')
            grid = [low] + [e for e in spectrum.get('energies', []) if low < e < high] + [high]
            density = evaluate_spectrum(spectrum, grid)
            maximum = max(density) or 1
            window = tk.Toplevel(self.root)
            window.title('Source photon spectral density')
            canvas = tk.Canvas(window, width=720, height=320, background='#101c2c')
            canvas.pack()
            points = [coordinate for e, v in zip(grid, density) for coordinate in (55 + 610*(e - low)/(high - low), 270 - 230*v/maximum)]
            canvas.create_line(*points, fill='#4bdae8')
            canvas.create_text(360, 300, text='Energy [{:g}, {:g}] keV'.format(low, high), fill='white')
            canvas.create_text(360, 18, text='Photon density per keV, normalized for preview', fill='white')
            ttk.Label(window, text='\n'.join(spectrum['comments']), wraplength=700).pack(padx=10, pady=8)
        except Exception as error:
            messagebox.showerror('Spectrum', str(error))

    def update_scan_view(self):
        """Show the calculated scan interval and enable only meaningful view controls."""
        grain_scan = self.result is not None and bool(self.result['matrices']) and not self.result.get('laue')
        self.integrated_box.configure(state='normal' if grain_scan else 'disabled')
        self.laue_labels_box.configure(state='normal' if self.result is not None and self.result.get('laue') else 'disabled')
        self.frame_box.configure(state='normal' if grain_scan and not self.integrated.get() else 'disabled')
        if not grain_scan:
            text = 'Calculate a grain-mode scan to browse frames. Powder rings have no individual frames.'
        else:
            settings = self.result['settings']
            start, end, step = [float(settings[key]) for key in ('start', 'end', 'step')]
            count = int(math.ceil((end - start)/step))
            if self.integrated.get():
                text = 'All {} frames over [{:g}, {:g}) deg are overlaid. Uncheck All frames to inspect one exposure.'.format(count, start, end)
            else:
                try:
                    frame = self.frame.get()
                except tk.TclError:
                    frame = 0
                frame = max(0, min(count - 1, frame))
                self.frame.set(frame)
                lower = start + frame*step
                text = 'Frame {} of 0-{}: omega [{:g}, {:g}) deg (end excluded). Use arrows or type a frame and press Enter.'.format(frame, count - 1, lower, min(lower + step, end))
        self.view_help.configure(text=text)
        if self.result is not None and self.result.get('laue'):
            self.view_help.configure(text='Fixed omega {}° | {}–{} keV | Exact harmonics summed per grain; energy colour uses their weight-averaged energy. No Laue intensity/response correction.'.format(self.result['settings']['omega'], self.result['settings']['energy_min'], self.result['settings']['energy_max']))

    def section(self, title):
        """Add a compact section heading."""
        ttk.Label(self.form, text=title, font=('Segoe UI', 9, 'bold')).grid(row=self.row, column=0, columnspan=3, sticky='w', pady=(10, 4))
        self.row += 1

    def field(self, title, key, choices=None, button=None):
        """Add one labelled entry or selection field."""
        ttk.Label(self.form, text=title).grid(row=self.row, column=0, sticky='w', pady=3)
        if choices is None:
            widget = ttk.Entry(self.form, textvariable=self.vars[key], width=16)
        else:
            widget = ttk.Combobox(self.form, textvariable=self.vars[key], values=choices, state='readonly', width=20)
        widget.grid(row=self.row, column=1, columnspan=1 if button else 2, sticky='ew', padx=4)
        if button:
            ttk.Button(self.form, text=button[0], command=button[1], width=7).grid(row=self.row, column=2)
        self.row += 1
        return widget

    def snapshot(self):
        """Capture user inputs without sharing Tk objects with worker threads."""
        return dict((key, variable.get()) for key, variable in self.vars.items())

    def browse_cif(self):
        """Select a CIF and populate its structure-block selector."""
        path = filedialog.askopenfilename(filetypes=[('Crystal files', '*.cif'), ('All files', '*.*')])
        if path:
            self.use_cif(path)

    def select_example(self, event=None):
        """Select a bundled/local example without changing its crystallographic data."""
        self.use_cif(os.path.join(ROOT, 'examples', self.example_name.get()))

    def use_cif(self, path):
        """Validate a selected CIF block list and update the active input atomically."""
        try:
            names = cif_blocks(path)
            if not names:
                raise ValueError('No unit-cell blocks found.')
            self.vars['cif'].set(path)
            self.blocks.configure(values=names)
            self.vars['block'].set(names[0])
            local = os.path.dirname(os.path.abspath(path)) == os.path.join(ROOT, 'examples')
            self.example_name.set(os.path.basename(path) if local else '')
        except Exception as error:
            self.example_name.set('')
            messagebox.showerror('CIF error', str(error))

    def browse_orientations(self):
        """Select an orientation file; contents are validated during calculation."""
        path = filedialog.askopenfilename(title='Bunge EulerZXZ angles in degrees', filetypes=[('Euler angle files', '*.csv *.txt'), ('All files', '*.*')])
        if path:
            self.vars['orientation_file'].set(path)
            self.vars['mode'].set('Imported orientations')

    def apply_preset(self, event=None):
        """Set exact image dimensions and pitches; centre the beam on the image."""
        preset = self.presets.get(self.vars['preset'].get())
        if preset:
            for key in ('columns', 'rows', 'pixel_y', 'pixel_z'):
                self.vars[key].set(str(preset[key]))
            self.vars['center_col'].set(str((preset['columns'] - 1)/2))
            self.vars['center_row'].set(str((preset['rows'] - 1)/2))

    def run(self):
        """Start a cancellable calculation and preserve its input snapshot."""
        if self.busy:
            return
        settings = self.snapshot()
        self.busy = True
        self.cancel_event.clear()
        self.run_button.configure(state='disabled')
        self.status.set('Calculating reflections and diffraction events...')
        self.worker = threading.Thread(target=self.work, args=(settings,), daemon=True)
        self.worker.start()

    def work(self, settings):
        """Run scientific calculations off the UI thread; report errors via queue."""
        try:
            progress = lambda done, total: self.messages.put(('progress', '{} / {} grains'.format(done, total)))
            self.messages.put(('result', calculate(settings, self.cancel_event, progress)))
        except Exception as error:
            self.messages.put(('error', str(error)))

    def poll(self):
        """Consume worker notifications on the Tk main thread."""
        try:
            while True:
                kind, value = self.messages.get_nowait()
                if kind == 'progress':
                    self.status.set(value)
                else:
                    self.busy = False
                    self.run_button.configure(state='normal')
                    if kind == 'error':
                        self.status.set(value)
                        if value != 'Calculation cancelled.':
                            messagebox.showerror('Calculation failed', value)
                    else:
                        self.show_result(value)
        except queue.Empty:
            pass
        self.root.after(100, self.poll)

    def show_result(self, result):
        """Display a completed calculation and refill the powder-line table."""
        self.result, self.selected_line = result, None
        self.zoom, self.view_center, self.pan_anchor = 1.0, (0.5, 0.5), None
        self.rectangle_mode.set(False)
        self.change_drag_mode()
        self.table.delete(*self.table.get_children())
        if result.get('laue'):
            columns = ('grain', 'hkl', 'energy', 'wavelength', 'omega', 'col', 'row', 'weight')
            self.table.configure(columns=columns)
            for key, title in zip(columns, ('Grain', 'HKL', 'Energy [keV]', 'Wavelength [Å]', 'Omega [deg]', 'Column', 'Row', 'Spectral weight')):
                self.table.heading(key, text=title)
                self.table.column(key, width=110, minwidth=85, anchor='center')
            for index, event in enumerate(result['events']):
                self.table.insert('', 'end', iid=str(index), values=(event['grain'], str(event['hkl'])) + tuple('{:.6g}'.format(event[key]) for key in ('energy_keV', 'wavelength_A', 'omega', 'col', 'row', 'spectral_weight')))
            self.frame.set(0)
            self.frame_box.configure(to=0)
            self.info.configure(text='{} | {} | Laue | fixed omega {}° | {}–{} keV'.format(result['settings']['block'], result['spacegroup'], result['settings']['omega'], result['settings']['energy_min'], result['settings']['energy_max']))
            self.table_note.configure(text='Each row is one signed HKL at its own energy. Harmonics share a plotted marker. Spectral weight S(E) × F² × P is approximate: no Laue Jacobian, absorption, anomalous scattering or detector response.')
            self.draw()
            return
        columns = ('index', 'hkl', 'd', 'inv_d', 'tth', 'mult', 'intensity', 'relative')
        self.table.configure(columns=columns)
        for key, title in zip(columns, ('Index', 'HKL representative', 'd [Å]', '1/d [Å⁻¹]', '2theta [deg]', 'Multiplicities', 'Intensity [a.u.]', 'Relative [%]')):
            self.table.heading(key, text=title)
            self.table.column(key, width=100, minwidth=80, anchor='center')
        self.table_note.configure(text='Powder table: azimuth-averaged P. Plot: relative F² × P × L. No absolute counts or detector response.')
        for i, line in enumerate(result['lines']):
            self.table.insert('', 'end', iid=str(i), values=(line['index'], str(line['hkl']), '{:.5f}'.format(line['d']), '{:.5f}'.format(1/line['d']), '{:.4f}'.format(line['two_theta']),
                              line['multiplicity'], '{:.3g}'.format(line['intensity']), '{:.3f}'.format(line['relative'])))
        settings = result['settings']
        nframes = 1 if not result['matrices'] else int(math.ceil((float(settings['end']) - float(settings['start']))/float(settings['step'])))
        self.frame.set(0)
        self.frame_box.configure(to=nframes - 1)
        self.info.configure(text='{} | {} | wavelength {:.5f} Å | {} | detector pixel view'.format(settings['block'], result['spacegroup'], result['wavelength'], settings['mode']))
        self.draw()

    def visible_data(self):
        """Get thresholded display data; full unthresholded results remain exportable."""
        if self.result is None:
            return [], []
        minimum = float(self.result['settings']['minimum'])
        if self.result.get('laue'):
            return [], [event for event in self.result['groups'] if event['intensity'] > 0 and 100*event['intensity']/self.result['intensity_max'] >= minimum]
        lines = [line for line in self.result['lines'] if line['relative'] >= minimum]
        maximum = self.result['intensity_max']
        try:
            frame = self.frame.get()
        except tk.TclError:
            frame = 0
        events = [event for event in self.result['events'] if 100*event['intensity']/maximum >= minimum and (self.integrated.get() or event['frame'] == frame)]
        return lines, events

    def draw(self):
        """Render schematic rings/spots; cap display only, never exported results."""
        self.update_scan_view()
        self.canvas.delete('all')
        self.drawn_spots = []
        self.hover_items = {}
        if self.result is None:
            self.canvas.create_text(30, 40, anchor='nw', fill='#b9c9dc', text='Load a crystal, define the detector, then Calculate.', font=('Segoe UI', 13))
            return
        d = self.result['detector']
        width, height = self.canvas.winfo_width(), self.canvas.winfo_height()
        scale, x0, y0, (left, top, right, bottom) = self.view_transform()
        self.zoom_info.set('{:.2f}x | Wheel: zoom at cursor; left-drag: {}'.format(self.zoom, 'select rectangle' if self.rectangle_mode.get() else 'pan'))
        x1, y1 = x0 + scale*d['columns'], y0 + scale*d['rows']
        cx, cy = x0 + (d['cx'] + 0.5)*scale, y0 + (d['cy'] + 0.5)*scale
        lines, events = self.visible_data()
        maximum = self.result['intensity_max']
        laue = self.result.get('laue', False)
        energy_colour = laue and self.result['settings']['colour'] == 'Energy'
        if energy_colour:
            low, high = float(self.result['settings']['energy_min']), float(self.result['settings']['energy_max'])
        if self.overlay.get() or not self.result['matrices']:
            for line in lines:
                for u0, v0, u1, v1, color in ring_segments(line, d, self.result['polarization'], maximum, bool(self.result['matrices'])):
                    coordinates = (x0 + (u0 + 0.5)*scale, y0 + (v0 + 0.5)*scale, x0 + (u1 + 0.5)*scale, y0 + (v1 + 0.5)*scale)
                    item = self.canvas.create_line(*coordinates, fill=color, width=3 if line is self.selected_line else 1)
                    self.hover_items[item] = ('ring', coordinates, line)
        limit = 20000
        stride = max(1, int(math.ceil(len(events)/limit)))
        labelled = []
        if laue and self.laue_labels.get():
            labelled = heapq.nlargest(10, (event for event in events if sum(member['spectral_weight'] for member in event['members']) > 0), key=lambda event: sum(member['spectral_weight'] for member in event['members']))
        plotted = events[::stride]
        plotted_ids = set(id(event) for event in plotted)
        plotted += [event for event in labelled if id(event) not in plotted_ids]
        for event in plotted:
            x, y = x0 + (event['col'] + 0.5)*scale, y0 + (event['row'] + 0.5)*scale
            color = intensity_color(event['intensity']/maximum)
            if energy_colour:
                low, high = float(self.result['settings']['energy_min']), float(self.result['settings']['energy_max'])
                color = intensity_color((event['energy_keV'] - low)/(high - low))
            selected = laue and self.selected_line is not None and self.selected_line in event.get('members', [])
            radius = 4 if selected else 2
            item = self.canvas.create_oval(x - radius, y - radius, x + radius, y + radius, fill=color, outline='white' if selected else '')
            self.hover_items[item] = ('spot', (x, y, x, y), event)
            if left <= x <= right and top <= y <= bottom:
                self.drawn_spots.append((x, y, event))
        # Mask arcs outside the detector's active rectangle (module gaps not modelled).
        if x0 <= cx <= x1 and y0 <= cy <= y1:
            self.canvas.create_line(cx - 6, cy, cx + 6, cy, fill='#f57676')
            self.canvas.create_line(cx, cy - 6, cx, cy + 6, fill='#f57676')
        clip_left, clip_top, clip_right, clip_bottom = max(left, x0), max(top, y0), min(right, x1), min(bottom, y1)
        for box in [(0, 0, width, clip_top), (0, clip_bottom, width, height), (0, clip_top, clip_left, clip_bottom), (clip_right, clip_top, width, clip_bottom)]:
            self.canvas.create_rectangle(*box, fill='#101c2c', outline='')
        # Draw the visible plot frame last so clipping masks cannot cover its edges.
        self.canvas.create_rectangle(clip_left, clip_top, clip_right, clip_bottom, outline='#9bb0c9', tags='plot_border')
        self.draw_laue_labels(labelled, scale, x0, y0, (left, top, right, bottom))
        title = '+Z up | +Y right | view along incident beam (+X)'
        if any(abs(value) > 1e-12 for value in d['tilts']):
            title = 'Detector pixel view | tilts X={:g}°, Y={:g}°, Z={:g}°'.format(*d['tilts'])
        self.canvas.create_text(left, top - 16, anchor='w', fill='#b9c9dc', text=title)
        self.canvas.create_text(left, bottom + 16, anchor='w', fill='#b9c9dc', text='{} x {} pixels | {:.2f} x {:.2f} mm'.format(d['columns'], d['rows'], d['columns']*d['py'], d['rows']*d['pz']))
        bar_x, bar_y, bar_height = width - 80, 60, max(60, height - 130)
        for index in range(100):
            self.canvas.create_rectangle(bar_x, bar_y + index*bar_height/100, bar_x + 16, bar_y + (index + 1)*bar_height/100,
                                         fill=intensity_color(1 - (index + 0.5)/100), outline='', tags='intensity_scale')
        for percent in (0, 25, 50, 75, 100):
            label = '{:.4g}'.format(low + (high - low)*percent/100) if energy_colour else str(percent)
            self.canvas.create_text(bar_x + 22, bar_y + (1 - percent/100)*bar_height, anchor='w', fill='#b9c9dc', text=label)
        scale_title = 'Energy\n[keV]' if energy_colour else ('Weight\n[%]' if laue else 'Intensity\n[%]')
        self.canvas.create_text(bar_x, bar_y - 25, fill='#b9c9dc', text=scale_title, tags='intensity_scale')
        text = '{} powder lines; {} grains; {} detected events; {} in current view.'.format(len(lines), len(self.result['matrices']), len(self.result['events']), len(events))
        if stride > 1:
            text += ' Preview decimated 1:{}; CSV contains all events.'.format(stride)
        if laue:
            text = '{} grains; {} HKL contributions; {} harmonic groups; {} visible. Weights are approximate.'.format(len(self.result['matrices']), len(self.result['events']), len(self.result['groups']), len(events))
        if self.result['matrices'] and not self.integrated.get() and not laue:
            settings = self.result['settings']
            try:
                start = float(settings['start']) + self.frame.get()*float(settings['step'])
                text += ' Omega [{:.4g}, {:.4g}) deg.'.format(start, min(start + float(settings['step']), float(settings['end'])))
            except tk.TclError:
                pass
        self.status.set(text)
        if self.hover_xy is not None:
            self.update_hover(*self.hover_xy)

    def draw_laue_labels(self, events, scale, x0, y0, bounds):
        """Label ranked harmonic groups using their strongest individual reflection.
        Args:
            events: Up to ten groups ranked by summed spectral weight.
            scale, x0, y0: Detector-to-canvas transform.
            bounds: Plot viewport in canvas coordinates, including detector margins.
        Returns:
            None. Draws labels and connecting lines on the canvas.
        """
        left, top, right, bottom = bounds
        occupied = []
        for event in events:
            x, y = x0 + (event['col'] + 0.5)*scale, y0 + (event['row'] + 0.5)*scale
            if not left <= x <= right or not top <= y <= bottom:
                continue
            member = max(event['members'], key=lambda entry: entry['spectral_weight'])
            text = 'Grain {} | ({}, {}, {})\n{:.2f} keV'.format(member['grain'], *member['hkl'], member['energy_keV'])
            if len(event['members']) > 1:
                text += ' (+{} hkl)'.format(len(event['members']) - 1)
            item = self.canvas.create_text(0, 0, text=text, anchor='nw', justify='left', fill='#ffffff', font=('Segoe UI', 9, 'bold'), tags='laue_spot_label')
            box = self.canvas.bbox(item)
            w, h = box[2] - box[0], box[3] - box[1]
            candidates = []
            for dx, dy in [(dx, dy) for dx in (10, -w - 10, w + 20, -2*w - 20) for dy in (-h - 8, 8, -2*h - 16, h + 16, -3*h - 24, 2*h + 24)]:
                tx, ty = max(left + 3, min(right - w - 3, x + dx)), max(top + 3, min(bottom - h - 3, y + dy))
                overlap = sum(max(0, min(tx + w + 3, b[2]) - max(tx - 3, b[0]))*max(0, min(ty + h + 3, b[3]) - max(ty - 3, b[1])) for b in occupied)
                candidates.append((overlap, (tx - x)**2 + (ty - y)**2, tx, ty))
            _, _, tx, ty = min(candidates)
            self.canvas.move(item, tx - box[0], ty - box[1])
            box = self.canvas.bbox(item)
            occupied.append(box)
            connector = self.canvas.create_line(x, y, max(box[0], min(box[2], x)), max(box[1], min(box[3], y)), fill='#c6d5e8', tags='laue_label_connector')
            background = self.canvas.create_rectangle(box[0] - 2, box[1] - 1, box[2] + 2, box[3] + 1, fill='#101c2c', outline='', tags='laue_label_background')
            self.canvas.tag_lower(connector, item)
            self.canvas.tag_lower(background, item)

    def clear_hover(self, event=None):
        """Clear the readout when the cursor leaves the plot canvas."""
        self.hover_xy = None
        self.cursor_info.set('Move over the detector to inspect pixel coordinates and relative signal intensity.')

    def hover_info(self, event):
        """Track cursor motion without recomputing diffraction or repainting the plot."""
        self.hover_xy = (event.x, event.y)
        self.update_hover(event.x, event.y)

    def update_hover(self, x, y):
        """Read zero-based detector coordinates and the nearest rendered signal within 6 screen pixels.
        Args:
            x, y: Cursor coordinates on the canvas, respecting zoom and pan.
        Returns:
            None. Updates the live readout; intensities are schematic feature weights,
            not rasterized pixel counts. Grey guide rings have no spot intensity.
        """
        if self.result is None:
            self.cursor_info.set('Calculate a diffraction pattern to inspect the detector.')
            return
        scale, x0, y0, (left, top, right, bottom) = self.view_transform()
        d = self.result['detector']
        col, row = (x - x0)/scale - 0.5, (y - y0)/scale - 0.5
        if not (left <= x <= right and top <= y <= bottom and -0.5 <= col < d['columns'] - 0.5 and -0.5 <= row < d['rows'] - 0.5):
            self.cursor_info.set('Outside detector — pixel coordinates and signal intensity unavailable.')
            return
        text = 'Pixel (col, row) = ({:.2f}, {:.2f})'.format(col, row)
        nearest = {}
        # Tk searches its spatial display list; only nearby rendered items need Python checks.
        for item in self.canvas.find_overlapping(x - 6, y - 6, x + 6, y + 6):
            feature = self.hover_items.get(item)
            if feature is None:
                continue
            kind, (ax, ay, bx, by), record = feature
            dx, dy = bx - ax, by - ay
            length2 = dx*dx + dy*dy
            fraction = max(0, min(1, ((x - ax)*dx + (y - ay)*dy)/length2)) if length2 else 0
            distance2 = (x - ax - fraction*dx)**2 + (y - ay - fraction*dy)**2
            if distance2 <= 36 and (kind not in nearest or distance2 < nearest[kind][0]):
                nearest[kind] = (distance2, record)
        if 'spot' in nearest:
            spot = nearest['spot'][1]
            intensity = spot['intensity']
            text += ' | Nearest spot I = {:.5g} a.u. ({:.2f}%) | Grain {} HKL {}'.format(intensity, 100*intensity/self.result['intensity_max'], spot['grain'], spot['hkl'])
            if self.result.get('laue'):
                text = 'Pixel ({:.2f}, {:.2f}) | Grain {} | Marker weight {:.5g} | Contributions: '.format(col, row, spot['grain'], intensity)
                text += '; '.join('{} at {:.5g} keV'.format(member['hkl'], member['energy_keV']) for member in spot['members'])
        elif 'ring' in nearest and not self.result['matrices']:
            line = nearest['ring'][1]
            point = pixel_lab(d, col, row)
            azimuth = math.atan2(point[2], point[1])
            intensity = line['sum_f2']*line['lorentz']*polarization_factor(line['two_theta'], azimuth, self.result['polarization'])
            text += ' | Nearest ring I = {:.5g} a.u. ({:.2f}%)'.format(intensity, 100*intensity/self.result['intensity_max'])
        else:
            text += ' | I = 0 (no nearby rendered signal)'
        if 'ring' in nearest:
            line = nearest['ring'][1]
            text += ' | Ring {} HKL {}'.format(line['index'], line['hkl'])
            if self.result['matrices']:
                text += ' (guide)'
        self.cursor_info.set(text)

    def highlight(self, event=None):
        """Highlight the ring selected in the reflection table."""
        selection = self.table.selection()
        collection = self.result['events'] if self.result and self.result.get('laue') else (self.result['lines'] if self.result else [])
        self.selected_line = collection[int(selection[0])] if selection and self.result else None
        self.draw()

    def spot_info(self, event):
        """Report the nearest rendered spot's grain, HKL and diffraction angle."""
        if self.drawn_spots:
            x, y, spot = min(self.drawn_spots, key=lambda item: (item[0] - event.x)**2 + (item[1] - event.y)**2)
            if math.hypot(x - event.x, y - event.y) < 12:
                if self.result.get('laue'):
                    messagebox.showinfo('Laue harmonic contributions', 'Grain {} | pixel ({:.5f}, {:.5f})\n'.format(spot['grain'], spot['col'], spot['row']) + '\n'.join(
                        'HKL {} | {:.7g} keV | weight {:.7g}'.format(member['hkl'], member['energy_keV'], member['spectral_weight']) for member in spot['members']))
                    return
                self.status.set('Grain {} | HKL {} | omega {:.5f} deg | 2theta {:.5f} deg | I={:.4g} a.u. | |F|²={:.4g}, P={:.4g}, L={:.4g}'.format(
                    spot['grain'], spot['hkl'], spot['omega'], spot['two_theta'], spot['intensity'], spot['f2'], spot['polarization'], spot['lorentz']))

    def save_settings(self):
        """Save editable inputs, not a potentially stale calculation result."""
        path = filedialog.asksaveasfilename(defaultextension='.json', filetypes=[('Settings', '*.json')])
        if path:
            try:
                with open(path, 'w', encoding='utf-8') as stream:
                    json.dump(settings_for_disk(self.snapshot()), stream, indent=2)
            except Exception as error:
                messagebox.showerror('Save failed', str(error))

    def load_settings(self):
        """Load input settings; user explicitly recalculates afterwards."""
        path = filedialog.askopenfilename(filetypes=[('Settings', '*.json')])
        if path:
            try:
                with open(path, encoding='utf-8') as stream:
                    settings = json.load(stream)
                if not isinstance(settings, dict):
                    raise ValueError('Expected a settings dictionary.')
                settings = settings_from_disk(settings)
                for key in ('beam_mode', 'energy_min', 'energy_max', 'omega', 'source', 'tube_voltage', 'spectrum_file', 'laue_polarization', 'colour', 'laue_weight', 'single_euler'):
                    settings.setdefault(key, DEFAULTS[key])
                # Files saved before tilt inputs existed describe an untilted detector.
                for key in ('tilt_x', 'tilt_y', 'tilt_z'):
                    settings.setdefault(key, '0')
                for key in DEFAULTS:
                    if key in settings:
                        self.vars[key].set(str(settings[key]))
                names = cif_blocks(self.vars['cif'].get())
                self.blocks.configure(values=names)
                path = os.path.abspath(self.vars['cif'].get())
                self.example_name.set(os.path.basename(path) if os.path.dirname(path) == os.path.join(ROOT, 'examples') else '')
                self.status.set('Settings loaded. Press Calculate to update the displayed result.')
            except Exception as error:
                messagebox.showerror('Load failed', str(error))

    def export_csv(self):
        """Export all calculated signed reflections, powder lines, spots and provenance."""
        if not self.result:
            return
        path = filedialog.asksaveasfilename(defaultextension='.csv', filetypes=[('CSV prefix', '*.csv')])
        if not path:
            return
        base = os.path.splitext(path)[0]
        if self.result.get('laue'):
            targets = [base + suffix for suffix in ('_laue.csv', '_laue_settings.json')]
            if any(os.path.exists(target) for target in targets) and not messagebox.askyesno('Overwrite exports?', 'Replace existing Laue exports?'):
                return
            try:
                export_laue(self.result, base)
                self.status.set('Saved Laue reflection table and settings including the source spectrum.')
            except Exception as error:
                messagebox.showerror('Export failed', str(error))
            return
        targets = [base + suffix for suffix in ('_reflections.csv', '_rings.csv', '_spots.csv', '_settings.json')]
        if any(os.path.exists(target) for target in targets) and not messagebox.askyesno('Overwrite exports?', 'One or more output files already exist. Replace them?'):
            return
        try:
            with open(targets[0], 'w', newline='', encoding='utf-8') as stream:
                writer = csv.writer(stream)
                writer.writerow(['h', 'k', 'l', 'd_A', 'inv_d_A_inverse', 'two_theta_deg', 'F_squared'])
                for record in self.result['reflections']:
                    writer.writerow(list(record['hkl']) + [record['d'], 1/record['d'], record['two_theta'], record['f2']])
            with open(targets[1], 'w', newline='', encoding='utf-8') as stream:
                writer = csv.writer(stream)
                writer.writerow(['index', 'hkl_representative', 'all_hkls', 'multiplicities', 'd_A', 'inv_d_A_inverse', 'two_theta_deg', 'sum_F_squared', 'polarization_mean', 'lorentz', 'intensity_au', 'relative_percent'])
                for line in self.result['lines']:
                    writer.writerow([line['index'], line['hkl'], ';'.join(str(hkl) for hkl in line['hkls']), line['multiplicity'], line['d'], 1/line['d'], line['two_theta'], line['sum_f2'], line['polarization'], line['lorentz'], line['intensity'], line['relative']])
            with open(targets[2], 'w', newline='', encoding='utf-8') as stream:
                writer = csv.writer(stream)
                writer.writerow(['grain_id', 'h', 'k', 'l', 'omega_deg', 'frame', 'column', 'row', 'F_squared', 'two_theta_deg', 'azimuth_rad', 'polarization', 'lorentz', 'intensity_au'])
                for spot in self.result['events']:
                    writer.writerow([spot['grain']] + list(spot['hkl']) + [spot[key] for key in ('omega', 'frame', 'col', 'row', 'f2', 'two_theta', 'azimuth', 'polarization', 'lorentz', 'intensity')])
            with open(targets[3], 'w', encoding='utf-8') as stream:
                json.dump(settings_for_disk(self.result['settings']), stream, indent=2)
            self.status.set('Saved {}_reflections.csv, _rings.csv, _spots.csv and _settings.json.'.format(base))
        except Exception as error:
            messagebox.showerror('Export failed', str(error))

    def export_orientations(self):
        """Save the displayed grain orientations as equivalent Bunge ZXZ Euler angles."""
        if not self.result or not self.result['matrices']:
            messagebox.showinfo('Orientations', 'Calculate a single, random or imported grain sample first.')
            return
        path = filedialog.asksaveasfilename(defaultextension='.csv', filetypes=[('Orientations', '*.csv')])
        if path:
            try:
                save_orientations(path, self.result['matrices'])
            except Exception as error:
                messagebox.showerror('Export failed', str(error))

    def export_png(self):
        """Save the current plot canvas at 300 DPI, preserving its displayed view."""
        if not self.result:
            return
        try:
            self.root.update_idletasks()
            # Capture the view before opening the dialog; no desktop screenshot is used.
            image = canvas_snapshot(self.canvas, dpi=300)
            path = filedialog.asksaveasfilename(defaultextension='.png', filetypes=[('PNG snapshot', '*.png')])
            if not path:
                return
            image.save(path, format='PNG', dpi=(300, 300))
            self.status.set('Saved PNG snapshot: {} | {} x {} pixels | 300 dpi.'.format(path, image.width, image.height))
        except Exception as error:
            messagebox.showerror('Export failed', str(error))

    def help(self):
        """Explain coordinate, orientation, exposure and intensity conventions."""
        messagebox.showinfo('Conventions & scope',
            'X: incident beam; Z: vertical up; Y: outward, right-handed.\n'
            'Laue mode: fixed omega, analytical wavelength for every signed HKL.\n'
            'Laue weight S(E) x F squared x P is not quantitative intensity.\n'
            'No rotation Lorentz factor is used in Laue mode. Harmonics overlap.\n'
            '+omega: counter-clockwise viewed from +Z toward the origin.\n'
            'U maps crystal Cartesian column vectors into the sample at omega=0.\n'
            'Crystal basis: a along x; b in xy; c has positive z; B=A^(-T), no 2*pi.\n'
            'Bunge EulerZXZ input/export: phi1_deg,Phi_deg,phi2_deg.\n'
            'U = Rz(phi1) Rx(Phi) Rz(phi2), crystal to sample at omega zero.\n'
            'Default Euler angles put cubic sample Z 27.37 degrees from [001].\n'
            'Zero tilt: columns +Y, rows -Z; zero-based pixel centres.\n'
            'Detector tilts: right-handed fixed lab X then Y then Z (Rz Ry Rx).\n'
            'Pivot: direct-beam point (distance,0,0); distance is along +X.\n\n'
            'All grains have equal weights and origin positions. Random rotations are uniform.\n'
            'Scan intervals are [start,end); step is frame exposure width. Exact Bragg events\n'
            'are assigned to frames; finite mosaicity/rocking width is not modelled.\n'
            'Intensity = |F|^2 x P x L; powder lines include multiplicity.\n'
            'Polarization: +1 horizontal, -1 vertical, 0 unpolarized; blank disables P.\n'
            'Powder table uses azimuth-averaged P; ring colours vary with azimuth.\n'
            'L: powder cos(theta)/sin(2theta)^2; rotation 1/|sin(2theta)cos(azimuth)|.\n'
            'Azimuth is measured from +Y toward +Z. Tangent rotation events are singular.\n'
            'Colour scale is linear: 100% = full-scan maximum spot intensity, or\n'
            'maximum full-ring intensity in powder mode; fixed across frames.\n'
            'CIF occupancies and ADPs are included. No anomalous scattering,\n'
            'absorption, texture weights, module gaps or detector response.\n'
            'Near absorption edges relative intensities may be inaccurate.\n'
            'Settings edits take effect only on Calculate. CSV exports the full last result.')

    def close(self):
        """Cancel pending work and close the window."""
        self.cancel_event.set()
        self.root.destroy()


def smoke_test(report_path):
    """Exercise bundled Gemmi, resources and Tk, writing a diagnostic JSON report."""
    report = dict(ok=False)
    root = None
    try:
        root = tk.Tk()
        root.withdraw()
        application = Application(root)
        if application.example_name.get() != 'Ni.cif' or application.vars['block'].get() != 'Ni' or application.crystal_values[-1].get() != '225':
            raise RuntimeError('Initial crystal metadata mismatch.')
        if application.crystal_info.get() != 'Crystal system: Cubic\nShort name: Fm3\u0305m':
            raise RuntimeError('Crystal system or short space-group name is missing.')
        application.use_cif(os.path.join(ROOT, 'examples', 'Al_fcc.cif'))
        application.update_crystal_info()
        if application.crystal_values[-1].get() != '225':
            raise RuntimeError('Selected crystal metadata mismatch.')
        application.vars['block'].set('missing')
        application.update_crystal_info()
        if 'unavailable' not in application.crystal_info.get() or application.crystal_values[-1].get() != '—':
            raise RuntimeError('Invalid block retained stale crystal metadata.')
        application.use_cif(DEFAULTS['cif'])
        if any(str(widget.cget('state')) != 'disabled' for widget in application.scan_fields):
            raise RuntimeError('Powder scan inputs should be disabled.')
        application.vars['mode'].set('Random grains')
        if any(str(widget.cget('state')) != 'normal' for widget in application.scan_fields):
            raise RuntimeError('Grain scan inputs should be enabled.')
        # Exercise every detector preset, including effective imaging-camera pitches.
        for name, preset in application.presets.items():
            application.vars['preset'].set(name)
            application.apply_preset()
            current = application.snapshot()
            if float(current['pixel_y']) != preset['pixel_y'] or float(current['pixel_z']) != preset['pixel_z']:
                raise RuntimeError('Preset pitch mismatch: {}.'.format(name))
            if float(current['center_col']) != (preset['columns'] - 1)/2 or float(current['center_row']) != (preset['rows'] - 1)/2:
                raise RuntimeError('Preset centre mismatch: {}.'.format(name))
        settings = dict(DEFAULTS)
        settings.update(mode='Random grains', number='20')
        result = calculate(settings)
        if not result['lines'] or not result['events']:
            raise RuntimeError('Smoke calculation produced no diffraction.')
        application.show_result(result)
        root.update()
        if str(application.table.item('0', 'values')[0]) != '0' or not application.canvas.find_withtag('intensity_scale'):
            raise RuntimeError('Missing zero-based table index or intensity scale.')
        if abs(float(application.table.item('0', 'values')[3]) - 1/result['lines'][0]['d']) > 0.0000051:
            raise RuntimeError('Reciprocal-spacing column mismatch.')
        if len(application.workspace_panes.panes()) != 2 or len(application.plot_panes.panes()) != 2:
            raise RuntimeError('Missing adjustable workspace/plot panes.')
        application.update_hover(0, 0)
        if 'Outside detector' not in application.cursor_info.get():
            raise RuntimeError('Outside-detector hover readout mismatch.')
        application.clear_hover()
        scale_max = application.result['intensity_max']
        application.canvas.configure(width=800, height=440)
        root.update_idletasks()
        # Verify that zoom preserves the detector point under its screen anchor.
        initial = application.view_transform()
        scale, x0, y0, (left, top, right, bottom) = initial
        px, py = left + 0.6*(right - left), top + 0.4*(bottom - top)
        original_point = ((px - x0)/scale, (py - y0)/scale)
        application.zoom_view(2, px, py)
        scale, x0, y0, viewport = application.view_transform()
        zoomed_point = ((px - x0)/scale, (py - y0)/scale)
        if application.zoom != 2 or max(abs(a - b) for a, b in zip(original_point, zoomed_point)) > 1e-8:
            raise RuntimeError('Zoom moved the detector point under its anchor.')
        application.zoom_view(0.8)
        if abs(application.zoom - 1.6) > 1e-12:
            raise RuntimeError('Zoom out failed.')
        mouse = type('Mouse', (), dict(x=px, y=py))()
        application.start_pan(mouse)
        old_center = application.view_center
        mouse.x += 0.1*(right - left)
        application.pan_view(mouse)
        if application.view_center == old_center:
            raise RuntimeError('Panning did not change the view centre.')
        application.restore_view()
        if application.view_transform() != initial or application.result['intensity_max'] != scale_max:
            raise RuntimeError('Restore view changed detector fit or intensity normalization.')
        if str(application.frame_box.cget('state')) != 'disabled':
            raise RuntimeError('Integrated view must disable frame selection.')
        application.integrated.set(False)
        application.frame.set(4)
        application.draw()
        if str(application.frame_box.cget('state')) != 'normal' or '[1, 1.25)' not in application.view_help.cget('text'):
            raise RuntimeError('Individual frame interval/control mismatch.')
        if application.result['intensity_max'] != scale_max:
            raise RuntimeError('Intensity normalization changed between frames.')
        # Map invisibly for real canvas dimensions, then exercise the bundled PNG encoder/fonts.
        root.attributes('-alpha', 0)
        root.deiconify()
        root.update()
        tilted = calculate(dict(settings, tilt_x='12', tilt_y='-18', tilt_z='23'))
        application.show_result(tilted)
        root.update()
        if application.result['detector']['tilts'] != (12, -18, 23):
            raise RuntimeError('Tilted detector workflow mismatch.')
        snapshot = canvas_snapshot(application.canvas)
        with tempfile.TemporaryDirectory() as folder:
            from PIL import Image
            path = os.path.join(folder, 'snapshot.png')
            snapshot.save(path, dpi=(300, 300))
            with Image.open(path) as saved:
                if abs(saved.info['dpi'][0] - 300) > 0.01 or saved.size != snapshot.size:
                    raise RuntimeError('PNG snapshot size/DPI mismatch.')
        report.update(ok=True, lines=len(result['lines']), spots=len(result['events']), presets=len(application.presets), examples=len(example_cifs()), tk=root.tk.call('info', 'patchlevel'))
        application.vars['beam_mode'].set('Polychromatic (Laue)')
        application.vars['mode'].set('Multiple orientations')
        if application.vars['energy_min'].get() != '5' or application.vars['energy_max'].get() != '80':
            raise RuntimeError('Laue default energy range is incorrect.')
        if application.vars['source'].get() != 'W' or application.vars['tube_voltage'].get() != '120':
            raise RuntimeError('Laue default source must be tungsten at 120 kV.')
        application.vars['single_euler'].set('10,20,30')
        if [variable.get() for variable in application.euler_values] != ['10', '20', '30']:
            raise RuntimeError('Euler entries did not follow settings.')
        application.euler_values[0].set('15')
        if application.vars['single_euler'].get() != '15,20,30':
            raise RuntimeError('Euler settings did not follow entries.')
        application.vars['single_euler'].set(DEFAULTS['single_euler'])
        application.euler_values[0].set('')
        application.add_orientation_button.invoke()
        if application.euler_rows[0][0].get() != '' or application.euler_rows[0][1].get() != DEFAULTS['single_euler'].split(',')[1]:
            raise RuntimeError('Adding a row must preserve unfinished angle edits.')
        application.vars['single_euler'].set(DEFAULTS['single_euler'])
        for index in range(12):
            application.add_orientation()
        if len(application.euler_rows) != 10 or str(application.add_orientation_button.cget('state')) != 'disabled':
            raise RuntimeError('Manual orientations must stop at ten rows.')
        application.euler_rows[1][0].set('42')
        restored = settings_from_disk(settings_for_disk(application.snapshot()))
        if len(multiple_orientations(restored['single_euler'])) != 10 or restored['single_euler'].split(';')[1].split(',')[0] != '42':
            raise RuntimeError('Multiple orientation edits/settings did not round trip.')
        application.remove_orientation(1)
        if len(application.euler_rows) != 9 or str(application.add_orientation_button.cget('state')) != 'normal':
            raise RuntimeError('Removing an orientation did not re-enable addition.')
        application.vars['single_euler'].set(DEFAULTS['single_euler'])
        application.remove_orientation(0)
        if len(application.euler_rows) != 1:
            raise RuntimeError('The last orientation must be retained.')
        report['multiple_orientation_controls'] = True
        if any(str(widget.cget('state')) != 'disabled' for widget in application.scan_fields):
            raise RuntimeError('Laue rotation-scan inputs must be disabled.')
        for target, voltages in TUBE_VOLTAGES.items():
            for voltage in voltages:
                source_spectrum(dict(source=target, tube_voltage=voltage))
        laue_result = calculate(dict(DEFAULTS, beam_mode='Polychromatic (Laue)', mode='Single orientation', source='W', colour='Energy', minimum='0'))
        application.show_result(laue_result)
        root.update()
        if not application.drawn_spots or application.table.heading('energy')['text'] != 'Energy [keV]':
            raise RuntimeError('Laue GUI did not display energy-resolved reflections.')
        application.laue_labels_box.invoke()
        label_items = application.canvas.find_withtag('laue_spot_label')
        expected_groups = sorted((group for group in laue_result['groups'] if sum(member['spectral_weight'] for member in group['members']) > 0), key=lambda group: sum(member['spectral_weight'] for member in group['members']), reverse=True)[:10]
        if len(label_items) != len(expected_groups):
            raise RuntimeError('Laue strongest-spot labels have the wrong count.')
        for item, group in zip(label_items, expected_groups):
            member = max(group['members'], key=lambda entry: entry['spectral_weight'])
            expected_text = 'Grain {} | ({}, {}, {})\n{:.2f} keV'.format(member['grain'], *member['hkl'], member['energy_keV'])
            if not application.canvas.itemcget(item, 'text').startswith(expected_text):
                raise RuntimeError('Laue label rank, HKL or energy is incorrect.')
        canvas_snapshot(application.canvas)
        application.laue_labels_box.invoke()
        if application.canvas.find_withtag('laue_spot_label'):
            raise RuntimeError('Laue labels did not switch off.')
        x, y, spot = application.drawn_spots[0]
        application.update_hover(x, y)
        if 'keV' not in application.cursor_info.get():
            raise RuntimeError('Laue hover omitted energy.')
        application.table.selection_set('0')
        application.highlight()
        canvas_snapshot(application.canvas)
        application.show_result(calculate(dict(DEFAULTS, beam_mode='Polychromatic (Laue)', mode='Single orientation', energy_min='1', energy_max='1.01', colour='Energy')))
        application.laue_labels_box.invoke()
        if application.canvas.find_withtag('laue_spot_label'):
            raise RuntimeError('Empty Laue result must not have labels.')
        root.update()
        application.show_result(result)
        root.update()
        if application.canvas.find_withtag('laue_spot_label') or str(application.laue_labels_box.cget('state')) != 'disabled':
            raise RuntimeError('Laue labels must be disabled for monochromatic results.')
        if application.table.heading('mult')['text'] != 'Multiplicities':
            raise RuntimeError('Returning to monochromatic table failed.')
        report.update(laue_spots=len(laue_result['events']), laue_groups=len(laue_result['groups']), source_spectra=sum(len(values) for values in TUBE_VOLTAGES.values()))
    except Exception:
        report['ok'] = False
        report['error'] = traceback.format_exc()
    finally:
        if root is not None:
            root.destroy()
        with open(report_path, 'w', encoding='utf-8') as stream:
            json.dump(report, stream, indent=2)
    return 0 if report['ok'] else 1


def main():
    """Launch the UI or perform an explicit bundled-runtime smoke test."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--smoke-test', metavar='REPORT_JSON')
    args = parser.parse_args()
    if args.smoke_test:
        return smoke_test(args.smoke_test)
    root = tk.Tk()
    Application(root)
    root.mainloop()
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
