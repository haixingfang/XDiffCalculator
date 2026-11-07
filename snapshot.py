# -*- coding: utf-8 -*-
# Author: Haixing Fang
# Email: haixingfang868@gmail.com
"""Render the current Tk plot into a high-resolution PNG without screen capture."""
import math
import os
from tkinter import font as tkfont
from PIL import Image, ImageDraw, ImageFont


def canvas_snapshot(canvas, dpi=300):
    """Render the displayed plot items, in drawing order, at the requested DPI.
    Args:
        canvas: Tk Canvas containing the current plot (including clipping masks).
        dpi: Output dots per inch; defaults to 300. Physical size follows Tk's DPI.
    Returns:
        image: Pillow RGB image ready to save with dpi=(dpi,dpi).
            Preserves zoom, pan, visible frame, preview decimation and legend.
    """
    if not math.isfinite(dpi) or dpi <= 0:
        raise ValueError('Snapshot DPI must be positive and finite.')
    scale = dpi/float(canvas.winfo_fpixels('1i'))
    width, height = canvas.winfo_width(), canvas.winfo_height()
    if width <= 1 or height <= 1:
        raise ValueError('The plot must have a visible size before export.')
    image = Image.new('RGB', (max(1, round(width*scale)), max(1, round(height*scale))), canvas.cget('background'))
    draw = ImageDraw.Draw(image)
    fonts = {}
    for item in canvas.find_all():
        if canvas.itemcget(item, 'state') == 'hidden' or 'zoom_rectangle' in canvas.gettags(item):
            continue
        kind = canvas.type(item)
        coordinates = [value*scale for value in canvas.coords(item)]
        fill = canvas.itemcget(item, 'fill') or None
        if kind in ('rectangle', 'oval'):
            x0, y0, x1, y1 = coordinates
            box = (min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1))
            outline = canvas.itemcget(item, 'outline') or None
            thickness = max(1, round(float(canvas.itemcget(item, 'width'))*scale))
            if kind == 'rectangle':
                draw.rectangle(box, fill=fill, outline=outline, width=thickness)
            else:
                draw.ellipse(box, fill=fill, outline=outline, width=thickness)
        elif kind == 'line':
            thickness = max(1, round(float(canvas.itemcget(item, 'width'))*scale))
            draw.line(coordinates, fill=fill, width=thickness)
        elif kind == 'text':
            specification = canvas.itemcget(item, 'font')
            if specification not in fonts:
                source = tkfont.Font(root=canvas, font=specification)
                actual = source.actual()
                size = actual['size']
                pixels = size*dpi/72 if size > 0 else -size*scale
                family = actual['family'].lower()
                regular, bold = ('segoeui.ttf', 'segoeuib.ttf')
                if family == 'arial':
                    regular, bold = 'arial.ttf', 'arialbd.ttf'
                elif family == 'tahoma':
                    regular, bold = 'tahoma.ttf', 'tahomabd.ttf'
                font_file = bold if actual['weight'] == 'bold' else regular
                path = os.path.join(os.environ.get('WINDIR', 'C:\\Windows'), 'Fonts', font_file)
                font = ImageFont.truetype(path, max(1, round(pixels)))
                fonts[specification] = (font, source.metrics('linespace')*scale)
            font, line_height = fonts[specification]
            lines = canvas.itemcget(item, 'text').split('\n')
            lengths = [draw.textlength(line, font=font) for line in lines]
            block_width, block_height = max(lengths or [0]), line_height*len(lines)
            anchor = canvas.itemcget(item, 'anchor')
            x, y = coordinates
            if anchor not in ('w', 'nw', 'sw'):
                x -= block_width if anchor in ('e', 'ne', 'se') else block_width/2
            if anchor not in ('n', 'nw', 'ne'):
                y -= block_height if anchor in ('s', 'sw', 'se') else block_height/2
            justify = canvas.itemcget(item, 'justify')
            ascent = font.getmetrics()[0]
            for index, (line, length) in enumerate(zip(lines, lengths)):
                offset = 0 if justify == 'left' else (block_width - length)/(2 if justify == 'center' else 1)
                draw.text((x + offset, y + index*line_height + ascent), line, font=font, fill=fill, anchor='ls')
        else:
            raise ValueError('Unsupported plot item for PNG snapshot: {}.'.format(kind))
    return image
