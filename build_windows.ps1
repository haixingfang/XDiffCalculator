# Author: Haixing Fang
# Email: haixingfang868@gmail.com
param([switch]$OneFile)
$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot
$pythonExe = Join-Path $PSScriptRoot '.venv\Scripts\python.exe'
if (-not (Test-Path -LiteralPath $pythonExe)) {
    throw 'Create .venv and install requirements-build.txt first; see README.md.'
}
# Keep -OneFile accepted for existing commands; every build now has one output.
$bundleMode = '--onefile'
$outputDirectory = 'dist'
& $pythonExe -m PyInstaller --noconfirm --clean --windowed $bundleMode --distpath $outputDirectory --name XDiffCalculator --copy-metadata gemmi --add-data 'detectors.json;.' --add-data 'examples;examples' app.py
if ($LASTEXITCODE -ne 0) { throw 'Executable build failed.' }
