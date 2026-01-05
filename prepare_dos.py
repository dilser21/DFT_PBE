# prepare_dos.py
dos_content = """&DOS
    prefix = 'cspbbr3',
    outdir = './tmp/',
    fildos = 'cspbbr3_dos.dat',
    Emin = -5.0,
    Emax = 10.0,
    DeltaE = 0.05
/
"""

pdos_content = """&PROJWFC
    prefix = 'cspbbr3',
    outdir = './tmp/',
    filpdos = 'cspbbr3_pdos.dat',
    Emin = -5.0,
    Emax = 10.0,
    DeltaE = 0.05
/
"""

with open("cspbbr3_dos.in", "w") as f:
    f.write(dos_content)

with open("cspbbr3_pdos.in", "w") as f:
    f.write(pdos_content)

print("Girdi dosyaları oluşturuldu: cspbbr3_dos.in ve cspbbr3_pdos.in")
