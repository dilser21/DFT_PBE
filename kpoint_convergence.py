import os
import subprocess

# --- AYARLAR ---
# Quantum Espresso komutun
qe_command = "q-e-qe-7.3.1/bin/pw.x"  
pseudo_dir = "./" 

# Sabitlenen ecutwfc değeri (Az önceki testten seçtik)
fixed_ecut = 80

# Denenecek K-point mesh değerleri (N x N x N)
k_grids = [2, 3, 4, 5, 6, 7, 8]

results = []

print(f"{'K-Grid':<10} | {'Total Energy (Ry)':<20} | {'Time (s)':<10}")
print("-" * 50)

for k in k_grids:
    # 1. Input dosyasını oluştur (ecutwfc sabit, K_POINTS değişken)
    input_content = f"""
&CONTROL
    calculation = 'scf'
    restart_mode = 'from_scratch'
    prefix = 'CsPbBr3'
    pseudo_dir = '{pseudo_dir}'
    outdir = './tmp/'
    tstress = .true.
    tprnfor = .true.
/
&SYSTEM
    ibrav = 0
    nat = 5
    ntyp = 3
    ecutwfc = {fixed_ecut}
    ecutrho = {fixed_ecut * 4}
/
&ELECTRONS
    conv_thr = 1.0d-8
    mixing_beta = 0.7
/
ATOMIC_SPECIES
 Cs  132.90545  Cs.nc.z_9.oncvpsp3.dojo.v4-str.UPF
 Pb  207.2      Pb.pbe-dn-kjpaw_psl.0.2.2.UPF
 Br  79.904     br_pbe_v1.4.uspp.F.UPF
ATOMIC_POSITIONS (crystal)
 Cs 0.000000 0.000000 0.000000
 Pb 0.500000 0.500000 0.500000
 Br 0.500000 0.500000 0.000000
 Br 0.500000 0.000000 0.500000
 Br 0.000000 0.500000 0.500000
CELL_PARAMETERS (angstrom)
   5.87 0.00 0.00
   0.00 5.87 0.00
   0.00 0.00 5.87
K_POINTS (automatic)
  {k} {k} {k} 0 0 0
"""

    filename = f"ktest_{k}.in"
    outfile = f"ktest_{k}.out"
    
    with open(filename, "w") as f:
        f.write(input_content)

    # 2. Çalıştır
    cmd = f"{qe_command} -in {filename} > {outfile}"
    subprocess.run(cmd, shell=True)

    # 3. Sonucu Oku
    energy = "N/A"
    time = "N/A"
    
    try:
        with open(outfile, "r") as f:
            lines = f.readlines()
            for line in lines:
                if "!    total energy" in line:
                    energy = float(line.split('=')[1].split('Ry')[0].strip())
                if "PWSCF        :" in line:
                    time = line.split('CPU')[0].split(':')[1].strip()
    except:
        energy = "Hata"

    results.append((k, energy))
    print(f"{k}x{k}x{k:<5} | {energy:<20} | {time:<10}")

    # os.remove(filename)
    # os.remove(outfile)

print("-" * 50)
print("K-Point testi tamamlandı!")
