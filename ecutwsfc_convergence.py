import os
import subprocess

# --- AYARLAR ---
# Quantum Espresso çalıştırma komutun (kendi bilgisayarındaki yola göre düzenle)
qe_command = "q-e-qe-7.3.1/bin/pw.x"
# Pseudo dosyalarının olduğu klasör (Nokta şu anki klasör demek)
pseudo_dir = "./" 

# Denenecek ecutwfc değerleri (Ry cinsinden)
cutoffs = [40, 50, 60, 70, 80, 90, 100]

# Sonuçları saklayacağımız liste
results = []

print(f"{'Cutoff (Ry)':<15} | {'Total Energy (Ry)':<20} | {'Time (s)':<10}")
print("-" * 50)

for ecut in cutoffs:
    # 1. Input dosyasını her döngüde yeniden oluşturuyoruz
    # (Senin cspbbr3_scf.in dosyanı baz aldım, ecutwfc kısmını {ecut} ile dinamik yaptım)
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
    ecutwfc = {ecut}   
    ecutrho = {ecut * 4} 
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
  4 4 4 0 0 0
"""

    # Dosyayı yaz
    filename = f"test_{ecut}.in"
    outfile = f"test_{ecut}.out"
    
    with open(filename, "w") as f:
        f.write(input_content)

    # 2. Quantum ESPRESSO'yu çalıştır
    # Komut: path/to/pw.x -in test_40.in > test_40.out
    cmd = f"{qe_command} -in {filename} > {outfile}"
    
    # Python üzerinden terminal komutu çalıştırıyoruz
    subprocess.run(cmd, shell=True)

    # 3. Sonucu Oku (Output dosyasından 'total energy' satırını bul)
    energy = "N/A"
    time = "N/A"
    
    try:
        with open(outfile, "r") as f:
            lines = f.readlines()
            for line in lines:
                if "!    total energy" in line:
                    energy = float(line.split('=')[1].split('Ry')[0].strip())
                if "PWSCF        :" in line: # Süreyi de alalım
                    time = line.split('CPU')[0].split(':')[1].strip()
    except:
        energy = "Hata"

    # Sonuçları listeye ekle ve ekrana bas
    results.append((ecut, energy))
    print(f"{ecut:<15} | {energy:<20} | {time:<10}")

    # Temizlik (İstersen yorum satırı yap, dosyalar kalsın)
    # os.remove(filename) 
    # os.remove(outfile)

print("-" * 50)
print("Test tamamlandı!")