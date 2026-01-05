import matplotlib.pyplot as plt
import numpy as np
import glob
import os

# --- AYARLAR ---
# Daha önce bulduğun Fermi enerjisi
e_fermi = 6.1234  

# DÜZELTİLEN KISIM: Senin dosya yapına tam uyan prefix
dosya_prefix = 'cspbbr3_pdos.dat.pdos_atm#' 

# Enerji aralığı (Band Gap çevresine odaklanalım)
xlim_min, xlim_max = -6, 6

# Veri saklama sözlükleri
data_pdos = {} 

print(f"'{dosya_prefix}*' ile başlayan dosyalar taranıyor...")

dosya_listesi = glob.glob(f"{dosya_prefix}*")
if not dosya_listesi:
    print("HATA: Hiç dosya bulunamadı! Lütfen dosya_prefix değişkenini kontrol et.")
    exit()
else:
    print(f"{len(dosya_listesi)} adet PDOS dosyası bulundu.")

# Klasördeki tüm pdos dosyalarını oku
for dosya in dosya_listesi:
    # Dosya adından atom ve orbital bilgisini çekme
    # Örnek: cspbbr3_pdos.dat.pdos_atm#1(Cs)_wfc#1(s)
    try:
        atom_tur = dosya.split('(')[1].split(')')[0]  # 'Cs', 'Pb', 'Br'
        orbital = dosya.split('(')[2].split(')')[0]   # 's', 'p', 'd'
    except IndexError:
        continue

    # Veriyi oku
    try:
        tmp_data = np.loadtxt(dosya)
        enerji = tmp_data[:, 0] - e_fermi
        ldos = tmp_data[:, 1] # Local DOS sütunu
    except:
        continue

    # Sözlüğe ekle veya topla
    if atom_tur not in data_pdos:
        data_pdos[atom_tur] = {}
    
    if orbital not in data_pdos[atom_tur]:
        data_pdos[atom_tur][orbital] = ldos
    else:
        if len(data_pdos[atom_tur][orbital]) == len(ldos):
            data_pdos[atom_tur][orbital] += ldos

# --- ÇİZİM ---
plt.figure(figsize=(10, 6))

print("Grafik oluşturuluyor...")

# Pb (Kurşun) - Genellikle İletim Bandını (Conduction Band) domine eder
if 'Pb' in data_pdos:
    if 's' in data_pdos['Pb']:
        plt.plot(enerji, data_pdos['Pb']['s'], color='black', label='Pb 6s', linewidth=1.5)
    if 'p' in data_pdos['Pb']:
        plt.plot(enerji, data_pdos['Pb']['p'], color='blue', label='Pb 6p', linewidth=1.5)

# Br (Brom) - Genellikle Valans Bandını (Valence Band) domine eder
if 'Br' in data_pdos:
    if 'p' in data_pdos['Br']:
        plt.plot(enerji, data_pdos['Br']['p'], color='red', label='Br 4p', linewidth=1.5)

plt.axvline(0, color='k', linestyle='--', alpha=0.5) # Fermi seviyesi
plt.xlabel(r'$E - E_F$ (eV)', fontsize=12)
plt.ylabel('Density of States (states/eV)', fontsize=12)
plt.title('CsPbBr$_3$ Partial DOS (Atomik Katkılar)', fontsize=14)
plt.xlim(xlim_min, xlim_max)
plt.ylim(0, 15) # Y ekseni ölçeği
plt.legend(loc='upper right')
plt.grid(True, alpha=0.3)
plt.fill_between(enerji, 0, 100, where=(enerji < 0), color='gray', alpha=0.1) 

plt.tight_layout()
plt.show()
