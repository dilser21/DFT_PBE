import numpy as np
import matplotlib.pyplot as plt

# --- AYARLAR ---
dosya_adi = 'cspbbr3_bands.dat.gnu'
fermi_enerjisi = 2.9468  # SCF çıktındaki "Fermi energy" değerini buraya yazman gerekebilir.
                      # Eğer grafik çok yukarıda veya aşağıda çıkarsa burayı güncelle.

# --- VERİYİ OKU ---
try:
    data = np.loadtxt(dosya_adi)
except OSError:
    print(f"Hata: {dosya_adi} bulunamadı. Python dosyasının veriyle aynı klasörde olduğundan emin ol.")
    exit()

k_path = data[:, 0]
energies = data[:, 1] - fermi_enerjisi

# --- GRAFİK ÇİZİMİ ---
plt.figure(figsize=(10, 6))

# Bantları çiz (Nokta nokta değil, çizgi olarak)
plt.plot(k_path, energies, 'b-', linewidth=1.5, alpha=0.8)

# --- SİMETRİ NOKTALARI (LOG DOSYANDAN ALINDI) ---
# Log çıktındaki koordinatlar:
# 0.0000 -> Gamma
# 0.5000 -> X
# 1.0000 -> M
# 1.7071 -> Gamma
# 2.5731 -> R

ozel_noktalar = [0.0, 0.5, 1.0, 1.7071, 2.5731]
etiketler = [r'$\Gamma$', 'X', 'M', r'$\Gamma$', 'R']

# Dikey çizgilerle simetri noktalarını belirt
for nokta in ozel_noktalar:
    plt.axvline(x=nokta, color='k', linestyle='--', linewidth=0.8)

plt.xticks(ozel_noktalar, etiketler, fontsize=12)
plt.ylabel(r'$E - E_F$ (eV)', fontsize=12)
plt.title('CsPbBr$_3$ Elektronik Bant Yapısı', fontsize=14)

# Fermi seviyesini çiz (0 çizgisi)
plt.axhline(0, color='r', linestyle=':', linewidth=1)

# Görünümü iyileştir (Enerji aralığını -4 ile +6 eV arasına odaklar, genelde bandgap buradadır)
plt.ylim(-5, 6)
plt.xlim(min(k_path), max(k_path))
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
