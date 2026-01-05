import numpy as np
import matplotlib.pyplot as plt

# --- AYARLAR ---
dosya_adi = 'cspbbr3_bands.dat.gnu'
fermi_enerjisi = 2.9468 

# --- VERİYİ OKU ---
try:
    data = np.loadtxt(dosya_adi)
except OSError:
    print(f"Hata: {dosya_adi} bulunamadı.")
    exit()

# --- DÜZELTME BURADA BAŞLIYOR ---
# Verinin yapısını analiz edelim:
# k_path sürekli artıp sıfırlanıyorsa, sıfırlandığı yeri bulup k-noktası sayısını (nk) çözeriz.
tum_k = data[:, 0]
tum_e = data[:, 1] - fermi_enerjisi

# k değerinin azaldığı (resetlendiği) yerleri bul
split_indices = np.where(np.diff(tum_k) < 0)[0] + 1

if len(split_indices) > 0:
    nk = split_indices[0] # Her bir banttaki k-noktası sayısı
else:
    nk = len(tum_k) # Tek bir bant varsa (nadirdir)

# Bant sayısı
nbnd = len(tum_e) // nk

# Veriyi (Bant Sayısı, k-noktası) şeklinde yeniden boyutlandır
bands = tum_e.reshape(nbnd, nk)
k_plot = tum_k[:nk] # X ekseni için sadece ilk seti almamız yeterli

# --- GRAFİK ÇİZİMİ ---
plt.figure(figsize=(10, 6))

# Her bandı ayrı ayrı çiz (Böylece sondan başa çizgi çekmez)
for band in bands:
    plt.plot(k_plot, band, 'b-', linewidth=1.5, alpha=0.8)

# --- SİMETRİ NOKTALARI ---
ozel_noktalar = [0.0, 0.5, 1.0, 1.7071, 2.5731]
etiketler = [r'$\Gamma$', 'X', 'M', r'$\Gamma$', 'R']

for nokta in ozel_noktalar:
    plt.axvline(x=nokta, color='k', linestyle='--', linewidth=0.8)

plt.xticks(ozel_noktalar, etiketler, fontsize=12)
plt.ylabel(r'$E - E_F$ (eV)', fontsize=12)
plt.title(r'CsPbBr$_3$ Elektronik Bant Yapısı', fontsize=14)

plt.axhline(0, color='r', linestyle=':', linewidth=1)

# Görünüm Ayarları
plt.ylim(-4, 6) # Aralığı biraz daralttım, gap daha net görünür
plt.xlim(min(k_plot), max(k_plot))
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
