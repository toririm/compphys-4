import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from pathlib import Path

cwd = Path(__file__).parent

# CSVファイルを読み込む
df = pd.read_csv(cwd / "LJ.csv")

# 時間(t)とMSDのデータを抽出
t = df["t"]
msd = df["MSD"]

# 線形回帰を実行
slope, intercept, r_value, p_value, std_err = stats.linregress(t, msd)

# 線形近似の結果を計算
line = slope * t + intercept

# プロットを作成
plt.figure(figsize=(10, 6))
plt.scatter(t, msd, color="blue", label="Data")
plt.plot(t, line, color="red", label=f"Linear fit (y = {slope:.4f}x + {intercept:.4f})")
plt.xlabel("Time (t)")
plt.ylabel("MSD")
plt.title("Linear Approximation of MSD")
plt.legend()
plt.grid(True)

# 結果を表示
print(f"Slope: {slope:.6f}")
print(f"Intercept: {intercept:.6f}")
print(f"R-squared: {r_value**2:.6f}")

# 拡散係数 D = slope / 6 を計算
print(f"Diffusion coefficient: {slope / 6:.6f}")

plt.show()
