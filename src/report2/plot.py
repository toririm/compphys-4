import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib import rcParams

rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = [
    "Hiragino Maru Gothic Pro",
]

cwd = Path(__file__).parent

# CSVファイルを読み込む
df = pd.read_csv(cwd / 'result.csv')

# プロットを作成
plt.figure(figsize=(10, 6))
plt.scatter(df['init_T'], df['D'], color='blue')
plt.xlabel('初期温度 T [K]')
plt.ylabel('拡散係数 D')
plt.title('初期温度と拡散係数の関係')
plt.grid(True)

# 原点を通る水平線と垂直線を追加
plt.axhline(y=0, color='r', linestyle='--')
plt.axvline(x=0, color='r', linestyle='--')

plt.savefig(cwd / 'plot.png')

# プロットを表示
plt.show()
