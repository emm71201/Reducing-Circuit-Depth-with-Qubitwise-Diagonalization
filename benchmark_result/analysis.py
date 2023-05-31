#%%
import pandas as pd
# %%
# with open("./benchmark_result/result_5.csv", "r") as f:
#     for line in f:
#         print(line)
ns = [5,10,15]
dfs = []
for n in ns:

    df = pd.read_csv(f"./benchmark_result/result_{n}.csv", )

    dfs.append(df)
# %%
dfs[0]
# %%
