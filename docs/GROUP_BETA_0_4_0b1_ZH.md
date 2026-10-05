# STEP_UPGRADE 组内研究 beta：v0.4.0b1

发布日期：2026-10-05。分发范围：组内研究试用，不是公开正式发布。
Git tag：`v0.4.0b1`；Python 包版本：`0.4.0b1`；安装包名：`step-upgrade`。

这是**降雨对象识别与谱系追踪** beta，不是已经验证的云体／对流识别工具，
也不是托管的云端服务。此次版本只明确安装、推荐入口和交付边界，不修改科学
算法、目录入口条件或冻结配置。授权和引用整理不作为本次组内试用的交付门槛；
这不构成新的授权声明。将来对外发布另行处理。

## 1. 推荐版本与使用原则

唯一推荐的科学配置：
`configs/rainfall_lineage_4km_1h_reference_v1.json`。

| 项目 | 本 beta 推荐值／语义 |
|---|---|
| 网格／时间间隔 | 4 km，1 小时 |
| 输入降雨 | 逐小时累计 RAINNC，mm；先差分得到 mm/h |
| 对象入口 | >=1 mm/h，保留完整目录，不增加强核心门槛 |
| grouping | radius=9 个格点的形态学连接；输出仍只保留实际湿像元 |
| tracking | 冻结的 overlap_ranked + coherent_adjacent 配置 |
| gap | 最多跳过一帧对象缺席；不能用于跨越缺失的原始累计文件 |
| 0.5／0.1 外扩 | 可作为独立测量层；不作为推荐 tracking 特征掩膜 |

9 格点半径是 grouping 参数，不是物理 storm 尺寸或最大相邻对象间距。
一个标签可能包含不连续的降雨片段，不等于一个独立对流核心。
请使用下面的 `run_reference_stream.py`，它会显式加载冻结配置。
直接调用裸 `identify()` / `track_with_graph()` 的默认值**不等于**本配置；
旧演算脚本也可能使用其他实验策略。

## 2. 安装环境

推荐 Python **3.13.5**；本 beta 安装最低要求 Python **3.12**。
锁定的 NumPy/SciPy 版本也要求 >=3.12。干净环境验证使用 macOS ARM64；
执行层使用 POSIX 的 fork、文件锁等机制，Linux/RCC 需在目标节点先通过自检。
其他 Python／操作系统组合尚未经过完整版本矩阵验证。Windows 原生不在当前
组内 beta 支持范围，WSL 也需要独立自检，不能沿用本机验证结论。

`requirements-beta.txt` 固定顶层依赖版本（不是所有传递依赖的完整 lock）：
NumPy 2.5.3、SciPy 1.18.1、scikit-image 0.26.0、xarray 2026.7.0、
netCDF4 1.7.4、Matplotlib 3.11.2 等。每个部署应保留 `pip freeze`；不要在一个
正在运行或准备恢复的环境中升级依赖。

### 新安装：macOS／Linux

在一个新的目录安装，避免覆盖旧 STEP。以下命令的 `python3.13` 必须已经可用：

```bash
git clone --branch v0.4.0b1 --single-branch \
  https://github.com/kjwtang/STEP_UPGRADE.git STEP_UPGRADE_beta_040b1
cd STEP_UPGRADE_beta_040b1
python3.13 -m venv .venv-step-beta
source .venv-step-beta/bin/activate
python -m pip install --only-binary=:all: -r requirements-beta.txt
python -m pip install --no-build-isolation -e .
python -m pip check
python -c "import step; from importlib.metadata import version; print(version('step-upgrade')); print(step.__file__)"
git describe --tags --exact-match
```

预期包版本为 `0.4.0b1`，Git tag 为 `v0.4.0b1`，`step.__file__` 指向本 checkout。
安装包名字不同，但原版与升级版都使用 `import step`，因此**不要把原版 STEP
和这个 beta 装进同一个虚拟环境**。当前安装方式需要保留源码 checkout，不能
只复制环境后删除源码。不要在运行过程中 `git pull` 或切换版本。

如果二进制安装失败，不要默默回退到随意编译或混合环境；先记录 Python、
CPU 架构和安装错误，检查是否有适合该平台的 wheel。

### RCC：沿用现有 module／conda 入口

先在 login node 准备环境，计算在 compute node 上执行：

```bash
module load python
source activate /project2/moyer/kjwtang/envs/wrfconda
python --version
```

先确认输出为 >=3.12；若现有 wrfconda 太旧，需要先选择较新的 Python 环境，
不要用旧 Python 强行安装这些固定版本。确认后，在新的 beta checkout 中：

```bash
python -m venv .venv-step-beta
source .venv-step-beta/bin/activate
python -m pip install --only-binary=:all: -r requirements-beta.txt
python -m pip install --no-build-isolation -e .
python -m pip check
```

这里的个人 conda 路径只适用于当前使用者；其他组员换成自己的基础环境。
每次进入 compute node 后，要进入这个 checkout 并重新激活 `.venv-step-beta`。

## 3. 首次运行：不需要下载任何研究数据

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python -u scripts/beta_smoke.py --output-dir results/beta_smoke_01 --workers 2
```

脚本生成小型合成累计降雨，跨过 1 月／2 月边界：12 个降雨区间需要 13 个
累计快照。它先跑完整控制，再暂停、换 chunk 大小／worker 后恢复，并验证
完整对象、ID、edges/scores、events、family roots、时间／网格与栅格哈希一致。

预期：`frames=12`、`nodes=11`、`edges=10`、`gap_edges=1`，八项 `checks`
全部为 `true`；只有完成全部检查才写入 `results/beta_smoke_01/SUCCESS`。
重跑请换新的输出目录。自检验证安装和执行一致性，不验证气象准确率。

完整回归测试：

```bash
python -m pytest -q
```

## 4. 真实数据的推荐入口

输入文件应为真实 NetCDF 文件，名字形如：
`cstm_d01_1996-07-01_00:00:00.nc` 或同样格式的 `.rain`。
目录会递归扫描，必须属于同一个 climate/member 序列，不能出现重复时间。

要求：

- `RAINNC` 单位为 `mm`，维度为 `(south_north, west_east)`，或额外带长度为 1
  的 `Time` 维度；只读取这个累计降雨变量，不加载其他三维气象变量。
- 时间使用文件名；不要求内部 `Times` 或 `XLAT/XLONG` 数组。
- 所有快照形状与投影属性一致，且 `DX=DY=4000` 米。需要全局属性
  `MAP_PROJ, CEN_LAT, CEN_LON, TRUELAT1, TRUELAT2, STAND_LON, MOAD_CEN_LAT,
  POLE_LAT, POLE_LON, DX, DY`；不要为了通过检查编造真实数据的投影属性。
- N 个逐小时降雨区间需要 N+1 个连续累计快照。缺文件、重复时间、累计下降／
  reset、网格变化会报错，不自动猜测、补零或修复。
- RAINNC 不自动包含其他降雨变量，例如存在 RAINC 时不会把它加进来。

先跑 24 小时的小规模组内试用：

```bash
python -u scripts/run_reference_stream.py /PATH/TO/ONE_SEQUENCE \
  --preset configs/rainfall_lineage_4km_1h_reference_v1.json \
  --start 0 --hours 24 --sequence-id member0_1996_beta_040b1 \
  --workers 4 --chunk-frames 6 --id-executor fresh \
  --output-dir results/beta_real24_01
```

`--start` 是按时间排序后的**累计前驱快照索引**，不是第一个降雨区间的索引。
默认全域，不进行裁剪。默认不保存大型标签数组；小案例需要画图时可在首次
运行时加 `--save-labels`，之后恢复必须保持该选项一致。

中断后，保持相同的数据计划、配置、sequence ID、源码和环境，重复原命令并
加 `--resume`。可以更改 workers、chunk 大小和已验证的执行器；不能通过修改
`--hours` 把一个旧计划悄悄扩展成新计划。不要相邻月份各自独立追踪再拼接编号；
应从一开始规划连续区间，跨月传递状态。长计划入口已支持，但试用先从 24 小时
开始，不把完整十年的运行可靠性当成已验证事实。

进度条表示已处理区间；`SUCCESS` 才是完成标记。多个进程不能写同一个输出目录。

## 5. 目前可以使用的功能

| 功能 | 当前状态／入口 |
|---|---|
| 2D 降雨对象识别 | 推荐；湿像元阈值、形态学 grouping、确定性标签 |
| 连续追踪、split／merge／复杂事件 | 推荐冻结配置；输出关系图，不是正确率认证 |
| 短时对象缺席的 gap 恢复 | 推荐冻结配置；带歧义／中间对象冲突保护 |
| 多 worker identification | 独立帧并行；tracking 状态按时间顺序推进 |
| 分块、跨月状态、checkpoint/resume | 推荐 runner；输出事务提交、完整性检查、单写者锁 |
| node／branch／family 编号 | 序列内不循环回收；split/merge 会结束旧 branch 并创建新 branch |
| 最终 family 目录 | `scripts/canonicalize_stream_families.py`；不改原始 branch/node |
| 保存标签、CSV、诊断信息 | runner 输出；大型标签由 `--save-labels` 显式启用 |
| 更换 chunk/worker 的精确审计 | `scripts/audit_reference_stream.py` |
| 0.5／0.1 seeded 外扩 | 可选测量／研究实验；`step.envelope.expand_seeded_envelope`，不作为推荐匹配输入 |
| 较强降雨阶段标签 | 可选离线诊断；不修改目录、tracking 或判定对流 |
| persistent worker／tile ID 模式 | 可选工程实验；先在目标机器做 exact/performance 对照 |

主输出在 `chunk_*/`：`objects.csv`、`edges.csv`、`events.csv`、
`family_map.csv`、`state.json`、`metadata.json`、`COMMITTED.json`；顶层有
`execution_contract.json`、`summary.json` 和完成后的 `SUCCESS`。
原始 CSV 的 family 是当时已知映射，未来 merge 可能改变最终 family root；
分析完整生命周期时使用最终 checkpoint 或 canonical export。

同一算法、输入与完整状态重放应保持 ID；不同科学配置重跑后的编号不能直接
当作同一对象来比较。过滤展示不能重新压缩编号。branch 生命周期不等于整个
family 生命周期，split/merge 后换 branch 编号是定义，不必然是错误。

## 6. 不应宣传为已经支持／验证的内容

没有验证云体识别、对流分类、CAPE／updraft 判据、云端托管服务或自动 reset
修复；没有观测真值准确率、跨模型／跨年泛化认证或十年运行保证。
不能把更少断裂／更多 links／更多降雨覆盖直接叫更准确。
0.1 掩膜直接参与 tracking 后事件/gap 大幅变化，仍为实验，不启用作 beta 默认。

推荐先用 4 workers、6–24 小时 chunk、小区间。全域资源随格点数、对象形状与
chunk 变化；已有数据不证明所有机器只需某个固定内存，也不证明 16/32 workers
能等比例加速。checkpoint 历史会随长序列增长，尚无有界历史清理方案。

## 7. 验证记录与反馈

科学与工程依据见 [冻结配置](TRACKING_REFERENCE_FREEZE.md)、
[连续执行与已有 JJA / 冬季测试](REFERENCE_STREAM_EXECUTION.md)、
[外扩对 tracking 的影响](TRACKING_GEOMETRY_SENSITIVITY_2026-10-05.md)。
不同重复执行不是新的独立气象样本。

2026-10-05 的独立新 venv 安装验证：Python 3.13.5、macOS 15.8.1 ARM64，
顶层依赖从包索引安装 binary wheels，使用现代 editable build 成功安装
`step-upgrade==0.4.0b1`，导入路径指向本 beta checkout。`pip check` 无损坏
依赖；一键自检八项 exact checks 全通过；完整回归 **174 passed, 1 skipped**。
跳过项需要另行固定的原版 STEP checkout，不表示已通过原版等价性测试。

干净环境仍出现 NumPy/NetCDF 的 binary-size 警告及 xarray/NumPy-2.5
弃用警告；安装和测试成功不等于消除了原生库兼容性风险。如目标机器出现
导入错误、segfault 或失败自检，先不要用于正式计算。RCC／Linux 环境仍需
组员自检，不能由 macOS 结果代替。

本机安装、自检和测试记录保存在忽略的 `results/beta_release_check_20261005/`，
不随 tag 上传。源码 tag 自身附带合成数据生成器，不需要复制本机输入或输出。

反馈请附：Git tag/commit、`python --version`、`pip freeze`、操作系统／节点、
输入变量单位和时间间隔、preset、完整命令、错误日志或少量问题案例。
不要把原始大数据、凭据或大型结果目录上传 GitHub。运行失败保留原输出，
不要删除证据，也不要在没有 `SUCCESS` 时把部分结果当完整目录。
