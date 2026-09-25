# STEP 原作者版 vs STEP_UPGRADE：当前 600×600 区域

此对照固定 RDCEP/STEP commit `ef75c083addcddbcabf6209dc95b1dff89055b0c`。
不把此前 STEP_UPGRADE 的旧 commit 当成原作者版，不修改/安装/覆盖旧 STEP 包。
原始模块按独立模块名加载。脚本核验原始仓库commit及两份源文件与commit内容完全相同。

## 一次性准备（login node）

在目前已经激活的新版独立环境中运行：

```bash
cd /project2/moyer/kjwtang/myproject/STEP_validation/STEP_UPGRADE
git pull --ff-only
python -m pip install -r requirements-validation.txt

git clone https://github.com/RDCEP/STEP.git ../STEP_original_reference
git -C ../STEP_original_reference checkout --detach ef75c083addcddbcabf6209dc95b1dff89055b0c
```

若 `STEP_original_reference` 已存在，先检查它的来源和状态，不覆盖、不强制checkout。
不要在原作者目录执行 pip install。用新版环境的 numpy/scipy/skimage 运行原始算法；
这是原始源代码对照，不声称重现2020年依赖环境。兼容性错误会记入对应run.log，不静默改原代码。

## 执行（compute node）

```bash
source /project2/moyer/kjwtang/myproject/STEP_validation/.venv/bin/activate
cd /project2/moyer/kjwtang/myproject/STEP_validation/STEP_UPGRADE

/usr/bin/time -v python -u scripts/compare_original_step.py \
  /project2/moyer/kjwtang/tmp/CIMP_1hr/2005.07 \
  --original-dir ../STEP_original_reference \
  --output-dir results/original_vs_upgrade_600 \
  --start 0 --hours 24 --crop 600 --threshold 1 --radius 9 \
  --workers 4 --timeout-seconds 900 --memory-gb 8
```

范围与刚才验证相同：中央600×600、前24个小时增量（读25个累计快照）、1mm/h阈值、disk半径9。
三组输入共享一次读取和差分后保存的同一降雨数组，SHA256、源文件清单和时间都记录。
各算法接收同样阈值截断后的降雨，min_size固定1。若存在缺测则停止对照，
因为原作者接口没有缺测屏障，不能偷偷填0再宣称同条件比较。
空间分辨率保持原样；脚本拒绝crop大于600或全域模式。

组合与目录：

| 组合 | identification | tracking | 目的 |
|---|---|---|---|
| A 原版 | original_id | original_track | 原作者端到端基线 |
| B 混合诊断 | new_id | new_id_original_track | 与C共用识别对象，隔离tracking差异 |
| C 新版 | new_id | new_track | 当前STEP_UPGRADE |

原版tracking不支持当前checkpoint协议，因此对24小时整体调用，不逐6小时重启。
新版在此也对同样24小时整体调用。识别共享计算一次，不把复用结果的时间当作0成本；
A核心总耗时=original_id+original_track，B=new_id+new_id_original_track，C=new_id+new_track。

## 参数与解释

原版默认 `tau=0.7, phi=0.003` 来自[原作者教程](https://github.com/relttira/STEP/wiki/Tutorial)，
只是启动对照的参考，不是本项目逐小时4km数据的已校准参数。其教程原始场景为3小时降雨。
`--original-km 30` 按120km/4km换算；原版还有方向判断回退，所以不是严格位移上限。
新版保留已跑实验的 `--new-tau 0.35 --new-displacement 20`，gap/event参数与现有验证一致。
两种算法的score及运动规则不同，相同tau数字不等价；本对照是明确配置的基线比较，不是参数已公平最优的结论。
如另有已经校准的原版参数，可以显式传入 `--original-tau / --original-phi / --original-km`，另建输出目录。
阈值和识别disk故意保持当前实验值，未照搬教程中的0.6和其他disk形状。

## 超时与内存

每个算法阶段单独子进程，默认最多900秒；共5阶段，串行执行，最坏约75分钟加I/O和出图。
8GiB限制包括Linux进程虚拟地址空间硬上限，以及0.2秒采样的进程树RSS限额。
进程树RSS可能重复计算共享页，是保守限制；不是精确独占物理内存。
识别4个worker也计入该阶段RSS。建议沿用4CPU/32GB并确保计算节点剩余时间足够。
外层 `/usr/bin/time` 的MaxRSS不是所有子进程独占内存总和。

原版超时/超内存/异常时停止该子进程，保留已完成结果，继续可独立执行的其他组。
不自动缩小区域、不改变参数、不将缺失结果记为0秒或成功。
原版像素两两矩阵可能很快触及内存限制；相应失败本身是可报告的可扩展性结果，
但不能据此声称新版识别/追踪更准确。`run.log` 可区分MemoryError和其他代码/依赖错误。
`SUCCESS`表示5阶段全部完成，`PARTIAL`表示存在缺失，程序返回码2。

## 输出检查

- `identification_000.png`至`005`：相同降雨、原版分组、新版分组。
- `tracking_000.png`至`005`：降雨、A/B/C的编号并排。列之间颜色/数字独立，不按相同数字判断是否同一系统。
- `identification_comparison.csv`：逐帧对象数、湿区IoU、共同雨像素分组ARI；全空/极少共同雨像素时ARI为空。
- `tracking_link_disagreements.json`：B与C在相同识别对象上的相邻帧同ID关联差异。
- `tracking_summary.json`：同ID关联数量、同一帧多个对象共享编号的次数、栅格ID数。
- `comparison_summary.png`、`stage_results.json`：完成阶段耗时，超时/失败标记，内存峰值。
- 各阶段目录：结果labels.npy、run.log、参数spec、性能记录；new_track还保存完整事件图表。
- `configuration.json`、`input_metadata.json`：版本、源码/输入hash、气象时间及输入文件来源。

首6帧预览与之前给出的6张图对应，统计覆盖24帧；`--plot-frames 24`可画完整时段。
不能只比较ID总数和ID持续时间：新版split/merge会开新branch，原版可能让多个对象共享ID，
旧编号也可能在以后再利用。因此目前不将栅格ID跨度称作真实风暴寿命或断轨率。
同ID关联差异表不包含新版split/merge跨branch边；应结合new_track/edges.csv人工检查。
先看右侧主要雨带在同一时刻是否分组不同，再查B/C是否对同一对象保留或丢失连接。
