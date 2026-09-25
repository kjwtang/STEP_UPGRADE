# RCC 高分辨率验证与性能测试

使用新版独立环境，不安装进旧 STEP 环境。以下命令在仓库根目录运行。
在 login node 安装依赖，在 compute node 跑计算：

```bash
python -m pip install -r requirements-validation.txt
```

## 1. 先确认数据

```bash
python scripts/validate_real_data.py /实际路径/降雨文件.nc --inspect
```

输入支持一个 NetCDF 文件或 `(time,y,x)` NPY 文件。目录不是输入文件。
从元数据确认变量名、维度、时间、单位。下方 `PRECIP`、`time y x` 是示例，
必须替换为检查到的名称。原始累计 RAINNC/RAINC 不可直接使用，需先处理累计量差分和模式重启。
已处理的逐小时 PRECIP(mm) 使用 `--rain-kind interval`，不重复差分。
`--units mm` 仅在你确认单位但文件元数据缺失时使用。
NPY 还必须提供 `--dt-hours 1 --units mm`，其时间连续性无法自动核验。

## 2. 小样本正确性验证

```bash
python -u scripts/validate_real_data.py /实际路径/降雨文件.nc \
  --variable PRECIP --dims time y x --rain-kind interval \
  --start 0 --hours 72 --crop 300 --workers 4 \
  --chunk-frames 24 --sequence-id present_2005_validation \
  --output-dir results/check_300_72h --gif
```

`--hours` 表示帧数，仅对逐小时数据等于小时数。默认 threshold=1 mm/h、radius=9 cells、
最大位移=20 cells/frame 均为待校准起点。不同分辨率之间比较必须调整网格参数，
保持相同物理尺度；这些命令扩大的是原分辨率下的空间范围，没有降采样。

检查 `summary.json` 的所有 checks=true：整段和独立分块识别、标签、最终 family 映射、
节点/边/事件一致；节点唯一；分支每帧唯一；边引用有效且时间前进。
图中数字为 branch ID；split/merge 后换 branch 是设计行为。颜色超过256个分支会复用。
`statistics.png` 给出对象数、面积、分支观测跨度和分支内部速度。
跨度含 gap 且受时段/裁剪边缘和 split/merge 截断，不等于真实风暴寿命。
`frame_*.png` 包括均匀抽样帧、首次各类型事件的前后帧及首次 gap 两端。
GIF 可选，绘图可能慢且占内存。所有坐标是裁剪后网格坐标，不是经纬度。

## 3. 高分辨率流式性能测试

先在 600×600 原始网格上跑24帧，每次只读6帧；再在独立进程中扩大到1000×1000，
最后 `--crop 0` 测全域。每次输出目录必须不同，已有目录不会覆盖。

```bash
/usr/bin/time -v python -u scripts/validate_real_data.py /实际路径/降雨文件.nc \
  --variable PRECIP --dims time y x --rain-kind interval \
  --start 0 --hours 24 --crop 600 --workers 4 --chunk-frames 6 \
  --sequence-id present_2005_benchmark --stream \
  --output-dir results/bench_600
```

流式模式按 chunk 读取、识别、追踪、保存并重新读取 checkpoint；不将全时段栅格放入内存。
边界时间会检查，不允许静默跳过缺帧。追踪在同一序列内串行，只有识别并行。
默认输出每块图表和 checkpoint；`--save-labels` 额外保存每块识别和追踪栅格。
预览图只画第一个 chunk，性能测试通常不加 GIF。全域使用同样机制；先依据小域的内存和对象数决定资源。
完整目录中 `SUCCESS` 文件表示脚本走完，失败目录可以留作诊断，不支持在同一目录自动重跑。

输出：

- `performance.csv`：每块读取、识别、追踪、保存、绘图耗时与对象数、内存峰值。
- `performance.png`：各阶段总耗时和峰值内存趋势。
- `summary.json`：每帧核心耗时、每秒处理百万网格数、对象/事件/gap 数量。
- `frame_statistics.csv`、`statistics.png`：逐帧对象数、区域平均/最大降雨、湿区和缺测比例。
- `chunk_*/`：节点、边、事件、family 表，时间/单位元数据和 checkpoint。

树状进程 RSS 每0.2秒采样，fork共享页面会重复计算，因此是近似观测，不是独占物理内存；
parent_peak_rss_mb 不含识别子进程。结合 `/usr/bin/time -v` 与作业结束后的
`sacct -j JOBID --format=JobID,Elapsed,TotalCPU,MaxRSS,State` 检查资源。
主进程 CPU 时间也不含识别子进程。总体 wall time包含数据输出与首块绘图，
但不包含依赖导入和最后生成性能汇总图的时间；整条命令的时间看 `/usr/bin/time -v`。
核心性能单独使用 identification+tracking。若进程内存查询被系统限制，
`memory_sampling_errors` 会记录原因；无法采样的值为 null，部分采样不能当作完整进程树峰值。

流式测试不自动验证整段等价性，也不是科学效果验收；小样本模式负责该检查。
目前追踪的对象提取、像素重叠和全局匹配仍可能随对象数量迅速变慢。
流式处理限制输入栅格占用，但不能保证任意高分辨率全域都能在32GB内运行。
观察最密集降雨场景，不能仅按网格数线性外推整个JJA时间。
跨月生产仍需外部清单/调度层；本脚本一次接受一个文件，不自动拼接原始小时文件。

逐步扩大范围时记录相同 commit、输入时间、阈值、半径、位移、workers 和 chunk 大小。
必要时额外对 workers=1/4、chunk=3/6 比较，其他条件不变。
代码只检查计算一致性；仍需人工核验相邻系统误连、断轨、分裂合并和gap的物理合理性。
