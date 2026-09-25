# Tracking revision 3：候选门限与匹配修复

本次不修改identification、降雨阈值、score权重、tau=.35、gap阈值或事件阈值。

候选集合为：当前质心附近、预测质心附近、以及原始/平移mask存在实质覆盖的对象的并集。
覆盖阈值沿用event_overlap；即使该阈值设为0，也要求非零像素重叠。
用包围盒中心KD-tree保守筛选潜在重叠，然后确认像素覆盖，不构建雨像素两两距离矩阵。
仅包围盒相交不能让远距离对象进入候选。评分距离仍是预测误差，以隔离此次候选规则变化。

低于tau的候选在分配前删除；有效候选按连通分量分配，每个parent提供不匹配选项。
目标为有效边总得分最大，不强制每个对象匹配。
split/merge仍在一对一分配之前判断；新增候选可能暴露之前漏掉的事件，需要查看图和边表。

状态tracking_config新增algorithm_revision=3。旧版本checkpoint拒绝续跑，必须从序列起点重跑新版。
旧诊断脚本若拿当前算法重放旧结果会停止，这是预期保护；旧诊断结论保留，不覆盖。

## RCC复跑（不重跑原作者版）

login node更新新版仓库：

```bash
cd /project2/moyer/kjwtang/myproject/STEP_validation/STEP_UPGRADE
git pull --ff-only
```

compute node激活新版环境，在仓库目录运行：

```bash
python -u scripts/recheck_tracking.py \
  results/original_vs_upgrade_600_mem24 \
  --output-dir results/tracking_revision3_600
```

此脚本只重新运行新版tracking；识别、小时降雨和两组原版tracking从旧目录只读复用。
复用使用符号链接，请保留源结果目录。追踪从零状态开始。
之后独立按6帧分块、经过JSON checkpoint重载，核对全部栅格及节点/边/事件/family表完全一致。
输出新的三组并排图及edge_changes.json；保存输入hash、源码版本和参数。
新旧branch ID可能整体变化，比较以(frame,local_label)标识对象，不直接按branch数字对齐。

```bash
cat results/tracking_revision3_600/edge_changes.json
```

重点看focus_pairs：旧36→76对应的对象对是否存在新版边（continue或事件边）；
旧12→26评分仍低于tau，本次不为接上它而修改评分阈值。
其他added/removed/changed_event边也必须检查，不能只看目标连接恢复。
若产生大量新增split/merge，需核查是否误把临近系统关联，而不是直接认为事件检测更好。

## 已验证的回归场景

- 运动反转/质心形变使预测失准，但当前位置附近对象仍可连。
- 两种质心门限都未命中，实质像素重叠仍可形成候选。
- 包围盒相交、mask不重叠且远离两质心时不产生候选。
- 原2×2反例选择.90而非.80；所有分数不合格时返回空匹配。
- split/merge、gap、跨chunk一致性和旧checkpoint拒绝。

真实RCC上的效果和性能需通过上述复跑确认；合成测试通过不等于科学验证完成。
