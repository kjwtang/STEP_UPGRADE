# 追踪断链诊断（保持算法不变）

版本注意：旧结果只能用产生它的代码重放。升级到0.3后，请先按
[revision 3指南](TRACKING_REVISION3.md)运行 `recheck_tracking.py`，不要直接用新版重放旧结果。
下述默认branch编号针对原先的诊断样本；新版编号可能变化。

输入是已完成的 `original_vs_upgrade_600_mem24` 输出目录。
复用保存的降雨、new_id标签、新版tracking标签、混合组tracking标签和参数，
不读取原始2GB文件、不重新跑原作者tracking。

```bash
python -u scripts/diagnose_tracking.py \
  results/original_vs_upgrade_600_mem24 \
  --output-dir results/diagnosis_600
```

脚本重放新版全部24帧，并逐帧校验与保存的结果完全一致；不一致则停止，
避免用不同代码/参数解释旧结果。默认对01时和05时的所有相邻对象对做诊断，
特别提取图中的 branch 12→26（00→01时）和36→76（04→05时）。
其他过渡可用 `--frames 1 2 3 4 5`，数字为从0开始的子帧序号。

输出 `diagnosis.json`、`REPORT.md` 和 `pair_decisions.csv`。
记录原始质心位移、预测误差、门限、原始/平移IoU、覆盖比例、强度比、评分、
是否进入候选、实际接受边以及关联端点上的其他边。
门限外评分仅为“如果计算会得到什么”的反事实诊断，不代表该pair实际进入过分配。
原因分类：

- `outside_all_candidate_gates`：所有候选入口均未通过（旧版本记为 `outside_prediction_gate`）。
- `below_score_threshold`：进入候选但评分低于tau。
- `endpoint_used_by_other_edge`：合格但端点分配给其他边，需看对应边，不直接称为bug。
- `eligible_but_not_assigned`：合格未被分配，需核查端点竞争和分配结果。
- `accepted_*`：实际边，包括continue/split/merge等，branch换号不一定断链。

原版同ID关系只用于对照，不是真值；原版多个对象共享ID时会形成多对多关系。
脚本另含一个2×2合成反例：tau=.35，分数矩阵[[.90,.34],[.80,0]]，
旧算法先匹配再过滤，选中.80而非可行的.90；revision 3先过滤，预期选中.90。
该反例证明旧分配存在问题，不能据此断言它就是实际样本两次断链的原因。
诊断脚本本身不修改阈值、运动门限或生产追踪代码。
