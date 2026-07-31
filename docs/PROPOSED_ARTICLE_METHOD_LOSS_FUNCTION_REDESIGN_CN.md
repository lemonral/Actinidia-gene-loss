# 旧文章标准基因丢失统计与功能展示重构方案

> 状态：**作者已批准并执行。23 单元证据/scaffold 统计、四面板丢失图、
> 逐套功能富集、具体 GO/KEGG 条目图及 39 个 scaffold 节点/末端功能分析
> 均已通过闭合检查。**
>
> 本方案的主分析固定使用旧文章标准：每套基因组中
> `loss = decayed + deleted`。23 套基因组、单倍型和亚基因组始终保留为
> 23 个独立分析单元，不计算物种平均值，也不合并亚基因组。

## 1. 先说明 `retained` 的含义

`retained` 不是丢失类型，也不会进入丢失基因前景集。它表示该参考基因在
目标基因组中存在精确的 SynOrths 锚定同源基因，即有明确的保留证据。

旧文章四种状态继续保持不变：

| 状态 | 旧文章规则 | 是否计入主丢失数 |
|---|---|---:|
| `retained` | 有精确 SynOrths 锚定同源基因 | 否 |
| `decayed` | 缺少精确锚定基因，但全基因组 tBLASTX 命中满足 identity ≥ 50%、bit score ≥ 50、e-value < 1e-5，且不设长度下限 | 是 |
| `deleted` | 属于旧文章缺失候选，但没有达到上述阈值的全基因组 tBLASTX 命中 | 是 |
| `not_called_loss` | 不在旧文章可判定候选范围内 | 否，也不作为确定保留 |

因此：

```text
每套基因组的丢失前景 = decayed + deleted
每套基因组的丢失总数 = decayed 数 + deleted 数
```

功能富集中的背景写成 `retained + decayed + deleted`，只是因为富集检验必须
有一个“本来有机会被判定为丢失的基因全集”。其中：

- 前景集只有 `decayed + deleted`；
- `retained` 只作为未丢失的可检测比较基因；
- `not_called_loss` 因为不可判定，从前景和背景同时排除。

以超几何检验为例，背景总数为 `N`，其中带某功能注释的基因为 `K`；丢失
基因为 `n`，其中带该注释的丢失基因为 `k`。检验的是 `k/n` 是否显著高于
`K/N`。如果背景也只放丢失基因，就会变成“丢失基因与自身比较”，无法得到
有意义的富集结果。**所以 `retained` 参与的是统计机会背景，不参与丢失基因
计数。**

## 2. 当前已经闭合的数据基础

当前旧文章矩阵是完整的 `23 × 35,547 = 817,581` 个单元–基因状态：

- `retained`: 633,957；
- `decayed`: 171,866；
- `deleted`: 7,961；
- `not_called_loss`: 3,797；
- 主丢失阳性 `decayed + deleted`: 179,827 个单元–基因记录；
- 3,616 个参考基因在全部 23 个单元中均为 `decayed` 或 `deleted`。

这些数字是“单元–基因记录”，同一个参考基因可以在多个目标基因组中分别
贡献一条记录。后续图中必须明确区分：

1. 唯一参考基因数；
2. 单元–基因记录数；
3. 树上的丢失事件数。

三者不能混写成同一种“基因数量”。

## 3. 主丢失统计图的重构

建议继续使用：

```text
results/figures/loss_evidence_classification/
```

作为主丢失证据图文件夹。实施时先在版本化临时目录生成并验证，全部检查
通过后再同步正式 PNG/PDF；现已按这一流程完成第一版。

### Panel A：23 套基因组的最粗略总体丢失

- 每一条横向柱代表一套原始基因组/单倍型/亚基因组；
- 橙色：`decayed`；
- 灰色：`deleted`；
- 柱末直接标注 `decayed + deleted` 的实际基因数量；
- 不显示比例作为主横坐标，主横坐标使用具体基因数；
- 23 套数据全部单独保留。

该面板回答：“按照旧文章标准，每套基因组共检测到多少丢失候选，其中多少
为 decayed、多少为 deleted？”

### Panel B：每套基因组中已经确定类型的 `decayed`

删除现有图片中汇总全部样本的旧 Panel B，不再用三根总体柱展示
frameshift/stop。新 Panel B 与 Panel A 使用完全相同的 23 行顺序，每套
基因组、单倍型或亚基因组各一行，不提供合并总体柱。

主图只展示目前具有明确编码破坏证据的三个互斥类别：

1. frameshift only；
2. in-frame premature stop only；
3. frameshift + premature stop。

每行用三段实色堆叠柱显示实际基因数量，柱末标注：

```text
已确定类型数 / 该套 decayed 总数（百分比）
```

因此这个面板回答的是“每套基因组分别有多少 decayed 已经分到哪一种确定
类型”，而不是把全部 23 套合并成一个总体，也不暗示所有 `decayed` 都已经
确定机制。主横坐标使用具体基因数，不用比例替代数量。

三个类别按原始 evidence flag 定义：

```text
frameshift only = frameshift=true 且 premature_stop=false
stop only       = frameshift=false 且 premature_stop=true
both            = frameshift=true 且 premature_stop=true
```

每套基因组必须满足：

```text
三个确定类型互斥；
三个确定类型之和 <= 该套基因组的 decayed 总数；
decayed - 三个确定类型之和 = 尚未得到上述明确编码破坏证据的 decayed。
```

`deleted` 不进入这个机制面板，因为该面板专门回答 `decayed` 中哪些已经
获得进一步的明确类型证据；`deleted` 继续在 Panel A 和 Panel C 中以灰色
显示。

现有的 N-terminal truncation candidate、C-terminal truncation candidate、
both-terminal truncation candidate、partial-alignment candidate 和 residual
sequence present/mechanism unresolved 不放入主图的“确定类型”堆叠柱。它们
作为候选/未解析证据进入补充表，仍按 23 套单元分别报告，不能汇总后冒充
确定机制。如果后续验证得到足够的起始密码子、剪接位点、外显子缺失或基因
裂解证据，再把通过同一证据门槛的新类别加入确定类型主图。

下列机制目前没有系统、同标准证据，因此不在主图中声称已经鉴定：

- 起始密码子丢失；
- 剪接供体/受体破坏；
- 外显子精确缺失；
- 基因融合/裂解；
- 转座子插入；
- 倒位、易位或结构变异断点；
- 启动子、调控区或表观遗传沉默。

如果以后取得独立的 TE、结构变异或转录证据，再增加正交证据层；不能把
“未解释的 decayed”直接改名为这些机制。

### Panel C：shared 与 non-shared

主图采用最直观且不需要物种聚合的定义：

- `shared`：同一参考基因在全部 23 套基因组中均为 `decayed` 或
  `deleted`；
- `non-shared`：该单元中为 `decayed` 或 `deleted`，但不属于上述全
  23 套共同丢失集合。

每套基因组显示四个可加和部分：

1. shared decayed；
2. shared deleted；
3. non-shared decayed；
4. non-shared deleted。

四部分之和必须严格等于 Panel A 的总丢失数。主图使用具体基因数量。比例和
分母保留在配套表中，不用比例替代数量。

### Panel D：树上内部节点与末端分枝

这里采用作者指定的 **23 端并列分析框架**，不再先把同一物种的亚基因组
合并成一个状态：

- 现有 13 谱系树只提供物种之间的骨架关系；
- 在每个物种的末端，把其亚基因组/单倍型展开成一个无内部先后次序的
  并列多分枝；
- *A. arguta* A–D、*A. chinensis* HY4A/HY4P、*A. deliciosa* A–F、
  *A. eriantha* HAP1/HAP2 分别作为独立末端；
- *A. × zhejiangensis* A/B 在本丢失展示框架中也作为两个并列末端；
- 其余单套基因组各保留一个末端；
- 总计正好 23 个末端，每个末端直接对应旧文章矩阵的一套原始数据。

这是为了丢失统计而构建的作者定义 scaffold，不是根据序列重新推断的
23“物种”系统树。图注和方法中会称它们为 23 个
`assembly-unit terminals`，不会声称它们是 23 个独立生物物种，也不会用
它替换正式的 13 谱系定年树。

在该框架上显示三类信息：

1. **单元末端事件**：某一套亚基因组/单倍型自身为
   `decayed/deleted`，直接计在该末端；
2. **共同上游节点事件**：一个并列组或更大单系分枝的所有后代末端都为
   `decayed/deleted` 时，可计在它们的共同上游节点；
3. **重复独立事件**：同一基因需要在两个或更多不相连的末端/节点解释。

因此，同一物种中“部分亚基因组丢失、部分 retained”不再压缩成一个
`partial/homeolog-specific` 物种状态；阳性的亚基因组直接作为各自的末端
丢失，retained 亚基因组保持未丢失。只有所有并列亚基因组均为阳性时，才把
最简共同事件放到该并列组的上游节点。

精确定位一个丢失发生在哪条分枝，必须使用其他分枝的 `retained` 证据作为
“该处未丢失”的信息；否则只能知道某处出现了丢失，不能判断祖先节点还是
多次独立事件。这里使用 `retained` 是为了定位事件，不会改变
`decayed + deleted` 的丢失定义。

图和表必须同时报告：

- 每条分枝的唯一丢失基因数；
- 每条分枝的事件数；
- 重复事件涉及的唯一基因数；
- 23 个独立单元末端各自的事件数；
- 因 `not_called_loss` 而无法精确定位的数。

这样不会把无法定位的基因静默删除，也不会把重复事件重复报告成唯一基因
总数。

### 已解决的 stop-only 统计范围差异

原严格机制总表的 “in-frame stop only” 为 3,266，而旧文章阳性范围内为
3,265。逐行集合差确认，唯一额外记录是
`CS09G1588.mRNA1 × act_arguta_c_legacy`：统一分析有 1 个提前终止密码子，
但旧文章把该行标为 `not_called_loss`，原因是
`not_called_outside_historical_scope`。这不是算法少算或重复，而是母集不同。

进一步按本 Panel B 要求只与旧文章 `decayed` 相交后，正确数目为：

```text
frameshift only = 11,559
in-frame stop only = 3,258
frameshift + stop = 5,071
```

3,266 个统一 stop-only 由 3,258 个旧 `decayed`、7 个旧 `deleted` 和
1 个旧 `not_called_loss` 组成。Panel B 只使用前 3,258 个；其余记录仍在
交叉表中透明保留。实施结果同时完成以下闭合：

```text
每套 decayed = 三个确定类型 + 候选类型 + 未解析余量
每套 positive loss = decayed + deleted
全部 positive unit-gene rows = 179,827
```

这里的第一条用于证据交叉表闭合；主图 Panel B 只画其中三个确定类型，并在
柱末报告其占全部 `decayed` 的比例。

## 4. 功能分析与更详细的功能图

功能分析继续直接使用上述每套基因组的 `decayed + deleted` 丢失基因。
不使用“必须只在一个物种丢失、其他 12 个物种全部 retained”的 1,167
基因严格集合作为主分析。该严格集合只保留为补充敏感性分析。

### 4.1 主功能富集

23 套基因组分别进行 GO/KEGG 富集：

```text
foreground(unit i) = unit i 的 decayed + deleted
background(unit i) = unit i 的 retained + decayed + deleted
```

- 不合并亚基因组；
- 不计算物种均值；
- 不要求其他单元 retained；
- GO biological process、molecular function、cellular component、KEGG
  orthology 和 KEGG pathway 分开校正；
- 完整显著结果和每个条目的贡献基因列表全部保留。

### 4.2 主功能图不再只显示“显著条目数”

目前“每套基因组显著条目数量”的热图只保留为补充 QC。新的出版主图建议：

#### Panel A：主要 GO biological process

- 纵轴：经过客观规则筛选的代表性 GO BP 名称；
- 横轴：23 套独立基因组；
- 点大小：该条目中的实际丢失基因数；
- 点颜色：`-log10(BH q-value)`；
- 配套表给出 fold enrichment、gene ratio 和具体丢失基因 ID。

#### Panel B：GO molecular function 与 cellular component

- MF 和 CC 分成上下两个子区；
- 使用与 Panel A 完全相同的点大小和颜色定义；
- 不把 MF、CC 与 BP 的 q-value 校正混在一起。

#### Panel C：主要 KEGG pathway

- 纵轴直接显示可阅读的 KEGG pathway 名称，不只显示 KO 编号；
- 点大小为丢失基因数；
- 点颜色为 `-log10(BH q-value)`；
- KEGG KO 的完整结果进入补充表；若有少数稳定、可解释的 KO，再在主图
  旁作为小面板展示。

#### Panel D：内部节点/末端分枝的代表性功能

- 使用主丢失统计图 Panel D 中已经确定的旧文章分枝丢失集合；
- 内部节点、末端分枝和重复独立事件分别分析；
- 每个分枝最多显示少量最稳定代表条目，避免文字遮挡；
- 各亚基因组末端分别报告功能；只有所有并列末端共同支持时，才报告其共同
  上游节点功能；
- 没有显著功能的分枝保留为空值，不补造条目。

这不是 PGLS。这里的问题是“某条分枝丢失了哪些功能”，而不是“某个连续
性状是否与丢失率存在系统发育校正后的相关性”。

### 4.3 代表性条目的客观选择规则

为避免人工挑选好看的功能条目，建议预先冻结以下规则：

1. 必须满足 BH `q <= 0.05`、foreground gene count ≥ 2、fold
   enrichment > 1；
2. 主图候选首先按“在多少个独立单元中显著”排序；
3. 再按中位 `-log10(q)`、中位 fold enrichment 和 term ID 排序；
4. GO 中具有祖先–后代关系且贡献基因集合 Jaccard 相似度 ≥ 0.70 的条目
   合并为一个代表簇；
5. 每个被合并条目、原始 q 值和贡献基因仍完整保留在补充表；
6. 建议主图上限：GO BP 10 个、GO MF 8 个、GO CC 6 个、KEGG pathway
   10 个；如果作者认为过多或过少，可在实施前调整。

这样主图会显示“具体富集了什么功能、多少丢失基因支持、显著性有多强”，
而不是只显示每套基因组有多少显著条目。

## 5. 已生成的配套表

已产生以下小型、可核验结果：

1. `unit_loss_summary.tsv`：23 套原始 decayed/deleted/总丢失；
2. `unit_shared_nonshared_summary.tsv`：每套 shared/non-shared 与证据类；
3. `branch_loss_summary.tsv`：内部节点、23 个单元末端、重复事件和无法
   精确定位的数量；
4. `loss_mechanism_crosswalk.tsv.gz`：每条 `decayed` 记录的全部原始机制
   evidence flag 和最终互斥展示类别；
5. `loss_mechanism_summary.tsv`：23 套单元 × 3 个主图确定类型，以及
   候选类型和未解析余量的分层计数；
6. `functional_enrichment_all.tsv.gz`：全部 GO/KEGG 统计；
7. `functional_top_terms.tsv`：主图代表条目和客观排名依据；
8. `functional_term_gene_membership.tsv.gz`：每个条目的具体丢失基因；
9. 路径无关的 manifest、SHA-256 和图形 validation 文件。

所有表都保留 `assembly_unit_id`。只有树事件表额外包含
`biological_lineage` 和 `branch_id`，不会用物种平均数替代原始单元。

## 6. 验证门槛

实施后必须同时满足：

- 23 套单元均存在且没有合并；
- 817,581 条原始状态唯一且完整；
- 每套 `decayed + deleted` 与旧文章矩阵逐行一致；
- Panel A 每套 `decayed + deleted` 与矩阵完全闭合；
- Panel B 的三个确定类型逐套互斥且其和不超过该套 `decayed`；
- Panel B 柱末的已定型数、`decayed` 总数和比例逐套精确一致；
- 补充交叉表中的确定类型、候选类型和未解析余量逐套闭合到 `decayed`；
- Panel C 的四个 shared/non-shared 证据类之和逐套等于 Panel A；
- shared 集合严格为 3,616 个全部 23 套阳性基因；
- 主图确定机制类型互斥；候选/未解析证据单列且不冒充确定机制；
- `deleted` 不进入机制面板；
- 分枝事件数与唯一基因数分开；
- 功能前景只含 decayed/deleted；
- 富集背景中的 retained 不出现在丢失基因列表；
- 所有功能点可追溯到 term–gene 明细；
- 拉丁学名斜体、单元后缀正体；
- 图中不出现 `article`、`method`、`threshold`、内部阶段代码或路径；
- PNG/PDF 目视检查、测试、隐私检查和 checksum 全部通过。

## 7. 已确认并实施的默认值

作者目前已经确认：

1. 主丢失标准采用旧文章的 `decayed + deleted`；
2. 23 套基因组/单倍型/亚基因组分别统计，不做物种聚合；
3. 树上统计采用 23 个单元末端的并列 scaffold，不把它称为重新推断的
   23 物种树；
4. 删除现有图片中汇总全部样本的旧 Panel B，改为逐套展示各自
   `decayed` 中已经确定的 frameshift、premature stop 和两者兼有类型。

已实施：

1. **shared 主定义**：全部 23 套均为 `decayed/deleted`；小分枝共同丢失
   在树面板 Panel D 中另外报告；
2. **功能主图条目上限**：BP 10、MF 8、CC 6、KEGG pathway 10，使用
   预先冻结的客观排序规则。完整 6,420 条显著结果和贡献基因仍保留在配套
   表中，计数热图降为 QC；
3. **scaffold 功能层**：39 个节点/末端的 56,602 个事件–基因记录使用
   33,998 个在全部 23 单元中可判定的基因作为共同背景，共保留 905 条
   显著 GO/KEGG 结果；没有显著结果的节点保留为真实零值。
