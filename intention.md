# 胃腺癌分型大作业：实施方案与报告骨架

> 定位：本科课程大作业，无 GPU、无服务器，全程笔记本 CPU + 公开数据即可完成。  
> 核心策略：**不做病理大图（WSI）深度学习，做 TCGA-STAD 多组学轻量分析**——模板要求的"数据科学探索、特征工程、多模态整合、鲁棒性评估、临床关联"全部能在表格数据上闭环。

---

## 一、总体降级思路：论文级 vs 本科大作业级

| 模板要求   | 论文级做法     | 本科降级做法（本方案）                                          |
| ------ | --------- | ---------------------------------------------------- |
| 数据科学探索 | 多中心百万级数据  | TCGA-STAD（约440例胃腺癌）+ 外部验证集 ACRG/GSE62254（约300例）      |
| 特征工程   | 深度特征/自监督  | 方差/MAD 过滤选高变异基因 + z-score 标准化                        |
| 多模态整合  | 跨注意力/基础模型 | mRNA + DNA甲基化 + 临床三张表，特征级拼接，与单模态做对比实验                |
| 分型方法   | 基础模型微调    | 无监督（层次聚类/一致聚类）+ 有监督验证（随机森林/SVM）                      |
| 鲁棒性评估  | 多中心前瞻验证   | 5折交叉验证 + 外部队列验证 + 子采样聚类稳定性 + bootstrap 置信区间          |
| 临床关联   | 前瞻性队列     | Kaplan-Meier 生存曲线 + log-rank 检验 + 卡方检验（亚型×Lauren/分期） |
| 可视化    | 交互式大屏     | PCA/UMAP、热图、ROC、KM 曲线、混淆矩阵（matplotlib/seaborn）       |

**算力评估**：440 样本 × 2万基因的表格，随机森林在笔记本 CPU 上秒级~分钟级完成。不需要任何 GPU。

---

## 二、实验流程总览

> 一句话概括：**通过 TCGA-STAD 的 mRNA 表达 + DNA 甲基化 + 临床变量（输入），预测每个病人的分子亚型（输出：四亚型 EBV/MSI/GS/CIN，或 MSI-H vs MSS 二分类），再用无监督结构、鲁棒性验证与临床关联，证明这个分型"可预测、稳健、有临床意义"。**

### 任务定义：通过什么预测什么

分型的本质是标准多分类任务：每个病人由多组学特征向量 x_i 描述（~1500 维），预测亚型标签 y_i（i = 1..约350）。两条实现路线：

|       | 形态 A：发现型（无监督）  | 形态 B：预测型（有监督，主线）        |
| ----- | -------------- | ----------------------- |
| 标签    | 没有，聚类自己"发现"亚型  | 用现成标签：四亚型 或 MSI 状态（已在手） |
| 回答的问题 | "数据里天然存在几类病人？" | "能否从组学数据预测病人已知亚型？"      |
| 难度    | 需自己解释簇的意义      | 标签现成，标准监督分类             |
| 定位    | 兜底 + 辅助验证      | 默认主线                    |

**主线选择（数据核对后确定，按优先级）**：

1. 首选：**MSI-H vs MSS 二分类**——标签就在 clinicalMatrix 的 `CDE_ID_3226963` 列，零外部依赖，且是 Kather 2019（Nat Med）的经典任务，文献背书充分；MSI-L 剔除或并入 MSS（报告中注明）
2. 升级：若从 cBioPortal 导出四亚型标注，做 EBV/MSI/GS/CIN 四分类（更贴"分型系统"主题）
3. 兜底：都拿不到 → 形态 A 聚类闭环，后续验证流程不变

**标签的临床价值**（为什么值得预测）：MSI-H → 高突变负荷，免疫检查点抑制剂响应好；EBV → PD-L1/L2 扩增，免疫治疗潜在获益；HER2 阳性 → 曲妥珠单抗靶向治疗；Lauren 弥漫型 → 预后差，化疗方案选择不同。

### 流程步骤表

| Step      | 做什么                                                                | 输入 → 输出                     | 产出图表                  |
| --------- | ------------------------------------------------------------------ | --------------------------- | --------------------- |
| 0 数据获取    | Xena 下四个文件 + cBioPortal 导出亚型标注 + GEO 下载 GSE62254                   | — → 原始数据包                   | —                     |
| 1 对齐与清洗   | 条形码截前 12 位做 join key；只留 01（原发肿瘤）样本；各表取交集；Ensembl→Symbol 映射         | 四张表 → 一张主表                  | 样本数交代表（Venn/UpSet）    |
| 2 EDA     | 表达分布、缺失值、各亚型样本量                                                    | 主表 → 数据质量结论                 | 缺失热图、箱线图、亚型条形图        |
| 3 特征工程    | MAD top 1500 基因 / top 变异 CpG；z-score；临床变量 one-hot                  | 主表 → X_rna / X_meth / X_cli | —                     |
| 4 无监督结构验证 | UMAP/PCA 按官方亚型着色；层次聚类 + 一致性矩阵                                      | X_rna → 二维散点 / 簇            | UMAP 散点图、聚类热图         |
| 5 主线建模    | 逻辑回归 / SVM / 随机森林三模型对比，选最优，粗调参，5 折 CV                              | X, y（亚型标签）→ 最优模型            | 模型对比表、CV AUC          |
| 6 多模态消融   | 仅 mRNA → +甲基化 → +临床，三档递进对比                                         | 各模态组合 → AUC 对比              | 消融对比表                 |
| 7 鲁棒性评估   | GSE62254 外部验证（交集基因、同标准化、不重训直接预测）；特征数敏感性（top 1000/2000/5000）        | 模型 → 外部 AUC                 | 外部验证表、敏感性折线图          |
| 8 解释与临床关联 | 混淆矩阵（哪两型易混+文献解释）；特征重要性 top 基因对照已知标志；KM + log-rank；亚型 × Lauren/分期卡方 | 模型 + 数据 → 结论                | 混淆矩阵、特征重要性图、KM 曲线、卡方表 |
| 9 报告成稿    | 按模板五部分组织以上产出                                                       | 全部产出 → 报告                   | —                     |

### 三个关键决策点

1. **D1 必须先确认亚型标注能拿到**：cBioPortal 导出 TCGA-STAD 的分子亚型列（或 Nature 2014 原文补充材料）。拿到 → 走"预测型"主线（本方案默认）；拿不到 → 切换"发现型"（Step 4 的聚类升为主线，自己定义亚型标签，Step 5-8 验证流程完全不变）。
2. **Step 6 是模板"多模态整合"要求的核心证据**：无论甲基化的加入是否提升 AUC 都要写——提升了说明表观层有增量信息，没提升说明转录层已充分，两种结果都是有效结论。
3. **每一步的产出图表直接进结果章节**：实验做完时结果章节的图也同步做完，不要实验和写作两张皮。

### 流程 ↔ 模板章节 ↔ 证据链对应

| 实验步骤 | 报告章节              | 证据链层级          |
| ---- | ----------------- | -------------- |
| 0-2  | 研究方法（1）数据科学探索     | 1 数据可信         |
| 3    | 研究方法（2）特征工程       | 1 数据可信         |
| 4    | 结果（4）可视化与比较       | 2 方法合理         |
| 5-6  | 研究方法（3）（4）+ 结果    | 2 方法合理         |
| 7    | 结果（1）分型鲁棒性评估      | 3 结果稳健         |
| 8    | 结果（2）（3）临床关联与结果解释 | 4 可解释 / 5 临床有用 |

### Step 0 实际产出：已下载数据盘点（2026-10-25 核对）

数据目录：`~/Downloads/生物信息学实验数据/`，9 个文件全部核对完毕。

| 文件                                                  | 内容           | 规模（行×列）                 | 值类型                 | 用途                        |
| --------------------------------------------------- | ------------ | ----------------------- | ------------------- | ------------------------- |
| TCGA.STAD.sampleMap_HiSeqV2.gz                      | mRNA 表达谱     | 20,530 基因 × 450 样本      | log2 归一化表达，约 0~15   | **模态一（主力）**，行是基因 Symbol   |
| TCGA.STAD.sampleMap_HumanMethylation450.gz          | DNA 甲基化      | 485,577 CpG 探针 × 398 样本 | β 值 0~1             | **模态二**，需 probeMap 聚合到基因级 |
| probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy | 甲基化探针注释      | 395,985 探针              | cg 探针 → 基因/坐标映射     | 甲基化特征工程用（首列列名带 #）         |
| mc3_gene_level_STAD_mc3_gene_level.txt.gz           | 体细胞突变（基因级计数） | 40,543 基因 × 439 样本      | 稀疏计数 0/1/2...       | 可选模态 + MSI sanity check   |
| Gistic2_all_thresholded.by_genes.gz                 | 拷贝数变异        | 24,776 基因 × 441 样本      | -2/-1/0/1/2 五档      | 可选模态                      |
| miRNA_HiSeq_gene.gz                                 | miRNA 表达     | 2,178 × 428             | log2，**大量 NA**      | 不建议用（缺失太多）                |
| RPPA.gz                                             | 蛋白/磷酸化蛋白     | 227 蛋白 × 357 样本         | 标准化信号               | 可选                        |
| STAD_clinicalMatrix                                 | 临床表          | 580 样本 × 108 变量         | 混合                  | 模态三 + 标签来源                |
| survival_STAD_survival.txt                          | 生存随访         | 511 样本                  | OS/DSS/DFI/PFI 四个终点 | 生存分析                      |

**核对后的关键发现：**

1. **MSI 标签已经在手**：clinicalMatrix 第 2 列 `CDE_ID_3226963` 就是 MSI 状态（MSS 352 / MSI-H 115 / MSI-L 75 / 缺失 38）。预测型主线可立即落地：**MSI-H vs MSS 二分类**（Kather 2019 的经典任务）；MSI-L 建议剔除或并入 MSS（需在报告中注明）。
2. **四亚型（EBV/MSI/GS/CIN）标签不在这些文件里**：仍需 cBioPortal 导出；若不想再折腾，"MSI 状态预测"本身就是标准且有文献背书的分型任务。
3. **各文件样本数不同**（450/398/441/428/357/439/580/511），且 HiSeqV2 表头混有 `-11` 正常样本——一律"截条形码前 12 位 + 只留 01 + 取交集"处理，预计入模 ~350 例。
4. **HiSeqV2 的行是基因 Symbol**（不是 Ensembl ID），与 GSE62254 对齐反而更省事。
5. clinicalMatrix 的 `histological_type` 含 Intestinal / Diffuse / Signet Ring 等，Lauren 分型可从这里近似提取（临床关联分析用）。
6. **甲基化 48.5 万探针不能直接与 mRNA 拼接**：先用 probeMap 把探针注释到基因（gene 列为 "." 的丢弃），按基因聚合（如取均值或变异最大的探针），降到基因级后再选 top 特征。

### 数据文件结构详解（2026-10-25 实测）

九个文件只有**三种结构**：

| 家族        | 结构             | 文件                                      | 对齐方式           |
| --------- | -------------- | --------------------------------------- | -------------- |
| A：特征×样本矩阵 | 行=特征，列=病人      | HiSeqV2、甲基化、Gistic2、mc3、miRNA、RPPA（6 个） | 按列名（条形码）join   |
| B：样本×变量表  | 行=病人，列=字段      | clinicalMatrix、survival                 | 按行索引 join      |
| C：注释映射表   | 行=探针，**无病人维度** | probeMap                                | 甲基化特征工程的"翻译字典" |

注意 A 与 B 方向相反：矩阵靠列名 join，临床/生存表靠行索引 join——写对齐代码时唯一要区分的事。

**家族 A 逐文件实测结构：**

| 文件      | 首列特征              | 样本数        | 值域（实测）                                   | 缺失    |
| ------- | ----------------- | ---------- | ---------------------------------------- | ----- |
| HiSeqV2 | 基因 Symbol（20,530） | 450（含 -11） | log2，~0-15                               | **0** |
| 甲基化450  | cg 探针 ID（485,577） | 398        | β∈[0,1]                                  | 少量    |
| Gistic2 | 基因 Symbol（24,776） | 441        | 整数档位 {-2,-1,0,1,2}                       | 无     |
| mc3     | 基因 Symbol（40,543） | 439        | **二值 {0,1}**（实测全表最大值=1，是"有无突变"指示矩阵，不是计数） | 无     |
| miRNA   | MIMAT ID（2,178）   | 428        | log2                                     | 大量 NA |
| RPPA    | 蛋白名（227）          | 357        | 标准化信号                                    | 少量    |

**家族 B：**

- clinicalMatrix：580 行 = **478 原发肿瘤 + 102 癌旁正常**（`sample_type` 列可直接过滤，比手解析条形码更稳），108 列。关键列：`CDE_ID_3226963`（MSI）、`histological_type`（Lauren 近似）、`pathologic_stage`/`pathologic_T`/`N`/`M`、`neoplasm_histologic_grade`（G1-G3）、`age_at_initial_pathologic_diagnosis`、`gender`、`h_pylori_infection`；另有 24 个 `_GENOMIC_ID_*` 管理列可忽略。缺失因列而异：分期缺 101、分级缺 69、MSI 缺 38（EDA 素材）。
- survival：511 行，四套终点（OS/DSS/DFI/PFI），每套"事件 0/1 + 天数"两列。主用 OS，DFI 缺失较多。

**家族 C：** probeMap 六列（id/gene/chrom/start/end/strand），实测 **344,298 / 395,985（87%）探针注释到基因**，gene 列为 "." 的丢弃。首列列名带 `#`，读入时处理。

**TCGA 条形码**：`TCGA-D7-6871-01` = 项目 - 医院代码 - 患者编号 - 样本类型 - …。**前 12 位（TCGA-D7-6871）是患者唯一键**；第 4 段 01=原发肿瘤、11=癌旁正常。同一患者的所有数据文件共享前 12 位。

---

## 三、五个部分逐项怎么写

### 1. 团队介绍（最简单，如实填写）

- 成员表：姓名、学号、专业
- 分工表（按人数裁剪，单人也可全包）：

| 角色        | 职责                       |
| --------- | ------------------------ |
| 数据获取与清洗   | Xena/cBioPortal 下载、缺失值处理 |
| 建模与分析     | 聚类、分类器、交叉验证              |
| 可视化与统计检验  | 画图、KM 曲线、卡方检验            |
| 报告撰写与文献调研 | 研究背景、结论、参考文献             |

- "应到/实到人数"：按课程要求如实填（一般是防代做签到的字段，写实际人数即可）。

### 2. 研究背景（重点，见第四节完整草稿）

五段式漏斗结构：**疾病负担 → 传统形态分型 → 分子分型 → 数据科学/AI 方法 → 现状不足与本文切入点**。每段 3-5 句，共约 600-800 字，配 10 篇左右参考文献。草稿见下文第四节，可直接改写使用。

### 3. 研究方法

按模板四个关键词逐一对应写，让阅卷人一眼看到全覆盖：

**（1）数据科学探索**

- 数据源：UCSC Xena（<https://xenabrowser.net）下载> TCGA-STAD 三个文件：
  - `TCGA-STAD.tpm.tsv`（mRNA 表达谱）
  - STAD 的 HumanMethylation450 甲基化谱（Xena 搜索 "STAD 450"）
  - `TCGA-STAD.GDC_phenotype.tsv` + `TCGA-STAD.survival.tsv`（临床与生存）
- EDA：样本分布、缺失值热图、表达分布箱线图、批次检查。这一节放 2-3 张图就是"数据科学探索"的实锤。

**（2）特征工程**

- MAD（中位绝对偏差）排序取 top 1500 高变异基因（甲基化同理取 top 变异 CpG 位点）
- 每个特征 z-score 标准化
- 说明为什么：胃腺癌样本仅 ~440 例，原始 2 万维直接进模型会过拟合，必须降维。这一句就是特征工程的动机论述。

**（3）多模态数据整合**

- 三模态：mRNA（转录组）、DNA 甲基化（表观组）、临床变量（分期/分级/年龄）
- 整合方式（本科够用）：特征级拼接后进模型
- **加分实验设计**：单模态（仅mRNA）vs 双模态 vs 三模态的性能对比表——这一张表就把"多模态整合与处理细节"写满了，而且结论通常是"加入甲基化略提升/不提升"，无论结果如何都有的写。

**（4）具体分型方法**

- 主线（有监督预测）：以 TCGA 官方分子亚型（EBV/MSI/GS/CIN）为标签——标注从 cBioPortal 的 TCGA-STAD 页面导出，或取 Nature 2014 原文补充材料——训练逻辑回归/SVM/随机森林做亚型预测，5 折分层交叉验证
- 辅助（无监督结构验证）：UMAP/PCA 按亚型着色 + 层次聚类，证明四亚型在数据空间天然可分（这是"可预测性"的结构性证据，也回答"分型是不是数据里真实存在的"）
- 兜底方案：若亚型标注获取困难 → 切换"发现型"路线（无监督聚类自己定义亚型标签，后续验证流程不变，报告里注明即可）
- 工具：`pandas + scikit-learn + lifelines + seaborn`，全部 pip 安装

### 4. 结果（对应模板四个要求）

**总原则：模型是"工具"，证据链才是主角。** 模型好坏不是报一个准确率，而是五级证据链：

> 1 数据可信（EDA）→ 2 方法合理（基线对比）→ 3 结果稳健（严格验证）→ 4 可解释（打开黑箱）→ 5 临床有用（KM/卡方）

结果章节证据清单（每张图表回答一个问题）：

| 产出                             | 回答的问题         | 模板要求   |
| ------------------------------ | ------------- | ------ |
| 样本数交代表 + EDA 图                 | 数据干净吗、样本够吗    | 数据科学探索 |
| 基线对比表（随机基线/仅临床/单模态/多模态）        | 比不动脑子的做法好吗    | 多模态整合  |
| 5 折 CV（均值±标准差）+ 外部 GSE62254 验证 | 换一批数据还准吗      | 分型鲁棒性  |
| 特征数敏感性（top 1000/2000/5000）     | 依赖随手选的参数吗     | 鲁棒性    |
| 混淆矩阵（哪两型易混 + 文献解释为什么）          | 模型的错误有规律吗     | 结果解释   |
| 特征重要性 top 基因 vs 已知文献标志         | 学到的是生物学还是噪声   | 结果解释   |
| KM 曲线 + log-rank + 卡方          | 分出的亚型临床上有什么区别 | 临床关联   |

**模型侧做到什么程度就够**：逻辑回归/SVM/随机森林三模型一张对比表（选型靠对比，不靠设计）；`n_estimators`/`max_depth` 各试 3-4 个值粗调即可，不搞 GridSearch 玄学；消融实验（单→双→三模态）必做。

**加分项——MSI sanity check**：mc3 突变矩阵按行求和 = 每个病人的"突变基因总数"。MSI-H 组的总突变数应显著高于 MSS 组（MSI 的定义就是突变累积）；模型预测出的 MSI-H 组突变负荷也应更高——模型输出与独立生物学证据互证，写进结果解释很出彩。

**（1）分型鲁棒性评估**

- 5 折分层交叉验证，报告 AUC 均值±标准差
- 外部验证：ACRG 队列 GSE62254（GEO 免费下载，约300例，带 Lauren 分型与生存数据），模型不重训直接预测
- 聚类稳定性：随机子采样 80% 样本重复聚类 100 次，计算共聚类比例（简化版 consensus clustering），画一致性矩阵热图

**（2）临床特征关联**

- Kaplan-Meier：各亚型总生存曲线 + log-rank 检验（`lifelines`）
- 卡方/Fisher 检验：亚型 × Lauren 分型 / TNM 分期 / 年龄 / 性别
- 若发现某亚型富集弥漫型或晚期，就是现成的"临床适用性"论据

**（3）结果解释**

- 对每个亚型取 top 差异表达基因，对照文献验证：MSI 样富集免疫相关基因、GS 样富集细胞黏附/ECM 基因——和 TCGA 原文对得上，解释就成立了

**（4）可视化与比较**

- PCA/UMAP 二维散点（按亚型着色）
- top 基因 × 样本 热图（按亚型分组）
- ROC 曲线、混淆矩阵、KM 曲线
- 与文献比较表：本作业结果 vs TCGA 四亚型 vs ACRG 四亚型的对应关系

### 5. 结论

写法公式（客观、不吹）：

> 本文基于 TCGA-STAD 公开多组学数据，用轻量数据科学流程复现了胃腺癌分子分型的"亚型发现—临床关联—鲁棒性验证"链条。结果显示……（填实际数字）。与 TCGA 原文相比……（一致/部分一致，差异可能源于样本量与特征选择策略）。局限性：样本量有限（~440例）、仅两个组学模态未纳入病理图像、存在平台批次效应、结论为关联性而非因果。后续可在更大队列和病理图像模态上扩展。

结论只陈述"做了什么、发现了什么数字、和文献对不对得上"，再老老实实列 3-4 条局限——这就是"客观、合理"。

---

## 四、研究背景草稿（可直接改写）

**P1 疾病负担与临床需求**  
胃癌是全球最常见的恶性肿瘤之一。GLOBOCAN 2022 估计全球年新发约 96.8 万例、死亡约 66.0 万例，均居恶性肿瘤第五位[1]。胃腺癌是其最主要的组织学类型。胃癌临床表现与治疗反应的高度异质性是疗效不佳的核心原因，精确分型因此成为预后判断与个体化治疗的基础。

**P2 传统形态学分型**  
经典分型体系包括 Borrmann 大体分型和 Lauren 组织学分型。Lauren（1965）将胃癌分为肠型、弥漫型与混合型[2]，肠型多见于老年男性、与环境因素相关，弥漫型发病年龄更轻、预后更差。WHO 消化系统肿瘤分类第5版（2019）进一步按形态学细分为管状腺癌、乳头状腺癌、黏液腺癌、印戒细胞癌等[3]。形态学分型高度依赖病理医师经验，观察者间一致性有限，且难以反映肿瘤内在的分子异质性。

**P3 分子分型**  
TCGA（2014）整合六种组学平台，将胃腺癌分为四个分子亚型：EBV 感染型（富集 PIK3CA 突变与 PD-L1/PD-L2 扩增）、微卫星不稳定型（MSI，高突变负荷、MLH1 沉默）、基因组稳定型（GS，富集弥漫型组织学、CDH1/RHOA 突变）和染色体不稳定型（CIN，显著非整倍体与受体酪氨酸激酶扩增）[4]。亚洲癌症研究组（ACRG）基于表达谱提出 MSI、MSS/TP53+、MSS/TP53−、MSS/EMT 四型，其中 EMT 型预后最差[5]。分子分型已开始指导临床决策，如 EBV 与 MSI 状态和免疫治疗响应相关。

**P4 数据科学与 AI 在分型中的应用**  
数据科学方法沿两条主线发展。组学方向：MoXGATE 用跨注意力机制融合 mRNA、DNA 甲基化与 miRNA 进行胃肠腺癌分子亚型分类，准确率达 95%[6]。病理图像方向：RMDL 提出重校准多示例学习框架用于全切片胃癌图像分类[7]；MKD 通过低秩多模态融合将病理报告文本知识蒸馏到图像模型，实现仅用图像的亚型推断[8]；GRACE 构建胃癌专科病理基础模型，可同时预测 Lauren 分型、组织学分级与 EBV 状态[9]。系统性基准 PathBench 显示各病理基础模型在胃癌 31 项分型相关任务上表现分化明显（HER2 预测最佳 AUC 仅约 0.62），提示任务仍有较大提升空间[10]。此外，Kather 等证明深度模型可从常规 H\&E 染色图像预测 MSI 状态并在包括胃癌在内的多癌种上验证[11]。


**P5 现状不足与本文切入点**  
现有研究存在三方面制约：其一，病理全切片图像为千兆像素级，深度模型训练依赖大量算力与精细标注；其二，多模态整合方法仍以特征拼接或重量级注意力为主，轻量方案的系统评估不足；其三，多数研究缺乏跨队列外部验证。本作业以公开 TCGA-STAD 多组学数据为对象，采用轻量数据科学方法复现"分子亚型发现—临床特征关联—鲁棒性验证"的完整链条，重点考察流程完整性与结果可解释性，而非追求新方法。

---

## 五、参考文献（编号与上文对应）

[1] Bray F, et al. Global cancer statistics 2022: GLOBOCAN estimates of incidence and mortality worldwide for 36 cancers in 185 countries. CA Cancer J Clin. 2024;74(3):229-263.

[2] Lauren P. The two histological main types of gastric carcinoma: diffuse and so-called intestinal-type carcinoma. Acta Pathol Microbiol Scand. 1965;64:31-49.

[3] WHO Classification of Tumours Editorial Board. WHO Classification of Tumours: Digestive System Tumours. 5th ed. Lyon: IARC; 2019.

[4] Cancer Genome Atlas Research Network. Comprehensive molecular characterization of gastric adenocarcinoma. Nature. 2014;513(7517):202-209.

[5] Cristescu R, et al. Molecular analysis of gastric cancer identifies subtypes associated with distinct clinical outcomes. Nat Med. 2015;21(5):449-456.

[6] Acharjee Dip S, et al. MoXGATE: Modality-aware cross-attention for multi-omic gastrointestinal cancer sub-type classification. arXiv:2506.06980. <https://arxiv.org/abs/2506.06980>

[7] Wang S, et al. RMDL: Recalibrated multi-instance deep learning for whole slide gastric image classification. Med Image Anal. 2021;70:101549. arXiv:2010.06440. <https://arxiv.org/abs/2010.06440>

[8] Dumbre S, Santra B. Multimodal knowledge distillation for gastric adenocarcinoma classification from whole-slide images. arXiv:2610.07913. <https://arxiv.org/abs/2610.07913>

[9] GRACE: A pathology foundation model for gastric cancer with real-world validation. arXiv:2606.04792. <https://arxiv.org/abs/2606.04792>

[10] PathBench: A comprehensive comparison benchmark for pathology foundation models towards precision oncology. arXiv:2505.20202. <https://arxiv.org/abs/2505.20202>

[11] Kather JN, et al. Deep learning can predict microsatellite instability directly from histology in gastrointestinal cancer. Nat Med. 2019;25(7):1054-1056.

> 注：arXiv 条目建议提交前逐一点开核对标题/作者；\[9\]\[10] 若需正式引用格式可到 arXiv 页面点 "Cite as" 复制。

---

## 六、关键技术路线（代码骨架）

```python
# 1. 数据准备（实际下载文件，2026-10-25 核对）
import pandas as pd
DATA = "/Users/foushen/Downloads/生物信息学实验数据/"
expr = pd.read_csv(DATA+"TCGA.STAD.sampleMap_HiSeqV2.gz", sep="\t", index_col=0).T  # 450样本×20530基因
clin = pd.read_csv(DATA+"TCGA.STAD.sampleMap_STAD_clinicalMatrix", sep="\t", index_col=0)
surv = pd.read_csv(DATA+"survival_STAD_survival.txt", sep="\t", index_col=0)
y_msi = clin["CDE_ID_3226963"]          # MSI 标签：MSS / MSI-H / MSI-L
lauren = clin["histological_type"]      # Lauren 分型近似来源

# 甲基化：先按 probeMap 聚合到基因级，再取 top 变异
meth = pd.read_csv(DATA+"TCGA.STAD.sampleMap_HumanMethylation450.gz", sep="\t", index_col=0)
pmap = pd.read_csv(DATA+"probeMap_illuminaMethyl450_hg19_GPL16304_TCGAlegacy",
                   sep="\t", comment="#", header=None,
                   names=["id","gene","chrom","start","end","strand"])

# 2. 特征工程：MAD top 1500
mad = expr.mad()
top_genes = mad.nlargest(1500).index
X = (expr[top_genes] - expr[top_genes].mean()) / expr[top_genes].std()

# 3. 辅助：无监督结构验证（UMAP + 层次聚类）
from sklearn.cluster import AgglomerativeClustering
labels = AgglomerativeClustering(n_clusters=4).fit_predict(X)

# 4. 主线：有监督亚型预测（随机森林 + 5折CV）
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score, StratifiedKFold
#  y 为官方分子亚型标注（cBioPortal 导出；拿不到则用第3步聚类标签兜底）
rf = RandomForestClassifier(n_estimators=500, n_jobs=-1)
auc = cross_val_score(rf, X, y, cv=StratifiedKFold(5), scoring="roc_auc_ovr")

# 5. 多模态对比实验
#    X_rna / X_meth 拼接前后分别跑第3-4步，整理成对比表

# 6. 生存分析
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
for k in sorted(set(labels)):
    mask = labels == k
    KaplanMeierFitter().fit(surv.OS[mask], surv.OS.time[mask], label=f"subtype{k}").plot()

# 7. 聚类稳定性（简化版 consensus）
import numpy as np
agree = np.zeros((len(X), len(X)))
for _ in range(100):
    idx = np.random.choice(len(X), int(len(X)*0.8), replace=False)
    sub = AgglomerativeClustering(n_clusters=4).fit_predict(X.iloc[idx])
    agree[np.ix_(idx, idx)] += (sub[:, None] == sub[None, :])
agree /= agree.max()   # 画热图即为一致性矩阵
```

外部验证：GEO 下载 GSE62254（ACRG 队列），取交集基因 → 同样标准化 → 用 TCGA 训练的 RF 直接预测 → 报告外部 AUC + KM 分层。

---

## 七、工作量与排期（单人 5-7 天）

| 天  | 任务（对应实验流程 Step）                          |
| -- | ---------------------------------------- |
| D1 | 文献调研 + 研究背景成稿；数据下载（Step 0）；**确认亚型标注可导出** |
| D2 | 数据对齐清洗 + EDA（Step 1-2）                   |
| D3 | 特征工程 + 无监督结构验证（Step 3-4）                 |
| D4 | 主线建模 + 5折CV + 多模态消融（Step 5-6）            |
| D5 | 外部验证 GSE62254 + 敏感性 + 临床关联（Step 7-8）     |
| D6 | 全部图表定稿                                   |
| D7 | 报告撰写与参考文献核对（Step 9）                      |

## 八、避坑清单

1. **别做 WSI 深度学习**：PatchGastric 等数据集单张全切片数 GB，训练需 GPU——课程若强制要求图像模态，退而求其次用"预训练 ResNet 冻结提特征 + sklearn 分类器"（300×300 小 patch CPU 可跑）。
2. **基因ID对齐**：本次下载的 HiSeqV2 行是基因 Symbol，与 GSE62254 直接可对齐；甲基化矩阵行是 cg 探针 ID，须经 probeMap 映射到基因。
3. **生存数据字段**：survival 表列为 `OS`（事件）与 `OS.time`（天数），别用错列；另有 DSS/DFI/PFI 备选终点。
4. **亚型标签**：MSI 状态在 clinicalMatrix 第 2 列 `CDE_ID_3226963`（MSS/MSI-H/MSI-L）；四亚型 EBV/GS/CIN 标签需 cBioPortal 导出，拿不到就以 MSI 预测为主线。
5. **别在结论里过度外推**：440 例样本 + 关联性分析，写成"提示/支持某结论"而不是"证明"。
6. **转置坑**：矩阵喂给 sklearn 前必须 `.T`（行=样本）。忘了转置不报错但结果全错——转置后立即 `print(X.shape)` 确认，应为约 `(350, 1500)` 而非 `(1500, 350)`。

---

## 九、写给计算机背景的概念手册

### 9.1 多组学是什么

中心法则：DNA → RNA → 蛋白质。"组学"（-ome = 全集）= 对某一层生物分子做**全量、无差别**的测量；"多组学" = 对同一批病人测多层，再联合分析。

| 组学    | 测什么           | 数据形态          | 计算机类比                  |
| ----- | ------------- | ------------- | ---------------------- |
| 基因组   | DNA 突变、拷贝数    | 突变列表 / CNV 矩阵 | 源代码                    |
| 表观基因组 | DNA 甲基化（调控开关） | CpG 位点 β 值矩阵  | 编译开关 / 宏定义：源码不变但决定哪段生效 |
| 转录组   | mRNA 表达量      | 基因×样本矩阵       | 运行时日志：哪些函数被调用、调用多少次    |
| 蛋白质组  | 蛋白丰度          | 蛋白×样本矩阵       | 实际运行的进程输出              |

为什么要多层联合：单层只看因果链一段——只看转录组发现 A 基因高表达，不知道是 DNA 扩增（基因组层）还是甲基化丢失（表观层）导致；TCGA 四亚型正是六平台联立才定义出来的。本作业的三模态 = 转录组（mRNA）+ 表观基因组（甲基化）+ 表型（临床）。

### 9.2 基因×样本矩阵

行=基因，列=病人，格=该基因在该病人肿瘤中的值。两种读法：**竖看一列** = 一个病人的完整分子画像（一个样本/观测）；**横看一行** = 一个基因在人群中的行为（差异分析、特征重要性的对象）。

矩阵不是随机数字——每一行都在讲生物学故事：CDH1（细胞黏附基因）低表达是弥漫型/GS 的特征，PIK3CA 高表达是 EBV 的特征，免疫相关基因高表达是 MSI 的特征。模型若真学到东西，它最倚重的特征必然是这些有已知故事的基因——这就是"特征重要性对照文献"一步的原理。

三模态 = 三张同构矩阵（行不同：基因 / CpG 位点 / 临床变量），列名（样本条形码）是同一批人，按列 join 后横向拼接——这就是"多模态对齐"的全部含义。

### 9.3 术语速查表

**疾病与病理**

| 术语            | 含义                                                |
| ------------- | ------------------------------------------------- |
| 胃腺癌           | 起源于胃内表面腺上皮的恶性肿瘤，占胃癌约 95%                          |
| 分化 / 分级 G1-G3 | 肿瘤细胞越像正常细胞=分化好（G1），越原始=分化差（G3），越凶险                |
| Lauren 分型     | 胃癌组织学二分法：肠型（腺体结构、老年多、预后相对好）vs 弥漫型（细胞散在浸润、年轻多、预后差） |
| 印戒细胞癌         | 弥漫型典型类型，细胞内黏液把核挤到边上，镜下像戒指                         |
| 转移            | 肿瘤扩散到其他器官                                         |
| H\&E 染色       | 病理标准染色（核蓝紫、胞质粉红），病理图像基本都用它                        |
| WSI           | 全切片病理图像，千兆像素级，本作业刻意避开                             |

**分子层**

| 术语                  | 含义                                      |
| ------------------- | --------------------------------------- |
| 基因                  | DNA 上编码一个功能产物的片段，人类约 2 万个蛋白编码基因         |
| 基因表达 / 表达量          | 基因被转录成 mRNA 的数量，反映当前活跃程度                |
| mRNA                | DNA 的临时工作拷贝，指导合成蛋白质                     |
| RNA-seq / TPM       | 测序测表达；TPM=归一化单位（每百万 reads 中的占比），跨样本可比   |
| DNA 甲基化 / CpG / β 值 | DNA 特定位置挂甲基基团通常压制表达；CpG=甲基化发生位置；β∈[0,1] |
| 突变                  | DNA 序列改变（碱基错字）；体细胞突变=后天发生在肿瘤细胞里         |
| CNV                 | 拷贝数变异：基因份数异常（正常 2 份）                    |
| 高突变负荷               | 突变数异常多（常见于 MSI），产生新抗原被免疫系统识别            |

**四个分子亚型（预测标签）**

| 亚型       | 核心特征                                      | 临床含义          |
| -------- | ----------------------------------------- | ------------- |
| EBV      | EB 病毒（疱疹病毒家族）相关，富集 PIK3CA 突变与 PD-L1/L2 扩增 | 免疫治疗潜在获益      |
| MSI      | 微卫星（基因组短重复序列）不稳定，错配修复 MLH1 沉默 → 突变累积      | 高突变负荷，免疫治疗响应好 |
| GS       | 基因组稳定型，富集弥漫型组织学，CDH1/RHOA（黏附基因）突变         | 预后偏差          |
| CIN      | 染色体不稳定型（最常见约一半），大片段丢失/获得、非整倍体             | 有靶向药机会        |
| HER2     | 生长因子受体，扩增则驱动生长                            | 曲妥珠单抗靶向治疗     |
| PD-L1/L2 | 肿瘤表面"免疫刹车"分子                              | 免疫检查点抑制剂适应症   |

**数据与临床统计**

| 术语                           | 含义                                             |
| ---------------------------- | ---------------------------------------------- |
| TCGA / STAD                  | 美国癌症基因组图谱计划 / 其中的胃腺癌队列（STomach ADenocarcinoma） |
| GEO / GSE62254               | NCBI 公共数据仓库 / ACRG 队列（韩国约 300 例，外部验证集）         |
| 队列 cohort                    | 按统一方案入组随访的一批患者                                 |
| Ensembl ID vs Symbol         | 基因的数据库唯一 ID（ENSG…，类似 GUID）vs 人类可读名（TP53）       |
| TNM 分期                       | T=侵犯深度、N=淋巴结转移、M=远处转移，组合成 I-IV 期               |
| OS / Kaplan-Meier / log-rank | 总生存时间 / 生存率随时间曲线（阶梯图）/ 曲线间差异显著性检验              |
| 卡方检验                         | 检验两个分类变量是否独立（亚型 × Lauren/分期）                   |
| 预后                           | 对疾病结局（复发/死亡风险）的预测                              |
