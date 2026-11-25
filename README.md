# HBSCM

HBSCM 将 Ma 等人（2026）的马蹄先验贝叶斯合成控制方法封装为 R 包。输入是一张长面板，每行对应一个单位在一个时期的记录。模型使用未受处理的观测拟合，并为受处理单位的完整观测历史生成反事实预测抽样。维护者为 Xiaohua Ma（mxh051218@163.com）。[论文](https://doi.org/10.1016/j.econmod.2026.107502) · [公开复现代码](https://doi.org/10.17632/4vx9dcbpvm.1)

## 安装

首次使用时，在 R 或 RStudio 中运行整段代码。它会从 GitHub 下载 v0.1.1，并按顺序安装模型依赖 `pblasso 1.0.8` 和 `HBSCM 0.1.1`。

```r
local({
  url <- "https://github.com/sun2901831024-cloud/HBSCM/releases/download/v0.1.1/HBSCM-0.1.1-release.zip"
  archive <- tempfile(fileext = ".zip")
  folder <- tempfile("HBSCM-")
  dir.create(folder)
  download.file(url, archive, mode = "wb")
  unzip(archive, exdir = folder)
  source(file.path(folder, "HBSCM", "install.R"), local = TRUE, encoding = "UTF-8")
})
library(HBSCM)
```

以后打开 R，只需运行：

```r
library(HBSCM)
```

Windows x64 配合 R 4.5.x 时使用预编译包。其他支持的 R 版本（R ≥ 4.1）和平台使用源码安装，需要与 R 匹配的 C++17、BLAS、LAPACK 和 Fortran 工具链。安装脚本会从 CRAN 补齐常规 R 依赖；指定版本的 `pblasso` 已随发布包提供。[原始依赖来源](https://doi.org/10.7910/DVN/B6SWA1)

## 从模拟面板估计 ATT

```r
library(HBSCM)
d <- simulateSA(N = 14, TT = 10, r = 2, tr.start = 7,
                p = 2, beta = c(1, 0), seed = 106)
set.seed(92)
fit <- HBSCM(d, index = c("id", "time"), Yname = "Y", Dname = "D",
             Xname = c("X1", "X2"), re = "both", r = 2,
             niter = 100, burn = 50)
effects <- heffSummary(fit, cumu = TRUE)
effects$est.avg
effects$est.eff
effPlot(effects, "eff")
```

`est.avg` 对处理后的单位–时期观测等权平均，报告后验均值以及 2.5% 和 97.5% 分位数。`est.eff` 给出逐期结果，相对时间 0 是处理前最后一期，1 是开始处理的第一期。用 `heffSummary(fit, usr.id = ...)` 查看指定单位；用 `rela.period = FALSE` 查看共同时间轴，`period_label` 给出原始时期。`fit$yct` 的每一列都是一次保留的后验预测抽样，包含观测噪声。示例用短链演示接口，正式分析应根据收敛和蒙特卡洛误差选择链长。

## 模型入口

| 函数 | 来源和用途 |
| --- | --- |
| `HBSCM()` | `function-code.R` 中的主模型，分别控制共同、单位、时间和因子部分的马蹄收缩。 |
| `bpCausalgroup()` | `e8.R` 与 `C++group1.cpp` 的分组扩展。用 `group = "列名"` 指定单位所属组；组系数和载荷沿用原 C++ 的条件均值更新。 |
| `bpCausaladp()` | `e7.R` 的自适应 Lasso 对照模型。 |
| `simulateSA()` | 生成带处理时点、真实处理效应和潜在因子的模拟面板。处理分配使用前两个载荷，生成时需 `r >= 2`。拟合时可取 `r = 0`。 |
| `heffSummary()`、`coefSummary()`、`effPlot()` | 汇总后验处理效应、系数与区间，并绘制反事实路径、处理效应或累积效应。 |

数据须有唯一的单位–时期键、数值型结局以及 0/1 处理指标。模型按指定列做完整观测筛选，剔除没有任何未处理观测的单位。保留的每一期都需要未处理观测；存在缺期的面板仍可使用。协变量按样本标准差缩放，系数对应缩放后的尺度。`re` 指定单位、时间、双向或无变截距；`Xname`、`Zname`、`Aname` 分别指定共同系数、单位变系数与时间变系数的协变量。更完整的说明可在 R 中运行 `?HBSCM`、`?heffSummary` 和 `citation("HBSCM")` 查看。

## 实证数据与示例

```r
source(system.file("examples", "smoking.R", package = "HBSCM"))
```

加州控烟示例用州名识别加州，以 1989 年为首个处理年，并用短链演示拟合流程。

| 包内原始文件 | 内容与读取方式 |
| --- | --- |
| `smoking.csv` | 加州控烟案例，使用 `read.csv()`。 |
| `repgermany.dta` | 德国统一案例，可用 `haven::read_dta()`。 |
| `data_main.RData` | 亚投行案例，包含数据框 `df.main`；可用 `load(path, envir = new.env())` 载入独立环境。 |

用 `system.file("extdata", "文件名", package = "HBSCM")` 获得包内数据的完整路径。

## 与原脚本的对应关系

- 为复现原链，主模型与分组模型保留了因子尺度更新中的原始写法。因子全局尺度更新的形状参数使用时间系数数目 `k3`，随后刷新时间块辅助量 `xxi`，因子块辅助量 `xg` 保持初始值。按理论公式重新实现此步会改变抽样链。
- 自适应 Lasso 的原脚本以自由变量 `T` 初始化因子，包内显式使用面板期数 `TT`，并按各系数块的实际长度初始化局部参数。主模型和分组模型现在也遵循调用者传入的 `ar1 = FALSE`。
- 修正了原效应摘要中反事实上下界列名的顺序，以及恰有两个处理后时期时第二个累计值缺失的问题。系数摘要读取模型实际返回的 `sigma2`，默认使用所有已保留抽样。单位和组系数汇总使用相应的原始标签。
- 固定种子的小面板核对包括主模型、分组模型及自适应 Lasso。每个模型的 15 项抽样结果与对应原脚本逐项相同。自适应 Lasso 的参照脚本显式设定 `T = TT`。该核对与实证短链示例检验软件封装；本次没有重新计算论文全部模拟表或评估链收敛。
