# 📁 项目结构 - 精简版

```
Time-Series-Data-Analysis/
│
├── 📄 Quarto 配置和样式
│   ├── _quarto.yml           ← 网站配置（关键）
│   ├── styles.css            ← 自定义样式
│   └── README.md             ← 项目说明
│
├── 📚 源代码文件（.qmd）
│   ├── index.qmd             ← 首页
│   ├── chap0_preface.qmd     ← 第 0 章：序言
│   ├── chap1.qmd             ← 第 1 章：数据初始化
│   ├── chap2.qmd             ← 第 2 章：频域分析
│   ├── chap3.qmd             ← 第 3 章：功率谱密度
│   ├── chap4.qmd             ← 第 4 章：FFT 深入
│   ├── chap5.qmd             ← 第 5 章：噪声
│   ├── chap6.qmd             ← 第 6 章：信号处理
│   ├── chap7.qmd             ← 第 7 章：时频谱
│   ├── chap8.qmd             ← 第 8 章：多通道信号
│   ├── chap10_appendix.qmd   ← 第 10 章：附录
│   └── chap_todo.qmd         ← 更新日志
│
├── 📂 资源文件夹
│   ├── Figure/               ← 所有图片
│   ├── Data/                 ← 数据文件
│   └── docs/                 ← 备份的原始文件
│
└── 📂 构建输出（自动生成）
    └── _build/html/          ← 最终的网站
```

## 🎯 工作流程

### 1. 编辑源文件
编辑 `chap*.qmd` 或 `index.qmd` 文件

### 2. 构建网站
```bash
quarto render
```

### 3. 预览
```bash
quarto preview
```

浏览器会自动打开 http://localhost:4200

### 4. 部署
推送到 GitHub，GitHub Pages 会自动构建发布

## 📝 添加新章节

1. 在根目录创建 `chapN.qmd`
2. 编辑 `_quarto.yml`，在 `sidebar.contents` 中添加：
```yaml
- text: "新章节名"
  file: chapN.qmd
```
3. 运行 `quarto render`

## 🔧 常用命令

```bash
# 完整清理和重建
rm -rf _build && quarto render

# 实时预览（自动检测文件变化）
quarto preview

# 只渲染特定文件
quarto render chap0_preface.qmd

# 检查环境
quarto check
```

## ✅ 关键配置

**_quarto.yml**:
- 定义网站标题、导航和侧边栏
- 配置主题为 `lumen`
- 指定输出目录为 `_build/html`
- 明确列出要渲染的文件

**styles.css**:
- 自定义网站颜色和样式
- 支持深色模式
- 响应式设计

**README.md**:
- 项目概述
- 快速开始指南

---

**就这么简单！** 专注于编写内容，Quarto 会处理其他的事。
